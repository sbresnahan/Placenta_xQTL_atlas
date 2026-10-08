#!/usr/bin/env python3
"""
52_gxe_scan.py — GxE cis-interaction scanner (Objective 2.1)

Model (per variant, mirroring tensorQTL's calculate_interaction_nominal):
    Y = b_g * G + b_e * E + b_gi * (G x E) + Z_cov * gamma + eps

Modes:
  cis-perm  — tier-1 discovery scan: adaptive phenotype-permutation scan
              (GTEx/tensorQTL convention) over one chromosome of an
              ancestry-specific BED. Writes one parquet row per phenotype with the
              top interaction variant, nominal/permutation/beta-approx
              p-values. Shard by chromosome (LSF array); merge with
              --mode merge.
  nominal   — full nominal interaction statistics for every cis variant x
              phenotype (optionally restricted to --phenotypes-file). Used
              for tier-2 ancestry-stratified scans and spot checks.
              Long-format bgzipped TSV output.
  merge     — concatenate per-chromosome cis-perm parquets for one
              modality x exposure, add Storey q-values on pval_beta via the
              05_qtl_mapping/compute_qvalues.R bridge (GTEx convention;
              QVALUE_RSCRIPT env var overrides the Rscript command), write
              the final parquet + sorted top table.

Genotypes are read per chromosome straight from the per-ancestry pgens with
tensorqtl.pgen.read_list (hardcalls; --dosages switches to dosages), variant
IDs are intersected across ancestries (IDs encode chr:pos:ref:alt, so the
intersection is allele-consistent), and ancestry blocks are stacked in the
manifest order — no pooled pgen is ever built. Missing hardcalls are imputed
to the within-ancestry-block variant mean. The MAF filter is applied to the
pooled (stacked) genotypes.

The exposure main effect is a model term, so the exposure must NOT also sit
in the covariate table: any covariate row named exactly like --exposure
(e.g. 'GA' when GA is the exposure) is dropped automatically.

Usage (tier 1, one ancestry/chromosome):
  python3 52_gxe_scan.py --mode cis-perm \
      --bed EAS_expression.bed.gz \
      --pgens EAS={QTL_DIR}/EAS_qtl \
      --manifest EAS_sample_manifest.tsv \
      --covariates EAS_covariates_expression.tsv \
      --exposures exposures.tsv --exposure GA \
      --chrom 21 --out EAS_expression_GA.chr21.gxe_cis.parquet

Usage (merge):
  python3 52_gxe_scan.py --mode merge \
      --chr-outputs 'pooled_expression_GA.chr*.gxe_cis.parquet' \
      --out pooled_expression_GA.gxe_cis.parquet
"""

import argparse
import glob
import os
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gxe_core


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------

def parse_pgens(pgen_args):
    """Parse --pgens ANC=prefix ... into an ordered dict."""
    out = {}
    for item in pgen_args:
        if '=' not in item:
            sys.exit(f"ERROR: --pgens entries must be ANC=prefix, got: {item}")
        anc, prefix = item.split('=', 1)
        out[anc] = prefix
    return out


def read_psam_iids(psam_path):
    """Sample IDs (IID) in pgen order. Handles #FID/IID and IID-only psams."""
    psam = pd.read_csv(psam_path, sep='\t', comment=None)
    psam.columns = [c.lstrip('#') for c in psam.columns]
    if 'IID' in psam.columns:
        return psam['IID'].astype(str).tolist()
    return psam.iloc[:, 0].astype(str).tolist()


def load_chromosome_genotypes(pgen_by_anc, samples_by_anc, chrom,
                              maf_threshold=0.01, dosages=False):
    """Read one chromosome from per-ancestry pgens and stack ancestry blocks.

    Returns
    -------
    G        : float32 ndarray [n_variants x n_samples_pooled], imputed
    var_ids  : ndarray of variant IDs (first-ancestry pvar order)
    var_pos  : ndarray of positions
    samples  : list of pooled sample IDs (block order = samples_by_anc order)
    """
    from tensorqtl import pgen as tq_pgen

    chrom_norm = str(chrom).replace('chr', '')
    reader = tq_pgen.read_dosages_list if dosages else tq_pgen.read_list

    per_anc = []
    for anc, prefix in pgen_by_anc.items():
        pvar = tq_pgen.read_pvar(f"{prefix}.pvar")
        pvar['chrom'] = pvar['chrom'].astype(str).str.replace('^chr', '', regex=True)
        pv = pvar[pvar['chrom'] == chrom_norm]
        pv['pvar_idx'] = pv.index.values  # positional index in the full pvar
        per_anc.append((anc, prefix, pv))
        log(f"    {anc}: {len(pv)} variants on chr{chrom_norm} ({os.path.basename(prefix)})")

    # Intersect variant IDs across ancestries (IDs encode chr:pos:ref:alt),
    # preserving the first ancestry's pvar order.
    anc0, _, pv0 = per_anc[0]
    keep = set(pv0['id'])
    for _, _, pv in per_anc[1:]:
        keep &= set(pv['id'])
    keep_ids = [i for i in pv0['id'] if i in keep]
    n_drop = len(pv0) - len(keep_ids)
    if n_drop:
        log(f"    variant-ID intersection dropped {n_drop} {anc0} variants "
            f"(not present in all ancestries)")
    if not keep_ids:
        sys.exit(f"ERROR: no shared variants on chr{chrom_norm}")

    blocks = []
    samples = []
    var_pos = None
    for anc, prefix, pv in per_anc:
        pv = pv.set_index('id')
        idx = pv.loc[keep_ids, 'pvar_idx'].values
        order = np.argsort(idx)  # pgenlib wants ascending variant indexes
        g = reader(f"{prefix}.pgen",
                   np.asarray(idx, dtype=np.uint32)[order])
        g = g[np.argsort(order)]  # back to keep_ids order
        g = g.astype(np.float32)
        if dosages:
            # missing dosages: anything outside [0, 2] (incl. -9) -> NaN
            g[(g < 0) | (g > 2)] = np.nan
        else:
            g[g == -9] = np.nan
        iids = read_psam_iids(f"{prefix}.psam")
        col = {s: i for i, s in enumerate(iids)}
        want = samples_by_anc[anc]
        miss = [s for s in want if s not in col]
        if miss:
            sys.exit(f"ERROR: {len(miss)} {anc} samples missing from pgen "
                     f"(e.g. {miss[:3]})")
        g = g[:, [col[s] for s in want]]
        # within-block variant-mean imputation (variants fully missing in
        # this block keep NaN and get the pooled-mean fallback below)
        if np.isnan(g).any():
            mu = np.nanmean(g, axis=1)
            ix = np.where(np.isnan(g) & np.isfinite(mu)[:, None])
            g[ix] = np.take(mu, ix[0])
        blocks.append(g)
        samples.extend(want)
        if var_pos is None:
            var_pos = pv.loc[keep_ids, 'pos'].values

    G = np.concatenate(blocks, axis=1)
    # variants fully missing in some block: pooled-mean fallback
    still = np.isnan(G).any(1)
    if still.any():
        Gsub = G[still]
        mu = np.nanmean(Gsub, axis=1)
        ix = np.where(np.isnan(Gsub) & np.isfinite(mu)[:, None])
        Gsub[ix] = np.take(mu, ix[0])
        G[still] = Gsub
    drop = ~np.isfinite(G).all(1)
    if drop.any():
        log(f"    dropping {int(drop.sum())} variants with no called genotypes "
            f"in any ancestry")
        G = G[~drop]
        keep_ids = list(np.asarray(keep_ids)[~drop])
        var_pos = var_pos[~drop]

    # pooled MAF + monomorphic filter (tensorQTL conventions)
    af = G.sum(1) / (2.0 * G.shape[1])
    maf = np.where(af > 0.5, 1.0 - af, af)
    mono = (G == G[:, [0]]).all(1)
    mask = (maf >= maf_threshold) & (~mono)
    n0 = G.shape[0]
    G = G[mask]
    keep_ids = list(np.asarray(keep_ids)[mask])
    var_pos = var_pos[mask]
    log(f"    analysis set: {G.shape[0]} variants x {G.shape[1]} samples after "
        f"MAF>={maf_threshold} + monomorphic filters ({n0 - G.shape[0]} dropped)")
    return G, np.asarray(keep_ids), np.asarray(var_pos), samples


def load_scan_inputs(args):
    """Load and align phenotypes, covariates, exposure and genotypes."""
    import tensorqtl

    manifest = pd.read_csv(args.manifest, sep='\t')
    if not {'array_id', 'ancestry'} <= set(manifest.columns):
        sys.exit("ERROR: manifest must have array_id and ancestry columns")
    anc_of = dict(zip(manifest['array_id'], manifest['ancestry']))

    exp_df = pd.read_csv(args.exposures, sep='\t', index_col=0)
    if args.exposure not in exp_df.index:
        sys.exit(f"ERROR: exposure '{args.exposure}' not in {args.exposures} "
                 f"(rows: {list(exp_df.index)})")
    exposure = exp_df.loc[args.exposure].dropna()
    log(f"  exposure '{args.exposure}': {len(exposure)} non-missing samples")

    cov = pd.read_csv(args.covariates, sep='\t', index_col=0)
    if args.exposure in cov.index:
        log(f"  dropping covariate row '{args.exposure}' (it is the exposure "
            f"main effect in the model)")
        cov = cov.drop(index=args.exposure)

    log(f"  loading phenotype BED: {args.bed}")
    pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.bed)
    pheno['__row__'] = np.arange(len(pheno))  # global row index (seed base)
    chrom_norm = str(args.chrom).replace('chr', '') if args.chrom else None
    pos_chrom = pheno_pos['chr'].astype(str).str.replace('^chr', '', regex=True)
    if chrom_norm is not None:
        m = (pos_chrom == chrom_norm).values
        pheno = pheno[m]
        pheno_pos = pheno_pos[m]
        log(f"  chr{chrom_norm}: {len(pheno)} phenotypes")
    if len(pheno) == 0:
        sys.exit(f"ERROR: no phenotypes on chr{chrom_norm}")

    # Analysis samples: BED columns with non-missing exposure and covariates,
    # grouped by ancestry. Primary Module-07 jobs pass exactly one ancestry.
    pgen_by_anc = parse_pgens(args.pgens)
    cov_samples = set(cov.columns)
    exp_samples = set(exposure.index)
    samples = []
    samples_by_anc = {}
    for anc in pgen_by_anc:
        block = [s for s in pheno.columns
                 if s in cov_samples and s in exp_samples
                 and anc_of.get(s) == anc]
        samples_by_anc[anc] = block
        samples.extend(block)
    if len(samples) < 20:
        sys.exit(f"ERROR: too few analysis samples ({len(samples)})")
    log(f"  analysis samples: {len(samples)} "
        f"({', '.join(f'{a}={len(samples_by_anc[a])}' for a in samples_by_anc)})")

    pheno = pheno[samples + ['__row__']]
    cov = cov[samples]
    if cov.isna().any().any():
        sys.exit("ERROR: NaN in covariate table after sample alignment")
    exposure = exposure[samples]

    G, var_ids, var_pos, g_samples = load_chromosome_genotypes(
        pgen_by_anc, samples_by_anc, args.chrom,
        maf_threshold=args.maf_threshold, dosages=args.dosages)
    assert g_samples == samples, "genotype/phenotype sample order mismatch"

    return dict(pheno=pheno, pheno_pos=pheno_pos, cov=cov, exposure=exposure,
                G=G, var_ids=var_ids, var_pos=var_pos, samples=samples)


# ---------------------------------------------------------------------------
# Scans
# ---------------------------------------------------------------------------

def cis_slice(var_pos, start, end, window):
    """[lo, hi) variant slice for [start-window, end+window] (sorted pos)."""
    lo = np.searchsorted(var_pos, start - window, side='left')
    hi = np.searchsorted(var_pos, end + window, side='right')
    return lo, hi


def run_cis_perm(args):
    data = load_scan_inputs(args)
    device = gxe_core.get_device()
    log(f"  device: {device}")

    G_t = torch_tensor(data['G'], device)
    Q_t, dof_base = gxe_core.make_residualizer(data['cov'].T.values, device)
    dof = dof_base - 2  # one exposure -> b_e and b_gi
    log(f"  n={len(data['samples'])}, k={data['cov'].shape[0]} covariates, "
        f"interaction dof={dof}")
    if dof < 10:
        sys.exit(f"ERROR: interaction dof too small ({dof})")

    e_t = torch_tensor(data['exposure'].values.astype(np.float32), device)
    g0_t, e0_t, ge0_t = gxe_core.prepare_interaction_terms(G_t, e_t, Q_t)
    af_t, ma_samples_t, ma_count_t = gxe_core.allele_stats(G_t)
    af = af_t.cpu().numpy()
    ma_samples = ma_samples_t.cpu().numpy()
    ma_count = ma_count_t.cpu().numpy()

    var_ids, var_pos = data['var_ids'], data['var_pos']
    blocks = tuple(int(b) for b in args.perm_blocks)
    log(f"  adaptive permutation blocks: {blocks} (cumulative "
        f"{np.cumsum(blocks).tolist()}), stop_p={args.stop_p}")

    rows = []
    pheno, pheno_pos = data['pheno'], data['pheno_pos']
    sample_cols = data['samples']
    for k, (pid, prow) in enumerate(pheno.iterrows(), 1):
        p_raw_t = torch_tensor(prow[sample_cols].values.astype(np.float32),
                               device)
        p0_t = gxe_core.prepare_phenotype(p_raw_t, Q_t)
        start = int(pheno_pos.loc[pid, 'start']) if 'start' in pheno_pos.columns \
            else int(pheno_pos.loc[pid, 'pos'])
        end = int(pheno_pos.loc[pid, 'end']) if 'end' in pheno_pos.columns else start
        lo, hi = cis_slice(var_pos, start, end, args.cis_window)
        if hi - lo < 1:
            continue
        sl = slice(lo, hi)
        Xt_w, Xinv_w, ok_w = gxe_core.build_window_design(
            g0_t[sl], e0_t, ge0_t[sl])
        b, b_se, tstat = gxe_core.window_fit(Xt_w, Xinv_w, p0_t, dof)
        t_gi = tstat[:, 2].cpu().numpy()
        r2 = gxe_core.t_to_r2(t_gi, dof)
        r2 = np.where(np.isfinite(r2), r2, -1.0)
        ix = int(np.argmax(r2))
        if r2[ix] < 0:
            continue  # every window variant masked (singular designs)
        r2_nom_max = float(r2[ix])

        rng = np.random.default_rng(args.seed + int(prow['__row__']))
        r2_perm, nperm, pval_perm = gxe_core.adaptive_permutation_scan(
            Xt_w, Xinv_w, p_raw_t, Q_t, dof, r2_nom_max, rng,
            blocks=blocks, stop_p=args.stop_p)
        beta = gxe_core.beta_approx_pval(r2_perm, r2_nom_max, dof)

        gv = lo + ix
        b_np = b[ix].cpu().numpy()
        bse_np = b_se[ix].cpu().numpy()
        t_np = tstat[ix].cpu().numpy()
        rows.append(dict(
            phenotype_id=pid,
            num_var=hi - lo,
            variant_id=var_ids[gv],
            start_distance=int(var_pos[gv] - start),
            end_distance=int(var_pos[gv] - end),
            af=float(af[gv]),
            ma_samples=int(ma_samples[gv]),
            ma_count=float(ma_count[gv]),
            b_g=float(b_np[0]), b_g_se=float(bse_np[0]),
            b_e=float(b_np[1]), b_e_se=float(bse_np[1]),
            b_gi=float(b_np[2]), b_gi_se=float(bse_np[2]),
            tstat_gi=float(t_np[2]),
            pval_g=float(gxe_core.r2_to_pval(gxe_core.t_to_r2(t_np[0], dof), dof)),
            pval_e=float(gxe_core.r2_to_pval(gxe_core.t_to_r2(t_np[1], dof), dof)),
            pval_nominal=float(gxe_core.r2_to_pval(r2_nom_max, dof)),
            dof=int(dof),
            nperm=int(nperm),
            pval_perm=float(pval_perm),
            **beta,
        ))
        if k % 2000 == 0:
            log(f"    {k}/{len(pheno)} phenotypes")

    out = pd.DataFrame(rows).set_index('phenotype_id')
    out.to_parquet(args.out)
    log(f"  wrote {args.out} ({len(out)} phenotypes)")
    sig = (out['pval_perm'] <= 0.05).sum() if len(out) else 0
    log(f"  pval_perm<=0.05: {sig}")


def run_nominal(args):
    data = load_scan_inputs(args)
    device = gxe_core.get_device()
    log(f"  device: {device}")

    keep_pheno = None
    if args.phenotypes_file:
        keep_pheno = set(pd.read_csv(args.phenotypes_file, header=None)[0])
        log(f"  restricting to {len(keep_pheno)} phenotypes from "
            f"{args.phenotypes_file}")

    G_t = torch_tensor(data['G'], device)
    Q_t, dof_base = gxe_core.make_residualizer(data['cov'].T.values, device)
    dof = dof_base - 2
    e_t = torch_tensor(data['exposure'].values.astype(np.float32), device)
    g0_t, e0_t, ge0_t = gxe_core.prepare_interaction_terms(G_t, e_t, Q_t)
    af_t, ma_samples_t, ma_count_t = gxe_core.allele_stats(G_t)
    af = af_t.cpu().numpy()
    ma_samples = ma_samples_t.cpu().numpy()
    ma_count = ma_count_t.cpu().numpy()
    var_ids, var_pos = data['var_ids'], data['var_pos']

    pheno, pheno_pos = data['pheno'], data['pheno_pos']
    sample_cols = data['samples']
    wrote_header = False
    n_rows = 0
    with open_out(args.out) as fh:
        for pid, prow in pheno.iterrows():
            if keep_pheno is not None and pid not in keep_pheno:
                continue
            p0_t = gxe_core.prepare_phenotype(
                torch_tensor(prow[sample_cols].values.astype(np.float32), device), Q_t)
            start = int(pheno_pos.loc[pid, 'start']) if 'start' in pheno_pos.columns \
                else int(pheno_pos.loc[pid, 'pos'])
            end = int(pheno_pos.loc[pid, 'end']) if 'end' in pheno_pos.columns else start
            lo, hi = cis_slice(var_pos, start, end, args.cis_window)
            if hi - lo < 1:
                continue
            sl = slice(lo, hi)
            Xt_w, Xinv_w, ok_w = gxe_core.build_window_design(
                g0_t[sl], e0_t, ge0_t[sl])
            b, b_se, tstat = gxe_core.window_fit(Xt_w, Xinv_w, p0_t, dof)
            t_gi = tstat[:, 2].cpu().numpy()
            pval = gxe_core.r2_to_pval(gxe_core.t_to_r2(t_gi, dof), dof)
            df = pd.DataFrame(dict(
                phenotype_id=pid,
                variant_id=var_ids[sl],
                tss_distance=(var_pos[sl] - start).astype(int),
                af=af[sl], ma_samples=ma_samples[sl].astype(int),
                ma_count=ma_count[sl],
                b_g=b[:, 0].cpu().numpy(), b_g_se=b_se[:, 0].cpu().numpy(),
                b_e=b[:, 1].cpu().numpy(), b_e_se=b_se[:, 1].cpu().numpy(),
                b_gi=b[:, 2].cpu().numpy(), b_gi_se=b_se[:, 2].cpu().numpy(),
                tstat_gi=t_gi, pval_nominal=pval, dof=int(dof),
            ))
            df = df[np.isfinite(df['pval_nominal'])]
            df.to_csv(fh, sep='\t', index=False, header=not wrote_header,
                      float_format='%.6g')
            wrote_header = True
            n_rows += len(df)
    log(f"  wrote {args.out} ({n_rows} variant x phenotype rows)")


def run_merge(args):
    """Concatenate per-chromosome cis-perm parquets + Storey q-values."""
    paths = sorted(glob.glob(args.chr_outputs))
    if not paths:
        sys.exit(f"ERROR: no chromosome outputs match: {args.chr_outputs}")
    log(f"  merging {len(paths)} chromosome outputs")
    df = pd.concat([pd.read_parquet(p) for p in paths])
    log(f"  {len(df)} phenotypes merged")

    # Storey q-values on pval_beta via the module-05 R bridge (GTEx
    # convention); non-finite pval_beta -> 1 for the pi0 fit, q forced to 1
    # (mirrors 27_run_tensorqtl.py).
    pval_beta = pd.to_numeric(df['pval_beta'], errors='coerce').to_numpy(float)
    nonfinite = ~np.isfinite(pval_beta)
    if nonfinite.any():
        log(f"  WARNING: {nonfinite.sum()} non-finite pval_beta; p=1 for the "
            f"Storey fit and qval=1 for those rows")
    q_in = pval_beta.copy()
    q_in[nonfinite] = 1.0
    rscript = os.environ.get('QVALUE_RSCRIPT', args.rscript)
    bridge = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..', '05_qtl_mapping', 'compute_qvalues.R')
    bridge = os.path.normpath(bridge)
    if not os.path.exists(bridge):
        sys.exit(f"ERROR: compute_qvalues.R bridge not found at {bridge}")
    with tempfile.TemporaryDirectory() as tmpd:
        in_tsv = os.path.join(tmpd, 'pval_beta.tsv')
        out_tsv = os.path.join(tmpd, 'qval.tsv')
        wrapped_bridge = os.path.join(tmpd, 'compute_qvalues.module07.R')
        pd.DataFrame({'pval_beta': q_in}).to_csv(in_tsv, sep='\t', index=False)

        # Module-07 contract: set the required package library from inside R,
        # not via R_LIBS_* in the launching shell/container. The temporary
        # bridge also masks ordinary .libPaths() calls so the required library
        # stays first if the Module-05 script changes its library paths.
        r_lib = args.r_package_lib.replace('\\', '\\\\').replace('"', '\\"')
        preamble = '\n'.join([
            f'.module07_r_lib <- "{r_lib}"',
            'if (!dir.exists(.module07_r_lib)) stop(',
            '    "Required Module-07 R package library does not exist: ",',
            '    .module07_r_lib',
            ')',
            '.module07_base_libPaths <- base::.libPaths',
            '.module07_base_libPaths(c(.module07_r_lib, '
            '.module07_base_libPaths()))',
            '.libPaths <- function(new) {',
            '    if (missing(new)) return(.module07_base_libPaths())',
            '    .module07_base_libPaths(c(.module07_r_lib, new))',
            '}',
            'if (normalizePath(.libPaths()[1], mustWork=TRUE) != '
            'normalizePath(.module07_r_lib, mustWork=TRUE)) stop(',
            '    "Failed to prepend required Module-07 R package library"',
            ')',
            '',
        ])
        with open(bridge) as src, open(wrapped_bridge, 'w') as dst:
            dst.write(preamble)
            dst.write(src.read())

        cmd = shlex.split(rscript) + [wrapped_bridge, in_tsv, out_tsv]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not os.path.exists(out_tsv):
            sys.exit("ERROR: Storey q-value computation failed.\n"
                     f"  stderr tail: {proc.stderr[-1500:]}")
        for line in proc.stdout.splitlines():
            if line.strip():
                log(f"  {line.strip()}")
        qdf = pd.read_csv(out_tsv, sep='\t')
        if len(qdf) != len(df):
            sys.exit("ERROR: compute_qvalues.R output row count mismatch")
        df['qval'] = qdf['qval'].values
        if nonfinite.any():
            df.loc[nonfinite, 'qval'] = 1.0

    df.to_parquet(args.out)
    log(f"  wrote {args.out} ({len(df)} phenotypes; "
        f"{(df['qval'] <= 0.05).sum()} at qval<=0.05)")
    top = df.reset_index().sort_values('pval_nominal')
    top_path = args.out.replace('.parquet', '_top.tsv')
    top.to_csv(top_path, sep='\t', index=False)
    log(f"  wrote {top_path}")


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def torch_tensor(x, device):
    import torch
    return torch.tensor(np.asarray(x), dtype=torch.float32, device=device)


class open_out:
    """Write plain or gzip-by-extension."""
    def __init__(self, path):
        self.path = path
    def __enter__(self):
        if self.path.endswith('.gz'):
            import gzip
            self.fh = gzip.open(self.path, 'wt')
        else:
            self.fh = open(self.path, 'w')
        return self.fh
    def __exit__(self, *exc):
        self.fh.close()
        return False


def main():
    p = argparse.ArgumentParser(description="GxE cis-interaction scanner")
    p.add_argument('--mode', required=True,
                   choices=['cis-perm', 'nominal', 'merge'])
    p.add_argument('--bed', help='pooled or per-ancestry phenotype BED(.gz)')
    p.add_argument('--pgens', nargs='+',
                   help='ANC=pgen_prefix entries, block order = this order')
    p.add_argument('--manifest', help='ancestry-specific sample manifest TSV')
    p.add_argument('--covariates', help='covariate TSV (rows = covariates)')
    p.add_argument('--exposures', help='exposures.tsv (rows = exposure_id)')
    p.add_argument('--exposure', help='exposure_id to scan')
    p.add_argument('--chrom', help='chromosome (required for cis-perm)')
    p.add_argument('--phenotypes-file',
                   help='optional one-column phenotype list (nominal mode)')
    p.add_argument('--cis-window', type=int, default=1000000)
    p.add_argument('--maf-threshold', type=float, default=0.01)
    p.add_argument('--perm-blocks', type=int, nargs='+',
                   default=[100, 400, 500, 9000],
                   help='adaptive permutation blocks (default: 100 400 500 '
                        '9000 = 100/500/1000/10000 cumulative)')
    p.add_argument('--stop-p', type=float, default=0.10,
                   help='early-stop threshold on the running empirical p '
                        '(default: 0.10)')
    p.add_argument('--seed', type=int, default=12345)
    p.add_argument('--dosages', action='store_true',
                   help='read dosages instead of hardcalls')
    p.add_argument('--chr-outputs', help='glob of per-chromosome parquets '
                                         '(merge mode)')
    p.add_argument('--rscript', default='Rscript',
                   help='Rscript command for the qvalue bridge (merge mode; '
                        'QVALUE_RSCRIPT env var wins)')
    p.add_argument(
        '--r-package-lib',
        default=os.environ.get(
            'R_PACKAGE_LIB',
            '/rsrch5/home/epi/bhattacharya_lab/software/'
            'R_package_library/ubuntu/4.3.1',
        ),
        help='R package library prepended inside R before the qvalue bridge '
             'runs (R_PACKAGE_LIB env var overrides the default)',
    )
    p.add_argument('--out', required=True)
    args = p.parse_args()

    if args.mode == 'merge':
        if not args.chr_outputs:
            sys.exit("ERROR: --chr-outputs required in merge mode")
        run_merge(args)
        return

    for req in ['bed', 'pgens', 'manifest', 'covariates', 'exposures',
                'exposure']:
        if getattr(args, req) is None:
            sys.exit(f"ERROR: --{req} required in {args.mode} mode")
    if not args.chrom:
        sys.exit("ERROR: --chrom required in cis-perm and nominal modes "
                 "(shard chromosomes across jobs; concatenate nominal TSVs "
                 "downstream)")
    if args.mode == 'cis-perm':
        run_cis_perm(args)
    else:
        run_nominal(args)


if __name__ == '__main__':
    main()
