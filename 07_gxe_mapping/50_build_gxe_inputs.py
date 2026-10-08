#!/usr/bin/env python3
"""
50_build_gxe_inputs.py — Build ancestry-specific GxE inputs (Objective 2.1)

Primary tier-1 discovery is ancestry-stratified.  This script prepares the
Module-05 products for those scans without forcing a cross-ancestry phenotype
intersection:

  1. inputs/{ANC}_{MOD}.bed.gz — all autosomal phenotypes available in each
     ancestry, z-scored within ancestry.  Splicing/IR ancestry-local suffixes
     are canonicalized to stable genomic event IDs so downstream results can be
     matched across ancestries without shrinking the discovery feature space.
  2. inputs/{ANC}_sample_manifest.tsv — one retained Module-05 individual per
     array_id, using replicate_collapsed/reports/ancestry_map_collapsed.tsv as
     the retained-run authority.
  3. inputs/{ANC}_metadata.tsv — rnaseq-level metadata for those final
     individuals, retaining same-cohort technical replicate rows only for HCP
     QC aggregation.
  4. inputs/exposures.tsv — one common exposure scale across all ancestries.
     Continuous exposures are transformed/range-filtered then z-scored across
     the combined retained individuals so ancestry-specific b_int estimates
     remain on the same exposure-SD scale for downstream comparison/meta.
  5. inputs/{ANC}_covariates_base.tsv — ancestry-specific genotype PCs, sex,
     cohort, GA, and harmonized cell-type proportions.  Maternal cell fraction
     is excluded and the dominant placental cell type is dropped as the
     compositional reference.  Correlation pruning is done within ancestry.
  6. --finalize-covariates --ancestry ANC --modality MOD --hcp-file F appends
     ancestry/modality-specific HCP factors and writes
     inputs/{ANC}_covariates_{MOD}.tsv.

No cross-ancestry phenotype or variant intersection is performed here.  The
primary scanner operates on each ancestry independently; cross-ancestry
matching/synthesis happens after the ancestry-specific scan merges.

Metadata provenance:
  Individual-level metadata are taken from the Module-05 final
  {ANC}_metadata.tsv files and restricted to the retained RNA run in
  replicate_collapsed/reports/ancestry_map_collapsed.tsv.  This keeps
  exposures/demographics aligned to the same representative run used by the
  replicate-collapse workflow.
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

MODALITIES = ("expression isoforms isoform_expression splicing "
              "intron_retention alt_TSS alt_polyA RNA_editing stability")
HERE = os.path.dirname(os.path.abspath(__file__))


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def read_bed(path):
    df = pd.read_csv(path, sep="\t", dtype={"#chr": str})
    return df


def canonicalize_phenotype_ids(df, modality, ancestry):
    """Replace ancestry-local splicing/IR IDs with stable genomic event IDs.

    Module 03 harmonizes splicing and intron retention within each ancestry.
    Consequently, LeafCutter meta-cluster numbers (clu_N) and MAJIQ IR
    unified indices are ancestry-local and cannot be intersected literally
    across ancestries.  The underlying junction/event coordinates embedded in
    phenotype_id are stable and are used here for the pooled Tier-1 identity.
    Other modalities already have stable IDs and are returned unchanged.
    """
    if modality == "splicing":
        # assemble_bed.py emits:
        #   {gene_id}__{chrom}_{start}_{end}_clu_{N}_{strand}
        # Drop only the ancestry-local cluster number.
        canonical = df["phenotype_id"].astype(str).str.replace(
            r"_clu_\d+_([+-])$", r"_\1", regex=True
        )
    elif modality == "intron_retention":
        # assemble_bed.py emits:
        #   {tss_gene_id}__IR_{gene_base}_{seqid}_{start}_{end}_{strand}_{idx}
        # Drop only the ancestry-local sequential unified index.
        canonical = df["phenotype_id"].astype(str).str.replace(
            r"_([+-])_\d+$", r"_\1", regex=True
        )
    else:
        return df

    changed = int((canonical != df["phenotype_id"].astype(str)).sum())
    if changed == 0 and len(df):
        sys.exit(
            f"ERROR: {modality} phenotype IDs for {ancestry} do not match "
            "the expected Module-03 ancestry-local ID format; refusing to "
            "canonicalize potentially non-equivalent phenotypes"
        )

    dup = canonical.duplicated(keep=False)
    if dup.any():
        examples = canonical[dup].head(5).tolist()
        sys.exit(
            f"ERROR: {modality} canonicalization creates duplicate phenotype "
            f"IDs for {ancestry}: {examples}. Investigate upstream "
            "harmonization before cross-ancestry matching."
        )

    df = df.copy()
    df["phenotype_id"] = canonical
    log(f"    {ancestry}: canonicalized {changed} {modality} phenotype IDs "
        "for cross-ancestry matching")
    return df


def zscore_rows(x):
    """Row-wise z-score; rows with sd < 1e-8 return NaN (dropped upstream)."""
    mu = x.mean(1)
    sd = x.std(1, ddof=0)
    sd = np.where(sd < 1e-8, np.nan, sd)
    return (x - mu[:, None]) / sd[:, None]


def build_ancestry_bed(modality, ancestry, qtl_dir, out_dir, meta_ids):
    """Build one ancestry-specific autosomal, canonicalized, z-scored BED."""
    path = os.path.join(qtl_dir, f"{ancestry}_{modality}.bed.gz")
    if not os.path.exists(path):
        sys.exit(f"ERROR: BED not found: {path}")
    df = read_bed(path)
    df["#chr"] = df["#chr"].astype(str)

    chrom_num = pd.to_numeric(
        df["#chr"].str.replace(r"^chr", "", case=False, regex=True),
        errors="coerce",
    )
    autosomal = chrom_num.between(1, 22) & (chrom_num % 1 == 0)
    n_nonauto = int((~autosomal).sum())
    if n_nonauto:
        log(f"    {ancestry}: dropping {n_nonauto} non-autosomal phenotypes")
    df = df.loc[autosomal].copy()
    df = canonicalize_phenotype_ids(df, modality, ancestry)

    keep_cols = [c for c in df.columns[4:] if c in meta_ids]
    n_drop = len(df.columns) - 4 - len(keep_cols)
    if n_drop:
        log(f"    {ancestry}: dropping {n_drop} BED samples without retained metadata")
    if not keep_cols:
        sys.exit(f"ERROR: no retained metadata samples in {path}")
    pos = df[["#chr", "start", "end", "phenotype_id"]].copy()
    vals = df[keep_cols].to_numpy(dtype=np.float64)
    z = zscore_rows(vals)
    zdf = pd.DataFrame(z, index=df.index, columns=keep_cols)
    bad = zdf.isna().any(axis=1)
    n_bad = int(bad.sum())
    if n_bad:
        log(f"    {ancestry}: dropping {n_bad} phenotypes with ~zero within-ancestry variance")
        pos = pos.loc[~bad].copy()
        zdf = zdf.loc[~bad].copy()

    out_df = pd.concat([pos.reset_index(drop=True), zdf.reset_index(drop=True)], axis=1)
    out_df = out_df.sort_values(
        ["#chr", "start"],
        key=lambda x: pd.to_numeric(
            x.astype(str).str.replace(r"^chr", "", case=False, regex=True),
            errors="raise",
        ) if x.name == "#chr" else x,
    )
    plain = os.path.join(out_dir, f"{ancestry}_{modality}.bed")
    out_df.to_csv(plain, sep="	", index=False, float_format="%.6g")
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-p", "bed", "-f", plain + ".gz"], check=True)
    log(f"    wrote {plain}.gz (+ .tbi): {len(out_df)} phenotypes x {len(keep_cols)} samples")
    return keep_cols

def load_module05_metadata(qtl_dir, ancestries, collapsed_ancestry_map,
                           exposure_config):
    """Load metadata with the same provenance as the Module-05 QTL inputs.

    Module 05 writes {ANC}_metadata.tsv after the final RNA/DNA intersection.
    Those files can retain multiple rnaseq_id rows for one array_id. For
    individual-level fields (exposures, ancestry, cohort, sex), select the
    retained RNA run from collapse_replicates.py's ancestry_map_collapsed.tsv.

    Return:
      meta_all: all final Module-05 rnaseq rows (used by HCP QC aggregation)
      meta_one: one retained-run row per array_id (used by GxE metadata)
    """
    meta_blocks = []
    for anc in ancestries:
        path = os.path.join(qtl_dir, f"{anc}_metadata.tsv")
        if not os.path.exists(path):
            sys.exit(f"ERROR: Module-05 metadata not found: {path}")
        df = pd.read_csv(path, sep="\t")
        missing = {"rnaseq_id", "array_id"} - set(df.columns)
        if missing:
            sys.exit(f"ERROR: {path} missing required columns: {sorted(missing)}")
        if "ancestry" in df.columns:
            bad = df[
                df["ancestry"].notna()
                & (df["ancestry"].astype(str) != str(anc))
            ]
            if len(bad):
                sys.exit(
                    f"ERROR: {path} contains {len(bad)} row(s) whose ancestry "
                    f"does not match file stratum {anc}"
                )
        df = df.copy()
        df["ancestry"] = anc
        meta_blocks.append(df)
        log(f"  Module-05 metadata {anc}: {len(df)} rnaseq rows, "
            f"{df['array_id'].nunique()} array_ids")

    meta_all = pd.concat(meta_blocks, ignore_index=True)
    if meta_all["rnaseq_id"].duplicated().any():
        bad = meta_all.loc[meta_all["rnaseq_id"].duplicated(keep=False),
                           "rnaseq_id"].astype(str).unique()[:5]
        sys.exit(
            "ERROR: rnaseq_id appears more than once across Module-05 metadata "
            f"files (examples: {list(bad)})"
        )
    id_counts = meta_all.groupby("rnaseq_id")["array_id"].nunique()
    if (id_counts > 1).any():
        bad = id_counts[id_counts > 1].index.astype(str).tolist()[:5]
        sys.exit(
            "ERROR: rnaseq_id maps to multiple array_ids in Module-05 metadata "
            f"(examples: {bad})"
        )

    if not os.path.exists(collapsed_ancestry_map):
        sys.exit(
            "ERROR: collapsed ancestry map not found: "
            f"{collapsed_ancestry_map}\n"
            "Run Module-05 collapse_replicates.py first or pass "
            "--collapsed-ancestry-map explicitly."
        )
    amap = pd.read_csv(collapsed_ancestry_map, sep="\t")
    required = {"sample_id", "assigned_ancestry", "cohort"}
    missing = required - set(amap.columns)
    if missing:
        sys.exit(
            f"ERROR: collapsed ancestry map missing columns {sorted(missing)}: "
            f"{collapsed_ancestry_map}"
        )
    if amap["sample_id"].duplicated().any():
        bad = amap.loc[amap["sample_id"].duplicated(keep=False),
                       "sample_id"].astype(str).unique()[:5]
        sys.exit(
            "ERROR: collapsed ancestry map has duplicate sample_id rows "
            f"(examples: {list(bad)})"
        )

    retained = amap[["sample_id", "assigned_ancestry", "cohort"]].rename(
        columns={"sample_id": "rnaseq_id",
                 "assigned_ancestry": "collapsed_ancestry",
                 "cohort": "collapsed_cohort"}
    )
    meta_retained = meta_all.merge(retained, on="rnaseq_id", how="inner",
                                   validate="one_to_one")
    if meta_retained.empty:
        sys.exit(
            "ERROR: no overlap between Module-05 metadata rnaseq_id values and "
            f"collapsed ancestry map: {collapsed_ancestry_map}"
        )

    # The collapsed map is the authority for the retained run's ancestry/cohort.
    anc_mismatch = (
        meta_retained["ancestry"].astype(str)
        != meta_retained["collapsed_ancestry"].astype(str)
    )
    if anc_mismatch.any():
        ex = meta_retained.loc[anc_mismatch,
                               ["rnaseq_id", "array_id", "ancestry",
                                "collapsed_ancestry"]].head(5)
        sys.exit(
            "ERROR: Module-05 metadata ancestry disagrees with collapsed "
            f"ancestry map for {int(anc_mismatch.sum())} retained run(s):\n"
            + ex.to_string(index=False)
        )
    meta_retained["ancestry"] = meta_retained["collapsed_ancestry"]

    if "cohort" in meta_retained.columns:
        cohort_mismatch = (
            meta_retained["cohort"].notna()
            & meta_retained["collapsed_cohort"].notna()
            & (meta_retained["cohort"].astype(str)
               != meta_retained["collapsed_cohort"].astype(str))
        )
        if cohort_mismatch.any():
            ex = meta_retained.loc[cohort_mismatch,
                                   ["rnaseq_id", "array_id", "cohort",
                                    "collapsed_cohort"]].head(5)
            sys.exit(
                "ERROR: Module-05 metadata cohort disagrees with collapsed "
                f"ancestry map for {int(cohort_mismatch.sum())} retained "
                "run(s):\n" + ex.to_string(index=False)
            )
    meta_retained["cohort"] = meta_retained["collapsed_cohort"]
    meta_retained = meta_retained.drop(
        columns=["collapsed_ancestry", "collapsed_cohort"]
    )

    # There must now be exactly one retained RNA run per individual.
    dup_array = meta_retained["array_id"].duplicated(keep=False)
    if dup_array.any():
        ex = meta_retained.loc[dup_array,
                               ["rnaseq_id", "array_id", "ancestry",
                                "cohort"]].head(10)
        sys.exit(
            "ERROR: collapsed ancestry map leaves >1 retained rnaseq_id for "
            "the same array_id; this violates the Module-05 replicate-collapse "
            "contract. Examples:\n" + ex.to_string(index=False)
        )

    # Exposure/demographic values come from the retained run. Still verify that
    # technical-replicate metadata agree for the individual-level fields Module
    # 07 may use; disagreement is a provenance error, not something to average.
    cfg = pd.read_csv(exposure_config, sep="\t")
    if "enabled" in cfg.columns:
        cfg = cfg[cfg["enabled"] == 1]
    exposure_cols = set(cfg["column"].dropna().astype(str)) \
        if "column" in cfg.columns else set()
    check_cols = [c for c in sorted(
        exposure_cols | {"ancestry", "sex", "GA", "ppBMI", "gdm", "ogtt"}
    ) if c in meta_all.columns]
    for c in check_cols:
        nuniq = meta_all.groupby("array_id")[c].nunique(dropna=True)
        bad = nuniq[nuniq > 1]
        if len(bad):
            sys.exit(
                f"ERROR: {len(bad)} final Module-05 array_id(s) have "
                f"discordant individual-level metadata '{c}' across rnaseq "
                "rows; refusing to choose/average values. Examples: "
                f"{bad.index.astype(str).tolist()[:5]}"
            )

    meta_one = meta_retained.set_index("array_id", drop=False)
    log(
        f"  retained-run metadata: {len(meta_one)} array_ids from "
        f"{collapsed_ancestry_map}"
    )
    return meta_all, meta_one


def build_exposures(config_path, meta_collapsed, samples, out_dir):
    cfg = pd.read_csv(config_path, sep="\t")
    cfg = cfg[cfg["enabled"] == 1]
    rows, summaries = {}, []
    for _, r in cfg.iterrows():
        eid, col = r["exposure_id"], r["column"]
        if col not in meta_collapsed.columns:
            log(f"  WARNING: exposure '{eid}' column '{col}' not in metadata "
                f"— skipping (enable only after the column is added)")
            continue
        v = pd.to_numeric(meta_collapsed[col], errors="coerce")
        lo, hi = r.get("min"), r.get("max")
        if pd.notna(lo):
            v = v.where(v >= float(lo))
        if pd.notna(hi):
            v = v.where(v <= float(hi))
        if str(r.get("transform", "none")).lower() == "log":
            v = np.log(v)
        v = v.reindex(samples)
        uniq = set(v.dropna().unique())
        binary = uniq <= {0.0, 1.0}
        if binary:
            out_v = v  # scanner centers internally
        else:
            out_v = (v - v.mean()) / v.std(ddof=0)  # pooled z-score
        rows[eid] = out_v
        summaries.append(dict(exposure_id=eid, column=col, binary=binary,
                              n=int(v.notna().sum()), n_missing=int(v.isna().sum()),
                              mean=float(v.mean()), sd=float(v.std(ddof=0))))
        log(f"  exposure '{eid}' ({col}): n={v.notna().sum()}, "
            f"{'binary' if binary else 'continuous, z-scored'}")
    if not rows:
        sys.exit("ERROR: no enabled exposures could be assembled")
    exp_df = pd.DataFrame(rows).T
    exp_df.index.name = "exposure_id"
    exp_df.to_csv(os.path.join(out_dir, "exposures.tsv"), sep="\t",
                  float_format="%.6g")
    pd.DataFrame(summaries).to_csv(
        os.path.join(out_dir, "exposure_manifest.tsv"), sep="\t", index=False)
    log(f"  wrote exposures.tsv ({len(rows)} exposures x {len(samples)} samples)")


def prune_correlated(cov, priority, threshold=0.9):
    """Module-05 iterative max-|r| pruning. cov: rows = covariates."""
    dropped = []
    cur = cov.copy()
    while cur.shape[0] > 1:
        c = np.abs(np.corrcoef(cur.values))
        np.fill_diagonal(c, 0.0)
        if not np.isfinite(c.max()) or c.max() <= threshold:
            break
        i, j = np.unravel_index(np.nanargmax(c), c.shape)
        n1, n2 = cur.index[i], cur.index[j]
        drop, keep = (n2, n1) if priority[n1] <= priority[n2] else (n1, n2)
        dropped.append(dict(covariate=drop, correlated_with=keep,
                            r_value=round(float(c[i, j]), 4)))
        cur = cur.drop(index=drop)
    return cur, dropped


def build_base_covariates(samples, meta_collapsed, qtl_dir, pcair_dir,
                           ancestry, out_dir, cor_threshold=0.9):
    """Build base covariates within one ancestry."""
    blocks, priority = [], {}

    sel = os.path.join(qtl_dir, f"{ancestry}_selected_pcs.txt")
    pcs_path = os.path.join(pcair_dir, f"{ancestry}_genotype_pcs.tsv")
    if not (os.path.exists(sel) and os.path.exists(pcs_path)):
        sys.exit(f"ERROR: genotype PCs missing for {ancestry}: {sel} / {pcs_path}")
    pc_names = [line.strip() for line in open(sel) if line.strip()]
    pcs = pd.read_csv(pcs_path, sep="	").set_index("sample_id")
    missing_pc_samples = [x for x in samples if x not in pcs.index]
    if missing_pc_samples:
        sys.exit(f"ERROR: {len(missing_pc_samples)} {ancestry} samples missing genotype PCs "
                 f"(examples: {missing_pc_samples[:5]})")
    if pc_names:
        block = pcs.loc[samples, pc_names].T
        block.index = [f"{ancestry}_{x}" for x in block.index]
        blocks.append(block)
        for j, name in enumerate(block.index):
            priority[name] = (1, j)
    log(f"  {ancestry} PCs: {len(pc_names)}")

    sex_map = {"M": 0, "F": 1, "m": 0, "f": 1, "Male": 0, "Female": 1}
    sex = meta_collapsed.loc[samples, "sex"].map(sex_map)
    if sex.isna().any():
        mode = sex.dropna().mode()
        if mode.empty:
            sys.exit(f"ERROR: no usable sex values for {ancestry}")
        sex = sex.fillna(float(mode.iloc[0]))
    blocks.append(pd.DataFrame({"sex": sex}, index=samples).T)
    priority["sex"] = (0, 100)

    co_vec = meta_collapsed.loc[samples, "cohort"].astype(str)
    co_counts = co_vec.value_counts()
    if len(co_counts) > 1:
        ref_co = co_counts.idxmax()
        for i, co in enumerate(sorted(c for c in co_counts.index if c != ref_co)):
            name = f"cohort_{co}"
            blocks.append(pd.DataFrame({name: (co_vec == co).astype(float)}, index=samples).T)
            priority[name] = (0, 200 + i)
        log(f"  {ancestry} cohort dummies: reference={ref_co}, {len(co_counts)-1} dummies")

    if "GA" in meta_collapsed.columns:
        ga = pd.to_numeric(meta_collapsed.loc[samples, "GA"], errors="coerce")
        sd = ga.std(ddof=0)
        ga = (ga - ga.mean()) / sd if np.isfinite(sd) and sd > 0 else ga * 0.0
        ga = ga.fillna(0.0)
        blocks.append(pd.DataFrame({"GA": ga}, index=samples).T)
        priority["GA"] = (0, 300)
        log(f"  {ancestry} GA covariate added (scanner drops it when GA is the exposure)")

    ct_path = os.path.join(qtl_dir, f"{ancestry}_deconvolution_harmonized.tsv")
    if os.path.exists(ct_path):
        ct = pd.read_csv(ct_path, sep="	").set_index("sample_id")
        ct = ct.drop(columns=["cohort"], errors="ignore").reindex(samples)
        ct = ct.apply(pd.to_numeric, errors="coerce")
        ct = ct.fillna(ct.mean())
        dominant = ct.mean().idxmax()
        maternal = next((c for c in ct.columns if str(c).strip().lower() == "maternal"), None)
        drop_ct = [dominant]
        if maternal is not None and maternal not in drop_ct:
            drop_ct.append(maternal)
        elif maternal is None:
            log(f"  WARNING: {ancestry} maternal cell fraction column not found")
        ct = ct.drop(columns=drop_ct)
        ct = np.arcsinh(ct)
        ct = ct - ct.mean()
        ct_t = ct.T
        ct_t.index = [f"ct_{c}" for c in ct_t.index]
        blocks.append(ct_t)
        for j, name in enumerate(ct_t.index):
            priority[name] = (2, j)
        log(f"  {ancestry} cell types: {ct_t.shape[0]} (dropped dominant: {dominant}; "
            f"excluded maternal fraction: {maternal if maternal is not None else 'not found'})")
    else:
        log(f"  WARNING: deconvolution not found: {ct_path}")

    if not blocks:
        sys.exit(f"ERROR: no base covariates assembled for {ancestry}")
    cov = pd.concat(blocks).fillna(0.0)
    log(f"  {ancestry}: assembled {cov.shape[0]} base covariates; pruning at |r|>{cor_threshold}")
    cov, dropped = prune_correlated(cov, priority, cor_threshold)
    cov.index.name = "covariate"
    out_path = os.path.join(out_dir, f"{ancestry}_covariates_base.tsv")
    cov.to_csv(out_path, sep="	", float_format="%.6g")
    log(f"  wrote {out_path}: {cov.shape[0]} covariates ({len(dropped)} pruned)")
    if dropped:
        pd.DataFrame(dropped).to_csv(
            os.path.join(out_dir, f"{ancestry}_covariate_pruning_base.tsv"),
            sep="	", index=False)
    return cov


def finalize_covariates(out_dir, ancestry, modality, hcp_file, cor_threshold=0.9):
    """Append ancestry-specific HCP rows to the ancestry base table."""
    base_path = os.path.join(out_dir, f"{ancestry}_covariates_base.tsv")
    base = pd.read_csv(base_path, sep="	", index_col=0)
    hcp = pd.read_csv(hcp_file, sep="	", index_col=0)
    missing_base = [c for c in hcp.columns if c not in base.columns]
    if missing_base:
        sys.exit(f"ERROR: HCP file {hcp_file} contains samples absent from the {ancestry} base covariates: {missing_base[:5]}")
    # Base covariates are built on the union of modality sample sets. HCPs are
    # modality-specific, so finalize only the samples present in this modality.
    base = base[hcp.columns]
    if hcp.isna().any().any():
        sys.exit(f"ERROR: NaN in HCP file {hcp_file}")
    priority = {name: (0 if name.startswith("cohort_") or name in ("sex", "GA")
                       else 1 if "_PC" in name else 2, i)
                for i, name in enumerate(base.index)}
    for j, name in enumerate(hcp.index):
        priority[name] = (3, j)
    cov = pd.concat([base, hcp])
    cov, dropped = prune_correlated(cov, priority, cor_threshold)
    cov.index.name = "covariate"
    out_path = os.path.join(out_dir, f"{ancestry}_covariates_{modality}.tsv")
    cov.to_csv(out_path, sep="	", float_format="%.6g")
    log(f"  wrote {out_path}: {cov.shape[0]} covariates ({hcp.shape[0]} HCPs in; "
        f"{sum(1 for d in dropped if str(d['covariate']).startswith('HCP_'))} HCPs pruned)")
    if dropped:
        pd.DataFrame(dropped).to_csv(
            os.path.join(out_dir, f"{ancestry}_covariate_pruning_{modality}.tsv"),
            sep="	", index=False)

def main():
    p = argparse.ArgumentParser(description="Build ancestry-specific GxE inputs")
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--collapsed-ancestry-map", default=None)
    p.add_argument("--metadata", default=None,
                   help="DEPRECATED: Module-07 reads Module-05 {ANC}_metadata.tsv")
    p.add_argument("--pcair-dir", default=None)
    p.add_argument("--ancestries", default="EAS EUR")
    p.add_argument("--modalities", default=MODALITIES)
    p.add_argument("--gxe-config", default=os.path.join(HERE, "gxe_config.tsv"))
    p.add_argument("--gxe-dir", default=None)
    p.add_argument("--cor-threshold", type=float, default=0.9)
    p.add_argument("--finalize-covariates", action="store_true")
    p.add_argument("--ancestry", default=None)
    p.add_argument("--modality", default=None)
    p.add_argument("--hcp-file", default=None)
    args = p.parse_args()

    out_dir = args.gxe_dir or os.path.join(args.results_dir, "gxe")
    inputs_dir = os.path.join(out_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    if args.finalize_covariates:
        if not (args.ancestry and args.modality and args.hcp_file):
            sys.exit("ERROR: --finalize-covariates needs --ancestry, --modality and --hcp-file")
        finalize_covariates(inputs_dir, args.ancestry, args.modality,
                            args.hcp_file, args.cor_threshold)
        return

    ancestries = args.ancestries.split()
    modalities = args.modalities.split()
    pcair_dir = args.pcair_dir or os.path.join(
        os.path.dirname(args.results_dir.rstrip("/")), "genotype_pcs")
    if args.metadata:
        log("WARNING: --metadata is deprecated and ignored")
    output_base = os.path.dirname(args.qtl_dir.rstrip("/"))
    collapsed_map = args.collapsed_ancestry_map or os.path.join(
        output_base, "replicate_collapsed", "reports", "ancestry_map_collapsed.tsv")
    meta_all, meta_one = load_module05_metadata(
        args.qtl_dir, ancestries, collapsed_map, args.gxe_config)

    all_manifest_samples = []
    for anc in ancestries:
        log(f"ancestry: {anc}")
        anc_meta = meta_one[meta_one["ancestry"].astype(str) == str(anc)].copy()
        if anc_meta.empty:
            sys.exit(f"ERROR: no retained Module-05 metadata for ancestry {anc}")
        meta_ids = set(anc_meta.index.astype(str))

        first_order = None
        modality_samples = {}
        for mod in modalities:
            log(f"  modality: {mod}")
            ms = build_ancestry_bed(mod, anc, args.qtl_dir, inputs_dir, meta_ids)
            modality_samples[mod] = ms
            if first_order is None:
                first_order = list(ms)
        all_samples = list(dict.fromkeys(
            (first_order or []) + [x for mod in modalities for x in modality_samples[mod]]))
        if not all_samples:
            sys.exit(f"ERROR: no GxE samples for {anc}")

        man = anc_meta.reindex(all_samples).dropna(subset=["ancestry"])
        manifest = pd.DataFrame({
            "array_id": man.index,
            "ancestry": man["ancestry"],
            "cohort": man.get("cohort", pd.Series(index=man.index, dtype=str)),
            "sex": man.get("sex", pd.Series(index=man.index, dtype=str)),
        })
        for col in ["GA", "ppBMI", "gdm", "ogtt"]:
            if col in anc_meta.columns:
                manifest[col] = pd.to_numeric(anc_meta.loc[man.index, col], errors="coerce").values
        manifest_path = os.path.join(inputs_dir, f"{anc}_sample_manifest.tsv")
        manifest.to_csv(manifest_path, sep="	", index=False)
        log(f"  wrote {manifest_path}: {len(manifest)} samples")
        all_manifest_samples.extend(manifest["array_id"].astype(str).tolist())

        pm = meta_all[
            (meta_all["ancestry"].astype(str) == str(anc))
            & meta_all["array_id"].astype(str).isin(set(manifest["array_id"].astype(str)))
        ].copy()
        if "cohort" in pm.columns:
            retained_cohort = anc_meta["cohort"].to_dict()
            expected = pm["array_id"].map(retained_cohort)
            cross_protocol = (pm["cohort"].notna() & expected.notna()
                              & (pm["cohort"].astype(str) != expected.astype(str)))
            if cross_protocol.any():
                log(f"  {anc} HCP metadata: dropping {int(cross_protocol.sum())} "
                    "non-primary cross-protocol rnaseq row(s)")
                pm = pm.loc[~cross_protocol].copy()
        pm_cols = [c for c in ["rnaseq_id", "array_id", "ancestry", "cohort", "sex"] if c in pm.columns]
        meta_path = os.path.join(inputs_dir, f"{anc}_metadata.tsv")
        pm[pm_cols].to_csv(meta_path, sep="	", index=False)
        log(f"  wrote {meta_path}: {len(pm)} rnaseq rows for {pm['array_id'].nunique()} array_ids")

        build_base_covariates(list(manifest["array_id"].astype(str)), meta_one,
                              args.qtl_dir, pcair_dir, anc, inputs_dir,
                              args.cor_threshold)

    # One common exposure scale, then each ancestry scan simply subsets columns.
    all_manifest_samples = list(dict.fromkeys(all_manifest_samples))
    build_exposures(args.gxe_config, meta_one, all_manifest_samples, inputs_dir)
    log("done. Next: 51_pooled_hcp.sh, now run per ancestry x modality.")


if __name__ == "__main__":
    main()
