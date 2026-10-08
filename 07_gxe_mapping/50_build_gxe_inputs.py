#!/usr/bin/env python3
"""
50_build_gxe_inputs.py — Build pooled multi-ancestry GxE inputs (Objective 2.1)

Tier-1 GxE discovery runs on pooled multi-ancestry samples. This script
builds the pooled inputs from the module-05 per-ancestry products:

  1. Pooled phenotype BEDs (per modality): stable phenotype intersection
     across ancestries (splicing/IR ancestry-local IDs are canonicalized to
     genomic event IDs first), z-scored WITHIN each ancestry (equalizes scale without
     removing between-ancestry mean structure into the phenotypes — the
     ancestry dummies in the covariates absorb that), then concatenated.
     Written bgzipped + tabix-indexed as inputs/pooled_{MOD}.bed.gz.
  2. inputs/pooled_sample_manifest.tsv — array_id, ancestry, cohort, sex and
     the raw exposure columns for every pooled sample.
  3. inputs/pooled_metadata.tsv — rnaseq_id/array_id/ancestry/cohort/sex for
     the pooled samples (input to 05_qtl_mapping/hcp_from_matrix.R, which
     averages technical replicates at the array_id level).
  4. inputs/exposures.tsv — rows = enabled exposures from gxe_config.tsv,
     columns = samples. Continuous exposures are z-scored pooled (beta_int
     is then per exposure SD); binary {0,1} exposures are left as-is
     (centering happens inside the scanner, as in tensorQTL). Values outside
     the configured [min, max] range are set to NA and those samples drop
     out of that exposure's scan. inputs/exposure_manifest.tsv summarizes
     per-exposure n/mean/sd/missingness.
  5. inputs/pooled_covariates_base.tsv — ancestry dummies (reference =
     largest ancestry), ancestry-specific genotype PCs (zero-filled outside
     each ancestry; the "ancestry PCs" of the aims' Zcov), sex, cohort
     dummies (reference = largest pooled cohort), GA (z-scored pooled,
     NA -> 0 = mean imputation; the scanner drops this row automatically
     when GA IS the exposure), and harmonized cell-type proportions
     (arcsinh, pooled mean-centered, dominant type dropped as the
     compositional reference) — the module-05 covariate schema, pooled.
     Iterative |r| > 0.9 pruning with the module-05 priority tiers
     (demographic > PC > cell type > HCP).

  6. --finalize-covariates --modality M --hcp-file F: append the pooled HCP
     factors estimated by 51_pooled_hcp.sh to the base table (HCPs lowest
     pruning priority) and write inputs/pooled_covariates_{M}.tsv — the
     per-modality covariate file the scanner reads.

Genotypes are NOT touched: the scanner reads per-ancestry pgens directly
(per-chromosome, variant-ID-intersected, ancestry-block stacked).

Metadata provenance:
  Individual-level metadata are taken from the Module-05 final
  {ANC}_metadata.tsv files and restricted to the retained RNA run in
  replicate_collapsed/reports/ancestry_map_collapsed.tsv. This keeps
  exposures/demographics aligned to the same representative run used by the
  replicate-collapse workflow. pooled_metadata.tsv retains same-cohort
  technical-replicate rnaseq_id rows for the final individuals so
  hcp_from_matrix.R can average their Picard metrics, but excludes non-primary
  cross-protocol rows using the retained run's cohort.

Usage:
  python3 50_build_gxe_inputs.py --qtl-dir <qtl_inputs> \
      --results-dir <qtl_results> \
      --collapsed-ancestry-map <replicate_collapsed/reports/ancestry_map_collapsed.tsv> \
      --pcair-dir <genotype_pcs> [--ancestries "EAS EUR"] [--modalities ...]

  python3 50_build_gxe_inputs.py --qtl-dir ... --results-dir ... \
      --finalize-covariates --modality expression \
      --hcp-file inputs/pooled_hcp_expression.tsv
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
            "pool potentially non-equivalent phenotypes"
        )

    dup = canonical.duplicated(keep=False)
    if dup.any():
        examples = canonical[dup].head(5).tolist()
        sys.exit(
            f"ERROR: {modality} canonicalization creates duplicate phenotype "
            f"IDs for {ancestry}: {examples}. Investigate upstream "
            "harmonization before pooling."
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


def build_pooled_bed(modality, ancestries, qtl_dir, out_dir, meta_ids):
    """Intersect phenotypes, z-score within ancestry, concatenate samples."""
    blocks = []
    pos_ref = None
    for anc in ancestries:
        path = os.path.join(qtl_dir, f"{anc}_{modality}.bed.gz")
        if not os.path.exists(path):
            sys.exit(f"ERROR: BED not found: {path}")
        df = read_bed(path)
        df["#chr"] = df["#chr"].astype(str)

        # Module 07 is autosomal only.  Upstream phenotype BEDs can contain
        # chrX/chrY/other contigs, so remove them before phenotype
        # intersection, positional consistency checks, and BED sorting.
        chrom_num = pd.to_numeric(
            df["#chr"].str.replace(r"^chr", "", case=False, regex=True),
            errors="coerce",
        )
        autosomal = chrom_num.between(1, 22) & (chrom_num % 1 == 0)
        n_nonauto = int((~autosomal).sum())
        if n_nonauto:
            log(f"    {anc}: dropping {n_nonauto} non-autosomal phenotypes")
        df = df.loc[autosomal].copy()

        # Splicing and IR are harmonized within ancestry upstream, so their
        # literal phenotype IDs contain ancestry-local numbering.  Convert
        # those IDs to stable genomic event IDs before ancestry intersection.
        df = canonicalize_phenotype_ids(df, modality, anc)

        # restrict to samples with metadata (exposures/covariates need it)
        keep_cols = [c for c in df.columns[4:] if c in meta_ids]
        n_drop = len(df.columns) - 4 - len(keep_cols)
        if n_drop:
            log(f"    {anc}: dropping {n_drop} BED samples without metadata")
        df = df[["#chr", "start", "end", "phenotype_id"] + keep_cols]
        blocks.append((anc, df))
        log(f"    {anc}: {len(df)} phenotypes x {len(keep_cols)} samples")

    # phenotype intersection + position consistency
    ids = set(blocks[0][1]["phenotype_id"])
    for _, df in blocks[1:]:
        ids &= set(df["phenotype_id"])
    if not ids:
        sys.exit(
            f"ERROR: no cross-ancestry phenotype overlap for {modality} "
            "after canonicalization; inspect upstream phenotype definitions"
        )
    log(f"    cross-ancestry phenotype intersection: {len(ids)}")
    base = blocks[0][1]
    base = base[base["phenotype_id"].isin(ids)]
    pos_ref = base[["phenotype_id", "#chr", "start", "end"]]
    value_blocks = []
    for anc, df in blocks:
        df = df.set_index("phenotype_id")
        sub = df.loc[pos_ref["phenotype_id"]]
        if not ((sub["#chr"].values == pos_ref["#chr"].values).all()
                and (sub["start"].values == pos_ref["start"].values).all()
                and (sub["end"].values == pos_ref["end"].values).all()):
            sys.exit(f"ERROR: phenotype positions differ across ancestries "
                     f"for {modality} ({anc}) — investigate before pooling")
        vals = sub.iloc[:, 3:].values.astype(np.float64)
        z = zscore_rows(vals)
        value_blocks.append(pd.DataFrame(
            z, index=sub.index, columns=sub.columns[3:]))
    pooled = pd.concat(value_blocks, axis=1)
    n_bad = int(pooled.isna().any(axis=1).sum())
    if n_bad:
        log(f"    dropping {n_bad} phenotypes with ~zero within-ancestry "
            f"variance")
        keep = ~pooled.isna().any(axis=1)
        pooled = pooled[keep]
        pos_ref = pos_ref[pos_ref["phenotype_id"].isin(pooled.index)]
    out_df = pd.concat([pos_ref.set_index("phenotype_id").loc[pooled.index]
                        .reset_index(), pooled.reset_index(drop=True)], axis=1)
    # Preserve the standard tensorQTL BED column order so bgzip/tabix see
    # chromosome/start/end in columns 1-3 and the header begins with #chr.
    sample_cols = list(pooled.columns)
    out_df = out_df[["#chr", "start", "end", "phenotype_id"] + sample_cols]
    out_df = out_df.sort_values(
        ["#chr", "start"],
        key=lambda s: pd.to_numeric(
            s.astype(str).str.replace(r"^chr", "", case=False, regex=True),
            errors="raise",
        ) if s.name == "#chr" else s,
    )
    plain = os.path.join(out_dir, f"pooled_{modality}.bed")
    out_df.to_csv(plain, sep="\t", index=False, float_format="%.6g")
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-p", "bed", "-f", plain + ".gz"], check=True)
    log(f"    wrote {plain}.gz (+ .tbi): {len(out_df)} phenotypes x "
        f"{len(out_df.columns) - 4} pooled samples")
    return list(out_df.columns[4:])


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
                          ancestries, out_dir, cor_threshold=0.9):
    blocks, priority = [], {}

    # ancestry dummies (reference = largest)
    anc_vec = meta_collapsed.loc[samples, "ancestry"].astype(str)
    counts = anc_vec.value_counts()
    ref = counts.idxmax()
    for i, a in enumerate(sorted(x for x in counts.index if x != ref)):
        blocks.append(pd.DataFrame({f"anc_{a}": (anc_vec == a).astype(float)},
                                   index=samples).T)
        priority[f"anc_{a}"] = (0, i)
    log(f"  ancestry dummies: reference={ref} (n={counts[ref]}), "
        f"{len(counts) - 1} dummies")

    # ancestry-specific genotype PCs (zero-filled outside the ancestry)
    for anc in ancestries:
        sel = os.path.join(qtl_dir, f"{anc}_selected_pcs.txt")
        pcs_path = os.path.join(pcair_dir, f"{anc}_genotype_pcs.tsv")
        if not (os.path.exists(sel) and os.path.exists(pcs_path)):
            log(f"  WARNING: PCs missing for {anc} ({sel} / {pcs_path})")
            continue
        pc_names = [l.strip() for l in open(sel) if l.strip()]
        pcs = pd.read_csv(pcs_path, sep="\t").set_index("sample_id")
        block = pd.DataFrame(0.0, index=samples,
                             columns=[f"{anc}_{p}" for p in pc_names])
        anc_samples = [s for s in samples if anc_vec[s] == anc]
        common = [s for s in anc_samples if s in pcs.index]
        block.loc[common, :] = pcs.loc[common, pc_names].values
        blocks.append(block.T)
        for j, name in enumerate(block.columns):
            priority[name] = (1, j)
        log(f"  {anc} PCs: {len(pc_names)} (zero-filled outside {anc})")

    # sex
    sex_map = {"M": 0, "F": 1, "m": 0, "f": 1, "Male": 0, "Female": 1}
    sex = meta_collapsed.loc[samples, "sex"].map(sex_map)
    if sex.isna().any():
        sex = sex.fillna(float(sex.mode()[0]))
    blocks.append(pd.DataFrame({"sex": sex}, index=samples).T)
    priority["sex"] = (0, 100)

    # cohort dummies (reference = largest pooled cohort)
    co_vec = meta_collapsed.loc[samples, "cohort"].astype(str)
    co_counts = co_vec.value_counts()
    if len(co_counts) > 1:
        ref_co = co_counts.idxmax()
        for i, co in enumerate(sorted(c for c in co_counts.index if c != ref_co)):
            blocks.append(pd.DataFrame({f"cohort_{co}": (co_vec == co).astype(float)},
                                       index=samples).T)
            priority[f"cohort_{co}"] = (0, 200 + i)
        log(f"  cohort dummies: reference={ref_co}, {len(co_counts) - 1} dummies")

    # GA (z-scored pooled, NA -> 0 = pooled-mean imputation)
    if "GA" in meta_collapsed.columns:
        ga = pd.to_numeric(meta_collapsed.loc[samples, "GA"], errors="coerce")
        ga = (ga - ga.mean()) / ga.std(ddof=0)
        ga = ga.fillna(0.0)
        blocks.append(pd.DataFrame({"GA": ga}, index=samples).T)
        priority["GA"] = (0, 300)
        log("  GA covariate added (scanner drops it when GA is the exposure)")

    # cell-type proportions (harmonized deconvolution; arcsinh, pooled
    # mean-centered, dominant type dropped as compositional reference)
    ct_blocks = []
    for anc in ancestries:
        path = os.path.join(qtl_dir, f"{anc}_deconvolution_harmonized.tsv")
        if not os.path.exists(path):
            log(f"  WARNING: deconvolution not found: {path}")
            continue
        d = pd.read_csv(path, sep="\t").set_index("sample_id")
        ct_blocks.append(d)
    if ct_blocks:
        common_cols = set(ct_blocks[0].columns)
        for d in ct_blocks[1:]:
            common_cols &= set(d.columns)
        common_cols -= {"cohort"}
        ct = pd.concat([d[sorted(common_cols)] for d in ct_blocks])
        ct = ct.reindex(samples)
        ct = ct.fillna(ct.mean())  # rare missing proportions -> pooled mean

        # Match Module 05's production covariate policy: the maternal cell
        # fraction is excluded, and the dominant placental cell type is also
        # dropped as the compositional reference.  Determine the dominant
        # reference from the full composition, then remove both terms before
        # transformation/pruning.
        dominant = ct.mean().idxmax()
        maternal = next(
            (c for c in ct.columns if str(c).strip().lower() == "maternal"),
            None,
        )
        drop_ct = [dominant]
        if maternal is not None and maternal not in drop_ct:
            drop_ct.append(maternal)
        elif maternal is None:
            log("  WARNING: maternal cell fraction column not found in "
                "harmonized deconvolution inputs")
        ct = ct.drop(columns=drop_ct)
        ct = np.arcsinh(ct)
        ct = ct - ct.mean()
        ct.index.name = "sample_id"
        ct_t = ct.T
        ct_t.index = [f"ct_{c}" for c in ct_t.index]
        blocks.append(ct_t)
        for j, name in enumerate(ct_t.index):
            priority[name] = (2, j)
        maternal_msg = maternal if maternal is not None else "not found"
        log(f"  cell types: {ct_t.shape[0]} (dropped dominant: {dominant}; "
            f"excluded maternal fraction: {maternal_msg})")

    cov = pd.concat(blocks)
    cov = cov.fillna(0.0)
    log(f"  assembled {cov.shape[0]} base covariates; pruning at |r|>{cor_threshold}")
    cov, dropped = prune_correlated(cov, priority, cor_threshold)
    log(f"  pruned {len(dropped)}; final base covariates: {cov.shape[0]}")
    cov.index.name = "covariate"
    cov.to_csv(os.path.join(out_dir, "pooled_covariates_base.tsv"), sep="\t",
               float_format="%.6g")
    if dropped:
        pd.DataFrame(dropped).to_csv(
            os.path.join(out_dir, "pooled_covariate_pruning_base.tsv"),
            sep="\t", index=False)
    return cov


def finalize_covariates(out_dir, modality, hcp_file, cor_threshold=0.9):
    """Append pooled HCP rows to the base table (lowest priority) -> per-modality file."""
    base_path = os.path.join(out_dir, "pooled_covariates_base.tsv")
    base = pd.read_csv(base_path, sep="\t", index_col=0)
    hcp = pd.read_csv(hcp_file, sep="\t", index_col=0)
    hcp = hcp[[c for c in base.columns if c in hcp.columns]]
    hcp = hcp.reindex(columns=base.columns)
    if hcp.isna().any().any():
        sys.exit(f"ERROR: HCP file {hcp_file} does not cover the base samples")
    priority = {name: (0 if n.startswith(("anc_", "cohort_")) or n in ("sex", "GA")
                       else 1 if "_PC" in n else 2, i)
                for i, n in enumerate(base.index)}
    for j, name in enumerate(hcp.index):
        priority[name] = (3, j)
    cov = pd.concat([base, hcp])
    cov, dropped = prune_correlated(cov, priority, cor_threshold)
    cov.index.name = "covariate"
    out_path = os.path.join(out_dir, f"pooled_covariates_{modality}.tsv")
    cov.to_csv(out_path, sep="\t", float_format="%.6g")
    log(f"  wrote {out_path}: {cov.shape[0]} covariates "
        f"({hcp.shape[0]} HCPs in, {sum(1 for d in dropped if d['covariate'].startswith('HCP_'))} pruned)")
    if dropped:
        pd.DataFrame(dropped).to_csv(
            os.path.join(out_dir, f"pooled_covariate_pruning_{modality}.tsv"),
            sep="\t", index=False)


def main():
    p = argparse.ArgumentParser(description="Build pooled GxE inputs")
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument(
        "--collapsed-ancestry-map", default=None,
        help=("Module-05 replicate-collapsed ancestry map. Default: "
              "{OUTPUT_BASE}/replicate_collapsed/reports/"
              "ancestry_map_collapsed.tsv")
    )
    p.add_argument(
        "--metadata", default=None,
        help=("DEPRECATED: Stage 0 now reads {ANC}_metadata.tsv from --qtl-dir "
              "and the collapsed ancestry map to match Module-05 provenance")
    )
    p.add_argument("--pcair-dir", default=None,
                   help="genotype PCs dir (default: {OUTPUT_BASE}/genotype_pcs)")
    p.add_argument("--ancestries", default="EAS EUR")
    p.add_argument("--modalities", default=MODALITIES)
    p.add_argument("--gxe-config", default=os.path.join(HERE, "gxe_config.tsv"))
    p.add_argument("--gxe-dir", default=None,
                   help="output dir (default: {results-dir}/gxe)")
    p.add_argument("--cor-threshold", type=float, default=0.9)
    p.add_argument("--finalize-covariates", action="store_true")
    p.add_argument("--modality", default=None)
    p.add_argument("--hcp-file", default=None)
    args = p.parse_args()

    out_dir = args.gxe_dir or os.path.join(args.results_dir, "gxe")
    inputs_dir = os.path.join(out_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)

    if args.finalize_covariates:
        if not (args.modality and args.hcp_file):
            sys.exit("ERROR: --finalize-covariates needs --modality and --hcp-file")
        finalize_covariates(inputs_dir, args.modality, args.hcp_file,
                            args.cor_threshold)
        return

    ancestries = args.ancestries.split()
    modalities = args.modalities.split()
    pcair_dir = args.pcair_dir or os.path.join(
        os.path.dirname(args.results_dir.rstrip("/")), "genotype_pcs")

    if args.metadata:
        log("WARNING: --metadata is deprecated and ignored; using Module-05 "
            "{ANC}_metadata.tsv plus the collapsed ancestry map")
    output_base = os.path.dirname(args.qtl_dir.rstrip("/"))
    collapsed_ancestry_map = args.collapsed_ancestry_map or os.path.join(
        output_base, "replicate_collapsed", "reports",
        "ancestry_map_collapsed.tsv"
    )
    meta, meta_c = load_module05_metadata(
        args.qtl_dir, ancestries, collapsed_ancestry_map, args.gxe_config
    )
    meta_ids = set(meta_c.index)

    # 1. pooled BEDs (sample order = ancestry blocks, BED order within)
    samples = None
    for mod in modalities:
        log(f"modality: {mod}")
        mod_samples = build_pooled_bed(mod, ancestries, args.qtl_dir,
                                       inputs_dir, meta_ids)
        if samples is None:
            samples = mod_samples
        elif set(mod_samples) != set(samples):
            log(f"  WARNING: {mod} sample set differs from the first "
                f"modality; the manifest keeps the union and scans intersect")
    # union of modality sample sets, first-modality order then extras
    all_samples = list(dict.fromkeys(
        (samples or []) + [s for mod in modalities for s in
                           pd.read_csv(os.path.join(inputs_dir, f"pooled_{mod}.bed.gz"),
                                       sep="\t", nrows=0).columns[4:]]))

    # 2. manifest
    man = meta_c.reindex(all_samples)
    man = man.dropna(subset=["ancestry"])
    manifest = pd.DataFrame({
        "array_id": man.index,
        "ancestry": man["ancestry"],
        "cohort": man.get("cohort", pd.Series(index=man.index, dtype=str)),
        "sex": man.get("sex", pd.Series(index=man.index, dtype=str)),
    })
    for col in ["GA", "ppBMI", "gdm", "ogtt"]:
        if col in meta_c.columns:
            manifest[col] = pd.to_numeric(meta_c.loc[man.index, col],
                                          errors="coerce").values
    manifest.to_csv(os.path.join(inputs_dir, "pooled_sample_manifest.tsv"),
                    sep="\t", index=False)
    log(f"wrote pooled_sample_manifest.tsv: {len(manifest)} samples")

    # 3. pooled metadata for hcp_from_matrix.R (rnaseq_id level)
    pm = meta[meta["array_id"].isin(set(manifest["array_id"]))].copy()
    # Same-cohort technical replicates contribute QC metrics to the collapsed
    # individual, matching hcp_from_matrix.R's averaging. Cross-protocol
    # non-primary runs are excluded: collapse_replicates.py's default policy
    # keeps the retained run/cohort as the individual's representation.
    if "cohort" in pm.columns:
        retained_cohort = meta_c["cohort"].to_dict()
        expected = pm["array_id"].map(retained_cohort)
        cross_protocol = (
            pm["cohort"].notna() & expected.notna()
            & (pm["cohort"].astype(str) != expected.astype(str))
        )
        n_cross_drop = int(cross_protocol.sum())
        if n_cross_drop:
            log(f"  pooled HCP metadata: dropping {n_cross_drop} non-primary "
                "cross-protocol rnaseq row(s)")
            pm = pm.loc[~cross_protocol].copy()
    pm_cols = [c for c in ["rnaseq_id", "array_id", "ancestry", "cohort", "sex"]
               if c in pm.columns]
    pm[pm_cols].to_csv(os.path.join(inputs_dir, "pooled_metadata.tsv"),
                       sep="\t", index=False)
    n_rep_extra = len(pm) - pm["array_id"].nunique()
    log(f"wrote pooled_metadata.tsv: {len(pm)} rnaseq rows for "
        f"{pm['array_id'].nunique()} array_ids ({n_rep_extra} same-cohort "
        "technical-replicate row(s) retained for HCP QC averaging)")

    # 4. exposures
    build_exposures(args.gxe_config, meta_c, list(manifest["array_id"]),
                    inputs_dir)

    # 5. base covariates
    build_base_covariates(list(manifest["array_id"]), meta_c, args.qtl_dir,
                          pcair_dir, ancestries, inputs_dir,
                          args.cor_threshold)
    log("done. Next: 51_pooled_hcp.sh (per-modality pooled HCPs + finalize).")


if __name__ == "__main__":
    main()
