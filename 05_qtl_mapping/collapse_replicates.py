#!/usr/bin/env python3
"""
collapse_replicates.py — Collapse technical-replicate RNA-seq runs to one
column per individual, upstream of within-ancestry pooling, normalization, and QTL
mapping.

Some individuals were sequenced as two runs (e.g. two placental quadrant
samples). Downstream QTL mapping requires one column per individual; this
module performs that collapse on the per-cohort pre-normalization products
and writes a parallel staging tree that the pooling / harmonization steps
consume in place of the original cohort directories.

Collapse policy (per pair of runs, per modality):

  * Count modalities (expression, isoform_expression):
    sum the per-run counts.
  * Ratio modalities (isoforms, alt_TSS, alt_polyA, splicing, stability,
    RNA_editing):
    sum the underlying numerator/denominator counts and recompute the
    ratio. Ratios themselves are never averaged.
  * intron_retention:
    MAJIQ PSI exposes no per-run counts at this stage, so the per-run PSI
    columns are averaged directly (same cohort -> identical feature set).
  * Cross-protocol pairs (runs from more than one cohort):
    default drop (--cross-protocol drop): a sample whose runs span
    cohorts indicates a labeling error, so the individual is removed
    entirely -- all runs dropped from every cohort file and from the
    collapsed ancestry map (action dropped_cross_protocol). With
    '--cross-protocol average' the count-exact collapse is applied across
    cohorts for the ratio modalities whose values are recomputed from
    per-run source files (isoforms, alt_TSS, alt_polyA, stability,
    RNA_editing); other modalities still drop the individual.
  * Concordance report (flag-only, never a veto):
    per pair x modality, Pearson and Spearman correlation between the two
    runs over the cohort's unnorm BED. Pairs with Spearman <
    --concordance-threshold (default 0.9) or too few features are flagged
    in the concordance report for manual review, but ARE still collapsed:
    technical replicates are always averaged. The only keep-primary case
    is single_run_present -- one run genuinely absent from the modality.

Primary-run rule: prefer the run present in the ancestry map, then
lexicographic order.

Staging tree contract (under --out-dir):
  * Touched cohorts (contain >=1 replicate pair) get a real directory with:
      output/unnorm/<modality>.bed
          for the seven intersection modalities (expression, isoforms,
          isoform_expression, alt_TSS, alt_polyA, RNA_editing, stability).
          Feature sets and non-replicated sample columns are unchanged.
      intermediate/splicing/leafcutter_perind_numers.counts.gz
          collapsed LeafCutter numerator counts (consumed by
          17_harmonize_within_ancestry.sh via --staging-dir).
      intermediate/intron_retention/retained_intron_psi.tsv.gz
          run-averaged MAJIQ PSI (consumed by the same step).
    Splicing and intron_retention per-cohort BEDs are NOT staged: those
    modalities flow through harmonize_within_ancestry.py from the staged
    intermediates, which re-derives the feature set.
  * Untouched cohorts are symlinked: <out-dir>/<cohort> -> <output_base>/<cohort>
    so downstream --cohort-dirs can point uniformly at the staging root.

Reports (under --out-dir/reports):
  replicate_collapses.tsv and {ANC}_replicate_collapses.tsv
      One row per individual x modality (x cohort for cross-protocol
      pairs): runs, primary run, action, kept column, stats.
  replicate_concordance.tsv and {ANC}_replicate_concordance.tsv
      One row per individual x modality: Pearson/Spearman + review flag.
  ancestry_map_collapsed.tsv
      Ancestry map with non-primary runs of replicate pairs removed and
      cross-protocol individuals removed entirely (drop-in replacement
      for downstream steps).

Usage:
  python3 collapse_replicates.py \
      --metadata placenta_QTL_cohort_metadata.tsv \
      --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
      --config config.yml \
      --ref-anno /path/to/normalized.gtf \
      --out-dir /path/to/collapse_staging

  # Inspect pairs and concordance without writing the staging tree:
  python3 collapse_replicates.py ... --dry-run
"""

import argparse
import gzip
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


BED_META_COLS = ["#chr", "start", "end", "phenotype_id"]

# Modalities with staged per-cohort unnorm BEDs (consumed by the
# intersection poolers via --cohort-dirs).
# The *_counts variants (tximport_counts.R outputs) are count-scale, so
# technical replicates are summed exactly like expression/isoform_expression.
COUNT_MODALITIES = {"expression", "isoform_expression",
                    "expression_counts", "isoform_expression_counts"}
RATIO_BED_MODALITIES = {"isoforms", "alt_TSS", "alt_polyA", "stability", "RNA_editing"}
STAGED_BED_MODALITIES = sorted(COUNT_MODALITIES | RATIO_BED_MODALITIES)

# Modalities staged at the intermediate level for 17_harmonize_within_ancestry.sh.
INTERMEDIATE_MODALITIES = {"splicing", "intron_retention"}

ALL_MODALITIES = sorted(set(STAGED_BED_MODALITIES) | INTERMEDIATE_MODALITIES)

# Per-cohort intermediate locations (relative to cohort dir), matching the
# aggregation scripts 10-15.
EXPR_QU_DIR = os.path.join("intermediate", "expression_qu")
ALT_DIR = os.path.join("intermediate", "alt_TSS_polyA")
SPLICING_NUMERS = os.path.join(
    "intermediate", "splicing", "leafcutter_perind_numers.counts.gz")
IR_PSI = os.path.join(
    "intermediate", "intron_retention", "retained_intron_psi.tsv.gz")
STAB_DIR = os.path.join("intermediate", "stability")
EDIT_MATRIX = os.path.join("intermediate", "RNA_editing", "edit_site_matrix.tsv")
EDIT_SITE_MAP = os.path.join("output", "RNA_editing.site_to_phenotype.tsv")

STAB_MIN_COUNT = 10  # per-VALUE floor, matching assemble_bed.load_featureCounts
# (values below the floor become NaN; this is per-observation, not a
# per-feature filter, so it is kept under union pooling — feature-level
# filtering happens once, on the pooled matrix, at the normalization stage)


# ---------------------------------------------------------------------------
# Config / metadata / pair discovery
# ---------------------------------------------------------------------------

def load_cohorts(config_path):
    """Extract cohort names from config.yml (the 'cohorts:' block)."""
    with open(config_path) as f:
        text = f.read()
    m = re.search(r"^cohorts:\s*$", text, re.MULTILINE)
    if not m:
        return []
    cohorts = []
    for line in text[m.end():].splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not line.startswith(" ") and not line.startswith("\t"):
            break
        if stripped.endswith(":") and ":" not in stripped[:-1]:
            cohorts.append(stripped[:-1])
    return cohorts


def load_output_base(config_path):
    with open(config_path) as f:
        text = f.read()
    m = re.search(r"^output_base:\s*(\S+)\s*$", text, re.MULTILINE)
    if not m:
        sys.stderr.write("ERROR: output_base not found in config.yml\n")
        sys.exit(1)
    return m.group(1)


def load_metadata(path):
    """Load pooled metadata; return DataFrame with rnaseq_id, array_id."""
    df = pd.read_csv(path, sep="\t")
    missing = {"rnaseq_id", "array_id"} - set(df.columns)
    if missing:
        sys.stderr.write(
            f"ERROR: metadata missing columns: {missing}\n"
            f"  Got: {list(df.columns)}\n")
        sys.exit(1)
    return df


def load_ancestry_map(path):
    df = pd.read_csv(path, sep="\t")
    missing = {"sample_id", "assigned_ancestry", "cohort"} - set(df.columns)
    if missing:
        sys.stderr.write(
            f"ERROR: ancestry map missing columns: {missing}\n"
            f"  Got: {list(df.columns)}\n")
        sys.exit(1)
    return df


def discover_pairs(meta_df):
    """Group runs by individual. Returns dict array_id -> sorted run list
    (only individuals with >1 distinct run)."""
    sub = meta_df[["rnaseq_id", "array_id"]].dropna().drop_duplicates()
    pairs = {}
    for array_id, grp in sub.groupby("array_id"):
        runs = sorted(grp["rnaseq_id"].unique())
        if len(runs) > 1:
            pairs[array_id] = runs
    return pairs


def map_runs_to_cohorts(output_base, cohorts):
    """Build run -> cohort from per-cohort samples.txt files."""
    run2cohort = {}
    for cohort in cohorts:
        samples_file = os.path.join(output_base, cohort, "samples.txt")
        if not os.path.exists(samples_file):
            sys.stderr.write(f"  WARN: {samples_file} not found, skipping cohort {cohort}\n")
            continue
        with open(samples_file) as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s in run2cohort and run2cohort[s] != cohort:
                    sys.stderr.write(
                        f"  WARN: run {s} listed in both {run2cohort[s]} and "
                        f"{cohort} samples.txt; keeping {run2cohort[s]}\n")
                    continue
                run2cohort[s] = cohort
    return run2cohort


def choose_primary(runs, ancestry_sample_ids):
    """Primary run: prefer a run present in the ancestry map, then
    lexicographic order."""
    in_map = [r for r in sorted(runs) if r in ancestry_sample_ids]
    if in_map:
        return in_map[0]
    return sorted(runs)[0]


def representative_run(info, runs_present):
    """The column that represents the individual in a given file.

    Same-cohort pair: the primary run when present, else the first present
    run (no duplicate risk — the file is per-cohort). Cross-protocol pair:
    always the primary run, so the individual stays single-column across
    cohorts after pooling.
    """
    if info["cross_protocol"]:
        return info["primary"]
    if info["primary"] in runs_present:
        return info["primary"]
    return runs_present[0]


# ---------------------------------------------------------------------------
# BED helpers
# ---------------------------------------------------------------------------

def load_bed(path):
    return pd.read_csv(
        path, sep="\t",
        dtype={"#chr": str, "start": int, "end": int, "phenotype_id": str})


def write_bed(df, path):
    df.to_csv(path, sep="\t", index=False, float_format="%g")


def concordance(bed, run1, run2):
    """Pearson/Spearman between two sample columns over pairwise-complete
    features. Returns (pearson, spearman, n_features)."""
    x = bed[run1]
    y = bed[run2]
    both = x.notna() & y.notna()
    n = int(both.sum())
    if n < 3:
        return np.nan, np.nan, n
    xr = x[both]
    yr = y[both]
    pearson = float(xr.corr(yr, method="pearson"))
    spearman = float(xr.corr(yr, method="spearman"))
    return pearson, spearman, n


# ---------------------------------------------------------------------------
# Per-modality collapsed-value handlers
#
# Each handler returns a numpy array aligned to the rows of the BED being
# rewritten, or None when the per-run count sources are missing (caller
# falls back to keep-primary). Ratio modalities recompute ratios from
# summed numerator/denominator counts — ratios are never averaged.
# ---------------------------------------------------------------------------

def _load_salmon_counts(quant_sf):
    """Load NumReads from one Salmon quant.sf, indexed by transcript Name."""
    df = pd.read_csv(quant_sf, sep="\t")
    return df.set_index("Name")["NumReads"].astype(float)


def _sum_salmon_runs(run_dirs, runs, subdir):
    """Sum NumReads across runs. subdir is the cohort-relative quant dir
    (e.g. expression_qu or alt_TSS_polyA/grp_1.upstream). Returns a Series
    indexed by target_id, or None if any run's quant.sf is missing."""
    summed = None
    for run in runs:
        q = os.path.join(run_dirs[run], subdir, run, "quant.sf")
        if not os.path.exists(q):
            return None
        s = _load_salmon_counts(q)
        summed = s if summed is None else summed.add(s, fill_value=0.0)
    return summed


def _load_tx2gene(ref_anno):
    """Transcript -> gene map from the normalized GTF, reusing
    assemble_bed.py (single source of truth for the mapping)."""
    scripts = Path(__file__).resolve().parent.parent / "03_phenotyping"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    from assemble_bed import transcript_to_gene_map  # lazy: needs gtfparse
    m = transcript_to_gene_map(Path(ref_anno))
    return dict(zip(m["transcript_id"], m["gene_id"]))


def collapse_counts(bed, runs_present):
    """expression / isoform_expression: sum the BED value columns.
    isoform_expression holds length-normalized TPM (assemble_bed.py
    tpm_from_counts); summed-TPM ranks equal TPM-of-pooled-counts ranks
    (the per-transcript length factor commutes with the run sum), so
    summing keeps collapsed samples on the same scale."""
    return bed[runs_present].sum(axis=1, min_count=1).to_numpy()


def collapse_isoforms(bed, runs_present, run_dirs, ref_anno, tx2gene_cache):
    """isoforms: sum QU-corrected transcript counts, recompute within-gene
    ratios. Ratios are molecule fractions: counts are length-normalized
    (count / EffectiveLength) before the within-gene division, matching the
    tpm_from_counts scale assemble_bed.py uses for single-run samples.
    The denominator uses ALL transcripts in the quants that map to
    the gene (same as assemble_bed.assemble_expression), not just the
    transcripts retained in the BED."""
    summed = _sum_salmon_runs(run_dirs, runs_present, EXPR_QU_DIR)
    if summed is None:
        return None
    # EffectiveLength from the first present run (constant across runs)
    q0 = os.path.join(run_dirs[runs_present[0]], EXPR_QU_DIR,
                      runs_present[0], "quant.sf")
    eff = pd.read_csv(q0, sep="\t").set_index("Name")["EffectiveLength"]
    eff = pd.to_numeric(eff, errors="coerce").reindex(summed.index)
    x = summed / eff
    x[(eff.isna()) | (eff <= 0)] = 0.0  # undefined length -> no abundance
    x = x.fillna(0.0)
    if tx2gene_cache.get("map") is None:
        tx2gene_cache["map"] = _load_tx2gene(ref_anno)
    tx2gene = tx2gene_cache["map"]
    genes = pd.Series(x.index.map(tx2gene), index=x.index)
    keep = genes.notna()
    x = x[keep]
    genes = genes[keep]
    totals = x.groupby(genes).sum()
    ratio = x / totals.loc[genes].to_numpy()
    pid = pd.Index(genes + "__" + x.index, name="phenotype_id")
    vec = pd.Series(ratio.to_numpy(), index=pid)
    return vec.reindex(bed["phenotype_id"]).to_numpy()


def collapse_alt(bed, runs_present, run_dirs, modality):
    """alt_TSS / alt_polyA: sum txrevise-group Salmon counts, recompute
    within-gene usage ratios per group (matching
    assemble_bed.assemble_alt_TSS_polyA)."""
    event = "upstream" if modality == "alt_TSS" else "downstream"
    parts = []
    for grp in ("grp_1", "grp_2"):
        subdir = os.path.join(ALT_DIR, f"{grp}.{event}")
        summed = _sum_salmon_runs(run_dirs, runs_present, subdir)
        if summed is None:
            return None
        genes = summed.index.str.split(".grp_").str[0]
        totals = summed.groupby(genes).sum()
        ratio = summed / totals.loc[genes].to_numpy()
        pid = summed.index.str.replace(".grp_", "__grp_", n=1, regex=False)
        pid = pid.str.replace(f".{event}.", f"_{event}_", n=1, regex=False)
        parts.append(pd.Series(ratio.to_numpy(),
                               index=pd.Index(pid, name="phenotype_id")))
    vec = pd.concat(parts)
    vec = vec[~vec.index.duplicated(keep="first")]
    return vec.reindex(bed["phenotype_id"]).to_numpy()


def parse_numers(path):
    """Parse leafcutter_perind_numers.counts.gz tolerating both header
    conventions (with or without a leading 'chrom' token).

    Returns (sample_ids, row_ids, count_matrix) where row_ids are the full
    junction strings (chrom:start:end:clu_N_strand) and count_matrix is a
    DataFrame (junctions x samples)."""
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split(" ")
        first = f.readline().strip().split(" ")
        # Samples-only header has one fewer token than data rows.
        if len(header) == len(first) - 1:
            sample_ids = header
        else:
            sample_ids = header[1:]
        row_ids = [first[0]]
        data = [[float(x) for x in first[1:]]]
        for line in f:
            parts = line.strip().split(" ")
            if len(parts) < 2:
                continue
            row_ids.append(parts[0])
            data.append([float(x) for x in parts[1:]])
    mat = pd.DataFrame(data, index=row_ids, columns=sample_ids)
    return sample_ids, row_ids, mat


def _parse_numers_row_id(row_id):
    """chrom:start:end:clu_N_strand -> (chrom, start, end, clu_n, strand)."""
    fields = row_id.split(":")
    cluster_field = fields[-1]
    end = int(fields[-2])
    start = int(fields[-3])
    chrom = ":".join(fields[:-3])
    sub = cluster_field.split("_")
    strand = sub[-1]
    clu_n = int(sub[-2])
    return chrom, start, end, clu_n, strand


def _parse_splicing_pid(pid):
    """Invert assemble_bed.assemble_splicing phenotype_id:
    gene_id + '__' + intron with ':' -> '_'. Returns
    (chrom, start, end, clu_n, strand) or None."""
    _, sep, rest = pid.partition("__")
    if not sep:
        return None
    f = rest.split("_")
    if len(f) < 6 or f[-3] != "clu":
        return None
    try:
        strand = f[-1]
        clu_n = int(f[-2])
        end = int(f[-4])
        start = int(f[-5])
        chrom = "_".join(f[:-5])
    except ValueError:
        return None
    return chrom, start, end, clu_n, strand


def collapse_splicing(bed, runs_present, numers_cache):
    """splicing: sum LeafCutter numerator counts, recompute within-cluster
    ratios. Cluster totals come from the full numers matrix (all junctions,
    including those filtered out of the BED)."""
    if numers_cache.get("mat") is None:
        return None
    mat = numers_cache["mat"]
    clu = numers_cache["clu"]
    summed = mat[runs_present].sum(axis=1)
    totals = summed.groupby(clu).sum()
    ratio = summed / totals.loc[clu].to_numpy()
    key_to_ratio = {}
    for row_id, r in zip(mat.index, ratio.to_numpy()):
        try:
            key_to_ratio[_parse_numers_row_id(row_id)] = r
        except (ValueError, IndexError):
            continue
    out = np.full(len(bed), np.nan)
    for i, pid in enumerate(bed["phenotype_id"]):
        key = _parse_splicing_pid(pid)
        if key is not None and key in key_to_ratio:
            out[i] = key_to_ratio[key]
    return out


def _load_featurecounts(path):
    """One sample's featureCounts output -> Series indexed by gene."""
    d = pd.read_csv(path, sep="\t", index_col="Geneid", skiprows=1)
    return d.iloc[:, 5].astype(float)


def collapse_stability(bed, runs_present, run_dirs):
    """stability: sum exonic and intronic featureCounts across runs,
    recompute the exon/intron ratio. The per-feature count floor (same
    threshold as assemble_bed.load_featureCounts) is applied to the SUMMED
    counts: pooling runs raises coverage, so a feature floored in a single
    run can still be reliably quantified from the pair."""
    exon = None
    intron = None
    for run in runs_present:
        ex_path = os.path.join(run_dirs[run], STAB_DIR, f"{run}.exonic.counts.txt")
        in_path = os.path.join(run_dirs[run], STAB_DIR, f"{run}.intronic.counts.txt")
        if not (os.path.exists(ex_path) and os.path.exists(in_path)):
            return None
        ex = _load_featurecounts(ex_path)
        itr = _load_featurecounts(in_path)
        exon = ex if exon is None else exon.add(ex, fill_value=0.0)
        intron = itr if intron is None else intron.add(itr, fill_value=0.0)
    exon[exon < STAB_MIN_COUNT] = np.nan
    intron[intron < STAB_MIN_COUNT] = np.nan
    ratio = exon / intron
    ratio.index.name = "phenotype_id"
    return ratio.reindex(bed["phenotype_id"]).to_numpy()


def _parse_edit_matrix(df):
    """'edit/cov' string matrix -> (num, den, fraction) DataFrames.
    Fraction = (num+0.5)/(den+0.5), NaN where den == 0 — matching
    prepare_rna_editing_phenotypes.ratios_to_fractions."""
    num = df.apply(lambda c: pd.to_numeric(c.str.split("/").str[0], errors="coerce"))
    den = df.apply(lambda c: pd.to_numeric(c.str.split("/").str[1], errors="coerce"))
    frac = (num + 0.5) / (den + 0.5)
    frac[den == 0] = np.nan
    return num, den, frac


def collapse_rna_editing(bed, runs_present, edit_cache):
    """RNA_editing: sum per-site edit/coverage counts across runs, recompute
    the smoothed per-site fraction, then average sites within each phenotype
    (cluster). Sites with zero coverage in both runs take the site's
    row-mean fraction across all samples — identical to the original
    imputation value (both runs were NaN for those sites)."""
    if edit_cache.get("num") is None:
        return None
    num = edit_cache["num"]
    den = edit_cache["den"]
    frac = edit_cache["frac"]
    site_map = edit_cache["site_map"]
    runs_present = [r for r in runs_present if r in num.columns]
    if not runs_present:
        return None
    e_sum = num[runs_present].sum(axis=1, min_count=1)
    c_sum = den[runs_present].sum(axis=1, min_count=1)
    site_frac = (e_sum + 0.5) / (c_sum + 0.5)
    missing = (c_sum == 0) | c_sum.isna()
    if missing.any():
        row_mean = frac.mean(axis=1, skipna=True)
        site_frac[missing] = row_mean[missing]
    pheno_of_site = site_map.groupby("site")["phenotype_id"].first()
    pheno = site_frac.index.map(pheno_of_site)
    keep = pheno.notna()
    vec = site_frac[keep].groupby(pheno[keep]).mean()
    vec.index.name = "phenotype_id"
    return vec.reindex(bed["phenotype_id"]).to_numpy()


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Collapse technical-replicate runs to one column per "
                    "individual in per-cohort unnorm BEDs + splicing/IR "
                    "intermediates (staging tree for pooling/harmonization).")
    p.add_argument("--metadata", required=True,
                   help="Pooled metadata TSV with rnaseq_id, array_id columns")
    p.add_argument("--ancestry-map", required=True,
                   help="TSV: sample_id, assigned_ancestry, cohort")
    p.add_argument("--config", required=True, help="config.yml")
    p.add_argument("--output-base", help="Override output_base from config")
    p.add_argument("--out-dir", required=True,
                   help="Staging root for collapsed outputs + reports/")
    p.add_argument("--ref-anno",
                   help="Normalized GTF (required when collapsing isoforms)")
    p.add_argument("--modalities", nargs="*", default=ALL_MODALITIES,
                   choices=ALL_MODALITIES, help="Default: all")
    p.add_argument("--cross-protocol", choices=["drop", "average"],
                   default="drop",
                   help="Policy for pairs whose runs come from >1 cohort "
                        "(default: drop -- the individual is removed "
                        "entirely; 'average' collapses ratio modalities "
                        "from per-run source files)")
    p.add_argument("--concordance-threshold", type=float, default=0.9,
                   help="Spearman below this -> flag pair for review (collapse "
                        "still proceeds; flag-only, never a veto)")
    p.add_argument("--concordance-min-features", type=int, default=100,
                   help="Minimum pairwise-complete features for concordance")
    p.add_argument("--dry-run", action="store_true",
                   help="Report pairs/actions only; write no staging files")
    args = p.parse_args()

    if "isoforms" in args.modalities and not args.ref_anno:
        sys.stderr.write("ERROR: --ref-anno is required to collapse isoforms "
                         "(transcript->gene map for ratio denominators)\n")
        sys.exit(1)

    output_base = args.output_base or load_output_base(args.config)
    cohorts = load_cohorts(args.config)
    if not cohorts:
        sys.stderr.write("ERROR: no cohorts found in config.yml\n")
        sys.exit(1)

    meta = load_metadata(args.metadata)
    amap = load_ancestry_map(args.ancestry_map)
    ancestry_ids = set(amap["sample_id"])
    run2ancestry = dict(zip(amap["sample_id"], amap["assigned_ancestry"]))

    pairs = discover_pairs(meta)
    run2cohort = map_runs_to_cohorts(output_base, cohorts)

    print(f"Cohorts: {cohorts}")
    print(f"Output base: {output_base}")
    print(f"Metadata: {len(meta)} rows; {len(pairs)} replicated individuals")

    pair_info = {}
    for array_id, runs in pairs.items():
        cohorts_of_pair = sorted({run2cohort.get(r, "unknown") for r in runs})
        primary = choose_primary(runs, ancestry_ids)
        ancestry = run2ancestry.get(primary)
        if ancestry is None:
            for r in runs:
                if r in run2ancestry:
                    ancestry = run2ancestry[r]
                    break
        pair_info[array_id] = {
            "runs": runs,
            "primary": primary,
            "cohorts": cohorts_of_pair,
            "cross_protocol": len(set(cohorts_of_pair)) > 1,
            "ancestry": ancestry if ancestry is not None else "unassigned",
        }
    n_cross = sum(1 for v in pair_info.values() if v["cross_protocol"])
    print(f"  Cross-protocol pairs (runs in >1 cohort): {n_cross} "
          f"(policy: {args.cross_protocol})")
    unknown_runs = [r for v in pair_info.values() for r in v["runs"]
                    if r not in run2cohort]
    if unknown_runs:
        sys.stderr.write(
            f"  WARN: {len(unknown_runs)} replicate runs not found in any "
            f"cohort samples.txt (e.g. {unknown_runs[:3]}); those runs are "
            f"ignored for cohort-level collapse\n")

    modalities = list(args.modalities)
    collapse_rows = []   # one row per pair x modality (x cohort if cross)
    concord_rows = []    # one row per pair x modality
    tx2gene_cache = {"map": None}

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    for cohort in cohorts:
        cohort_dir = os.path.join(output_base, cohort)
        same_pairs = {a: v for a, v in pair_info.items()
                      if not v["cross_protocol"]
                      and sum(1 for r in v["runs"] if run2cohort.get(r) == cohort) >= 2}
        cross_here = {a: v for a, v in pair_info.items()
                      if v["cross_protocol"]
                      and any(run2cohort.get(r) == cohort for r in v["runs"])}
        touched = bool(same_pairs or cross_here)

        staging_cohort = os.path.join(out_dir, cohort)
        if not touched:
            if not args.dry_run and not os.path.exists(staging_cohort):
                os.symlink(cohort_dir, staging_cohort)
            continue

        print(f"\n{'='*60}\nCohort {cohort}: {len(same_pairs)} same-cohort pair(s), "
              f"{len(cross_here)} cross-protocol pair(s)")

        if not args.dry_run:
            os.makedirs(os.path.join(staging_cohort, "output", "unnorm"),
                        exist_ok=True)

        # Per-pair x modality actions, consulted by the intermediate writers
        actions = {}

        # ---- BED-level pass (all modalities: collapse or concordance) ----
        for modality in modalities:
            bed_path = os.path.join(cohort_dir, "output", "unnorm", f"{modality}.bed")
            if not os.path.exists(bed_path):
                sys.stderr.write(f"  WARN: {bed_path} not found, skipping {modality}\n")
                continue
            bed = load_bed(bed_path)
            sample_cols = [c for c in bed.columns if c not in BED_META_COLS]
            stage_bed = modality in STAGED_BED_MODALITIES

            numers_cache = {"mat": None, "clu": None}
            edit_cache = {"num": None, "den": None, "frac": None, "site_map": None}
            if modality == "splicing":
                numers_path = os.path.join(cohort_dir, SPLICING_NUMERS)
                if os.path.exists(numers_path):
                    _, _, mat = parse_numers(numers_path)
                    clu = []
                    for r in mat.index:
                        try:
                            clu.append(_parse_numers_row_id(r)[3])
                        except (ValueError, IndexError):
                            clu.append(None)
                    numers_cache["mat"] = mat
                    numers_cache["clu"] = clu
                else:
                    sys.stderr.write(f"  WARN: {numers_path} not found; "
                                     f"splicing pairs will be keep-primary\n")
            if modality == "RNA_editing":
                em_path = os.path.join(cohort_dir, EDIT_MATRIX)
                sm_path = os.path.join(cohort_dir, EDIT_SITE_MAP)
                if os.path.exists(em_path) and os.path.exists(sm_path):
                    em = pd.read_csv(em_path, sep="\t", index_col=0, dtype=str)
                    num, den, frac = _parse_edit_matrix(em)
                    edit_cache.update(num=num, den=den, frac=frac,
                                      site_map=pd.read_csv(sm_path, sep="\t"))
                else:
                    sys.stderr.write(
                        f"  WARN: RNA_editing intermediates not found under "
                        f"{cohort_dir}; pairs will be keep-primary\n")

            modified = False
            for array_id, info in {**same_pairs, **cross_here}.items():
                runs_present = [r for r in info["runs"] if r in sample_cols]
                base_row = {
                    "array_id": array_id,
                    "ancestry": info["ancestry"],
                    "cohorts": ";".join(info["cohorts"]),
                    "runs": ";".join(info["runs"]),
                    "primary_run": info["primary"],
                    "modality": modality,
                }
                if not runs_present:
                    continue

                cross = info["cross_protocol"]
                rep = representative_run(info, runs_present)
                note = ""
                if not cross and info["primary"] not in runs_present:
                    note = f"primary {info['primary']} absent from BED; kept {rep}"

                # ---- Decide action ----
                pearson = spearman = np.nan
                n_feat = 0
                handler_runs = runs_present
                if cross:
                    # Cross-cohort individuals are labeling errors by
                    # policy: the individual is removed entirely (all runs,
                    # every cohort, and the collapsed ancestry map). The
                    # '--cross-protocol average' escape hatch collapses only
                    # the ratio modalities whose values are recomputed from
                    # per-run source files (Salmon quants, featureCounts,
                    # edit matrix); expression / isoform_expression are
                    # summed from BED columns, only comparable within a
                    # cohort's feature set.
                    can_average = (args.cross_protocol == "average"
                                   and modality in RATIO_BED_MODALITIES
                                   and rep in runs_present)
                    action = "collapsed" if can_average else "dropped_cross_protocol"
                    if can_average:
                        handler_runs = [r for r in info["runs"] if r in run2cohort]
                    else:
                        note = (note + "; " if note else "") + \
                            "individual dropped: runs span >1 cohort"
                    if rep in runs_present:
                        # Record the cross-protocol pair once (from the
                        # representative's cohort) in the concordance report.
                        concord_rows.append({
                            "array_id": array_id, "ancestry": info["ancestry"],
                            "cohort": cohort, "run1": info["runs"][0],
                            "run2": info["runs"][1], "modality": modality,
                            "n_features": 0, "pearson": np.nan,
                            "spearman": np.nan, "flag": "cross_protocol",
                        })
                elif len(runs_present) == 1:
                    action = "single_run_present"
                else:
                    pearson, spearman, n_feat = concordance(
                        bed, runs_present[0], runs_present[1])
                    if n_feat < args.concordance_min_features:
                        flag = "insufficient_data"
                    elif not np.isnan(spearman) and \
                            spearman < args.concordance_threshold:
                        flag = "discordant"
                    else:
                        flag = "ok"
                    concord_rows.append({
                        "array_id": array_id, "ancestry": info["ancestry"],
                        "cohort": cohort, "run1": runs_present[0],
                        "run2": runs_present[1], "modality": modality,
                        "n_features": n_feat, "pearson": pearson,
                        "spearman": spearman, "flag": flag,
                    })
                    # Technical replicates are always averaged;
                    # concordance flags are for review in the report,
                    # never a veto.
                    action = "collapsed"
                    if flag != "ok":
                        note = (note + "; " if note else "") + \
                            f"concordance flag: {flag}"

                # ---- Execute ----
                if action == "collapsed" and stage_bed:
                    run_dirs = {r: os.path.join(output_base, run2cohort[r])
                                for r in handler_runs if r in run2cohort}
                    if modality in COUNT_MODALITIES:
                        vec = collapse_counts(bed, handler_runs)
                    elif modality == "isoforms":
                        vec = collapse_isoforms(bed, handler_runs, run_dirs,
                                                args.ref_anno, tx2gene_cache)
                    elif modality in ("alt_TSS", "alt_polyA"):
                        vec = collapse_alt(bed, handler_runs, run_dirs, modality)
                    elif modality == "stability":
                        vec = collapse_stability(bed, handler_runs, run_dirs)
                    elif modality == "RNA_editing":
                        vec = collapse_rna_editing(bed, handler_runs, edit_cache)
                    else:
                        vec = None
                    if vec is None:
                        # Per-run source files missing: average the BED
                        # columns directly (same cohort -> identical
                        # feature set) rather than discarding a run.
                        vec = bed[runs_present].mean(axis=1).to_numpy()
                        note = (note + "; " if note else "") + \
                            "per-run sources missing; averaged BED columns"
                    bed[rep] = vec
                    modified = True

                # splicing collapse is realized at the numers level and
                # intron_retention at the PSI level (sections below); for
                # staged-BED modalities the non-representative columns are
                # dropped here.
                actions[(array_id, modality)] = action
                if stage_bed:
                    if action == "dropped_cross_protocol":
                        drop_cols = list(runs_present)  # individual removed
                    else:
                        drop_cols = [r for r in runs_present if r != rep]
                    if drop_cols:
                        bed = bed.drop(columns=drop_cols)
                        sample_cols = [c for c in sample_cols if c not in drop_cols]
                        modified = True
                kept = rep if (rep in runs_present
                               and action != "dropped_cross_protocol") else ""
                collapse_rows.append({**base_row, "action": action,
                                      "kept_column": kept,
                                      "spearman": spearman, "note": note})

            if stage_bed and not args.dry_run:
                out_bed = os.path.join(staging_cohort, "output", "unnorm",
                                       f"{modality}.bed")
                write_bed(bed, out_bed)
                tag = "collapsed" if modified else "unchanged"
                print(f"  {modality}: wrote staged BED ({tag}; "
                      f"{bed.shape[0]} features x "
                      f"{bed.shape[1] - len(BED_META_COLS)} samples)")

        # ---- Splicing numers (consumed by 17_harmonize_within_ancestry.sh) ----
        if "splicing" in modalities:
            numers_path = os.path.join(cohort_dir, SPLICING_NUMERS)
            if os.path.exists(numers_path) and not args.dry_run:
                sample_ids, row_ids, mat = parse_numers(numers_path)
                mat_out = mat.copy()
                drop = []
                for array_id, info in {**same_pairs, **cross_here}.items():
                    runs_present = [r for r in info["runs"] if r in mat_out.columns]
                    if not runs_present:
                        continue
                    rep = representative_run(info, runs_present)
                    action = actions.get((array_id, "splicing"), "collapsed")
                    if action == "dropped_cross_protocol":
                        drop.extend(runs_present)  # individual removed
                        continue
                    if action == "collapsed" and not info["cross_protocol"]:
                        same_cohort_runs = [r for r in runs_present
                                            if run2cohort.get(r) == cohort]
                        if len(same_cohort_runs) >= 2:
                            mat_out[rep] = mat[same_cohort_runs].sum(axis=1)
                    drop.extend([r for r in runs_present if r != rep])
                mat_out = mat_out.drop(columns=list(dict.fromkeys(drop)))
                out_num = os.path.join(staging_cohort, SPLICING_NUMERS)
                os.makedirs(os.path.dirname(out_num), exist_ok=True)
                with gzip.open(out_num, "wt") as f:
                    # Samples-only header (no 'chrom' token), matching the
                    # per-cohort numers convention expected by assemble_bed.py.
                    f.write(" ".join(mat_out.columns) + "\n")
                    for row_id, vals in zip(mat_out.index, mat_out.to_numpy()):
                        f.write(row_id + " " +
                                " ".join(str(int(v)) for v in vals) + "\n")
                print(f"  splicing numers: wrote staged file "
                      f"({mat_out.shape[1]} samples)")

        # ---- Intron retention PSI (averaged: MAJIQ exposes no per-run ----
        # ---- counts, so same-cohort run PSI columns are averaged) ----
        if "intron_retention" in modalities:
            psi_path = os.path.join(cohort_dir, IR_PSI)
            if os.path.exists(psi_path) and not args.dry_run:
                psi = pd.read_csv(psi_path, sep="\t", dtype={"seqid": str})
                drop = []
                for array_id, info in {**same_pairs, **cross_here}.items():
                    runs_present = [r for r in info["runs"] if r in psi.columns]
                    if not runs_present:
                        continue
                    rep = representative_run(info, runs_present)
                    action = actions.get((array_id, "intron_retention"),
                                         "collapsed")
                    if action == "dropped_cross_protocol":
                        drop.extend(runs_present)  # individual removed
                        continue
                    if action == "collapsed" and not info["cross_protocol"]:
                        same_cohort_runs = [r for r in runs_present
                                            if run2cohort.get(r) == cohort]
                        if len(same_cohort_runs) >= 2:
                            psi[rep] = psi[same_cohort_runs].mean(axis=1)
                    drop.extend([r for r in runs_present if r != rep])
                psi = psi.drop(columns=list(dict.fromkeys(drop)))
                out_psi = os.path.join(staging_cohort, IR_PSI)
                os.makedirs(os.path.dirname(out_psi), exist_ok=True)
                psi.to_csv(out_psi, sep="\t", index=False,
                           float_format="%g", compression="gzip")
                print(f"  intron_retention PSI: wrote staged file "
                      f"(dropped {len(set(drop))} run column(s))")

    # ---- Reports ----
    reports_dir = os.path.join(out_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    coll_df = pd.DataFrame(collapse_rows)
    conc_df = pd.DataFrame(concord_rows)

    if not coll_df.empty:
        coll_df = coll_df.sort_values(["ancestry", "array_id", "modality"])
        coll_df.to_csv(os.path.join(reports_dir, "replicate_collapses.tsv"),
                       sep="\t", index=False)
        for anc, grp in coll_df.groupby("ancestry"):
            grp.to_csv(os.path.join(reports_dir, f"{anc}_replicate_collapses.tsv"),
                       sep="\t", index=False)
    if not conc_df.empty:
        conc_df = conc_df.sort_values(["ancestry", "array_id", "modality"])
        conc_df.to_csv(os.path.join(reports_dir, "replicate_concordance.tsv"),
                       sep="\t", index=False)
        for anc, grp in conc_df.groupby("ancestry"):
            grp.to_csv(os.path.join(reports_dir, f"{anc}_replicate_concordance.tsv"),
                       sep="\t", index=False)

    # Collapsed ancestry map: drop non-primary runs of same-cohort pairs;
    # cross-protocol individuals are removed entirely under the drop policy.
    drop_runs = set()
    for info in pair_info.values():
        if info["cross_protocol"] and args.cross_protocol == "drop":
            drop_runs.update(info["runs"])
        else:
            drop_runs.update(r for r in info["runs"] if r != info["primary"])
    amap_out = amap[~amap["sample_id"].isin(drop_runs)].copy()
    amap_path = os.path.join(reports_dir, "ancestry_map_collapsed.tsv")
    amap_out.to_csv(amap_path, sep="\t", index=False)

    # Duplicate-representation check: an individual should reduce to one
    # kept column per modality across all cohorts.
    if not coll_df.empty:
        kept = coll_df[coll_df["kept_column"].astype(str) != ""]
        dup = kept.groupby(["array_id", "modality"])["kept_column"].nunique()
        dup = dup[dup > 1]
        if len(dup):
            sys.stderr.write(
                f"  WARN: {len(dup)} individual x modality combinations keep "
                f">1 column across cohorts; review replicate_collapses.tsv: "
                f"{list(dup.index[:5])}\n")

    # ---- Summary ----
    print(f"\n{'='*60}\nCOLLAPSE SUMMARY")
    print(f"  Replicated individuals: {len(pair_info)} "
          f"({n_cross} cross-protocol)")
    if not coll_df.empty:
        for action, n in coll_df["action"].value_counts().items():
            print(f"    {action}: {n}")
    if not conc_df.empty:
        n_disc = int((conc_df["flag"] == "discordant").sum())
        n_ins = int((conc_df["flag"] == "insufficient_data").sum())
        print(f"  Concordance: {len(conc_df)} pair x modality tests; "
              f"{n_disc} discordant (Spearman < {args.concordance_threshold}), "
              f"{n_ins} insufficient-data")
    print(f"  Reports: {reports_dir}")
    print(f"  Collapsed ancestry map: {amap_path} "
          f"({len(amap)} -> {len(amap_out)} samples)")
    if args.dry_run:
        print("  DRY RUN — no staging files written")


if __name__ == "__main__":
    main()
