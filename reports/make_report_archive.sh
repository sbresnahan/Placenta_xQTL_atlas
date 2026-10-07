#!/bin/bash
# make_report_archive.sh — assemble every file report_placenta_xqtl.Rmd needs
# into a single tarball, in the exact relative layout the Rmd expects
# (data/results, data/qc/..., gene_map.tsv at the root).
#
# Run on seadragon AFTER modules 28 (grouped, ungrouped, grouped-independent,
# ungrouped-independent) and 29 have completed for every ancestry:
#
#   export CONFIG=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml
#   bash reports/make_report_archive.sh
#
# Then copy placenta_xqtl_report_inputs_<date>.tar.gz off the cluster and
# unpack it into reports/ before rendering:
#   tar xzf placenta_xqtl_report_inputs_<date>.tar.gz -C reports/
#
# Optional env vars:
#   SCRIPTS_DIR  — explicit scripts dir (default: self-located ../05_qtl_mapping)
#   OUTPUT_BASE  — default: from config.yml
#   QTL_DIR      — default: ${OUTPUT_BASE}/qtl_inputs
#   RESULTS_DIR  — default: ${OUTPUT_BASE}/qtl_results
#   PC_DIR       — default: ${OUTPUT_BASE}/genotype_pcs
#   ANCESTRIES   — space-separated labels (default: "EAS EUR")
#   ANCESTRY_MAP — default: ${OUTPUT_BASE%/*}/pooled/pooled_sample_ancestry_RNAseq.tsv
#   METADATA_TSV — default: ${OUTPUT_BASE%/*}/pooled/placenta_QTL_cohort_metadata.tsv
#   NORM_DIR     — default: ${OUTPUT_BASE}/normalized_modalities/int
#   GENO_BASE     — module-01 data root (default: parent of OUTPUT_BASE;
#                   e.g. /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL)
#   GENO_COHORTS  — cohort names used by module 01
#                   (default: "GUSTO NIEHS_RICHS SNUH MALI_G3A NIGMS")
#   POOLED_GENO_DIR — pooled genotype directory
#                   (default: ${GENO_BASE}/pooled/genotypes)
#   OUT          — output tarball path (default: ./placenta_xqtl_report_inputs_<date>.tar.gz)
#   KEEP_STAGING — 1 to keep the staging directory after archiving (default: 0)
#
# Exit status: 0 on success (COMPLETE or INCOMPLETE manifest); 1 only when a
# core input (per-ancestry covariates table, expression grouped top table,
# or the normalized GTF) is missing — a report cannot be rendered without
# those. All other gaps are archived as-is and listed loudly in MANIFEST.txt.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
SCRIPTS_DIR="${SCRIPTS_DIR:-$REPO_DIR/05_qtl_mapping}"
CONFIG="${CONFIG:?ERROR: CONFIG env var required (path to config.yml)}"

CONFIG_GET="$SCRIPTS_DIR/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="$REPO_DIR/03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $REPO_DIR/03_phenotyping" >&2
    exit 1
fi
eval "$(python3 "$CONFIG_GET" "$CONFIG")"

OUTPUT_BASE="${OUTPUT_BASE:?ERROR: OUTPUT_BASE not set (env or config.yml)}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
PC_DIR="${PC_DIR:-${OUTPUT_BASE}/genotype_pcs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
GTF="${NORMALIZED_GTF:?ERROR: normalized_gtf not found in config.yml}"
DATE_TAG="$(date +%Y%m%d)"
STAGING="placenta_xqtl_report_inputs_${DATE_TAG}"
OUT="${OUT:-${PWD}/placenta_xqtl_report_inputs_${DATE_TAG}.tar.gz}"
KEEP_STAGING="${KEEP_STAGING:-0}"

echo "== make_report_archive.sh =="
echo "  CONFIG:      $CONFIG"
echo "  OUTPUT_BASE: $OUTPUT_BASE"
echo "  RESULTS_DIR: $RESULTS_DIR"
echo "  QTL_DIR:     $QTL_DIR"
echo "  PC_DIR:      $PC_DIR"
echo "  GTF:         $GTF"
echo "  ANCESTRIES:  $ANCESTRIES"
echo "  Staging:     $STAGING"
echo "  Output:      $OUT"

rm -rf "$STAGING"
mkdir -p "$STAGING/data/results" "$STAGING/data/qc/qtl_inputs" \
         "$STAGING/data/qc/hcp_optimization" \
         "$STAGING/data/qc/hcp_optimization_modalities" \
         "$STAGING/data/qc/genotype_pcs" "$STAGING/data/extras"

MISSING_CORE=()
MISSING_OTHER=()

copy_req() { # copy_req <src> <dst_dir> <core|other> <label>
    local src="$1" dst="$2" level="$3" label="$4"
    if [ -f "$src" ]; then
        cp "$src" "$dst/"
    elif [ "$level" = core ]; then
        MISSING_CORE+=("$label")
    else
        MISSING_OTHER+=("$label")
    fi
}

# ---- 1. Result top tables (27/29 outputs) -------------------------------
# Grouped (9 modalities), ungrouped (8), grouped-independent (9),
# ungrouped-independent (8), and combined top TSVs per ancestry — the glob
# copies every top table present; the audit below verifies expected layers.
n_top=0
for f in "$RESULTS_DIR"/*_top.tsv; do
    [ -f "$f" ] || continue
    cp "$f" "$STAGING/data/results/"
    n_top=$((n_top + 1))
done
echo "  results: staged $n_top top tables"

# Core: per-ancestry covariates + expression grouped top table (report
# cannot render without these).
for ANC in $ANCESTRIES; do
    copy_req "$RESULTS_DIR/${ANC}_expression_cisqtl_top.tsv" \
             "$STAGING/data/results" core "${ANC}_expression_cisqtl_top.tsv"
done

# Expected-completeness audit (non-core): grouped 9 + ungrouped 8 +
# grouped-independent 9 + ungrouped-independent 8 + combined 1 +
# combined-independent 1 = 36 top tables per ancestry.
# modality file stems as written by 27/28/29 (display-case: alt_TSS, alt_polyA, RNA_editing)
MODS_GROUPED="expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability"
MODS_UNGROUPED="isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability"
for ANC in $ANCESTRIES; do
    for M in $MODS_GROUPED; do
        copy_req "$RESULTS_DIR/${ANC}_${M}_cisqtl_top.tsv" \
                 "$STAGING/data/results" other "grouped:${ANC}_${M}"
        copy_req "$RESULTS_DIR/${ANC}_${M}_cisqtl_independent_top.tsv" \
                 "$STAGING/data/results" other "independent:${ANC}_${M}"
    done
    for M in $MODS_UNGROUPED; do
        copy_req "$RESULTS_DIR/${ANC}_${M}_ungrouped_cisqtl_top.tsv" \
                 "$STAGING/data/results" other "ungrouped:${ANC}_${M}"
        copy_req "$RESULTS_DIR/${ANC}_${M}_ungrouped_cisqtl_independent_top.tsv" \
                 "$STAGING/data/results" other "ungrouped-independent:${ANC}_${M}"
    done
    copy_req "$RESULTS_DIR/${ANC}_combined_cisqtl_top.tsv" \
             "$STAGING/data/results" other "combined:${ANC}"
    copy_req "$RESULTS_DIR/${ANC}_combined_cisqtl_independent_top.tsv" \
             "$STAGING/data/results" other "independent:${ANC}_combined"
done

# ---- 1.5 Report extras (35_extract_report_extras.py outputs) --------------
# Targeted extractions powering the DevBrain-style figure panels: cross-
# ancestry and cross-modality lead-pair lookups, showcase-locus regional
# scans, lead-variant annotations, gnomAD constraint, PC-AiR PCs. All
# non-core: the report renders without them (panels skip with a note).
EXTRAS_DIR="$RESULTS_DIR/report_extras"
if [ -d "$EXTRAS_DIR" ]; then
    n_extras=0
    for f in "$EXTRAS_DIR"/*.tsv "$EXTRAS_DIR"/*.tsv.gz; do
        [ -f "$f" ] || continue
        cp "$f" "$STAGING/data/extras/"
        n_extras=$((n_extras + 1))
    done
    echo "  extras: staged $n_extras files from report_extras/"
    for label in cross_ancestry_lookup.tsv.gz cross_modality_lookup.tsv.gz \
                 lead_variant_annotations.tsv.gz gene_constraint.tsv; do
        copy_req "$EXTRAS_DIR/$label" "$STAGING/data/extras" other \
                 "extras/$label (run 35_extract_report_extras.py)"
    done
    # SuSHiE aggregation summaries (module 33) — power the conditional
    # fine-mapping panels; present only after the fine-mapping run completes.
    copy_req "$RESULTS_DIR/finemap/aggregated/finemap_locus_summary.tsv" \
             "$STAGING/data/extras" other \
             "extras/finemap_locus_summary.tsv (run 33_aggregate_finemap.py)"
    copy_req "$RESULTS_DIR/finemap/aggregated/finemap_credible_sets.tsv.gz" \
             "$STAGING/data/extras" other \
             "extras/finemap_credible_sets.tsv.gz (run 33_aggregate_finemap.py)"
else
    echo "  extras: $EXTRAS_DIR not found — run 35_extract_report_extras.py"
    echo "          (DevBrain-style panels will be skipped in the report)"
    MISSING_OTHER+=("extras/ (35_extract_report_extras.py not run)")
fi

# ---- 2. QC inputs ---------------------------------------------------------
for ANC in $ANCESTRIES; do
    copy_req "$QTL_DIR/${ANC}_covariates.tsv" \
             "$STAGING/data/qc/qtl_inputs" core "${ANC}_covariates.tsv"
    copy_req "$QTL_DIR/${ANC}_covariate_pruning.tsv" \
             "$STAGING/data/qc/qtl_inputs" other "${ANC}_covariate_pruning.tsv"
    copy_req "$QTL_DIR/${ANC}_outliers.tsv" \
             "$STAGING/data/qc/qtl_inputs" other "${ANC}_outliers.tsv"
    copy_req "$QTL_DIR/${ANC}_covariate_correlation.png" \
             "$STAGING/data/qc/qtl_inputs" other "${ANC}_covariate_correlation.png"
    # Per-modality HCP optimization (module 25b): one table + curve per
    # ancestry x modality group (9 modalities + combined). The older
    # per-ancestry hcp_optimization/{ANC}_optimal_hcp.{tsv,png} outputs are
    # superseded by these and are no longer staged.
    for f in "$QTL_DIR/hcp_optimization_modalities/${ANC}_"*"_optimal_hcp.tsv"; do
        [ -f "$f" ] || continue
        copy_req "$f" "$STAGING/data/qc/hcp_optimization_modalities" other \
                 "hcp_optimization_modalities/$(basename "$f")"
    done
    for f in "$QTL_DIR/hcp_optimization_modalities/${ANC}_"*"_optimal_hcp.png"; do
        [ -f "$f" ] || continue
        copy_req "$f" "$STAGING/data/qc/hcp_optimization_modalities" other \
                 "hcp_optimization_modalities/$(basename "$f")"
    done
    copy_req "$PC_DIR/${ANC}_genotype_pcs_scree.png" \
             "$STAGING/data/qc/genotype_pcs" other "${ANC}_genotype_pcs_scree.png"
    copy_req "$PC_DIR/${ANC}_genotype_pcs.tsv" \
             "$STAGING/data/qc/genotype_pcs" other "${ANC}_genotype_pcs.tsv"
done

# ---- 3. Gene map + gene bodies from the normalized GTF --------------------
# Regenerated in full at archive time: every gene in the GTF,
# version-stripped IDs matching the pipeline's BED files.
if [ ! -f "$GTF" ]; then
    MISSING_CORE+=("normalized_gtf ($GTF)")
else
    python3 - "$GTF" "$STAGING/gene_map.tsv" "$STAGING/data/gene_bodies.tsv" <<'PYEOF'
import re
import sys

gtf_path, map_out, bodies_out = sys.argv[1:4]
attr_re = re.compile(r'(\S+) "([^"]*)"')

seen = {}
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(attr_re.findall(f[8]))
        gid = attrs.get("gene_id", "").split(".")[0]
        if not gid or gid in seen:
            continue
        seen[gid] = (
            attrs.get("gene_name", ""),
            attrs.get("gene_type", attrs.get("gene_biotype", "")),
            f[0], f[3], f[4], f[6],
        )

with open(map_out, "w") as out:
    out.write("gene_id\tsymbol\tbiotype\tdescription\n")
    for gid, (sym, bio, _c, _s, _e, _t) in sorted(seen.items()):
        out.write(f"{gid}\t{sym}\t{bio}\tNA\n")

with open(bodies_out, "w") as out:
    out.write("gene_id\tchr\tstart\tend\tstrand\n")
    for gid, (_sym, _bio, chrom, start, end, strand) in sorted(seen.items()):
        out.write(f"{gid}\t{chrom}\t{start}\t{end}\t{1 if strand == '+' else -1}\n")

print(f"  gene_map.tsv / gene_bodies.tsv: {len(seen)} genes from {gtf_path}")
PYEOF
fi

# ---- 3.5 Sample-attrition log (raw FASTQ -> QTL mapping) -------------------
# Per-sample matrix + stage-level summary covering every stage where samples
# can drop, in chronological pipeline order: raw FASTQ staging, Salmon
# quant, QU correction, assembled BEDs, within-cohort INT (per ancestry x modality;
# runs upstream of the genotype intersection on one run per individual),
# genotype intersection + outlier exclusion ({ANC}_metadata.tsv from
# scripts 23+24), the final covariate mapping set, and QTL inputs (per
# modality). The summary counts every stage at two levels: run (per
# rnaseq_id) and individual (unique array_id). BED column namespaces from
# upstream phenotyping stages are stripped before membership tests. All
# inputs are non-core: gaps produce warnings and zero-filled stage columns,
# never a failed archive.
ATTRITION_DIR="$STAGING/data/qc/attrition"
mkdir -p "$ATTRITION_DIR"
NORM_DIR="${NORM_DIR:-${OUTPUT_BASE}/normalized_modalities/int}"
ANCESTRY_MAP="${ANCESTRY_MAP:-${OUTPUT_BASE%/*}/pooled/pooled_sample_ancestry_RNAseq.tsv}"
METADATA_TSV="${METADATA_TSV:-${OUTPUT_BASE%/*}/pooled/placenta_QTL_cohort_metadata.tsv}"

echo "  ANCESTRY_MAP: $ANCESTRY_MAP"
echo "  METADATA_TSV: $METADATA_TSV"
echo "  NORM_DIR:     $NORM_DIR"

if ! python3 - "$ATTRITION_DIR" "$CONFIG" "$OUTPUT_BASE" "$QTL_DIR" \
         "$NORM_DIR" "$ANCESTRY_MAP" "$METADATA_TSV" "$ANCESTRIES" <<'PYEOF'
import csv
import gzip
import os
import re
import sys
attrition_dir, config_path, output_base, qtl_dir, norm_dir, \
    ancestry_map, metadata_tsv, ancestries = sys.argv[1:9]
ancestries = ancestries.split()
MODS = ["expression", "isoforms", "isoform_expression", "splicing",
        "intron_retention", "alt_TSS", "alt_polyA", "RNA_editing", "stability"]
EXPR_MODS = ["expression", "isoforms", "isoform_expression"]


def warn(msg):
    print(f"  ATTRITION WARN: {msg}", file=sys.stderr)


# ---- minimal config.yml parser (cohorts section only; config_get.py style) --
def parse_cohorts(path):
    cohorts, cur, in_cohorts = {}, None, False
    for raw in open(path):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        key = line.strip()
        if indent == 0:
            in_cohorts = key == "cohorts:"
            continue
        if not in_cohorts:
            continue
        if indent == 2 and key.endswith(":"):
            cur = key[:-1]
            cohorts[cur] = {}
        elif indent >= 4 and cur and ":" in key:
            k, v = key.split(":", 1)
            cohorts[cur][k.strip()] = v.strip()
    return cohorts


def read_samples(path):
    try:
        with open(path) as fh:
            return [l.strip() for l in fh if l.strip()]
    except OSError:
        return []


def tsv_columns(path, skip):
    """Header fields of a BED/TSV after the first `skip` columns."""
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        return [c for c in header[skip:] if c]
    except OSError:
        return None


def read_table(path):
    try:
        with open(path) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return None


def strip_namespace(col, id_set):
    """Try stripping a leading '{prefix}_' namespace from a sample column;
    return the matching base sample ID, or None. Tries every underscore
    split point (shortest prefix first) so multi-underscore prefixes such
    as 'NIEHS_RICHS_SAMPLE1' are handled. Ported from
    26_harmonize_modalities.py."""
    parts = col.split('_')
    for i in range(1, len(parts)):
        cand = '_'.join(parts[i:])
        if cand in id_set:
            return cand
    return None


def normalize_cols(cols, id_set):
    """Map BED sample columns to base IDs, stripping stage/cohort namespace
    prefixes where needed. Returns (normalized_set, n_stripped)."""
    out, n_stripped = set(), 0
    for c in cols:
        if c in id_set:
            out.add(c)
            continue
        base = strip_namespace(c, id_set)
        if base is not None:
            out.add(base)
            n_stripped += 1
        else:
            out.add(c)  # unmatched columns simply never match a sample
    return out, n_stripped


cohorts = parse_cohorts(config_path)
if not cohorts:
    sys.exit("no cohorts parsed from config.yml")

# ---- per-cohort raw/quant stage evidence ------------------------------------
# fastq_map: tab-delimited, 2 fields (R1, sample_id) or 3 (R1, R2, sample_id);
# paths relative to the cohort's fastq_dir; no header (script 02 convention).
fastq = {}      # sample_id -> (r1_present, r2_present or None)
salmon = {}     # sample_id -> quant.sf exists
qu = {}         # sample_id -> QU-adjusted quant.sf exists
cohort_of_samples = {}
for coh, cfg in sorted(cohorts.items()):
    fq_dir = cfg.get("fastq_dir", "")
    fmap = cfg.get("fastq_map", "")
    if fmap and os.path.isfile(fmap):
        for line in open(fmap):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            sid = fields[-1].strip()
            if sid in ("sample_id",) or fields[0].strip() in ("R1_path", "R1"):
                continue  # header-ish line
            r1 = os.path.isfile(os.path.join(fq_dir, fields[0].strip())) \
                if fields[0].strip() else False
            r2 = os.path.isfile(os.path.join(fq_dir, fields[1].strip())) \
                if len(fields) >= 3 and fields[1].strip() else None
            fastq[sid] = (r1, r2)
    else:
        warn(f"fastq_map not found for {coh}: {fmap}")
    sfile = cfg.get("samples_file", "")
    sids = read_samples(sfile)
    if not sfile:
        warn(f"samples_file missing for {coh}")
    for sid in sids:
        cohort_of_samples.setdefault(sid, coh)
    expr_dir = os.path.join(output_base, coh, "intermediate", "expression")
    qu_dir = os.path.join(output_base, coh, "intermediate", "expression_qu")
    for sid in sids:
        salmon[sid] = os.path.isfile(os.path.join(expr_dir, sid, "quant.sf"))
        qu[sid] = os.path.isfile(os.path.join(qu_dir, sid, "quant.sf"))

# ---- assembled per-cohort BED columns (rnaseq_id space) ----------------------
bed_cols = {}
for mod in EXPR_MODS:
    for coh in cohorts:
        cols = tsv_columns(
            os.path.join(output_base, coh, "output", f"{mod}.bed.gz"), 4)
        if cols is None:
            warn(f"assembled BED not found: {coh}/{mod}.bed.gz")
        bed_cols[(coh, mod)] = set(cols or [])
# Union across cohorts: a sample appears in exactly one cohort's BED, and the
# pooled metadata's cohort column uses study names (not config keys), so
# per-cohort lookup via the sample's cohort is unreliable.
bed_union = {m: set().union(*(cols for (c, mm), cols in bed_cols.items()
                              if mm == m))
             for m in EXPR_MODS}

# ---- ancestry map ------------------------------------------------------------
ancestry_of = {}
rows = read_table(ancestry_map) if os.path.isfile(ancestry_map) else None
if rows:
    sample_col = next((c for c in ("sample_id", "rnaseq_id")
                       if c in rows[0]), list(rows[0].keys())[0])
    anc_col = next((c for c in rows[0] if "ancestry" in c.lower()), None)
    for r in rows:
        ancestry_of[r[sample_col]] = (r.get(anc_col) or "NA") if anc_col else "NA"
else:
    warn(f"ancestry map not found: {ancestry_map}")

# ---- pooled metadata (rnaseq_id -> array_id, cohort) -------------------------
array_of, meta_cohort_of = {}, {}
rows = read_table(metadata_tsv) if os.path.isfile(metadata_tsv) else None
if rows:
    coh_col = next((c for c in rows[0] if c.lower() in ("cohort", "cohort_id")),
                   None)
    for r in rows:
        rid = r.get("rnaseq_id")
        if rid and rid not in ("rnaseq_id", "sample_id"):
            array_of[rid] = r.get("array_id", "")
            if coh_col:
                meta_cohort_of[rid] = r.get(coh_col, "")
else:
    warn(f"pooled metadata not found: {metadata_tsv}")

# ---- per-ancestry final sets -------------------------------------------------
meta_rnaseq = {}    # ANC -> set of rnaseq_id (post-intersection, post-outlier)
meta_array = {}     # ANC -> {rnaseq_id: array_id}
cov_cols = {}       # ANC -> set of array_id (final mapping set)
int_cols = {}       # (ANC, mod) -> set of array_id
int_imputed = []    # modalities whose INT stage was carried forward
qtl_cols = {}       # (ANC, mod) -> set of array_id
for anc in ancestries:
    mrows = read_table(os.path.join(qtl_dir, f"{anc}_metadata.tsv"))
    if mrows and "rnaseq_id" in mrows[0]:
        meta_rnaseq[anc] = {r["rnaseq_id"] for r in mrows}
        meta_array[anc] = {r["rnaseq_id"]: r.get("array_id", "")
                           for r in mrows if r.get("array_id")}
    else:
        warn(f"{anc}_metadata.tsv missing or lacks rnaseq_id (run 23+24)")
    cc = tsv_columns(os.path.join(qtl_dir, f"{anc}_covariates.tsv"), 1)
    if cc is None:
        warn(f"{anc}_covariates.tsv not found")
    cov_cols[anc] = set(cc or [])
    for mod in MODS:
        # TMM->VST schema: isoform_expression normalizes to _vst.bed; the
        # ratio modalities keep _int.bed. Try VST first (canonical when both
        # exist in a transitional tree).
        for pat, skip in ((f"{anc}_{mod}_vst.bed", 4),
                          (f"{anc}_{mod}_vst.bed.gz", 4),
                          (f"{anc}_{mod}_int.bed", 4),
                          (f"{anc}_{mod}_int.bed.gz", 4)):
            p = os.path.join(norm_dir, pat)
            if os.path.isfile(p):
                int_cols[(anc, mod)] = set(tsv_columns(p, skip) or [])
                break
        else:
            # Normalized BED missing (e.g. intermediates cleaned): carry forward
            # the previous stage — assume every post-intersection/post-outlier
            # sample was retained through normalization (stated in the report).
            warn(f"normalized BED not found: {anc}_{mod} — carrying forward "
                 f"metadata-stage membership")
            int_cols[(anc, mod)] = set(meta_array.get(anc, {}).values())
            int_imputed.append(f"{anc}_{mod}")
        for pat, skip in ((f"{anc}_{mod}_harmonized.bed", 4),
                          (f"{anc}_{mod}.bed.gz", 4)):
            p = os.path.join(qtl_dir, pat)
            if os.path.isfile(p):
                qtl_cols[(anc, mod)] = set(tsv_columns(p, skip) or [])
                break
        else:
            warn(f"QTL-input BED not found: {anc}_{mod}")

# ---- normalize BED column ID spaces ------------------------------------------
# INT and QTL-input BED columns can carry stage/cohort namespace
# prefixes (e.g. 'NIEHS_RICHS_SRR...' written by upstream phenotyping
# stages). Strip them against the known ID universe so membership tests
# match; without this, namespaced stages silently score 0.
known_ids = set(array_of) | {v for v in array_of.values() if v}
for anc in ancestries:
    known_ids |= meta_rnaseq.get(anc, set())
    known_ids |= {v for v in meta_array.get(anc, {}).values() if v}
    known_ids |= cov_cols.get(anc, set())
ns_notes = []
for key, cols in list(int_cols.items()):
    int_cols[key], n_s = normalize_cols(cols, known_ids)
    if n_s:
        ns_notes.append(f"int {key[0]}_{key[1]} ({n_s} cols)")
for key, cols in list(qtl_cols.items()):
    qtl_cols[key], n_s = normalize_cols(cols, known_ids)
    if n_s:
        ns_notes.append(f"qtl {key[0]}_{key[1]} ({n_s} cols)")

# ---- per-sample matrix -------------------------------------------------------
all_ids = []
seen = set()


def add_id(sid):
    if sid and sid not in seen:
        seen.add(sid)
        all_ids.append(sid)


for src in (fastq, salmon, qu, ancestry_of, array_of):
    for sid in src:
        add_id(sid)
for (coh, _mod), cols in bed_cols.items():
    for sid in cols:
        add_id(sid)
for anc in ancestries:
    for sid in meta_rnaseq.get(anc, set()):
        add_id(sid)

# Chronological pipeline order: normalization (INT) runs upstream of the
# genotype intersection (one run per individual), so in_int_* precedes
# in_metadata_*; covariates are the final per-individual mapping set.
stage_cols = (["fastq_r1_present", "fastq_r2_present", "has_quant_sf",
               "has_qu_quant_sf"] + [f"in_{m}_bed" for m in EXPR_MODS]
              + [f"in_int_{a}_{m}" for a in ancestries for m in MODS]
              + [f"in_metadata_{a}" for a in ancestries]
              + [f"in_covariates_{a}" for a in ancestries]
              + [f"in_qtl_{a}_{m}" for a in ancestries for m in MODS])

matrix = []
for sid in all_ids:
    coh = meta_cohort_of.get(sid) or cohort_of_samples.get(sid, "")
    arr = array_of.get(sid) or next(
        (meta_array[a][sid] for a in ancestries
         if sid in meta_array.get(a, {})), "")
    anc = ancestry_of.get(sid) or next(
        (a for a in ancestries if sid in meta_rnaseq.get(a, set())), "NA")
    fq1, fq2 = fastq.get(sid, (None, None))
    row = {
        "rnaseq_id": sid, "cohort": coh, "array_id": arr, "ancestry": anc,
        "fastq_r1_present": "" if fq1 is None else int(fq1),
        "fastq_r2_present": "" if fq2 is None else int(fq2),
        "has_quant_sf": int(salmon.get(sid, False)),
        "has_qu_quant_sf": int(qu.get(sid, False)),
    }
    for m in EXPR_MODS:
        row[f"in_{m}_bed"] = int(sid in bed_union[m])
    for a in ancestries:
        row[f"in_metadata_{a}"] = int(sid in meta_rnaseq.get(a, set()))
        row[f"in_covariates_{a}"] = int(
            (arr or sid) in cov_cols.get(a, set()))
    for a in ancestries:
        for m in MODS:
            key = arr or sid
            # INT BEDs are pre-harmonization (rnaseq_id columns);
            # accept either ID space
            ccols = int_cols.get((a, m), set())
            row[f"in_int_{a}_{m}"] = int(sid in ccols or key in ccols)
            row[f"in_qtl_{a}_{m}"] = int(key in qtl_cols.get((a, m), set()))
    matrix.append(row)

samples_path = os.path.join(attrition_dir, "attrition_samples.tsv")
with open(samples_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["rnaseq_id", "cohort", "array_id",
                                       "ancestry"] + stage_cols,
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(matrix)

# ---- stage-level summary (long format) ---------------------------------------
def strata_for(col):
    """(stratum_type, stratum) pairs meaningful for a stage column."""
    out = [("total", "all")]
    m = re.fullmatch(r"in_(?:metadata|covariates|int|qtl)_(\w+?)(?:_(.+))?",
                     col)
    anc_match = next((a for a in ancestries if f"_{a}_" in col
                      or col.endswith(f"_{a}")), None)
    if anc_match:
        out.append(("ancestry", anc_match))
    else:
        for c in sorted({r["cohort"] for r in matrix if r["cohort"]}):
            out.append(("cohort", c))
        for a in ancestries:
            out.append(("ancestry", a))
    return out


# Every stage is counted at two levels: run (one row per RNA-seq run /
# rnaseq_id) and individual (unique array_id; technical replicates share an
# array_id). Run-level counts can exceed individual-level counts without
# any sample gain — both runs of a replicate pair map to one individual.
summary = []
for col in stage_cols:
    vals = [(r[col], r["ancestry"], r["cohort"], r["array_id"] or r["rnaseq_id"])
            for r in matrix]
    for stype, sname in strata_for(col):
        if stype == "total":
            sub = vals
        elif stype == "cohort":
            sub = [v for v in vals if v[2] == sname]
        else:
            sub = [v for v in vals if v[1] == sname]
        summary.append({"stage": col, "level": "run",
                        "stratum_type": stype, "stratum": sname,
                        "n_present": sum(1 for v in sub if v[0] == 1),
                        "n_total": sum(1 for v in sub if v[0] != "")})
        summary.append({"stage": col, "level": "individual",
                        "stratum_type": stype, "stratum": sname,
                        "n_present": len({v[3] for v in sub if v[0] == 1}),
                        "n_total": len({v[3] for v in sub if v[0] != ""})})

summary_path = os.path.join(attrition_dir, "attrition_summary.tsv")
with open(summary_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["stage", "level", "stratum_type",
                                       "stratum", "n_present", "n_total"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary)

print(f"  attrition log: {len(matrix)} samples x {len(stage_cols)} stages")
if int_imputed:
    print(f"  NOTE: INT stage carried forward from the post-outlier "
          f"metadata stage for {len(int_imputed)} ancestry x modality "
          f"cell(s) (INT BEDs not on disk): {', '.join(int_imputed)}")
if ns_notes:
    print(f"  NOTE: stripped stage/cohort namespace prefixes from BED "
          f"columns: {', '.join(ns_notes)}")
for row in summary:
    if row["stratum_type"] == "total":
        print(f"    {row['stage']:<38} {row['level']:<10} "
              f"{row['n_present']}/{row['n_total']}")
PYEOF
then
    MISSING_OTHER+=("attrition log (python section failed — see warnings)")
else
    :
fi

# ---- 3.6 Genotype QC summaries (module 01 outputs) -------------------------
# Module 01 does not place all report inputs in one directory:
#   cohort_imputationQC.R writes
#     ${GENO_BASE}/${COHORT}/genotypes/imputed/qc/report/
#       ${COHORT}_rsq_pass_per_chr.tsv
#   mega_imputationQC.R writes pooled report tables under
#     ${GENO_BASE}/pooled/genotypes/report/
#   and reads the canonical genotype ancestry map from
#     ${GENO_BASE}/pooled/genotypes/pooled_sample_ancestry.tsv
# Resolve those repository-defined locations directly instead of requiring a
# synthetic GENO_QC_DIR. All are non-core: absence only blanks genotype-QC
# report panels and is recorded in MANIFEST.txt.
GENO_BASE="${GENO_BASE:-$(dirname "$OUTPUT_BASE")}"
GENO_COHORTS="${GENO_COHORTS:-GUSTO NIEHS_RICHS SNUH MALI_G3A NIGMS}"
POOLED_GENO_DIR="${POOLED_GENO_DIR:-${GENO_BASE}/pooled/genotypes}"
POOLED_GENO_REPORT="${POOLED_GENO_DIR}/report"

mkdir -p "$STAGING/data/qc/genotype_qc"

echo "  GENO_BASE:          $GENO_BASE"
echo "  GENO_COHORTS:       $GENO_COHORTS"
echo "  POOLED_GENO_DIR:    $POOLED_GENO_DIR"
echo "  POOLED_GENO_REPORT: $POOLED_GENO_REPORT"

# Per-cohort Rsq summaries produced by 01_genotype_imputation/cohort_imputationQC.R.
for COHORT in $GENO_COHORTS; do
    f="${GENO_BASE}/${COHORT}/genotypes/imputed/qc/report/${COHORT}_rsq_pass_per_chr.tsv"
    copy_req "$f" "$STAGING/data/qc/genotype_qc" other \
             "genotype_qc/${COHORT}_rsq_pass_per_chr.tsv"
done

# Pooled summaries produced by 01_genotype_imputation/mega_imputationQC.R.
for name in pooled_sample_counts.tsv pooled_variant_summary.tsv pooled_variants_per_chr.tsv; do
    copy_req "${POOLED_GENO_REPORT}/${name}" "$STAGING/data/qc/genotype_qc" other \
             "genotype_qc/${name}"
done

# Canonical pooled genotype ancestry table consumed by mega_imputationQC.R.
copy_req "${POOLED_GENO_DIR}/pooled_sample_ancestry.tsv" \
         "$STAGING/data/qc/genotype_qc" other \
         "genotype_qc/pooled_sample_ancestry.tsv"

# ---- 4. Manifest -----------------------------------------------------------
{
    echo "placenta_xqtl_report_inputs_${DATE_TAG}"
    echo "generated: $(date -Iseconds) on $(hostname)"
    echo "config:    $CONFIG"
    echo "ancestries: $ANCESTRIES"
    echo
    if [ ${#MISSING_CORE[@]} -eq 0 ] && [ ${#MISSING_OTHER[@]} -eq 0 ]; then
        echo "STATUS: COMPLETE — all expected report inputs present."
    elif [ ${#MISSING_CORE[@]} -eq 0 ]; then
        echo "STATUS: INCOMPLETE — core inputs present; missing non-core files:"
        printf '  %s\n' "${MISSING_OTHER[@]}"
    else
        echo "STATUS: INCOMPLETE — MISSING CORE INPUTS (report will not render):"
        printf '  %s\n' "${MISSING_CORE[@]}"
        echo "missing non-core files:"
        printf '  %s\n' "${MISSING_OTHER[@]:-none}"
    fi
    echo
    echo "contents (file, bytes, lines):"
    (cd "$STAGING" && find . -type f ! -name MANIFEST.txt | sort | while read -r f; do
        printf '%s\t%s\t%s\n' "$f" "$(stat -c %s "$f")" "$(wc -l < "$f")"
    done)
} > "$STAGING/MANIFEST.txt"

if [ ${#MISSING_CORE[@]} -gt 0 ]; then
    echo "ERROR: missing core inputs — see $STAGING/MANIFEST.txt" >&2
    printf '  %s\n' "${MISSING_CORE[@]}" >&2
    exit 1
fi
if [ ${#MISSING_OTHER[@]} -gt 0 ]; then
    echo "WARNING: archive INCOMPLETE — ${#MISSING_OTHER[@]} non-core file(s) missing (see MANIFEST.txt):"
    printf '  %s\n' "${MISSING_OTHER[@]}"
fi

tar czf "$OUT" -C "$(dirname "$STAGING")" "$(basename "$STAGING")"
echo "Wrote: $OUT ($(du -h "$OUT" | cut -f1))"
if [ "$KEEP_STAGING" != "1" ]; then
    rm -rf "$STAGING"
fi
echo "Done."
