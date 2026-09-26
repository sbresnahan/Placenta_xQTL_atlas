#!/bin/bash
# make_report_archive.sh — assemble every file report_placenta_xqtl.Rmd needs
# into a single tarball, in the exact relative layout the Rmd expects
# (data/results, data/qc/..., gene_map.tsv at the root).
#
# Run on seadragon AFTER modules 28 (all three layers) and 29 have completed
# for every ancestry:
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
         "$STAGING/data/qc/hcp_optimization" "$STAGING/data/qc/genotype_pcs"

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
# Grouped (8 modalities), ungrouped (7), combined, and independent top TSVs
# per ancestry — glob covers every layer 29 writes.
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

# Expected-completeness audit (non-core): grouped 8 + ungrouped 7 +
# combined 1 + independent 9 = 25 top tables per ancestry.
# modality file stems as written by 27/28/29 (display-case: alt_TSS, alt_polyA, RNA_editing)
MODS_GROUPED="expression isoforms splicing intron_retention alt_TSS alt_polyA RNA_editing stability"
MODS_UNGROUPED="isoforms splicing intron_retention alt_TSS alt_polyA RNA_editing stability"
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
    done
    copy_req "$RESULTS_DIR/${ANC}_combined_cisqtl_top.tsv" \
             "$STAGING/data/results" other "combined:${ANC}"
    copy_req "$RESULTS_DIR/${ANC}_combined_cisqtl_independent_top.tsv" \
             "$STAGING/data/results" other "independent:${ANC}_combined"
done

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
    copy_req "$QTL_DIR/hcp_optimization/${ANC}_optimal_hcp.tsv" \
             "$STAGING/data/qc/hcp_optimization" other "${ANC}_optimal_hcp.tsv"
    copy_req "$QTL_DIR/hcp_optimization/${ANC}_optimal_hcp.png" \
             "$STAGING/data/qc/hcp_optimization" other "${ANC}_optimal_hcp.png"
    copy_req "$PC_DIR/${ANC}_genotype_pcs_scree.png" \
             "$STAGING/data/qc/genotype_pcs" other "${ANC}_genotype_pcs_scree.png"
done

# ---- 3. Gene map + gene bodies from the normalized GTF --------------------
# Regenerated in full each round (replaces the ad-hoc round-3 map): every
# gene in the GTF, version-stripped IDs matching the pipeline's BED files.
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
