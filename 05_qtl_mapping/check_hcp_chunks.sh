#!/bin/bash
# =============================================================================
# check_chunks.sh — verify sharded Picard QC chunk outputs from
#                   19a_picard_sharded.sh
#
# A complete chunk file has all 6 metric groups in its header:
#   AlignMetrics.* (5)  InsertSize.* (3)  RnaMetrics.* (8)
#   GcBias.* (3)        DupMetrics.* (3)  SubjectBias.* (2)
# plus the leading sample column = 25 columns total.
#
# Usage:
#   bash check_chunks.sh                 # uses default QC_DIR below
#   QC_DIR=/path/to/qc_metrics bash check_chunks.sh
#
# Download this file rather than copy-pasting (paste corrupts tab/quote
# escapes, which is what produced the earlier bogus BAD/UNREADABLE flags).
# =============================================================================

QC_DIR="${QC_DIR:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/hcp/qc_metrics}"

if [ ! -d "$QC_DIR" ]; then
    echo "ERROR: directory not found: $QC_DIR" >&2
    exit 1
fi

shopt -s nullglob
files=( "$QC_DIR"/*.chunk*.qcmetrics.tsv )
echo "Found ${#files[@]} chunk files in $QC_DIR"
if [ ${#files[@]} -eq 0 ]; then
    echo "(jobs still running, or wrong directory?)"
    exit 0
fi

# NB: do not name this variable GROUPS — that is a special bash variable
# (process group IDs); assignments to it are silently discarded.
METRIC_GROUPS="AlignMetrics InsertSize RnaMetrics GcBias DupMetrics SubjectBias"
nfull=0
for f in "${files[@]}"; do
    base=$(basename "$f")
    if [ ! -s "$f" ]; then
        echo "EMPTY (job still running or died before writing): $base"
        continue
    fi
    header=$(head -n 1 "$f")
    missing=""
    for g in $METRIC_GROUPS; do
        case "$header" in
            *"$g".*) : ;;                       # group present
            *) missing="$missing $g" ;;         # group absent
        esac
    done
    if [ -n "$missing" ]; then
        echo "INCOMPLETE:${missing}  <- $base"
    else
        nfull=$((nfull+1))
    fi
done
echo "----"
echo "Complete (all 6 metric groups): $nfull / ${#files[@]}"

