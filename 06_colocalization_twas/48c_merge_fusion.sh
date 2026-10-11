#!/bin/bash
# =============================================================================
# 48c_merge_fusion.sh — merge completed FUSION feature-shard outputs
# =============================================================================
# Submitted after the entire 48a global shard array ends. A shard is complete
# when its process produced an output and .done sentinel; individual phenotype
# failures inside that shard are tolerated and recorded by 48b_fusion_assoc.R.
#
# The merge publishes:
#   <weight_set>_<trait>.twas.tsv      successful phenotypes only, with the
#                                      original FUSION schema
#   <weight_set>_<trait>.failures.tsv failed/skipped phenotypes with the error
#                                      recorded by the resilient FUSION wrapper
#
# Legacy shard outputs produced before FUSION.STATUS/FUSION.ERROR existed are
# restart-compatible: rows with a nonmissing TWAS.P are retained, while legacy
# NA/blank TWAS.P rows are moved to the failure audit with a generic reason.
#
# Env: ISOTWAS_DIR, WEIGHT_SET, TRAIT_ID, SHARD_MANIFEST
# =============================================================================
set -eo pipefail

SHARD_MANIFEST="${SHARD_MANIFEST:?ERROR: SHARD_MANIFEST not set}"
if [ ! -s "$SHARD_MANIFEST" ]; then
    echo "ERROR: FUSION shard manifest missing/empty: $SHARD_MANIFEST" >&2
    exit 1
fi

PLAN_ID="$(basename "$(dirname "$SHARD_MANIFEST")")"
PER_SHARD_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}/per_shard/${TRAIT_ID}/${PLAN_ID}"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}"
FINAL="${OUT_DIR}/${WEIGHT_SET}_${TRAIT_ID}.twas.tsv"
FAILURES="${OUT_DIR}/${WEIGHT_SET}_${TRAIT_ID}.failures.tsv"
mkdir -p "$OUT_DIR"

TMP_FINAL="${FINAL}.tmp.${LSB_JOBID:-$$}"
TMP_FAILURES="${FAILURES}.tmp.${LSB_JOBID:-$$}"
trap 'rm -f "$TMP_FINAL" "$TMP_FAILURES"' EXIT
rm -f "$TMP_FINAL" "$TMP_FAILURES"

# Fail only for shard-level/infrastructure incompleteness. Phenotype-level
# failures are represented inside completed shard outputs and are filtered
# below instead of aborting the merge.
N_EXPECTED=0
while IFS=$'\t' read -r ARRAY_INDEX CHR SHARD_INDEX N_FEATURES SHARD_POS; do
    [ "$ARRAY_INDEX" = "array_index" ] && continue
    N_EXPECTED=$((N_EXPECTED + 1))
    printf -v SHARD_STEM 'chr%02d.shard%04d' "$CHR" "$SHARD_INDEX"
    SHARD_OUT="${PER_SHARD_DIR}/${SHARD_STEM}.twas.tsv"
    DONE="${SHARD_OUT}.done"
    if [ ! -f "$DONE" ] || [ ! -s "$SHARD_OUT" ]; then
        echo "ERROR: ${SHARD_STEM} is incomplete: expected $SHARD_OUT and $DONE" >&2
        exit 1
    fi
done < "$SHARD_MANIFEST"

if [ "$N_EXPECTED" -eq 0 ]; then
    echo "ERROR: no feature shards in $SHARD_MANIFEST" >&2
    exit 1
fi

# Normalize and merge with Python's stdlib CSV reader so old shard files
# (without FUSION.STATUS/FUSION.ERROR) can coexist with newly resilient shards.
python3 - "$SHARD_MANIFEST" "$PER_SHARD_DIR" "$WEIGHT_SET" "$TRAIT_ID" "$TMP_FINAL" "$TMP_FAILURES" <<'PY'
import csv
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
per_shard_dir = Path(sys.argv[2])
weight_set = sys.argv[3]
trait_id = sys.argv[4]
final_path = Path(sys.argv[5])
fail_path = Path(sys.argv[6])

status_col = "FUSION.STATUS"
error_col = "FUSION.ERROR"
base_fields = None
success_rows = []
failure_rows = []

with manifest_path.open(newline="", encoding="utf-8") as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))

for mrow in manifest:
    chrom = int(mrow["chromosome"])
    shard = int(mrow["shard_index"])
    stem = f"chr{chrom:02d}.shard{shard:04d}"
    path = per_shard_dir / f"{stem}.twas.tsv"

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise SystemExit(f"ERROR: shard output has no header: {path}")

        this_base = [x for x in fields if x not in (status_col, error_col)]
        if base_fields is None:
            base_fields = this_base
        elif this_base != base_fields:
            raise SystemExit(
                "ERROR: incompatible FUSION columns across shards: "
                f"{path}\nexpected={base_fields}\nobserved={this_base}"
            )

        has_status = status_col in fields
        for row in reader:
            if has_status:
                status = (row.get(status_col) or "failed").strip().lower()
                error = (row.get(error_col) or "").strip()
            else:
                # Pre-patch shards cannot identify runtime exceptions, but an
                # NA/blank TWAS.P is FUSION's legacy representation of a model
                # that was skipped/failed. Preserve valid legacy rows while
                # moving those failed models into the failure audit table.
                twas_p = (row.get("TWAS.P") or "").strip().lower()
                status = "failed" if twas_p in {"", "na", "nan"} else "ok"
                error = (
                    "Legacy FUSION output has missing TWAS.P; detailed failure "
                    "reason was not recorded" if status != "ok" else ""
                )

            if status == "ok":
                success_rows.append({k: row.get(k, "") for k in base_fields})
                continue

            failure_rows.append({
                "weight_set": weight_set,
                "trait_id": trait_id,
                "shard": stem,
                "ID": row.get("ID", ""),
                "FILE": row.get("FILE", ""),
                "CHR": row.get("CHR", ""),
                "P0": row.get("P0", ""),
                "P1": row.get("P1", ""),
                "FUSION.STATUS": status or "failed",
                "FUSION.ERROR": error or "FUSION phenotype failed/skipped",
            })

if base_fields is None:
    raise SystemExit("ERROR: no FUSION shard headers found")

with final_path.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=base_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(success_rows)

failure_fields = [
    "weight_set", "trait_id", "shard", "ID", "FILE", "CHR", "P0", "P1",
    "FUSION.STATUS", "FUSION.ERROR",
]
with fail_path.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=failure_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(failure_rows)

print(f"merge summary: {len(success_rows)} successful phenotypes, {len(failure_rows)} failed/skipped phenotypes")
PY

if [ ! -s "$TMP_FINAL" ] || [ ! -s "$TMP_FAILURES" ]; then
    echo "ERROR: merge did not produce expected temporary outputs" >&2
    exit 1
fi

mv -f "$TMP_FINAL" "$FINAL"
mv -f "$TMP_FAILURES" "$FAILURES"
trap - EXIT

N_SUCCESS=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$FINAL")
N_FAILED=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$FAILURES")
echo "wrote $FINAL (${N_SUCCESS} successful phenotype rows from ${N_EXPECTED} chromosome x feature shards)"
echo "wrote $FAILURES (${N_FAILED} failed/skipped phenotype rows)"
