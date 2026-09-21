COHORT="SNUH"
OUTDIR=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/${COHORT}/genotypes/imputed/qc
REPORT_DIR="${OUTDIR}/report"

ancestries=$(awk 'NR>1 {print $2}' "${OUTDIR}/${COHORT}_ancestry_assignments.tsv" | sort -u)

for anc in $ancestries; do
    rm -f "${OUTDIR}/${COHORT}_${anc}_samples.txt"
    rm -f "${OUTDIR}/${COHORT}_${anc}_passing.snplist"
    rm -f "${OUTDIR}/${COHORT}_${anc}_passing.log"
    rm -f "${OUTDIR}/${COHORT}_passing_variants_${anc}_Rsq"*".tsv"
done

rm -f "${OUTDIR}/${COHORT}_cis_window_variant_counts.tsv"
rm -f "${OUTDIR}/${COHORT}_QC_summary_per_chr.tsv"
rm -f "${OUTDIR}/${COHORT}_passing_variants_union_Rsq"*".tsv"
rm -f "${OUTDIR}/${COHORT}_passing_variants_intersect_Rsq"*".tsv"
rm -f "${OUTDIR}/${COHORT}_passing_variants_ancestry_annot.tsv"

rm -f "${REPORT_DIR}/${COHORT}_snps_per_mb_per_chr.png"
rm -f "${REPORT_DIR}/${COHORT}_cis_window_hist.png"
rm -f "${REPORT_DIR}/${COHORT}_cis_window_per_chr.png"
rm -f "${REPORT_DIR}/${COHORT}_cis_low_coverage_summary.tsv"
rm -f "${REPORT_DIR}/${COHORT}_cis_median_summary.tsv"
rm -f "${REPORT_DIR}/${COHORT}_cis_window_variant_counts.tsv"
rm -f "${REPORT_DIR}/${COHORT}_QC_summary_per_chr.tsv"