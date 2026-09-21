#!/bin/sh

DIR=/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/SNUH/RNA
for FILE in ${DIR}/*_1.fastq; do
  LIBID=$(basename ${FILE} _1.fastq)
  bsub -env LIBID=${LIBID} < run_salmon_SNUH.lsf
done