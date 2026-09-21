#!/bin/sh

DIR=/rsrch5/home/epi/bhattacharya_lab/data/mapqtl/NIEHS_RICHS

for FILE in ${DIR}/*.fastq.gz; do
  LIBID=$(basename ${FILE} .fastq.gz)
  bsub -env LIBID=${LIBID} < run_salmon_NIEHS_RICHS.lsf
done