#!/bin/bash

mkdir -p results/aligned

for r1 in results/trimmed/*_r1.trimmed.fastq.gz; do
    r2=${r1/_r1.trimmed.fastq.gz/_r2.trimmed.fastq.gz}
    sample=$(basename "$r1" _r1.trimmed.fastq.gz)

    echo "Aligning sample: $sample"

    STAR \
      --runMode alignReads \
      --genomeDir data/reference/star_index \
      --readFilesIn "$r1" "$r2" \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "results/aligned/${sample}_" \
      --outTmpDir ~/star_tmp \
      --runThreadN 2

done

echo "All samples aligned."
