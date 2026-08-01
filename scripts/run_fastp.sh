#!/bin/bash

mkdir -p results/trimmed

for r1 in data/*_r1.fastq.gz; do
    r2=${r1/_r1.fastq.gz/_r2.fastq.gz}
    sample=$(basename "$r1" _r1.fastq.gz)

    echo "Processing sample: $sample"

    fastp \
      -i "$r1" \
      -I "$r2" \
      -o "results/trimmed/${sample}_r1.trimmed.fastq.gz" \
      -O "results/trimmed/${sample}_r2.trimmed.fastq.gz" \
      -h "results/trimmed/${sample}_fastp.html" \
      -j "results/trimmed/${sample}_fastp.json"
done

echo "All samples processed."
