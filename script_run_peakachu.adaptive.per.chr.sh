#!/bin/bash

# Input BAMs
TREAT="P13.bam"
CTRL="P13_In.bam"

# Output base name
OUTBASE="P13_unique.peakachu.adaptive.chr."

# Get list of chromosomes from BAM header
chroms=$(samtools idxstats "$TREAT" | cut -f1 | grep -v '*' | grep -v '^$')

for chr in $chroms; do
    echo "Processing $chr..."

    # Create chromosome-specific BAMs
    samtools view -b "$TREAT" "$chr" > "${chr}_treat.bam"
    samtools view -b "$CTRL" "$chr" > "${chr}_ctrl.bam"

    # Index them
    samtools index "${chr}_treat.bam"
    samtools index "${chr}_ctrl.bam"

    f=1.5 # fold change
    Q=0.1 # q value

  peakachu adaptive \
    -t "${chr}_treat.bam" \
    -c "${chr}_ctrl.bam" \
    -o "${OUTBASE}_$chr.unique.peakachu.adaptive" \
    -p 1 \
    -M 200 \
    --norm_method deseq \
    --gff_folder gencode.v48.basic.annotation/ \
    --features "transcript,exon,five_prime_UTR,three_prime_UTR,start_codon,stop_codon" \
    -f ${f} \
    -Q ${Q}

    # Optional: clean up
    rm "${chr}_treat.bam" "${chr}_treat.bam.bai"
    rm "${chr}_ctrl.bam" "${chr}_ctrl.bam.bai"

    echo "Done with $chr."
done