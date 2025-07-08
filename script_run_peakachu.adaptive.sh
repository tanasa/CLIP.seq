#!/bin/bash

# PEAKachu peak caller parameters :

f=1.5    # -f FC_CUTOFF, --fc_cutoff FC_CUTOFF
Q=1      # -Q PADJ_THRESHOLD, --padj_threshold PADJ_THRESHOLD

peakachu adaptive \
-t P13.bam \
-c P13_In.bam \
-o "P13_unique.peakachu.adaptive" \
-p 2 \
-M 200 \
--norm_method deseq \
--gff_folder gencode.v48.basic.annotation/ \
--features "transcript,exon,five_prime_UTR,three_prime_UTR,start_codon,stop_codon" \
-f ${f} \
-Q ${Q}
# -m MAD_MULTIPLIER, --mad_multiplief=2

# --features "transcript,exon,five_prime_UTR,three_prime_UTR,start_codon,stop_codon" 
# cut -f3 gencode.v48.basic.annotation.gff3 | grep -v "^#" | sort | uniq
# CDS
# exon
# five_prime_UTR
# gene
# start_codon
# stop_codon
# stop_codon_redefined_as_selenocysteine
# three_prime_UTR
# transcript
# --features "exon,CDS,UTR,transcript,gene,start_codon,stop_codon" 