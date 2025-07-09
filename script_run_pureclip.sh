#!/bin/bash

pureclip \
-i A27_13.Aligned.sortedByCoord.out.bam \
-bai A27_13.Aligned.sortedByCoord.out.bam.bai \
-ibam A27_13_In.Aligned.sortedByCoord.out.bam \
-ibai A27_13_In.Aligned.sortedByCoord.out.bam.bai \
-g hg38.ucsc.fa \
-nt 25 \
-bw 50 \
-o A27_13_PureCLIP.crosslink_sites.cov_inputSignal.bed \
-or A27_13_PureCLIP.crosslink_sites.cov_inputSignal.all.sites.bed \
-p A27_13_PureCLIP.crosslink_sites.cov_inputSignal.learnt.param.bed \
-iv 'chr1;chr2;chr3;'


# pureclip \
#            -i {input.bam_file} \
#            -bai {input.bam_index} \
#            -g {params.genome_fasta} \
#            -o {output.crosslink_site_file} \
#            -or {output.binding_region_file} \
#            -nt {threads}

#    pureclip [OPTIONS] <-i BAM FILE> <-bai BAI FILE> <-g GENOME FILE> <-o OUTPUT BED FILE>

# DESCRIPTION

#    Protein-RNA interaction site detection using a non-homogeneous HMM.
#    -h, --help
#          Display the help message.
#    --version
#          Display version information.
#    -i, --in BAM
#          Target bam files. Valid filetype is: .bam.
#    -bai, --bai BAI
#          Target bam index files. Valid filetype is: .bai.
#    -g, --genome FILE
#          Genome reference file. Valid filetypes are: .fa, .fasta, .fa.gz, and .fasta.gz.
#    -o, --out FILE 
           Output file to write crosslink sites. Valid filetype is: .bed.
#    -or, --or FILE
#          Output file to write binding regions. Valid filetype is: .bed.
#    -p, --par FILE
#          Output file to write learned parameters.

#  Options:

#    -ctr, --ctr
#          Assign crosslink sites to read start positions. Note: depends on RT enzyme, buffer conditions and likel>
#          protein. Default: assign crosslink sites to positions upstream of read starts.
#    -st, --st NUM
#          Scoring scheme. Default: 0 -> score_UC (log posterior probability ratio of most likely and second most
#          likely state). In range [0..3].
#    -iv, --inter STR
#          Genomic chromosomes to learn HMM parameters, e.g. 'chr1;chr2;chr3'. Contigs have to be in the same orde>
#          in BAM file. Useful to reduce runtime and memory consumption. Default: all contigs from reference file >
#          used (useful when applying to transcript-wise alignments or poor data).
#    -chr, --chr STR
#          Contigs to apply HMM, e.g. 'chr1;chr2;chr3;'. Contigs have to be in the same order as in BAM file.
#    -bc, --bc NUM
#          Flag to set parameters according to binding characteristics of protein: see description in section belo>
#          range [0..1].
#    -bw, --bdw NUM
#          Bandwidth for kernel density estimation used to access enrichment. NOTE: Increasing the bandwidth incre>
#          runtime and memory consumption. Default: 50. In range [1..500].
#    -bwn, --bdwn NUM
#          Bandwidth for kernel density estimation used to estimate n for binomial distributions. For proteins tha>
#          rather sliding along the RNA or showing long crosslink clusters this should be increased, e.g. to 100
#          (should be <= 4*bdw). Default: same as bdw. In range [1..500].
#    -dm, --dm NUM
#          Distance used to merge individual crosslink sites to binding regions. Default: 8
#    -ld, --ld
#          Use higher precision to store emission probabilities, state poster posterior probabilities etc. (i.e. l>
#          double). Should not be necessary anymore, due to computations in log-space. Note: increases memory
#          consumption. Default: double.
#    -ts, --ts NUM
#          Size of look-up table for log-sum-exp values. Default: 600000
#    -tmv, --tmv NUM
#          Minimum value in look-up table for log-sum-exp values. Default: -2000
#    -ur, --ur NUM
#          Flag to define which read should be selected for the analysis: 1->R1, 2->R2. Note: PureCLIP uses read s>
#          corresponding to 3' cDNA ends. Thus if providing paired-end data, only the corresponding read should be
#          selected (e.g. eCLIP->R2, iCLIP->R1). If applicable, used for input BAM file as well. Default: uses read
#          starts of all provided reads assuming single-end or pre-filtered data. In range [1..2].


#  Options for incorporating covariates:

#    -is, --is FILE
#          Covariates file: position-wise values, e.g. smoothed reads start counts (KDEs) from input data. Valid
#          filetype is: .bed.
#    -ibam, --ibam FILE
#          File containing mapped reads from control experiment, e.g. eCLIP input. Valid filetype is: .bam.
#    -ibai, --ibai FILE
#          File containing BAM index corresponding to mapped reads from control experiment Valid filetype is: .bai.
#    -fis, --fis FILE
#          Fimo input motif score covariates file. Valid filetype is: .bed.
#    -nim, --nim NUM
#          Max. motif ID to use. Default: Only covariates with motif ID 1 are used.

#  Advanced user options:

#    -upe, --upe
#          Use (n dependent) pseudo emission probabilities for crosslink state.
#    -m, --mibr NUM
#          Maximum number of iterations within BRENT algorithm. In range [1..1000].
#    -w, --mibw NUM
#          Maximum number of iterations within Baum-Welch algorithm. In range [0..500].
#    -g1kmin, --g1kmin NUM
#          Minimum shape k of 'non-enriched' gamma distribution (g1.k). In range [1.5..inf].
#    -g1kmax, --g1kmax NUM
#          Maximum shape k of 'non-enriched' gamma distribution (g1.k).
#    -g2kmin, --g2kmin NUM
#          Minimum shape k of 'enriched' gamma distribution (g2.k).
#    -g2kmax, --g2kmax NUM
#          Maximum shape k of 'enriched' gamma distribution (g2.k).
#    -fk, --fk
#          When incorporating input signal, do not constrain 'non-enriched' shape parameter k <= 'enriched' gamma
#          parameter k.
#    -mkn, --mkn NUM
#          Max. k/N ratio (read start sites/N) used to learn truncation probabilities for 'non-crosslink' and
#          'crosslink' emission probabilities (high ratios might originate from mapping artifacts that can disturb
#          parameter learning). Default: 1.0 In range [0.5..1.5].
#    -b1p, --b1p NUM
#          Initial value for binomial probability parameter of 'non-crosslink' state. Default: 0.01.
#    -b2p, --b2p NUM
#          Initial value for binomial probability parameter of 'crosslink' state. Default: 0.15.
#    -mtp, --mtp NUM
#          Min. transition probability from state '2' to '3' (helpful for poor data, where no clear distinction be>
#          'enriched' and 'non-enriched' is possible). Default: 0.0001.
#    -mk, --mkde NUM
#          Minimum KDE value used for fitting left-truncated gamma distributions. Default: corresponding to single>
#          read start.
#    -ntp, --ntp NUM
#          Only sites with n >= ntp are used to learn binomial probability parameters (bin1.p, bin2.p). Default: 10
#    -ntp2, --ntp2 NUM
#          Only sites with n >= ntp2 are used to learn probability of transition from state '2' to '2' or '3'. Use>
#          for data with low truncation rate at crosslink sites or in general high fraction of non-coinciding read
#          starts. Default: 0
#    -antp, --antp
#          Automatically choose n threshold (-ntp, -ntp2) to estimate parameters linked to crosslink states based >
#          expected read start count at crosslink sites.
#    -pa, --pat NUM
#          Length threshold for internal poly-X stretches to get excluded.
#    -ea1, --epal
#          Exclude intervals containing poly-A stretches from learning.
#    -ea2, --epaa
#          Exclude intervals containing poly-A stretches from analysis.
#    -et1, --eptl
#          Exclude intervals containing poly-U stretches from learning.
#    -et2, --epta
#          Exclude intervals containing poly-U stretches from analysis.
#    -mrtf, --mrtf NUM
#          Fit gamma shape k only for positions with min. covariate value.
#    -mtc, --mtc NUM
#          Maximum number of read starts at one position used for learning. For sites with counts above threshold >
#          whole covered regions will be ignored for learning! Default: 500. In range [50..50000].
#    -mtc2, --mtc2 NUM
#          Maximum number of read starts at one position stored. For sites with counts above threshold the count w>
#          be truncated. Influences k and n. Default: 65000. In range [5000..65000].
#    -pet, --pet NUM
#          Prior enrichment threshold: a KDE threshold corresponding to 7 read start counts at one position will be
#          used for initial classification of 'non-enriched' and 'enriched' site. Default: 7 In range [2..50].

#  General user options:
#    -nt, --nt NUM
#          Number of threads used for learning.
#    -nta, --nta NUM
#          Number of threads used for applying learned parameters. Increases memory usage, if greater than number >
#          chromosomes used for learning, since HMM will be build for multiple chromosomes in parallel. Default:
#          min(nt, no. of chromosomes/transcripts used for learning).
#    -oa, --oa
#          Outputs all sites with at least one read start in extended output format.
#    -q, --quiet
#          Set verbosity to a minimum.
#    -v, --verbose
#          Enable verbose output.
#    -vv, --very-verbose
#          Enable very verbose output.

# PARAMETER SETTINGS FOR PROTEINS WITH DIFFERENT BINDING CHARACTERISTICS

#    By default, the parameters are set to values optimized for proteins binding to short defined binding regions,>
#    proteins binding to short specific motifs such as PUM2 and RBFOX2. With the -bc option this behaviour can be
#    changed:

#    0
#          Short defined. Default. Equivalent to: -bdwn 50 -ntp 10 -ntp2 0 -b1p 0.01 -b2p 0.15.
#    1
#          Larger clusters. Proteins causing larger crosslink clusters with relatively lower read start counts, e.>
#          proteins binding to low complexity motifs. Equivalent to: -bdwn 100 -antp -b2p 0.01 -b2p 0.1.

#    In case of different binding characteristics adjust parameters -bdw, -bdwn, -b1p, -b2p, -antp or see
#    http://pureclip.readthedocs.io/en/latest/PureCLIPTutorial/userOptions.html for more information.

# EXAMPLES

#    pureclip -i target.bam -bai target.bai -g ref.fasta -o called_crosslinksites.bed -nt 10  -iv '1;2;3;'
#          Learn HMM parameters on chromosomes 1-3, use 10 threads for learning and otherwise default parameters.
#    pureclip -i target.rep1.bam -bai target.rep1.bai -i target.rep2.bam -bai target.rep2.bai -g ref.fasta -o call>
#          Include individual replicates (currently only supported for two), while learning parameters on whole
#          datasets.
#    pureclip -i target.bam -bai target.bai -g ref.fasta -o called_crosslinksites.bed -nt 10  -iv '1;2;3;' -bc 1 
#          Use parameter settings for proteins causing larger crosslink clusters
#    pureclip -i target.bam -bai target.bai -g ref.fasta -o called_crosslinksites.bed -nt 10  -iv '1;2;3;' -bc 1 ->
#          Use parameter settings for proteins causing larger crosslink clusters and decrease initial probability
#          parameter for 'crosslink' state for data with high fraction of non-coinciding read starts.
#    pureclip -i target.bam -bai target.bai -g ref.fasta -o called_crosslinksites.bed -nt 10  -iv '1;2;3;' -bdw 25 
#          Use decreased bandwidth of 25 bp to access enrichment.
