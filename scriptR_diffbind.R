############################################################ 13P vs 13P + Cisplatin
############################################################ 16 hr CLIP

library(DiffBind)
library(rtracklayer)

############################################################
############################################################

# load sample sheet
samples <- read.csv("samples_PLAT-seq.csv", stringsAsFactors = FALSE)

# construct the DBA object
dba_obj <- dba(sampleSheet = samples)

############################################################
############################################################ raw counts 

dba_obj  <- dba.count(
  dba_obj,
  summits = 0,    
  score = DBA_SCORE_READS,   # ensures integer counts      
  bParallel = TRUE      
)

counts_raw <- dba.peakset(dba_obj, bRetrieve = TRUE)
head(counts_raw)

count_data_raw <- as.data.frame(counts_raw)
rownames(count_data_raw) <- paste0(count_data_raw$seqnames, ":", count_data_raw$start, "-", count_data_raw$end)
head(count_data_raw)

#                   seqnames  start    end width strand X16h_A27_13_1
# chr1:16732-17135       chr1  16732  17135   404      *             5
# chr1:189974-190113     chr1 189974 190113   140      *            10
# chr1:630844-631242     chr1 630844 631242   399      *            23
# chr1:923950-925283     chr1 923950 925283  1334      *           121
# chr1:956832-957420     chr1 956832 957420   589      *          1484
# chr1:960945-961761     chr1 960945 961761   817      *            40
#                   X16h_A27_13_2 X16h_A27_Cis13_1 X16h_A27_Cis13_2
# chr1:16732-17135               7               82               67
# chr1:189974-190113             6               59               64
# chr1:630844-631242            23              482              158
# chr1:923950-925283           116              361              416
# chr1:956832-957420          1552             7250             5458
# chr1:960945-961761            23              278              209

write.csv(count_data_raw, file = "16h_DiffBind_raw_counts.csv", row.names = TRUE)

# The counts stored in dba_obj at this point are indeed raw counts — 
# i.e., unnormalized read/fragment counts overlapping each peak region for each sample.

# These counts are:
# Directly computed from the BAM files and the consensus peak set
# Not normalized, not scaled, and not log-transformed
# Suitable for exporting or feeding into downstream tools like DESeq2 or edgeR.

x_cols_raw <- grep("^X", names(count_data_raw), value = TRUE)

boxplot(count_data_raw[, x_cols_raw],
        las = 2,
        main = "Boxplot of Raw Read Counts per Sample",
        ylab = "Read Counts",
        col = "lightgray",
        outline = FALSE )


png("16h_boxplot_raw_counts.png", width = 800, height = 600)
boxplot(count_data_raw[, x_cols_raw],
        las = 2,
        main = "Boxplot of Raw Read Counts per Sample",
        ylab = "Read Counts",
        col = "lightgray",
        outline = FALSE)
dev.off()

boxplot(log2(count_data_raw[, x_cols_raw] + 1),
        las = 2,
        main = "Log2-Transformed Read Counts per Sample",
        ylab = "Log2(Read Counts + 1)",
        col = "skyblue",
        outline = FALSE)


png("16h_boxplot_log2_counts.png", width = 800, height = 600)
boxplot(log2(count_data_raw[, x_cols_raw] + 1),
        las = 2,
        main = "Log2-Transformed Read Counts per Sample",
        ylab = "Log2(Read Counts + 1)",
        col = "skyblue",
        outline = FALSE)
dev.off()

############################################################
############################################################ official documentation

# which score to use in the binding affinity matrix : options in dba.count()

#            DBA_SCORE_NORMALIZED     normalized reads, as set by ‘dba.normalize’                                           

#            DBA_SCORE_READS          raw read count for interval using only reads from ChIP                                
#            DBA_SCORE_CONTROL_READS  raw read count for interval using only reads from Control                             

#            DBA_SCORE_READS_FOLD     raw read count for interval from ChIP divided by read count for interval from control 
#            DBA_SCORE_READS_MINUS    raw read count for interval from ChIP minus read count for interval from control      

#            DBA_SCORE_RPKM           RPKM for interval using only reads from ChIP                                          
#            DBA_SCORE_RPKM_FOLD      RPKM for interval from ChIP divided by RPKM for interval from control                 
#            DBA_SCORE_RPKM_MINUS     RPKM for interval from ChIP minus RPKM for interval from control                      
#            DBA_SCORE_SUMMIT         summit height (maximum read pileup value)                                             
#            DBA_SCORE_SUMMIT_ADJ     summit height (maximum read pileup value), normalized to relative library size        
#            DBA_SCORE_SUMMIT_POS     summit position (location of maximum read pileup)     ''

############################################################
############################################################

# save a copy :  
# dba_obj_cp = dba_obj 

dba.plotHeatmap(dba_obj, correlations = FALSE, scale = "row")
dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)

png("16h_heatmap_dba_obj_raw.png", width = 1000, height = 800)
dba.plotHeatmap(dba_obj, correlations = FALSE, scale = "row")
dev.off()

png("16h_pca_dba_obj_raw.png", width = 1000, height = 800)
dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
dev.off()

############################################################
############################################################ dba.normalize()

# dba_obj <- dba.normalize(dba_obj)

# Enables normalization of datasets using a variety of methods,
# including background, spike-in, and parallel factor normalization.
# Alternatively, it allows a user to specify library sizes and
# normalization factors directly, or retrieve computed ones.

# method: Underlying method, or vector of methods, for which to normalize.
#          Supported methods:
#  • ‘DBA_EDGER’ use ‘edgeR’ package for analysis
#  • ‘DBA_DESEQ2’ use ‘DESeq2’ package for analysis
#  • ‘DBA_ALL_METHODS’ normalize for both both ‘edgeR’ and ‘DESeq2’

# normalize: Either user-supplied normalization factors in a numeric
#          vector, or a specification of a method to use to calculate
#          normalization factors.
#
#          Methods can be specified using one of the following:
#
#            • ‘DBA_NORM_RLE’ ("RLE") RLE normalization (native to
#              ‘DBA_DESEQ2’, and available for ‘DBA_EDGER’).
#
#            • ‘DBA_NORM_TMM’ ("TMM") TMM normalization (native to
#              ‘DBA_EDGER’, and available for ‘DBA_DESEQ2’).
#
#            • ‘DBA_NORM_NATIVE’ ("native") Use native method based on
#              ‘method’: ‘DBA_NORM_RLE’ for ‘DBA_DESEQ2’ or
#              ‘DBA_NORM_TMM’ for ‘DBA_EDGER’.
#
#            • ‘DBA_NORM_LIB’ ("lib") Normalize by library size only.
#              Library sizes can be specified using the ‘library’
#              parameter.  Normalization factors will be calculated to
#              give each equal weight in a manner appropriate for the
#              analysis ‘method’.  See also the ‘libFun’ parameter,
#              which can be used to scale the normalization factors for
#              ‘DESeq2.’

# library: Either user-supplied library sizes in a numeric vector, or a
#          specification of a method to use to calculate library sizes.
#
#          Library sizes can be based on one of the following:
#
#            • ‘DBA_LIBSIZE_FULL’ ("full") Use the full library size
#              (total number of reads in BAM/SAM/BED file)
#
#            • ‘DBA_LIBSIZE_PEAKREADS’ ("RiP") Use the number of reads
#              that overlap consensus peaks.


# The default normalization parameters are as follows:
#
#        • ‘normalize=DBA_NORM_LIB’
#
#        • ‘library=DBA_LIBSIZE_FULL’
#
#        • ‘background=FALSE’

# Basic Normalization Methods

# 1. Default normalization (DESeq2 with RLE)
# norm_default <- dba.normalize(dba_obj, bRetrieve=TRUE) 

# 2. Library size normalization
# norm_lib <- dba.normalize(dba_obj, 
#                         normalize=DBA_NORM_LIB, 
#                         bRetrieve=TRUE)

# 3. TMM normalization (edgeR method)
# norm_tmm <- dba.normalize(dba_obj, 
#                         normalize=DBA_NORM_TMM, 
#                         bRetrieve=TRUE)

# 4. RLE normalization (DESeq2 method)
# norm_rle <- dba.normalize(dba_obj, 
#                         normalize=DBA_NORM_RLE, 
#                         bRetrieve=TRUE)

############################################################
############################################################

# In DiffBind, method and normalize serve different purposes:

# method Parameter
# Controls the statistical framework/algorithm used for analysis

# rmethod=DBA_DESEQ2    # Uses DESeq2 statistical methods
# method=DBA_EDGER     # Uses edgeR statistical methods  
# method=DBA_ALL_METHODS # Uses both DESeq2 and edgeR

# 👉 This affects how differential binding is computed — for example, how dispersion is estimated, 
# hypothesis testing is done, and fold changes are calculated.

# normalize Parameter
# Controls how the raw counts are normalized before statistical analysis

# rnormalize=DBA_NORM_RLE   # Relative Log Expression (DESeq2's method)
# normalize=DBA_NORM_TMM    # Trimmed Mean of M-values (edgeR's method)
# normalize=DBA_NORM_LIB    # Simple library size normalization
# normalize=DBA_NORM_NONE   # No normalization

# 👉 This affects the input matrix of normalized counts that goes into the modeling. 
# It’s especially important when comparing samples with very different sequencing depths or biases.

# Option 1: Default DESeq2

norm1 <- dba.normalize(dba_obj, 
                      method=DBA_DESEQ2,
                      normalize=DBA_NORM_RLE)

# Option 2: edgeR with TMM
norm2 <- dba.normalize(dba_obj, 
                      method=DBA_EDGER,
                      normalize=DBA_NORM_TMM)

# Option 3: Library size only
norm3 <- dba.normalize(dba_obj, 
                      method=DBA_DESEQ2,
                      normalize=DBA_NORM_LIB)

# Here are several ways to view the normalization factors after running dba.normalize():

# View RLE normalization factors (DESeq2)
norm1$norm$DESeq2$norm.facs

# View TMM normalization factors (edgeR)
norm2$norm$edgeR$norm.facs

# View library size factors
norm3$norm$DESeq2$norm.facs

norm1$norm$DESeq2$norm.facs
#   16h_A27_13_1    16h_A27_13_2 16h_A27_Cis13_1 16h_A27_Cis13_2 
#      0.4717863       0.4855319       2.3659001       2.0088122 
norm2$norm$edgeR$norm.facs
# [1] 0.9586388 0.9314120 1.0907274 1.0268025
norm3$norm$DESeq2$norm.facs
# [1] 0.3204215 0.3490627 1.7557752 1.5747406

# 🔍 Why the Differences?

# 1. Library Size Normalization (norm3)
# Simple scaling by total number of reads per sample.
# Can be heavily biased by highly expressed regions or global differences in read depth.
# This approach assumes all samples have the same composition — which is rarely true in ChIP/CLIP/PLAT-seq.
# It does not account for differences in signal vs background.
# ✅ Gives the most extreme values.

# 2. DESeq2 RLE Normalization (norm1)
# Uses Relative Log Expression: assumes most features are not differentially bound, and centers log ratios around zero.
# Strongly penalizes global shifts in signal — espeially effective if there's a systematic enrichment in one condition (e.g., Cis13 has stronger signal).
# Sensitive to changes in composition, not just depth.
# ✅ This is why your Cis13 samples get high factors — they appear globally "heavier" and need to be scaled down.

# 3. edgeR TMM Normalization (norm2)
# Trimmed Mean of M-values: like RLE but robust to outliers and large fold-changes.
# Excludes extreme values when computing scaling factors.
# Generally more moderate and stable across datasets with some imbalance.
# ✅ This gives more balanced correction, often less aggressive than DESeq2.

############################################################
############################################################

# Option 1: Default DESeq2
norm1 <- dba.normalize(dba_obj, 
                      method=DBA_DESEQ2,
                      normalize=DBA_NORM_RLE)

# Option 2: edgeR with TMM
norm2 <- dba.normalize(dba_obj, 
                      method=DBA_EDGER,
                      normalize=DBA_NORM_TMM)

# Option 3: Library size only
norm3 <- dba.normalize(dba_obj, 
                      method=DBA_DESEQ2,
                      normalize=DBA_NORM_LIB)

############################################################
############################################################

# If a design formula is specified, it must be composed from the following allowable factors:
#
#            • Tissue
#
#            • Factor
#
#            • Condition
#
#            • Treatment
#
#            • Replicate
#
#            • Caller

# Categories: when automatically generating contrasts, attribute or vector of attributes to base contrasts on:
#
#            • DBA_ID
#
#            • DBA_TISSUE
#
#            • DBA_FACTOR
#
#            • DBA_CONDITION
#
#            • DBA_TREATMENT
#
#            • DBA_REPLICATE
#
#            • DBA_CALLER


# Yes, there is a 1-to-1 correspondence between the design formula factors and the categories constants:
# Case sensitive: Condition ≠ condition
# Exact spelling: Must match exactly
# Sample sheet columns: Must use the design factor names as column headers
# Categories constants: Use DBA_ prefix in R code

# If you have :
# dba_contrast <- dba.contrast(dba_obj, categories=DBA_CONDITION)
# You should make sure your design formula is compatible:
# dba_obj <- dba(dba_obj, design = ~Condition)

############################################################
############################################################ define the contrast 13P vs Combo treatment

# You don't have to run dba.normalize() before dba.contrast(), but the order affects how normalization is applied.

# Workflow 1: Normalize First (Recommended)
# 1. Normalize first
# dba_obj <- dba.normalize(dba_obj, 
#                        method=DBA_DESEQ2,
#                        normalize=DBA_NORM_RLE)

# 2. Set up contrasts
# dba_obj <- dba.contrast(dba_obj, categories=DBA_CONDITION)

# 3. Analyze
# dba_obj <- dba.analyze(dba_obj)

# Workflow 2: Normalize During Analysis
# 1. Set up contrasts first
# dba_obj <- dba.contrast(dba_obj, categories=DBA_CONDITION)

# 2. Normalize and analyze together
# dba_obj <- dba.analyze(dba_obj, 
#                      method=DBA_DESEQ2,
#                      normalize=DBA_NORM_RLE)

############################################################

# a. Normalize first
dba_obj <- dba.normalize(dba_obj, 
                        method=DBA_DESEQ2,
                        normalize=DBA_NORM_RLE)

# b. Set up contrasts
dba_obj <- dba.contrast(dba_obj, 
                        categories=DBA_CONDITION,
                        minMembers = 2) 

# c. Analyze
dba_obj <- dba.analyze(dba_obj, 
                       bParallel = TRUE,
                       bBlacklist = FALSE,
                       bGreylist = FALSE)  

# d. Check normalization effects
dba.plotPCA(dba_obj, label=DBA_CONDITION)

# e. Show
dba.show(dba_obj, bContrasts=TRUE)

#     Factor Group Samples Group2 Samples2 DB.DESeq2
# 1 Condition Cis13       2    P13        2      2594

############################################################
############################################################

# retrieve significant sites (FDR < 0.05)

db_results <- dba.report(
  dba_obj,
  method  = DBA_DESEQ2,
  th      = 0.05,      
  bUsePval = FALSE     
)

# head(db_results,2)
# GRanges object with 2 ranges and 6 metadata columns:
#      seqnames              ranges strand |      Conc Conc_Cis13  Conc_P13
#         <Rle>           <IRanges>  <Rle> | <numeric>  <numeric> <numeric>
#  347     chr1   58783054-58783876      * |   11.1666    12.0809   8.05175
#  590     chr1 211577949-211579384      * |   12.4559    13.1557  11.04362
#           Fold      p-value          FDR
#      <numeric>    <numeric>    <numeric>
#  347   4.01963 5.03350e-164 3.60902e-160
#  590   2.09269 1.65827e-116 5.94489e-113
#  -------
#  seqinfo: 25 sequences from an unspecified genome; no seqlengths

# tail(db_results,2)
# GRanges object with 2 ranges and 6 metadata columns:
#       seqnames            ranges strand |      Conc Conc_Cis13  Conc_P13
#          <Rle>         <IRanges>  <Rle> | <numeric>  <numeric> <numeric>
#  4820    chr22 20116873-20117256      * |   9.01599    8.82372   9.18562
#  1777    chr14 23047272-23048116      * |   7.96311    8.18123   7.70605
#            Fold   p-value       FDR
#       <numeric> <numeric> <numeric>
#  4820 -0.349808 0.0180740 0.0499646
#  1777  0.448027 0.0180765 0.0499646
#  -------
#  seqinfo: 25 sequences from an unspecified genome; no seqlengths

length(db_results)
sum(db_results$Fold>0)
sum(db_results$Fold<0)

# [1] 2594
# [1] 1280
# [1] 1314

############################################################
############################################################ write to bed

rtracklayer::export(db_results, "16h_Diffbind_13P_vs_combo_diffSites.bed")

results_df <- as.data.frame(db_results)
write.csv(results_df, file = "16h_DiffBind_13P_vs_combo_diffSites_results.csv", row.names = FALSE)

dba.plotMA(dba_obj)         # MA-plot
dba.plotVolcano(dba_obj)    # volcano plot

png("16h_Diffbind_13P_vs_combo_diffSites_dba_MAplot.png", width=800, height=800)
dba.plotMA(dba_obj)
dev.off()

png("16h_Diffbind_13P_vs_combo_diffSites_dba_VolcanoPlot.png", width=800, height=800)
dba.plotVolcano(dba_obj)
dev.off()

############################################################
############################################################

dba.show(dba_obj, bContrast = TRUE)

fold <- mcols(db_results)$Fold          # extract log2 Fold column
gain_in_comp <- db_results[ fold > 0 ]  # enriched in competition Combo  
gain_in_13P  <- db_results[ fold < 0 ]  # enriched in 13P experiment   

length(gain_in_comp)   
length(gain_in_13P)    

rtracklayer::export(gain_in_comp, "16h_Diffbind_13P_vs_combo_diffSites_combo_gainSites.bed")
rtracklayer::export(gain_in_13P,  "16h_Diffbind_13P_vs_combo_diffSites_13P_gainSites.bed")

############################################################  
############################################################  display the HEATMAP of DB loci

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
dba.plotHeatmap(dba_obj,  contrast=1, correlations=FALSE, scale="row", colScheme=hmap)

############################################################
############################################################
############################################################
############################################################ extracting NORMALIZED COUNTS

samples <- read.csv("samples_PLAT-seq.csv")

dba_obj2 <- dba(sampleSheet = samples)

# Count with normalized scores
dba_obj2 <- dba.count(
  dba_obj2,
  summits = 0,
  score = DBA_SCORE_NORMALIZED,
  bParallel = TRUE
)

# Extract normalized counts
counts_norm2 <- dba.peakset(dba_obj2, bRetrieve = TRUE)

# Convert to data frame
count_data_norm2 <- as.data.frame(counts_norm2)

# Set rownames to genomic coordinates
rownames(count_data_norm2) <- paste0(count_data_norm2$seqnames, ":", count_data_norm2$start, "-", count_data_norm2$end)

# Preview
head(count_data_norm2)
tail(count_data_norm2)

write.csv(count_data_norm2, file = "16_DiffBind_normalized_counts_dba_score.csv", row.names = TRUE)

############################################################
############################################################

# 📊 Boxplot: normalized counts

# Select columns with sample counts
x_cols_norm2 <- grep("^X", names(count_data_norm2), value = TRUE)

# Boxplot
boxplot(count_data_norm2[, x_cols_norm2],
        las = 2,
        main = "Boxplot of Normalized Read Counts per Sample",
        ylab = "Normalized Read Counts",
        col = "lightgreen",
        outline = FALSE)

# Save as PNG
png("boxplot_normalized_counts2.png", width = 800, height = 600)
boxplot(count_data_norm2[, x_cols_norm2],
        las = 2,
        main = "Boxplot of Normalized Read Counts per Sample",
        ylab = "Normalized Read Counts",
        col = "lightgreen",
        outline = FALSE)
dev.off()

# 📊 Boxplot: log2-transformed normalized counts

boxplot(log2(count_data_norm2[, x_cols_norm2] + 1),
        las = 2,
        main = "Log2-Transformed Normalized Read Counts",
        ylab = "Log2(Normalized Read Counts + 1)",
        col = "skyblue",
        outline = FALSE)

# Save as PNG
png("boxplot_log2_normalized_counts2.png", width = 800, height = 600)
boxplot(log2(count_data_norm2[, x_cols_norm2] + 1),
        las = 2,
        main = "Log2-Transformed Normalized Read Counts",
        ylab = "Log2(Normalized Read Counts + 1)",
        col = "skyblue",
        outline = FALSE)
dev.off()

############################################################
############################################################

dba.plotHeatmap(dba_obj2, correlations = FALSE, scale = "row")
dba.plotPCA(dba_obj2, DBA_CONDITION, label = DBA_ID)


png("heatmap_dba_obj2_norm.png", width = 1000, height = 800)
dba.plotHeatmap(dba_obj2, correlations = FALSE, scale = "row")
dev.off()

# 📊 Save PCA plot for dba_obj2, colored by condition and labeled by ID:

png("pca_dba_obj2_norm.png", width = 1000, height = 800)
dba.plotPCA(dba_obj2, DBA_CONDITION, label = DBA_ID)
dev.off()

############################################################
############################################################ Why Use a Blocking Factor?

# In experiments where samples are not independent, a blocking factor lets you:
# Model and remove variation due to a known nuisance variable (like individual differences, batches, etc.)
# Avoid confounding, so the effect you care about (like treatment) isn't masked or inflated by these other sources.

# 🧪 Common Examples:
# Scenario  Blocking Factor
# Paired samples (before vs. after per patient) PatientID
# Technical replicates from the same library  ReplicateID or Batch

# A blocking factor is a known grouping variable that captures non-interesting, structured variability in your data. 
# You include it in your model so that it doesn’t interfere with your ability 
# to detect real differences from the factor you do care about (like treatment or condition).

############################################################
############################################################

# Blocking Methods (For Complex Designs)

# . DBA_DESEQ2_BLOCK
# rdba_obj <- dba.analyze(dba_obj, method=DBA_DESEQ2_BLOCK)

# . DBA_EDGER_BLOCK
# dba_obj <- dba.analyze(dba_obj, method=DBA_EDGER_BLOCK)


# for example : a Complete Workflow Context
# Set up contrast with blocking design
# dba_obj <- dba.contrast(dba_obj, 
#                       design="~Replicate + Condition",
#                       block=DBA_REPLICATE)

# Option 1: DESeq2 with blocking
# dba_obj <- dba.analyze(dba_obj, method=DBA_DESEQ2_BLOCK)

# Option 2: edgeR with blocking (alternative)
# dba_obj <- dba.analyze(dba_obj, method=DBA_EDGER_BLOCK)

############################################################  
############################################################  Other normalization methods :

# *. Spike-in normalization

# *. Parallel factor normalization
# We can perform a ChIP on the same sample for an alternative factor that is known not to change 
# its binding patterns under the conditions of our experiment. 

# *. "housekeeping" normalization 
# In THOR tool, where the "housekeeping" normalization is based on focusing on histone marks
# such as H3K4me3 that should be consistently bound in the promoter regions associated with
# housekeeping genes

############################################################ Quantile normalization
############################################################ Other normalization methods :

# Quantile normalization is a statistical technique used to make the distribution of values in different samples identical 
# in statistical properties, especially useful for high-throughput data like gene expression or ChIP-seq signal.

# 🔹 What It Does (Intuition)
# It forces all samples to have the same distribution of values — typically the same quantiles (percentiles).
# Imagine sorting the expression values in each sample. Quantile normalization makes sure that:
# The lowest value across all samples becomes the same.
# The 2nd lowest value across all samples becomes the same.
# And so on...
# This is done by replacing each value with the average of the values at that rank across all samples.

############################################################ Loess normalization
############################################################ Other normalization methods :

# LOESS normalization (Locally Estimated Scatterplot Smoothing) is a method that corrects systematic biases 
# in high-throughput experiments, especially two-color microarrays, and occasionally ChIP-seq or RNA-seq log-ratios 
# (after transformation). It’s a type of local regression used to remove intensity-dependent or spatial biases.

# 🔍 What It Does
# LOESS normalization adjusts values locally based on a smooth curve, rather than applying a global scaling factor.

# Typically, it's used to correct:
# Intensity-dependent bias (e.g., brighter spots are consistently biased high or low)

# M vs A plot bias, where:
# M = log2(sample1/sample2) → log-ratio
# A = average log-intensity (log2((sample1 + sample2)/2))

# LOESS fits a curve through the M vs A plot to remove trends — the idea is that the M values should center around 0 
# (no systematic bias), unless true differential expression/binding is present.

# 📈 Visual Explanation

# In an MA-plot:
# If most points are above or below the M=0 line, that’s a systematic bias.
# LOESS normalization fits a smooth curve through the cloud of points and subtracts this curve, flattening the bias.

# ⚠️ Not for Raw Counts
# LOESS requires continuous data, like log-transformed intensities, so it should not be applied to raw RNA-seq counts. 
# Instead, it's used after transformation (e.g., on log ratios or normalized scores).

############################################################  
############################################################  