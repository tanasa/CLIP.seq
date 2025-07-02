#####################################################################
#####################################################################

# BAM files :

# 16h_A27_13_1.Aligned.sortedByCoord.out.bam
# 16h_A27_13_2.Aligned.sortedByCoord.out.bam
# 16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam
# 16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam

# PEAK files :

# 31002 16h_A27_13_1_16h_A27_13_2.pval.0.01.narrowPeaks
# 37574 16h_A27_Cis13_1_16h_A27_Cis13_2.pval.0.01.narrowPeaks

# 10998 16h_A27_13_1.pval.0.01.narrowPeaks
# 12920 16h_A27_13_2.pval.0.01.narrowPeaks
# 12400 16h_A27_Cis13_1.pval.0.01.narrowPeaks
# 9030 16h_A27_Cis13_2.pval.0.01.narrowPeaks

#####################################################################
##################################################################### methods to count the reads on the peaks in R
#####################################################################

# For production pipelines, summarizeOverlaps or featureCounts (which are more optimized for I/O and large datasets) 
# are generally preferred. 

#####################################################################
#####################################################################
##################################################################### 

# findOverlaps() tells you which features overlap which events
# countOverlaps() tells you how many events overlap each feature

     ## S4 method for signature 'GenomicRanges,GenomicRanges'
     # findOverlaps(query, subject,
     #    maxgap=-1L, minoverlap=0L,
     #    type=c("any", "start", "end", "within", "equal"),
     #    select=c("all", "first", "last", "arbitrary"),
     #    ignore.strand=FALSE)

     ## S4 method for signature 'GenomicRanges,GenomicRanges'
     # countOverlaps(query, subject,
     #    maxgap=-1L, minoverlap=0L,
     #    type=c("any", "start", "end", "within", "equal"),
     #    ignore.strand=FALSE)

library(GenomicAlignments)
library(rtracklayer)
library(GenomicRanges)

# Read peak file
peaks <- import("16h_A27_13_1_16h_A27_13_2.pval.0.01.narrowPeaks", 
                format = "narrowPeak")

# Keep the canonical chromosomes 
peaks_canonical <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
export(peaks_canonical, "16h_A27_13_1_16h_A27_13_2.pval.0.01.canonical.narrowPeaks", 
                        format = "narrowPeak")

peaks = peaks_canonical


# Remove chrM and scaffold sequences
# chroms_to_remove <- c("chrM", 
#                      grep("^GL", seqlevels(peaks_canonical), value = TRUE),
#                      grep("^KI", seqlevels(peaks_canonical), value = TRUE))

# peaks_canonical <- dropSeqlevels(peaks_canonical, chroms_to_remove, pruning.mode = "coarse")

# Check result
table(seqnames(peaks))

# BAM file paths
bam_files <- c("16h_A27_13_1.Aligned.sortedByCoord.out.bam",
               "16h_A27_13_2.Aligned.sortedByCoord.out.bam", 
               "16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam",
               "16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam")

# Count reads in peaks for each BAM file
count_matrix <- matrix(nrow = length(peaks), ncol = length(bam_files))
colnames(count_matrix) <- gsub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files))

#####################################################################
#####################################################################
##################################################################### METHOD1 : overlaps() function
##################################################################### library(GenomicRanges)
#####################################################################
#####################################################################

for(i in 1:length(bam_files)) {
    reads <- readGAlignments(bam_files[i])
    overlaps <- countOverlaps(peaks, reads)
    count_matrix[,i] <- overlaps
}

# Create results data frame
results_method1 <- data.frame(
    chr = as.character(seqnames(peaks)), # Chromosome name
    start = start(peaks),               # Start position
    end = end(peaks),                   # End position
    strand = as.character(strand(peaks)), # Strand (if applicable, though '*' in your example)
    name = mcols(peaks)$name,           # Peak name (e.g., peak_0, peak_1)
    score = mcols(peaks)$score,         # Score
    signalValue = mcols(peaks)$signalValue, # Signal Value
    pValue = mcols(peaks)$pValue,       # pValue
    qValue = mcols(peaks)$qValue,       # qValue
    peak = mcols(peaks)$peak,           # Peak summit offset
    count_matrix                        # Your read count matrix
)

# You can then check the structure and head of the new data frame
str(results_method1)
head(results_method1, 2)
tail(results_method1, 2)

write.csv(results_method1, "the.matrix.peaks_and_counts.method1.csv", row.names = FALSE)

#####################################################################
#####################################################################
############################# Write summary statistics : compute FRIP
#####################################################################
#####################################################################

sample_names <- colnames(count_matrix)
summary_stats <- data.frame(
    Sample = sample_names,
    Total_Reads = colSums(count_matrix),
    Peaks_with_Reads = colSums(count_matrix > 0),
    Mean_Reads_per_Peak = round(colMeans(count_matrix), 2)
)
summary_stats$FrIP = summary_stats$Peaks_with_Reads / summary_stats$Total_Reads * 100 

write.table(summary_stats, 
            file = "the.matrix.peaks_and_counts.method1.summary_stats.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("GenomicAlignments results written to:\n")
cat("- genomicalignments_peak_counts.txt\n")
cat("- genomicalignments_peak_counts.csv\n")
cat("- peaks_with_counts.bed\n")
cat("- genomicalignments_summary_stats.txt\n")

#####################################################################
#####################################################################
##################################################################### library(GenomicAlignments)
##################################################################### summarizeOverlaps
#####################################################################
#####################################################################

#     ‘summarizeOverlaps’ extends ‘findOverlaps’ by providing options to
#     resolve reads that overlap multiple features.

# Count reads
counts2 <- summarizeOverlaps(peaks, bam_files, mode = "Union")
count_matrix2 <- assay(counts2)
colnames(count_matrix2) <- gsub(".Aligned.sortedByCoord.out.bam", "", basename(bam_files))

# Create results data frame
results_method2 <- data.frame(
    chr = as.character(seqnames(peaks)), # Chromosome name
    start = start(peaks),               # Start position
    end = end(peaks),                   # End position
    strand = as.character(strand(peaks)), # Strand (if applicable, though '*' in your example)
    name = mcols(peaks)$name,           # Peak name (e.g., peak_0, peak_1)
    score = mcols(peaks)$score,         # Score
    signalValue = mcols(peaks)$signalValue, # Signal Value
    pValue = mcols(peaks)$pValue,       # pValue
    qValue = mcols(peaks)$qValue,       # qValue
    peak = mcols(peaks)$peak,           # Peak summit offset
    count_matrix2                        # Your read count matrix
)

# You can then check the structure and head of the new data frame
str(results_method2)
head(results_method2, 2)
tail(results_method2, 2)

# Write results

write.csv(results_method2, "the.matrix.peaks_and_counts.method2.summarizeOverlaps.csv",  
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#####################################################################
#####################################################################
##################################################################### featureCounts
##################################################################### library(Rsubread)
#####################################################################
#####################################################################

library(Rsubread)
library(GenomicAlignments)
library(rtracklayer)
library(GenomicRanges)
library(SummarizedExperiment)

# Why SAF is Required:
# featureCounts() expects specific input formats:
# GTF/GFF files (gene annotation)
# SAF format (Simple Annotation Format)
# ✅ SAF format must have 5 required columns:

# GeneID – unique identifier per feature (e.g., "peak_1")
# Chr – chromosome (e.g., "chr1")
# Start – 1-based start coordinate
# End – end coordinate
# Strand – "+", "-" or "." (for unknown)

# Read and convert narrowPeak to SAF format
peaks_df <- data.frame(peaks)

# Create SAF format
saf <- data.frame(
    GeneID = paste0("peak_", 1:nrow(peaks_df)),
    Chr = peaks_df$seqnames,
    Start = peaks_df$start + 1,  # Convert to 1-based
    End = peaks_df$end,
    Strand = "."
)

# Count reads
counts <- featureCounts(
    files = c("16h_A27_13_1.Aligned.sortedByCoord.out.bam",
              "16h_A27_13_2.Aligned.sortedByCoord.out.bam", 
              "16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam",
              "16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam"),
    annot.ext = saf,
    isPairedEnd = FALSE,
    nthreads = 4
)

# Check what's actually in the counts object
str(counts)
# Check if it's really a list
is.list(counts)
# See what elements are in the list
names(counts)
# Check the length
length(counts)

# str(counts)
# List of 4
# $ counts    : int [1:31002, 1:4] 8 7 18 10 14 4 6 7 49 15 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:31002] "peak_1" "peak_2" "peak_3" "peak_4" ...
#  .. ..$ : chr [1:4] "16h_A27_13_1.Aligned.sortedByCoord.out.bam" "16h_A27_13_2.Aligned.sortedByCoord.out.bam" "16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam" "16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam"
# $ annotation:'data.frame': 31002 obs. of  6 variables:
#  ..$ GeneID: chr [1:31002] "peak_1" "peak_2" "peak_3" "peak_4" ...
#  ..$ Chr   : chr [1:31002] "chr1" "chr1" "chr1" "chr1" ...
#  ..$ Start : int [1:31002] 16736 189979 631077 632840 634242 812808 821013 822058 904855 905846 ...
# ..$ End   : int [1:31002] 17109 190109 631244 632991 634459 813023 821272 822556 905637 906073 ...
#  ..$ Strand: chr [1:31002] "." "." "." "." ...
#  ..$ Length: int [1:31002] 374 131 168 152 218 216 260 499 783 228 ...
# $ targets   : chr [1:4] "16h_A27_13_1.Aligned.sortedByCoord.out.bam" "16h_A27_13_2.Aligned.sortedByCoord.out.bam" "16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam" "16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam"
# $ stat      :'data.frame': 14 obs. of  5 variables:
#  ..$ Status                                       : chr [1:14] "Assigned" "Unassigned_Unmapped" "Unassigned_Read_Type" "Unassigned_Singleton" ...
#  ..$ 16h_A27_13_1.Aligned.sortedByCoord.out.bam   : int [1:14] 4875919 0 0 0 0 0 0 0 0 0 ...
#  ..$ 16h_A27_13_2.Aligned.sortedByCoord.out.bam   : int [1:14] 5253593 0 0 0 0 0 0 0 0 0 ...
#  ..$ 16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam: int [1:14] 22330965 0 0 0 0 0 0 0 0 0 ...
#  ..$ 16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam: int [1:14] 20616403 0 0 0 0 0 0 0 0 0 ...

peak_count_table <- cbind(counts$annotation, counts$counts)

# Clean column names
colnames(peak_count_table) <- gsub(".Aligned.sortedByCoord.out.bam", "", 
                              basename(colnames(peak_count_table)))

# Check the cleaned names
colnames(peak_count_table)

write.csv(peak_count_table, 
          "the.matrix.peaks_and_counts.method3.featureCounts.csv", 
          row.names = FALSE)

#####################################################################
#####################################################################
#####################################################################
##################################################################### diffbind

library(DiffBind)

# Step 1: Define samples
samples <- data.frame(
  SampleID = c("13_1", "13_2", "Cis13_1", "Cis13_2"),
  Tissue = "A27",
  Factor = "Drug",                   # Not used here, just placeholder
  Condition = c("P13", "P13", "Cis13", "Cis13"),
  Replicate = c(1, 2, 1, 2),
  bamReads = c(
    "16h_A27_13_1.Aligned.sortedByCoord.out.bam",
    "16h_A27_13_2.Aligned.sortedByCoord.out.bam",
    "16h_A27_Cis13_1.Aligned.sortedByCoord.out.bam",
    "16h_A27_Cis13_2.Aligned.sortedByCoord.out.bam"
  ),
  Peaks = rep("16h_A27_13_1_16h_A27_13_2.pval.0.01.canonical.narrowPeaks", 4),
  PeakCaller = rep("narrow", 4)
)

# Step 2: Save sample sheet
write.csv(samples, "samples_diffbind.csv", row.names = FALSE)

# Step 3: Load into DiffBind
dba_obj <- dba(sampleSheet = "samples_diffbind.csv")

# DBA_SCORE_READS will give raw counts.
# If you have paired-end reads, set fragmentSize to 0 or appropriate length.
# If you don't care about strand, set mapRate = 0.

counts_dba <- dba.count(
  DBA = dba_obj,
  summits = FALSE,         # Set to an integer (e.g., 150) to extend around peak summits
  score = DBA_SCORE_READS, # Report raw read counts
  bUseSummarizeOverlaps = TRUE, # Use Rsubread's summarizeOverlaps via GenomicAlignments, which is efficient
  # fragmentSize = 0, # Auto-detect fragment size for paired-end, or provide fixed size for single-end.
                      # Or set to 0 for paired-end data to use fragment length.
  minOverlap = 1,     # Minimum overlap (in base pairs) between read and peak to be counted
  filter = 0          # Do not filter peaks based on minimum reads at this stage
)

# Step 5: Extract and save count matrix
counts <- dba.peakset(counts_dba, bRetrieve = TRUE)

# Format into data frame
count_matrix <- as.data.frame(counts)
write.csv(count_matrix, 
          "the.matrix.peaks_and_counts.method4.diffBind.csv", 
           row.names = FALSE)

# Optional: View summary
dba.show(counts_dba)

#####################################################################
#####################################################################
##################################################################### 

# dba.count(DBA)
#            DBA_SCORE_NORMALIZED     normalized reads, as set by ‘dba.normalize’                                           
#            DBA_SCORE_READS          raw read count for interval using only reads from ChIP                                
#            DBA_SCORE_CONTROL_READS  raw read count for interval using only reads from Control                             
#            DBA_SCORE_READS_FOLD     raw read count for interval from ChIP divided by read count for interval from control 
#            DBA_SCORE_READS_MINUS    raw read count for interval from ChIP minus read count for interval from control 

#####################################################################
#####################################################################