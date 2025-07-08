library(DESeq2)
library(ggplot2)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) # gene annotation 
library(cowplot)
library(tidyverse)

# Read the file
peak_data <- read.table("the.matrix.peaks_and_counts.method1.csv", 
                 header = TRUE, 
                 sep = ",",    # Use sep = "\t" if it's a tab-separated file
                 stringsAsFactors = FALSE)

# View structure
str(peak_data)
head(peak_data, 2)
tail(peak_data, 2)

# Peak Annotation and Assignment

peaks_signal200 <- subset(peak_data, signalValue >= 200)
peaks_signal300 <- subset(peak_data, signalValue >= 300)
peaks_signal400 <- subset(peak_data, signalValue >= 400)

# Check the number of peaks in each subset
cat("Peaks with signal ≥ 200:", nrow(peaks_signal200), "\n")
cat("Peaks with signal ≥ 300:", nrow(peaks_signal300), "\n")
cat("Peaks with signal ≥ 400:", nrow(peaks_signal400), "\n")

peaks = peaks_signal200 # use all the peaks called by Generich
head(peaks_signal200)
dim(peaks_signal200)

library(ggplot2)

# Select the relevant columns
summary_data <- peaks[, c("score", "signalValue", "pValue", "qValue", "peak")]
# Display summary statistics
summary(summary_data)

# Plot: Signal Value vs pValue
ggplot(peaks, aes(x = signalValue, y = pValue)) +
  geom_point(color = "darkblue") +
  labs(x = "Signal Value", y = "-log10(p-value)", title = "Peak Signal vs Significance") +
  theme_minimal() + 
  xlim(0, 5000) +
  ylim(0, 20)



# library(ChIPseeker)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# c(-1000, 100)	Core promoter	Strict promoter-only binding

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) # For gene annotation 

peaks_gr <- GRanges(
  seqnames = peaks$chr,
  ranges = IRanges(start = peaks$start, end = peaks$end),
  strand = peaks$strand,
  mcols = peaks[, c("name", "score", "signalValue", "pValue", "qValue", "peak", 
                    "X16h_A27_13_1", "X16h_A27_13_2",
                    "X16h_A27_Cis13_1", "X16h_A27_Cis13_2")] # metadata columns
)

# Keep peaks on the canonical chromosomes
canonical_chrs <- paste0("chr", c(1:19, "X", "Y"))
peaks_gr_canonical <- peaks_gr[seqnames(peaks_gr) %in% canonical_chrs]

peaks_gr = peaks_gr_canonical
seqnames(peaks_gr)

# Introducing a second variable
gr_peaks = peaks_gr
seqnames(gr_peaks)

head(peaks, 2)

# Performing differential analysis : setting up DESeq2 analysis

# Extract the count columns only
count_columns <- grep("^X16h", colnames(peaks))  # assuming all sample names start with X16h
counts <- as.matrix(peaks[, count_columns])

# Assign rownames using peak name (or coordinates)
rownames(counts) <- peaks$name  # or paste0(peak_df$chr, ":", peak_df$start, "-", peak_df$end)
head(counts)

colnames(counts)
colnames(counts) <- c("A27_13_1", "A27_13_2", "A27_Cis13_1", "A27_Cis13_2")

col_data <- data.frame(
  row.names = colnames(counts),
  condition = c("A27_13", "A27_13", "A27_Cis13", "A27_Cis13")
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ condition)

col_data$condition <- relevel(factor(col_data$condition), ref = "A27_13")

# QC checks : Filter low count peaks
cat("Peaks retained before filtering:", nrow(dds), "\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Peaks retained after filtering:", nrow(dds), "\n")

# Check library sizes
cat("Library sizes:\n")
print(colSums(counts(dds)))

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "A27_Cis13", "A27_13"))

# The order of the contrast means:
# log2FoldChange = log2(expression in A27_Cis13 / expression in A27_13)
# Interpretation:
# A positive log2FoldChange means: expression is higher in A27_Cis13
# A negative log2FoldChange means: expression is higher in A27_13

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# View top 6 results
head(res_ordered)

cat("\n=== Size Factors ===\n")
size_factors <- sizeFactors(dds)
print(size_factors)

# Extract normalized counts from DESeqDataSet
norm_counts <- counts(dds, normalized = TRUE)

# Convert normalized counts matrix to a data frame and add peak_id as a new column
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$peak_id <- rownames(norm_counts_df)

# Preview the first few rows of normalized counts
cat("Preview of normalized counts:\n")
head(norm_counts_df, 2)

write.csv(norm_counts_df, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_norm_counts.csv", row.names = FALSE)

print("plotting dispersion")

# Combine side-by-side
options(repr.plot.width = 6, repr.plot.height = 6)

plotDispEsts(dds)
dispersionFunction(dds)

# Total number of peaks analyzed
cat("Total peaks analyzed:", nrow(res), "\n")

# Handle NA values (optional, but safer)
res_clean <- res[!is.na(res$pvalue) & !is.na(res$padj), ]
print(nrow(res_clean)) 

# --- By raw p-value ---
cat("Peaks with p-value < 0.05:", sum(res_clean$pvalue < 0.05), "\n")
cat("Peaks with p-value < 0.1 :", sum(res_clean$pvalue < 0.1), "\n")

# --- By adjusted p-value (FDR) ---
cat("Peaks with padj < 0.05:", sum(res_clean$padj < 0.05), "\n")
cat("Peaks with padj < 0.1 :", sum(res_clean$padj < 0.1), "\n")

# --- By log2FC (upregulated) ---
print("")
cat("Peaks with log2FC > 1 :", sum(res_clean$log2FoldChange > 1), "\n")
cat("Peaks with log2FC > 2 :", sum(res_clean$log2FoldChange > 2), "\n")

# --- Combined filters: log2FC > 1/2 and p-value ---
cat("log2FC > 1 & p-value < 0.05 :", sum(res_clean$log2FoldChange > 1 & res_clean$pvalue < 0.05), "\n")
cat("log2FC > 2 & p-value < 0.05 :", sum(res_clean$log2FoldChange > 2 & res_clean$pvalue < 0.05), "\n")
cat("log2FC > 1 & p-value < 0.1  :", sum(res_clean$log2FoldChange > 1 & res_clean$pvalue < 0.1), "\n")

# --- Combined filters: log2FC > 1/2 and padj ---
print("")
cat("log2FC > 1 & padj < 0.05 :", sum(res_clean$log2FoldChange > 1 & res_clean$padj < 0.05), "\n")
cat("log2FC > 2 & padj < 0.05 :", sum(res_clean$log2FoldChange > 2 & res_clean$padj < 0.05), "\n")
cat("log2FC > 1 & padj < 0.1  :", sum(res_clean$log2FoldChange > 1 & res_clean$padj < 0.1), "\n")

# === DOWNREGULATED peaks ===
print("")
cat("Peaks with log2FC < -1 :", sum(res_clean$log2FoldChange < -1), "\n")
cat("Peaks with log2FC < -2 :", sum(res_clean$log2FoldChange < -2), "\n")

# --- Combined filters: log2FC < -1/-2 and p-value ---
cat("log2FC < -1 & p-value < 0.05 :", sum(res_clean$log2FoldChange < -1 & res_clean$pvalue < 0.05), "\n")
cat("log2FC < -2 & p-value < 0.05 :", sum(res_clean$log2FoldChange < -2 & res_clean$pvalue < 0.05), "\n")
cat("log2FC < -1 & p-value < 0.1  :", sum(res_clean$log2FoldChange < -1 & res_clean$pvalue < 0.1), "\n")

# --- Combined filters: log2FC < -1/-2 and padj ---
print("")
cat("log2FC < -1 & padj < 0.05 :", sum(res_clean$log2FoldChange < -1 & res_clean$padj < 0.05), "\n")
cat("log2FC < -2 & padj < 0.05 :", sum(res_clean$log2FoldChange < -2 & res_clean$padj < 0.05), "\n")
cat("log2FC < -1 & padj < 0.1  :", sum(res_clean$log2FoldChange < -1 & res_clean$padj < 0.1), "\n")


# Merging these two dataframes

head(res, 2)
head(peaks, 2)

res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

combined <- merge(peaks, res_df, by.x = "name", by.y = "peak_id")
head(combined, 2)

head(combined, 2)
write.csv(combined, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.csv", row.names = FALSE)

# Combine with Norm Counts : convert normalized counts to a dataframe and add peak ID
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$peak_id <- rownames(norm_counts_df)

# Merge normalized counts with 'combined' data
combined_final <- merge(combined, norm_counts_df, by.x = "name", by.y = "peak_id")

# Preview the result
cat("Preview of final combined table:\n")
print(colnames(combined_final, 2))
print(head(combined_final, 2))
write.csv(combined_final, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.and.norm.csv", row.names = FALSE)





# Clean up: remove rows with NA pvalue or padj
combined_clean <- combined[!is.na(combined$pvalue) & !is.na(combined$padj), ]

# Total number of peaks
cat("Total peaks:", nrow(combined_clean), "\n\n")

# === UPREGULATED ===
cat("UPREGULATED: log2FoldChange > 1\n")
cat("  pvalue < 0.05 :", sum(combined_clean$log2FoldChange > 1 & combined_clean$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(combined_clean$log2FoldChange > 1 & combined_clean$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(combined_clean$log2FoldChange > 1 & combined_clean$padj < 0.1), "\n\n")

# === DOWNREGULATED ===
cat("DOWNREGULATED: log2FoldChange < -1\n")
cat("  pvalue < 0.05 :", sum(combined_clean$log2FoldChange < -1 & combined_clean$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(combined_clean$log2FoldChange < -1 & combined_clean$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(combined_clean$log2FoldChange < -1 & combined_clean$padj < 0.1), "\n")

# Clean the data
combined_clean <- combined[!is.na(combined$padj) & !is.na(combined$log2FoldChange), ]

# === Filter 1: UPREGULATED ===
upregulated <- combined_clean[combined_clean$padj < 0.1 & combined_clean$log2FoldChange > 1, ]
cat("Number of upregulated peaks:", nrow(upregulated), "\n")
# Save to file
# write.csv(upregulated, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.upregulated_padj_lt_0.1_log2FC_gt_plus1.csv", row.names = FALSE)

# === Filter 2: DOWNREGULATED ===
downregulated <- combined_clean[combined_clean$padj < 0.1 & combined_clean$log2FoldChange < -1, ]
cat("Number of downregulated peaks:", nrow(downregulated), "\n")
# Save to file
# write.csv(downregulated, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.downregulated_padj_lt_0.1_log2FC_lt_minus1.csv", row.names = FALSE)

# --- Save UPREGULATED peaks to BED ---
bed_up <- upregulated[, c("chr", "start", "end", "name")]
# write.table(bed_up, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.upregulated.peaks.padj_0.1_log2FC_plus1.bed",
#            sep = "\t", 
#            quote = FALSE, 
#            row.names = FALSE, 
#            col.names = FALSE)

# --- Save DOWNREGULATED peaks to BED ---
bed_down <- downregulated[, c("chr", "start", "end", "name")]
# write.table(bed_down, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.downregulated.peaks.padj_0.1_log2FC_minus1.bed",
#            sep = "\t", 
#            quote = FALSE, 
#            row.names = FALSE, 
#            col.names = FALSE)

head(upregulated, 2)
head(downregulated, 2)

cat("New thresholds : log2FoldChange > 0.5, padj < 0.1")

# Other threholds : log2FoldChange > 0.5, padj < 0.1
dim(combined_clean)

# Filter out rows with NA in padj or log2FC
# Define thresholds

upregulated_peaks <- combined_clean %>%
  filter(log2FoldChange > 0.5, padj < 0.1)

downregulated_peaks <- combined_clean %>%
  filter(log2FoldChange < -0.5, padj < 0.1)

# Print summaries
cat("Upregulated peaks (log2FC > 0.5, padj < 0.1):", nrow(upregulated_peaks), "\n")
cat("Downregulated peaks (log2FC < -0.5, padj < 0.1):", nrow(downregulated_peaks), "\n")

# Print upregulated and downregulated tables
cat("\n=== Upregulated Peaks ===\n")
head(upregulated_peaks, 2)

cat("\n=== Downregulated Peaks ===\n")
head(downregulated_peaks, 2)

# Extract BED format: chr, start, end, name, score, strand

upregulated_peaks_bed <- upregulated_peaks %>% select(chr, start, end, name, score, strand)
downregulated_peaks_bed <- downregulated_peaks %>% select(chr, start, end, name, score, strand)

# Save BED files
write_tsv(upregulated_peaks_bed, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.upregulated_peaks.log2FC.0.5.padj.0.1.bed", col_names = FALSE)
write_tsv(downregulated_peaks_bed, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.downregulated_peaks.log2FC.0.5.padj.0.1.bed", col_names = FALSE)

# Load required libraries
library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)

print("Boxplot of Raw vs log2 Normalized Counts")

# Prepare raw counts
raw_counts <- as.data.frame(counts(dds, normalized = FALSE))
raw_counts$peak_id <- rownames(raw_counts)
raw_long <- pivot_longer(raw_counts, -peak_id, names_to = "Sample", values_to = "Count")
raw_long$log2_count <- log2(raw_long$Count + 1)

# Prepare normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$peak_id <- rownames(norm_counts)
norm_long <- pivot_longer(norm_counts, -peak_id, names_to = "Sample", values_to = "Count")
norm_long$log2_count <- log2(norm_long$Count + 1)

# Create consistent sample color palette
sample_list <- unique(raw_long$Sample)
sample_colors <- setNames(viridis(length(sample_list), option = "D"), sample_list)

# Plot p1: Raw counts
p1 <- ggplot(raw_long, aes(x = Sample, y = log2_count, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7) +
  scale_fill_manual(values = sample_colors, name = "Sample") +
  theme_minimal(base_size = 12) +
  labs(title = "Raw Counts (log2 scale)",
       y = "log2(Counts + 1)", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot p2: Normalized counts
p2 <- ggplot(norm_long, aes(x = Sample, y = log2_count, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7) +
  scale_fill_manual(values = sample_colors, name = "Sample") +
  theme_minimal(base_size = 12) +
  labs(title = "Normalized Counts (log2 scale)",
       y = "log2(Counts + 1)", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots side-by-side
options(repr.plot.width = 12, repr.plot.height = 6)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = 'h')

# Load required libraries
library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)

# === Raw counts (log10 transformed) ===
raw_counts <- counts(dds, normalized = FALSE)
raw_log_counts <- log10(raw_counts + 1)

raw_df <- as.data.frame(raw_log_counts)
raw_df$peak_id <- rownames(raw_df)
raw_long <- pivot_longer(raw_df, -peak_id, names_to = "Sample", values_to = "log10_count")

# === Normalized counts (log10 transformed) ===
norm_counts <- counts(dds, normalized = TRUE)
norm_log_counts <- log10(norm_counts + 1)

norm_df <- as.data.frame(norm_log_counts)
norm_df$peak_id <- rownames(norm_df)
norm_long <- pivot_longer(norm_df, -peak_id, names_to = "Sample", values_to = "log10_count")

# Optional: Define color palette
sample_list <- unique(raw_long$Sample)
sample_colors <- setNames(viridis(length(sample_list), option = "D"), sample_list)

# === Plot 1: Raw counts ===
p1 <- ggplot(raw_long, aes(x = log10_count, color = Sample)) +
  geom_density(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = sample_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density Plot of log10 Raw Counts (Before Normalization)",
    x = "log10(Counts + 1)",
    y = "Density"
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# === Plot 2: Normalized counts ===
p2 <- ggplot(norm_long, aes(x = log10_count, color = Sample)) +
  geom_density(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = sample_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Density Plot of log10 Normalized Counts",
    x = "log10(Normalized Counts + 1)",
    y = "Density"
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# === Set figure size and display ===
options(repr.plot.width = 14, repr.plot.height = 6)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, align = 'h')

# Normalization works as intended —  the right-hand plot shows the samples now have more similar distributions, 
# which is the goal of DESeq2 normalization.



# Apply All Shrinkage Methods + MA & Volcano Plots

# Load required libraries
library(DESeq2)
library(ggplot2)
library(apeglm)
library(ashr)

# Apply LFC shrinkage with different methods

# Run DESeq2
dds2 <- DESeq(dds)

# Check available coefficients
resultsNames(dds2)
# For example, you might see: "Intercept", "condition_A27_Cis13_vs_A27_13"

# Apply shrinkage

# 1. Normal (can use contrast)
res_normal <- lfcShrink(dds2, contrast = c("condition", "A27_Cis13", "A27_13"), type = "normal")

# 2. apeglm (must use coef)
res_apeglm <- lfcShrink(dds2, coef = "condition_A27_Cis13_vs_A27_13", type = "apeglm")

# 3. ashr (must use coef)
res_ashr <- lfcShrink(dds2, coef = "condition_A27_Cis13_vs_A27_13", type = "ashr")

head(res_normal, 2)
head(res_apeglm, 2)
head(res_ashr, 2)

# Step 1: Convert shrinkage results to data frames and add peak_id
res_normal_df  <- as.data.frame(res_normal)
res_apeglm_df  <- as.data.frame(res_apeglm)
res_ashr_df    <- as.data.frame(res_ashr)

res_normal_df$peak_id <- rownames(res_normal_df)
res_apeglm_df$peak_id <- rownames(res_apeglm_df)
res_ashr_df$peak_id   <- rownames(res_ashr_df)

# Step 2: Add padj values from unshrunken results (if not already present)
res_full <- results(dds2, contrast = c("condition", "A27_Cis13", "A27_13"))
padj_vals <- data.frame(
  peak_id = rownames(res_full),
  padj = res_full$padj
)

# Merge padj values into each shrinkage result
res_normal_df  <- merge(res_normal_df, padj_vals, by = "peak_id", suffixes = c("", "_padj"))
res_apeglm_df  <- merge(res_apeglm_df, padj_vals, by = "peak_id", suffixes = c("", "_padj"))
res_ashr_df    <- merge(res_ashr_df, padj_vals, by = "peak_id", suffixes = c("", "_padj"))

# Step 3: Subset relevant columns and rename log2FC
res_normal_df <- res_normal_df[, c("peak_id", "log2FoldChange", "padj")]
res_apeglm_df <- res_apeglm_df[, c("peak_id", "log2FoldChange", "padj")]
res_ashr_df   <- res_ashr_df[, c("peak_id", "log2FoldChange", "padj")]

colnames(res_normal_df)[2] <- "log2FC_normal"
colnames(res_apeglm_df)[2] <- "log2FC_apeglm"
colnames(res_ashr_df)[2]   <- "log2FC_ashr"

# Step 4: Merge all shrinkage results into one
lfc_all <- Reduce(function(x, y) merge(x, y, by = c("peak_id", "padj")),
                  list(res_normal_df, res_apeglm_df, res_ashr_df))

# Step 5: Merge with peaks dataframe
peaks$peak_id <- peaks$name
combined_all <- merge(peaks, lfc_all, by = "peak_id")

# Step 6: Save result
write_tsv(
  combined_all,
  "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.log2FCshrinkage.tsv"
)

# Step 7: Count significantly differential peaks
cat("Number of differential bound peaks after shrinkage\n")

count_sig <- function(df, lfc_col, padj_thresh = 0.1, lfc_thresh = 0.5) {
  up <- sum(df[[lfc_col]] >  lfc_thresh & df$padj < padj_thresh, na.rm = TRUE)
  down <- sum(df[[lfc_col]] < -lfc_thresh & df$padj < padj_thresh, na.rm = TRUE)
  data.frame(Up = up, Down = down, Total = up + down)
}

sig_normal <- count_sig(lfc_all, "log2FC_normal")
sig_apeglm <- count_sig(lfc_all, "log2FC_apeglm")
sig_ashr   <- count_sig(lfc_all, "log2FC_ashr")

sig_summary <- rbind(Normal = sig_normal, Apeglm = sig_apeglm, Ashr = sig_ashr)
cat("")
print(sig_summary)





print("MA plots")

ggplot_MA <- function(res, title, lfc_thresh = 0.5, padj_thresh = 0.1) {
  res_df <- as.data.frame(res)
  res_df$baseMean <- res_df$baseMean + 1e-8  # avoid log(0)
  res_df$peak_id <- rownames(res_df)
  
  # Assign significance categories
  res_df <- res_df %>%
    mutate(
      category = case_when(
        padj < padj_thresh & log2FoldChange >  lfc_thresh ~ "Upregulated",
        padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      log10_baseMean = log10(baseMean)
    )
  
  ggplot(res_df, aes(x = log10_baseMean, y = log2FoldChange, color = category)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_hline(yintercept = 0, color = "black") +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "gray"
    )) +
    coord_cartesian(ylim = c(-5, 5)) +  # Fix y-axis range
    labs(
      title = title,
      x = "log10(Mean of normalized counts)",
      y = "Shrunken log2 Fold Change",
      color = "Significance"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
}

library(cowplot)

# Create individual plots
p1 <- ggplot_MA(res_normal, "MA Plot (Normal)")
p2 <- ggplot_MA(res_apeglm, "MA Plot (apeglm)")
p3 <- ggplot_MA(res_ashr, "MA Plot (ashr)")

# Display side-by-side with shared y-axis scale
options(repr.plot.width = 18, repr.plot.height = 5)
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3)

# ✅ Possible Reasons Why the Plots Look Similar
# 1. High Expression Peaks Dominate the Signal
# All shrinkage methods tend to agree on well-expressed features with strong signal. 
# These are less affected by shrinkage, so the log2 fold changes remain similar across methods.

# 2. Low Count Features Are Already Suppressed
# If DESeq2's prefiltering removed many low-count features, or if your dataset has relatively well-behaved count distributions, 
# there won’t be much room for shrinkage to do its job. All methods then yield similar results.

# 3. Not Enough Variability Between Methods
# In some datasets, especially with:

# Moderate-to-large sample sizes
# Clean replicates
# Few extreme outliers
# → All shrinkage methods produce comparable results.

# 4. Plot Limits Mask Subtle Differences
# The axis ranges in the plots (e.g. ylim = c(-3, 3) or xlim = c(-5, 5)) may compress visible differences. The methods may differ more in the extreme values, which get visually flattened due to shared limits.

# print("Volcano plots")

# log2FoldChange = 0.5
# padj = 0.1

# === Volcano Plot Function ===
# plot_volcano <- function(res_obj, title) {
#  res_df <- as.data.frame(res_obj)
#  res_df$peak_id <- rownames(res_df)
#  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
#    geom_point(alpha = 0.5, color = "grey40") +
#    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
#    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
#    theme_minimal(base_size = 14) +
#    labs(title = title, x = "Shrunken log2 Fold Change", y = "-log10 Adjusted p-value")
# }

# === Show volcano plots ===
# library(cowplot)
# volcano_normal <- plot_volcano(res_normal, "Volcano Plot (Normal)")
# volcano_apeglm <- plot_volcano(res_apeglm, "Volcano Plot (apeglm)")
# volcano_ashr <- plot_volcano(res_ashr, "Volcano Plot (ashr)")

# options(repr.plot.width = 18, repr.plot.height = 5)
# plot_grid(volcano_normal, volcano_apeglm, volcano_ashr, labels = c("A", "B", "C"), ncol = 3)

# Volcano plot
library(EnhancedVolcano)

# Set thresholds
pval_cutoff <- 0.1
lfc_cutoff <- 0.5
# Optional: convert to data frame for label inspection
res_clean_df <- as.data.frame(res_clean)

options(repr.plot.width = 12, repr.plot.height = 10)
EnhancedVolcano(res_clean_df,
    lab = rownames(res_clean_df),
    x = 'log2FoldChange',
    y = 'padj',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'adjusted p-value'),
    title = 'A27_Cis13 vs A27_13 (DESeq2)',
    pCutoff = pval_cutoff,
    FCcutoff = lfc_cutoff,
    pointSize = 2.0,
    labSize = 3.0,
    boxedLabels = FALSE,
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 3.0,
    colAlpha = 0.9,
    ylim = c(0, 100),  # <-- Extend Y-axis here
    col = c("grey60", "forestgreen", "royalblue", "red3")
)

print("The pheatmap of the most down-regulated peaks :")
library(pheatmap)

# No Z-score standardization

# Step 1: Sort by most negative log2FoldChange
# res_down <- combined_final[order(combined_final$log2FoldChange), ]

# Step 2: Select top 300 most downregulated peaks
# top300_down <- head(res_down, 300)

# Step 3: Extract normalized expression matrix
# sample_cols <- c("A27_13_1", "A27_13_2", "A27_Cis13_1", "A27_Cis13_2")
# mat_down <- as.matrix(top300_down[, sample_cols])

# Step 4: Column annotation
# annotation_col <- data.frame(
#  Condition = c("A27_13", "A27_13", "A27_Cis13", "A27_Cis13")
# )
# rownames(annotation_col) <- sample_cols

# Step 5: Plot heatmap (without scaling)
# pheatmap(mat_down,
#         cluster_rows = FALSE,
#         cluster_cols = TRUE,
#         annotation_col = annotation_col,
#         show_rownames = FALSE,
#         fontsize_row = 6,
#         fontsize_col = 10,
#         main = "Top 300 Downregulated Peaks (Raw Normalized Counts)")


# Step 1: Sort by most negative log2FoldChange
res_down <- combined_final[order(combined_final$log2FoldChange), ]

# Step 2: Take top 300 most downregulated peaks
top300_down <- head(res_down, 300)

# Step 3: Extract the 4 normalized columns
sample_cols <- c("A27_13_1", "A27_13_2", "A27_Cis13_1", "A27_Cis13_2")
mat_down <- as.matrix(top300_down[, sample_cols])

# Step 4: Row-wise Z-score normalization (for better heatmap contrast)
mat_down_scaled <- t(scale(t(mat_down)))

# Step 5: Optional column annotation
annotation_col <- data.frame(
  Condition = c("A27_13", "A27_13", "A27_Cis13", "A27_Cis13")
)
rownames(annotation_col) <- sample_cols

# Step 6: Plot heatmap
options(repr.plot.width = 10, repr.plot.height = 10)
pheatmap(mat_down_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Top 300 Downregulated Loci (Standardized Normalized Counts)")



print("Working on RLD and VST transformations:")
head(res, 2)
tail(res, 2)



print("RLD and VST transformations")

# Effects of transformations on the variance
rld <- rlog(dds, blind = FALSE)  
vsd <- vst(dds, blind = FALSE) 
ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(rld))
# meanSdPlot(assay(vsd))

print("RLD and VST transformations")

library(DESeq2)
library(vsn)
library(gridExtra)
library(grid)

# Transformations
rld <- rlog(dds, blind = FALSE)
vsd <- vst(dds, blind = FALSE)
ntd <- normTransform(dds)

# Capture the grid plots into objects
p1 <- grid.grabExpr(meanSdPlot(assay(ntd), main = "NTD"))
p2 <- grid.grabExpr(meanSdPlot(assay(rld), main = "RLD"))
p3 <- grid.grabExpr(meanSdPlot(assay(vsd), main = "VST"))

# Set plot size for notebook
options(repr.plot.width = 12, repr.plot.height = 4)

# Arrange side by side
grid.arrange(p1, p2, p3, ncol = 3)




library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# ==== 1. rlog transformation and PCA ====
rld <- rlog(dds, blind = FALSE)
pca_rld <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar_rld <- round(100 * attr(pca_rld, "percentVar"))

pca_rld_plot <- ggplot(pca_rld, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Dark2") +
  xlab(paste0("PC1: ", percentVar_rld[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_rld[2], "% variance")) +
  ggtitle("PCA (rlog)") +
  theme_minimal(base_size = 14)

# ==== 2. rlog MDS ====
dists_rld <- dist(t(assay(rld)))
mds_rld <- cmdscale(as.matrix(dists_rld))
mds_rld_df <- data.frame(MDS1 = mds_rld[, 1], MDS2 = mds_rld[, 2],
                         condition = colData(dds)$condition)

mds_rld_plot <- ggplot(mds_rld_df, aes(MDS1, MDS2, color = condition)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("MDS (rlog)") +
  theme_minimal(base_size = 14)

# ==== 3. vst transformation and PCA ====
vsd <- vst(dds, blind = FALSE)
pca_vsd <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar_vsd <- round(100 * attr(pca_vsd, "percentVar"))

pca_vsd_plot <- ggplot(pca_vsd, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Dark2") +
  xlab(paste0("PC1: ", percentVar_vsd[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_vsd[2], "% variance")) +
  ggtitle("PCA (vst)") +
  theme_minimal(base_size = 14)

# ==== 4. vst MDS ====
dists_vsd <- dist(t(assay(vsd)))
mds_vsd <- cmdscale(as.matrix(dists_vsd))
mds_vsd_df <- data.frame(MDS1 = mds_vsd[, 1], MDS2 = mds_vsd[, 2],
                         condition = colData(dds)$condition)

mds_vsd_plot <- ggplot(mds_vsd_df, aes(MDS1, MDS2, color = condition)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("MDS (vst)") +
  theme_minimal(base_size = 14)

# ==== 5. Combine all four plots in a 2x2 grid ====
options(repr.plot.width = 14, repr.plot.height = 10)

plot_grid(
  pca_rld_plot, mds_rld_plot,
  pca_vsd_plot, mds_vsd_plot,
  labels = c("A", "B", "C", "D"),
  ncol = 2, align = "hv"
)


library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(cowplot)

# === Color Palette ===
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# === rlog Heatmap ===
rlog_matrix <- assay(rld)
sampleDists_rlog <- dist(t(rlog_matrix))
sampleDistMatrix_rlog <- as.matrix(sampleDists_rlog)

p1 <- pheatmap(sampleDistMatrix_rlog,
               clustering_distance_rows = sampleDists_rlog,
               clustering_distance_cols = sampleDists_rlog,
               col = colors,
               fontsize_row = 10,
               fontsize_col = 10,
               cellwidth = 40,
               cellheight = 40,
               angle_col = 45,
               main = "Sample Distance Heatmap (rlog)",
               silent = TRUE)

# === vst Heatmap ===
vsd_matrix <- assay(vsd)
sampleDists_vsd <- dist(t(vsd_matrix))
sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

p2 <- pheatmap(sampleDistMatrix_vsd,
               clustering_distance_rows = sampleDists_vsd,
               clustering_distance_cols = sampleDists_vsd,
               col = colors,
               fontsize_row = 10,
               fontsize_col = 10,
               cellwidth = 40,
               cellheight = 40,
               angle_col = 45,
               main = "Sample Distance Heatmap (vst)",
               silent = TRUE)

# === Plot side-by-side ===
options(repr.plot.width = 16, repr.plot.height = 8)
grob1 <- p1$gtable
grob2 <- p2$gtable

plot_grid(grob1, grob2, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1))

print("Methods to use : GLM-PCA for PCA and PoissonDistance to calculate the sample distances")
# Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package.
# This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating
# the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows
# instead of columns, so we need to transpose the counts in dds.

library(glmpca)
library("PoiClaClu")
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(grid)

print("PCA by using GLMPCA library. RLOG and VSD transformations are more suitable than scale().")

# === GLM-PCA plot ===
gpca <- glmpca(counts(dds), L = 2, fam = "poi")  # Use counts, not assay(dds)

gpca.dat <- as.data.frame(gpca$factors)
gpca.dat$sample <- colnames(dds)
gpca.dat$condition <- colData(dds)$condition

p_gpca <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = condition)) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample), size = 4, box.padding = 0.4, max.overlaps = 8) +
  coord_fixed() +
  theme_minimal(base_size = 16) +
  labs(title = "GLM-PCA", x = "GLM-PC1", y = "GLM-PC2") +
  scale_color_brewer(palette = "Set2") +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

# === Poisson distance heatmap ===
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)

sample_names <- colnames(dds)
rownames(samplePoisDistMatrix) <- sample_names
colnames(samplePoisDistMatrix) <- sample_names

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheat <- pheatmap(samplePoisDistMatrix,
                  clustering_distance_rows = poisd$dd,
                  clustering_distance_cols = poisd$dd,
                  col = colors,
                  fontsize_row = 10,
                  fontsize_col = 10,
                  cellwidth = 50,
                  cellheight = 50,
                  angle_col = 45,
                  main = "Poisson Distance Heatmap",
                  silent = TRUE)

# Convert to grob for use with cowplot
g_poisson <- ggdraw(grobTree(pheat$gtable)) + theme(plot.margin = margin(5, 5, 5, 5))

# Set display size
options(repr.plot.width = 16, repr.plot.height = 8)

# Combine plots side-by-side
combined_plot <- plot_grid(
  g_poisson, p_gpca,
  labels = c("A", "B"),
  label_size = 16,
  nrow = 1,
  rel_widths = c(1.2, 1)
)

# Show plot
print(combined_plot)









print("Performing Surrogate Variable Analysis")
print("SVA analysis")
# SV1, SV2, ... are surrogate variables — latent (hidden) factors estimated from the data that capture unwanted variation 
# (like batch effects, technical noise, or hidden biological subtypes).
# You can think of them as "virtual covariates" — constructed purely from the structure of your data — 
# that explain sources of variation not included in your model (like treatment or condition).
table(dds$condition)

# Load necessary libraries
library(DESeq2)
library(sva)

# === Step 1: Normalize counts and extract matrix ===
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)

# === Step 2: Filter low-expressed genes ===
dat <- dat[rowMeans(dat) > 1, ]

# === Step 3: Remove rows with NA or zero variance ===
dat <- dat[complete.cases(dat), ]
dat <- dat[apply(dat, 1, var) > 1e-6, ]  # Remove invariant genes

# === Step 4: Check 'condition' variable ===
if (!"condition" %in% colnames(colData(dds))) {
  stop("The 'condition' column does not exist in colData(dds).")
}
dds$condition <- as.factor(dds$condition)
print(table(dds$condition))

# === Step 5: Create model matrices ===
mod <- model.matrix(~ 0 + condition, data = colData(dds))
colnames(mod) <- levels(dds$condition)

mod0 <- model.matrix(~ 1, data = colData(dds))

# === Step 6: Check model sanity ===
if (any(is.na(dat))) stop("Data matrix still contains NA values.")
if (any(is.na(mod)) || any(is.na(mod0))) stop("Model matrices contain NA values.")
if (qr(mod)$rank != ncol(mod)) stop("Design matrix 'mod' is not full rank.")

# === Step 7: Estimate number of surrogate variables (optional) ===
# n.sv <- num.sv(dat, mod, method = "leek")  # Alternative method
n.sv <- 1  # If known/estimated manually

# === Step 8: Run SVA ===
svseq <- sva(dat, mod, mod0, n.sv = n.sv)

# === Step 9: Check surrogate variable matrix ===
print(dim(svseq$sv))
print(head(svseq$sv))

# === Step 10: Visualize surrogate variables ===
par(mfrow = c(1, ncol(svseq$sv)))
for (i in 1:ncol(svseq$sv)) {
  stripchart(
    svseq$sv[, i] ~ dds$condition,
    vertical = TRUE,
    method = "jitter",
    pch = 21,
    bg = "steelblue",
    col = "black",
    frame.plot = FALSE,
    ylim = range(svseq$sv[, i]) * 1.2,
    main = paste0("SV", i),
    ylab = "Surrogate Variable",
    xlab = "Condition"
  )
  abline(h = 0, col = "gray50", lty = 2)
}

# === Step 11: Add surrogate variables to colData(dds) ===
for (i in 1:ncol(svseq$sv)) {
  colData(dds)[[paste0("SV", i)]] <- svseq$sv[, i]
}

# Clone original DESeqDataSet object
ddssva <- dds

# Add surrogate variables SV1 and SV2 from svseq
ddssva$SV1 <- svseq$sv[, 1]
# ddssva$SV2 <- svseq$sv[, 2]

# Update the design formula to include surrogate variables
design(ddssva) <- ~ SV1 + condition
# design(ddssva) <- ~ SV1 + SV2 + condition

# Run DESeq2 with SVA-adjusted design
ddssva <- DESeq(ddssva)

# Check available results
resultsNames(ddssva)

# === Get results for comparisons between conditions ===
# Update these to match your actual levels in `dds$condition`

# View your levels to confirm names:
levels(ddssva$condition)

# Replace below with your real condition names
# Example: "A27_Cis13" vs "A27_13"
res_ddsva_Cis13_vs_13 <- results(ddssva, contrast = c("condition", "A27_Cis13", "A27_13"))
res_ddsva_Cis13_vs_13_df = data.frame(res_ddsva_Cis13_vs_13)

# Save the results
# write.csv(res_ddsva_Cis13_vs_13_df, file = "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.csv")

# combining res_ddsva_Cis13_vs_13_df with peaks and normalized counts
head(res_ddsva_Cis13_vs_13_df, 2)
res_ddsva_Cis13_vs_13_df$peak_id <- rownames(res_ddsva_Cis13_vs_13_df)

head(peaks, 2)

# Merge with peaks by 'name' and 'peak_id'
combined2 <- merge(peaks, res_ddsva_Cis13_vs_13_df, by.x = "name", by.y = "peak_id")

# Preview the merged data
head(combined, 2)

# Write the merged results to CSV
write.csv(combined2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.csv", row.names = FALSE)

# Combine with normalized counts: convert normalized counts to a data frame and add peak ID

norm_counts2 = counts(ddssva, normalized = TRUE)
norm_counts2_df <- as.data.frame(norm_counts2)
norm_counts2_df$peak_id <- rownames(norm_counts2_df)

# Merge normalized counts with 'combined' data
combined_final2 <- merge(combined2, norm_counts2_df, by.x = "name", by.y = "peak_id")

# Preview the result
cat("Preview of final combined table:\n")
print(colnames(combined_final2), 2)
print(head(combined_final2, 2))

# Write the final combined results with normalized counts to CSV
write.csv(combined_final2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.and.norm.csv", row.names = FALSE)

# === Summary of DE results ===
cat("Differential binding: A27_Cis13 vs A27_13\n")
cat("Genes with p-value < 0.05:", nrow(subset(res_ddsva_Cis13_vs_13, pvalue < 0.05)), "\n")
cat("Genes with padj < 0.1:", nrow(subset(res_ddsva_Cis13_vs_13, padj < 0.1)), "\n")

# Clean up: remove rows with NA pvalue or padj
# res_clean2 <- res_ddsva_Cis13_vs_13[!is.na(res_ddsva_Cis13_vs_13$pvalue) & !is.na(res_ddsva_Cis13_vs_13$padj), ]
# res_clean2_df = data.frame(res_clean2)
# res_clean2_df$peak_id = rownames(res_clean2_df)
# head(res_clean2_df, 2)

# Clean up: remove rows with NA pvalue or padj
cat("before cleaning the dataframe\n")
print(nrow(combined_final2))
combined_final2 <- combined_final2[!is.na(combined_final2$pvalue) & !is.na(combined_final2$padj), ]
cat("after cleaning the dataframe\n")
print(nrow(combined_final2))

# Total number of peaks
cat("Total peaks analyzed (non-NA):", nrow(combined_final2), "\n\n")

# === UPREGULATED ===
cat("UPREGULATED: log2FoldChange > 1\n")
cat("  pvalue < 0.05 :", sum(combined_final2$log2FoldChange > 1 & combined_final2$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(combined_final2$log2FoldChange > 1 & combined_final2$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(combined_final2$log2FoldChange > 1 & combined_final2$padj < 0.1), "\n\n")

# === DOWNREGULATED ===
cat("DOWNREGULATED: log2FoldChange < -1\n")
cat("  pvalue < 0.05 :", sum(combined_final2$log2FoldChange < -1 & combined_final2$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(combined_final2$log2FoldChange < -1 & combined_final2$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(combined_final2$log2FoldChange < -1 & combined_final2$padj < 0.1), "\n")

# === Filtered tables ===

# UPREGULATED
cat("\n")
upregulated2 <- combined_final2[combined_final2$padj < 0.1 & combined_final2$log2FoldChange > 0.5, ]
cat("Number of upregulated peaks: padj < 0.1 & log2FoldChange > 0.5 ", nrow(upregulated2), "\n")
# write.csv(upregulated2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.upregulated.padj0.1_log2FC.plus.0.5.csv", 
#          row.names = FALSE)

# DOWNREGULATED
cat("\n")
downregulated2 <- combined_final2[combined_final2$padj < 0.1 & combined_final2$log2FoldChange < -0.5, ]
cat("Number of downregulated peaks: padj < 0.1 & log2FoldChange < -0.5 ", nrow(downregulated2), "\n")
# write.csv(downregulated2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.downregulated.padj0.1_log2FC.minus.0.5.csv", 
#          row.names = FALSE)

# === Optional BED output ===
# Only if combined_final2 has columns: chr, start, end, name

cat("\nNumber of up-regulated peaks:", nrow(upregulated2), "\n")
head(upregulated2, 2)
cat("\nNumber of down-regulated peaks:", nrow(downregulated2), "\n")
head(downregulated2, 2)

bed_up2 <- upregulated2[, c("chr", "start", "end", "name")]
# write.table(bed_up2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.upregulated.padj0.1_log2FC.plus.0.5.bed",
#              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed_down2 <- downregulated2[, c("chr", "start", "end", "name")]
# write.table(bed_down2, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.SVA.downregulated.padj0.1_log2FC.minus.0.5.bed",
#              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Transform count data
vsd2 <- vst(ddssva, blind = TRUE)
rld2 <- rlog(ddssva, blind = TRUE)

# Get PCA data
pca_vsd <- plotPCA(vsd2, intgroup = "condition", returnData = TRUE)
pca_rld <- plotPCA(rld2, intgroup = "condition", returnData = TRUE)

# Variance explained
percentVar_vsd <- round(100 * attr(pca_vsd, "percentVar"))
percentVar_rld <- round(100 * attr(pca_rld, "percentVar"))

# PCA plot for VST
p1 <- ggplot(pca_vsd, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA after SVA (VST)",
    x = paste0("PC1 (", percentVar_vsd[1], "%)"),
    y = paste0("PC2 (", percentVar_vsd[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "right")

# PCA plot for RLOG
p2 <- ggplot(pca_rld, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA after SVA (RLOG)",
    x = paste0("PC1 (", percentVar_rld[1], "%)"),
    y = paste0("PC2 (", percentVar_rld[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "right")

# Show both plots side by side with legends
options(repr.plot.width = 14, repr.plot.height = 6)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)























print("RUVseq analysis")

# Load required libraries
library(RUVSeq)
library(DESeq2)

# === Step 1: Create SeqExpressionSet from DESeq2 object ===
set <- newSeqExpressionSet(as.matrix(counts(dds, normalized = FALSE)),
                           phenoData = data.frame(condition = dds$condition, 
                                                  row.names = colnames(dds)))

# === Step 2: Filter low-expressed peaks ===
keep <- rowSums(counts(set) > 4) >= 4
set <- set[keep, ]

# === Step 3: Normalize between lanes ===
set <- betweenLaneNormalization(set, which = "upper")

# === Step 4: Estimate empirical control genes (non-DE) ===
dds_temp <- dds[keep, ]
dds_temp <- DESeq(dds_temp)
res_temp <- results(dds_temp)

# Define empirical control genes as those with high p-value (non-DE)
not.sig <- rownames(res_temp)[which(res_temp$pvalue > 0.1)]
empirical <- rownames(set)[rownames(set) %in% not.sig]

# Apply RUVg with k = 2 unwanted factors
set <- RUVg(set, empirical, k = 2)

# === Step 5: Add unwanted factors to DESeq2 object ===
ddsruv <- dds[keep, ]
ddsruv$W1 <- pData(set)$W_1
# ddsruv$W2 <- pData(set)$W_2

# Specify model design
design(ddsruv) <- ~ W1 + condition
# design(ddsruv) <- ~ W1 + W2 + condition

# === Step 6: Run DESeq2 with RUV-adjusted model ===
ddsruv <- DESeq(ddsruv)

# === Step 7: Inspect model variables ===
resultsNames(ddsruv)

# Set up single plot
par(
  mfrow = c(1, 1),
  mar = c(4, 4, 3, 1),   # margins: bottom, left, top, right
  cex.main = 1.2,        # title size
  cex.axis = 1.0,        # axis tick label size
  cex.lab = 1.1,         # axis title size
  las = 1                # horizontal y-axis labels
)

# Plot W_1 grouped by condition
stripchart(
  pData(set)$W_1 ~ pData(set)$condition,
  vertical = TRUE,
  method = "jitter",
  pch = 21,
  bg = "steelblue",
  col = "black",
  frame.plot = FALSE,
  main = "RUV Factor W1",
  ylab = "W1 Value",
  xlab = "Condition",
  cex = 1.1
)
abline(h = 0, lty = 2, col = "gray60", lwd = 1)

# View your levels to confirm names:
levels(ddsruv$condition)

# Replace below with your real condition names
# Example: "A27_Cis13" vs "A27_13"
res_ddsruv_Cis13_vs_13 <- results(ddsruv, contrast = c("condition", "A27_Cis13", "A27_13"))
res_ddsruv_Cis13_vs_13

# Save the results
write.csv(as.data.frame(res_ddsruv_Cis13_vs_13), file = "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.RUV.csv")

# === Summary of DE results ===
cat("Differential binding (RUV-adjusted): A27_Cis13 vs A27_13\n")
cat("Genes with p-value < 0.05:", nrow(subset(res_ddsruv_Cis13_vs_13, pvalue < 0.05)), "\n")
cat("Genes with padj < 0.1:", nrow(subset(res_ddsruv_Cis13_vs_13, padj < 0.1)), "\n")

# Clean up: remove rows with NA pvalue or padj
res_clean <- res_ddsruv_Cis13_vs_13[!is.na(res_ddsruv_Cis13_vs_13$pvalue) & !is.na(res_ddsruv_Cis13_vs_13$padj), ]

# Total number of peaks
cat("Total peaks analyzed (non-NA):", nrow(res_clean), "\n\n")

# === UPREGULATED ===
cat("UPREGULATED: log2FoldChange > 1\n")
cat("  pvalue < 0.05 :", sum(res_clean$log2FoldChange > 1 & res_clean$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(res_clean$log2FoldChange > 1 & res_clean$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(res_clean$log2FoldChange > 1 & res_clean$padj < 0.1), "\n\n")

# === DOWNREGULATED ===
cat("DOWNREGULATED: log2FoldChange < -1\n")
cat("  pvalue < 0.05 :", sum(res_clean$log2FoldChange < -1 & res_clean$pvalue < 0.05), "\n")
cat("  padj   < 0.05 :", sum(res_clean$log2FoldChange < -1 & res_clean$padj < 0.05), "\n")
cat("  padj   < 0.1  :", sum(res_clean$log2FoldChange < -1 & res_clean$padj < 0.1), "\n")

# === Filtered tables ===
cat("\n")
upregulated <- res_clean[res_clean$padj < 0.1 & res_clean$log2FoldChange > 0.5, ]
cat("Number of upregulated peaks: padj < 0.1 & log2FC > 0.5 :", nrow(upregulated), "\n")
write.csv(upregulated, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.RUV.upregulated.padj0.1_log2FC.plus.0.5.csv", 
          row.names = FALSE)

cat("\n")
downregulated <- res_clean[res_clean$padj < 0.1 & res_clean$log2FoldChange < -0.5, ]
cat("Number of downregulated peaks: padj < 0.1 & log2FC < -0.5 :", nrow(downregulated), "\n")
write.csv(downregulated, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.RUV.downregulated.padj0.1_log2FC.minus.0.5.csv", 
          row.names = FALSE)

# === Optional BED output ===
if (all(c("chr", "start", "end", "name") %in% colnames(res_clean))) {
  
  bed_up <- upregulated[, c("chr", "start", "end", "name")]
  write.table(bed_up, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.RUV.upregulated.padj0.1_log2FC.plus.0.5.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  bed_down <- downregulated[, c("chr", "start", "end", "name")]
  write.table(bed_down, "the.matrix.peaks_and_counts.method1.combined_peak_DESeq2_results.RUV.downregulated.padj0.1_log2FC.minus.0.5.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Preview
cat("\n")
head(upregulated, 2)
head(downregulated, 2)


# Transform count data from ddsruv
vsd3 <- vst(ddsruv, blind = TRUE)
rld3 <- rlog(ddsruv, blind = TRUE)

# Get PCA data
pca_vsd <- plotPCA(vsd3, intgroup = "condition", returnData = TRUE)
pca_rld <- plotPCA(rld3, intgroup = "condition", returnData = TRUE)

# Variance explained
percentVar_vsd <- round(100 * attr(pca_vsd, "percentVar"))
percentVar_rld <- round(100 * attr(pca_rld, "percentVar"))

# PCA plot for VST
p1 <- ggplot(pca_vsd, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA after RUV (VST)",
    x = paste0("PC1 (", percentVar_vsd[1], "%)"),
    y = paste0("PC2 (", percentVar_vsd[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "right")

# PCA plot for RLOG
p2 <- ggplot(pca_rld, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA after RUV (RLOG)",
    x = paste0("PC1 (", percentVar_rld[1], "%)"),
    y = paste0("PC2 (", percentVar_rld[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "right")

# Show both plots side by side with legends
options(repr.plot.width = 14, repr.plot.height = 6)
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)


