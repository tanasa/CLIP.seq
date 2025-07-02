####################################################
####################################################

library(eulerr)
library(rtracklayer)
library(GenomicRanges)
library(UpSetR)
library(VennDiagram)

####################################################
####################################################

# Read each peak file separately
peaks_13_1 <- import("16h_A27_13_1.pval.0.01.narrowPeaks", format = "narrowPeak")
peaks_13_2 <- import("16h_A27_13_2.pval.0.01.narrowPeaks", format = "narrowPeak")
peaks_Cis13_1 <- import("16h_A27_Cis13_1.pval.0.01.narrowPeaks", format = "narrowPeak")
peaks_Cis13_2 <- import("16h_A27_Cis13_2.pval.0.01.narrowPeaks", format = "narrowPeak")

cat("peaks_13_1:", length(peaks_13_1), "\n")
cat("peaks_13_2:", length(peaks_13_2), "\n")
cat("peaks_Cis13_1:", length(peaks_Cis13_1), "\n")
cat("peaks_Cis13_2:", length(peaks_Cis13_2), "\n")

####################################################
####################################################

# Define allowed chromosomes: chr1 to chr22 + chrX, chrY, chrM
allowed_chroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# Filter each peak set
peaks_13_1 <- subset(peaks_13_1, as.character(seqnames(peaks_13_1)) %in% allowed_chroms)
peaks_13_2 <- subset(peaks_13_2, as.character(seqnames(peaks_13_2)) %in% allowed_chroms)
peaks_Cis13_1 <- subset(peaks_Cis13_1, as.character(seqnames(peaks_Cis13_1)) %in% allowed_chroms)
peaks_Cis13_2 <- subset(peaks_Cis13_2, as.character(seqnames(peaks_Cis13_2)) %in% allowed_chroms)

cat("peaks_13_1:", length(peaks_13_1), "\n")
cat("peaks_13_2:", length(peaks_13_2), "\n")
cat("peaks_Cis13_1:", length(peaks_Cis13_1), "\n")
cat("peaks_Cis13_2:", length(peaks_Cis13_2), "\n")

####################################################
####################################################

# 🔁 Intersect: 13_1 with 13_2

hits_13 <- findOverlaps(peaks_13_1, peaks_13_2)
intersect_13 <- pintersect(peaks_13_1[queryHits(hits_13)],
                           peaks_13_2[subjectHits(hits_13)])
length(intersect_13)  # number of overlapping regions

# 🔁  Intersect: Cis13_1 with Cis13_2

hits_Cis13 <- findOverlaps(peaks_Cis13_1, peaks_Cis13_2)
intersect_Cis13 <- pintersect(peaks_Cis13_1[queryHits(hits_Cis13)],
                              peaks_Cis13_2[subjectHits(hits_Cis13)])
length(intersect_Cis13)  # number of overlapping regions

#################################################### export the results
####################################################

export(intersect_13, "16h_peaks_intersect_replicates_13_1_13_2.bed", format = "BED")
export(intersect_Cis13, "16h_peaks_intersect_replicates_Cis13_1_Cis13_2.bed", format = "BED")

####################################################
#################################################### compute multi-way intersection

# Put all peak sets in a list
peak_sets <- list(peaks_13_1, peaks_13_2, peaks_Cis13_1, peaks_Cis13_2)

# Use Reduce to find overlapping ranges step by step
intersect_all <- Reduce(function(x, y) {
  hits <- findOverlaps(x, y)
  pintersect(x[queryHits(hits)], y[subjectHits(hits)])
}, peak_sets)

# Number of shared peaks
length(intersect_all)

# Show number of shared peaks
cat("Number of peaks in the 4-way intersection:", length(intersect_all), "\n")

# Display first few intersected peaks on screen
print(intersect_all[1:min(5, length(intersect_all))])

# Save the result to BED
export(intersect_all, "16h_peaks_intersection_all_4_sets.bed", format = "BED")
cat("Saved to: 16h_peaks_intersection_all_4_sets.bed\n")

#################################################### 
#################################################### library(VennDiagram)

# Lengths
n_13_1 <- length(peaks_13_1)
n_13_2 <- length(peaks_13_2)
n_Cis13_1 <- length(peaks_Cis13_1)
n_Cis13_2 <- length(peaks_Cis13_2)

# Overlaps
overlap_13 <- length(unique(queryHits(findOverlaps(peaks_13_1, peaks_13_2))))
overlap_Cis13 <- length(unique(queryHits(findOverlaps(peaks_Cis13_1, peaks_Cis13_2))))

venn.diagram(
  x = list(`13_1` = 1:n_13_1, `13_2` = (n_13_1 - overlap_13 + 1):(n_13_1 + n_13_2 - overlap_13)),
  category.names = c("13_1", "13_2"),
  filename = "16h_peaks_venn_13.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  col = "black",
  fill = c("skyblue", "lightpink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

venn.diagram(
  x = list(`Cis13_1` = 1:n_Cis13_1, `Cis13_2` = (n_Cis13_1 - overlap_Cis13 + 1):(n_Cis13_1 + n_Cis13_2 - overlap_Cis13)),
  category.names = c("Cis13_1", "Cis13_2"),
  filename = "16h_peaks_venn_Cis13.png",
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 500,
  col = "black",
  fill = c("lightgreen", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)

######################################################
###################################################### library(eulerr)

# 13 peaks
plot(euler(c("13_1" = n_13_1, "13_2" = n_13_2, "13_1&13_2" = overlap_13)), 
     fills = c("skyblue", "lightpink"), main = "13 Overlap")

# Cis13 peaks
plot(euler(c("Cis13_1" = n_Cis13_1, "Cis13_2" = n_Cis13_2, "Cis13_1&Cis13_2" = overlap_Cis13)), 
     fills = c("lightgreen", "orange"), main = "Cis13 Overlap")

######################################################
######################################################

fit <- euler(c(
  "13_1" = length(peaks_13_1),
  "13_2" = length(peaks_13_2),
  "Cis13_1" = length(peaks_Cis13_1),
  "Cis13_2" = length(peaks_Cis13_2),
  
  "13_1&13_2" = length(findOverlaps(peaks_13_1, peaks_13_2)),
  "13_1&Cis13_1" = length(findOverlaps(peaks_13_1, peaks_Cis13_1)),
  "13_1&Cis13_2" = length(findOverlaps(peaks_13_1, peaks_Cis13_2)),
  "13_2&Cis13_1" = length(findOverlaps(peaks_13_2, peaks_Cis13_1)),
  "13_2&Cis13_2" = length(findOverlaps(peaks_13_2, peaks_Cis13_2)),
  "Cis13_1&Cis13_2" = length(findOverlaps(peaks_Cis13_1, peaks_Cis13_2)),

  "13_1&13_2&Cis13_1&Cis13_2" = length(intersect_all)  # 4-way overlap
))

plot(fit,
     fills = list(fill = c("skyblue", "lightpink", "lightgreen", "orange"), alpha = 0.6),
     quantities = TRUE,
     main = "Euler Diagram of Peak Set Intersections")


png("16h_peaks_4intersection_peak_set_euler_diagram.png", width = 1600, height = 1600, res = 300)
plot(fit,
     fills = list(fill = c("skyblue", "lightpink", "lightgreen", "orange"), alpha = 0.6),
     quantities = TRUE,
     main = "Euler Diagram of Peak Set Intersections")
dev.off()

######################################################
###################################################### display UpsetR plot

# Combine and reduce all peaks into one universe
all_peaks <- reduce(c(peaks_13_1, peaks_13_2, peaks_Cis13_1, peaks_Cis13_2))

# Build membership matrix (TRUE/FALSE for each peak's presence in each set)
membership_matrix <- data.frame(
  "13_1" = overlapsAny(all_peaks, peaks_13_1),
  "13_2" = overlapsAny(all_peaks, peaks_13_2),
  "Cis13_1" = overlapsAny(all_peaks, peaks_Cis13_1),
  "Cis13_2" = overlapsAny(all_peaks, peaks_Cis13_2)
)

# Convert to 0/1
membership_matrix[] <- lapply(membership_matrix, as.integer)

# Optional: remove rows with no membership (i.e., 0000)
membership_matrix <- membership_matrix[rowSums(membership_matrix) > 0, ]

# Update the sets argument to match column names

upset(
  membership_matrix,
  sets = c("X13_1", "X13_2", "Cis13_1", "Cis13_2"),
  order.by = "freq",
  mainbar.y.label = "Peak Set Intersections",
  sets.x.label = "Peaks per Set"
)

# A cleaner version if you want to stick to original set names:

colnames(membership_matrix) <- sub("^X", "", colnames(membership_matrix))  # remove leading "X"

# Now the sets match original names
png("16h_peaks_intersection_peak_set_upset_plot.png", width = 1600, height = 1200, res = 200)

upset(
  membership_matrix,
  sets = c("13_1", "13_2", "Cis13_1", "Cis13_2"),
  order.by = "freq",
  mainbar.y.label = "Peak Set Intersections",
  sets.x.label = "Peaks per Set"
)

dev.off()

######################################################
######################################################
######################################################