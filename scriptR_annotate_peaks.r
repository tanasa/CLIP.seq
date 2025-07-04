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

# Peak Annotation

peaks_signal200 <- subset(peak_data, signalValue >= 200)
peaks_signal300 <- subset(peak_data, signalValue >= 300)
peaks_signal400 <- subset(peak_data, signalValue >= 400)

# Check the number of peaks in each subset
cat("Peaks with signal ≥ 200:", nrow(peaks_signal200), "\n")
cat("Peaks with signal ≥ 300:", nrow(peaks_signal300), "\n")
cat("Peaks with signal ≥ 400:", nrow(peaks_signal400), "\n")


peaks = peaks_signal400 
head(peaks_signal400)
dim(peaks_signal400)

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



library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) # For gene annotation 
# c(-1000, 100)	Core promoter	Strict promoter-only binding

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Define the TSS region (e.g., -3kb to +3kb around the TSS)
# Adjust these values based on your biological question

# tssRegion <- c(-3000, 100)
tssRegion <- c(-3000, 0)

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

# Annotate peaks
tssRegion <- c(-3000, 0)

# Custom order (example: prioritize Exon over UTRs)
my_priority <- c("Promoter", "Exon", "5UTR", "3UTR", "Intron", "Downstream", "Intergenic")

peak_anno <- annotatePeak(
    peak = peaks_gr,
    tssRegion = tssRegion,
    TxDb = txdb,
    # Use OrgDb for gene symbol mapping if desired
    annoDb = "org.Hs.eg.db",
    genomicAnnotationPriority = my_priority,
    verbose = TRUE
  )

print("the annotation of the peaks:")
plotAnnoPie(peak_anno)

# write the peak_annotation into an external file
write.csv(peak_anno, "the.matrix.peaks_and_counts.method1.peak.annotations.ChIPSeekR.csv", row.names = FALSE)

# print("Promoters regions located within 1 kilobase (kb) of a Transcription Start Site (TSS)")
# "Promoter (1-2kb)" means peaks falling between 1kb and 2kb from the TSS.
# "Promoter (2-3kb)" means peaks falling between 2kb and 3kb from the TSS.

# Load libraries
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Define genomic feature regions
promoters_gr <- promoters(txdb, upstream = 3000, downstream = 0)
fiveUTR_gr   <- unlist(fiveUTRsByTranscript(txdb, use.names = TRUE))
threeUTR_gr  <- unlist(threeUTRsByTranscript(txdb, use.names = TRUE))
exons_gr     <- unlist(exonsBy(txdb, by = "tx"))
introns_gr   <- unlist(intronsByTranscript(txdb))

# Initialize annotation vector
categories <- rep("", length(peaks_gr))

# Helper function to assign categories with priority
assign_category <- function(query, target_gr, label) {
  hits <- findOverlaps(query, target_gr)
  idx <- queryHits(hits)
  categories[idx[categories[idx] == ""]] <<- label
}

# Assign annotations by priority
assign_category(peaks_gr, promoters_gr, "promoter")
assign_category(peaks_gr, fiveUTR_gr, "5UTR")
assign_category(peaks_gr, threeUTR_gr, "3UTR")
assign_category(peaks_gr, exons_gr, "exon")
assign_category(peaks_gr, introns_gr, "intron")
categories[categories == ""] <- "intergenic"

# Create final annotation data frame
annotated_df <- data.frame(
  peak_name = if (!is.null(mcols(peaks_gr)$name)) mcols(peaks_gr)$name else paste0("peak_", seq_along(peaks_gr)),
  annotation = categories
)

# Print annotation counts
print(table(annotated_df$annotation))

# Create summary for pie chart
annotation_summary <- as.data.frame(table(annotated_df$annotation))
colnames(annotation_summary) <- c("Region", "Count")

# Plot pie chart
ggplot(annotation_summary, aes(x = "", y = Count, fill = Region)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Genomic Annotation of Peaks") +
  theme_void() +
  theme(legend.title = element_blank())


library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)

library(EnsDb.Hsapiens.v86)
ensembl.hs86.transcripts <- transcripts(EnsDb.Hsapiens.v86)

# 1. Convert gene features to GRanges
# annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, feature = "gene")
# 2. Annotate peaks
# annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
# head(annotatedPeaks, 2)

# Convert your PEAKS to NCBI style (not the database)
seqlevelsStyle(peaks_gr) <- "NCBI"

# Remove problematic sequences
peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")

peaks_ensembl <- annotatePeakInBatch(peaks_gr, 
                                     AnnotationData = ensembl.hs86.transcripts)

write.table(as.data.frame(peaks_ensembl), 
            file = "the.matrix.peaks_and_counts.method1.peak.annotations.ChIPpeakAnno.csv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

head(peaks_ensembl, n = 2)

# https://bioconductor.org/packages/devel/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html

table(peaks_ensembl$insideFeature) 

# library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# annotation_data <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
# binOverFeature(peaks_gr, 
#               featureSite = c("FeatureStart", "FeatureEnd", "bothEnd"),
#               nbins = 20,
#               annotationData = annotation_data,
#              xlab = "peak distance from TSS (bp)", 
#              ylab = "peak count", 
#               main = "Distribution of aggregated peak numbers around TSS")

chromosome_region <- assignChromosomeRegion(peaks_gr,
                                            TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                            nucleotideLevel = FALSE,
                                            precedence=c("Promoters",
                                                         "immediateDownstream", 
                                                         "fiveUTRs", 
                                                         "threeUTRs",
                                                         "Exons", 
                                                         "Introns"))

# optional helper function to rotate x-axis labels for barplot(): 
# ref: https://stackoverflow.com/questions/10286473/rotating-x-axis-labels-in-r-for-barplot

rotate_x <- function(data, rot_angle) {
  plt <- barplot(data, xaxt = "n")
  text(plt, par("usr")[3], 
       labels = names(data), 
       srt = rot_angle, adj = c(1.1,1.1), 
       xpd = TRUE, cex = 0.6)
}

rotate_x(chromosome_region[["percentage"]], 45)

genomicElementDistribution(peaks_gr, 
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                           promoterRegion = c(upstream=3000, downstream=0)
                           )

# library(UpSetR)
# library(repr)
# library(reactome.db)

# res <- genomicElementUpSetR(peaks_gr,
#                            TxDb.Hsapiens.UCSC.hg38.knownGene)

# options(repr.plot.width = 18, repr.plot.height = 12)

# upset(res[["plotData"]], 
#      nsets = length(colnames(res$plotData)), 
#      nintersects = NA)

# anno_data <- genes(EnsDb.Hsapiens.v86)
# annotated_peaks <- annotatePeakInBatch(peaks_gr,
#                                       AnnotationData = anno_data,
#                                       output = "both")

# enriched_go <- getEnrichedGO(annotated_peaks, 
#                             orgAnn = "org.Hs.eg.db", 
#                             feature_id_type = "ensembl_gene_id")

# enriched_path <- getEnrichedPATH(annotated_peaks,
#                                 orgAnn = "org.Hs.eg.db",
#                                 feature_id_type = "ensembl_gene_id",
#                                 # pathAnn = "reactome.db",
#                                 pathAnn= "KEGGREST")

# enrichmentPlot(enriched_go)
# enrichmentPlot(enriched_path)

# ensembl.hs86.transcript <- transcripts(EnsDb.Hsapiens.v86)
# peaks_ensembl <- annotatePeakInBatch(peaks_gr, 
#                                    AnnotationData = ensembl.hs86.transcript, 
#                                    output = "both")

# head(peaks_ensembl, n = 2)



# library(annotatr)

# Load built-in annotations
# annotations <- c("hg38_genes_promoters", 
#                 "hg38_genes_5UTRs", 
#                 "hg38_genes_3UTRs", 
#                 "hg38_genes_exons", 
#                 "hg38_genes_introns")

# annot_regions <- build_annotations(genome = 'hg38', annotations = annotations)
# cat("Using EnsDb.Hsapiens.v86 directly")

# Extract different genomic features from EnsDb
# edb <- EnsDb.Hsapiens.v86

# Get transcripts, exons, and other features
# transcripts_ensdb <- transcripts(edb)
# exons_ensdb <- exons(edb)
# genes_ensdb <- genes(edb)



library(EnsDb.Hsapiens.v86)

edb <- EnsDb.Hsapiens.v86 # Get your EnsDb object

# You can directly provide the EnsDb object to ChIPseeker
# ChIPseeker's annotatePeak function is designed to handle this.
# Ensure your peaks GRanges object (e.g., macs_peak_gr) is loaded.

# Annotate peaks using the EnsDb object
peak_annotation_cs <- annotatePeak(
  peak = peaks_gr,
  tssRegion = c(-3000, 0), # Define promoter region
  TxDb = edb,                 # Provide the EnsDb object here
  annoDb = "org.Hs.eg.db",    # For gene symbols, if installed
  verbose = FALSE
)

# View summary and convert to data.frame
print(peak_annotation_cs)

annotated_df_cs <- as.data.frame(peak_annotation_cs)
head(annotated_df_cs, 2)

write.table(as.data.frame(peaks_ensembl), 
            file = "the.matrix.peaks_and_counts.method1.peak.annotations.ChIPSeekR.ensdb.csv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)



# Annotate with a gtf file

gtf_file <- "gencode.v48.basic.annotation.gtf"
txdb_gtf <- makeTxDbFromGFF(gtf_file, format = "gtf")

promoters_gr <- promoters(txdb_gtf, upstream = 3000, downstream = 0)
genes_gr     <- genes(txdb_gtf)

exons_gr     <- exons(txdb_gtf)
introns_gr   <- unlist(intronsByTranscript(txdb_gtf))

# 5' UTRs
fiveUTR_gr <- unlist(fiveUTRsByTranscript(txdb_gtf, use.names = TRUE))

# 3' UTRs
threeUTR_gr <- unlist(threeUTRsByTranscript(txdb_gtf, use.names = TRUE))

# Annotations

annoData <- genes_gr  # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))

annoData <- promoters_gr  # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))

annoData <- exons_gr   # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))

annoData <- introns_gr  # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))

annoData <- fiveUTR_gr  # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))

annoData <- threeUTR_gr  # or exons_gr, introns_gr, etc.
annotatedPeaks <- annotatePeakInBatch(peaks_gr, AnnotationData = annoData)
head(data.frame(annotatedPeaks, 2))





# other packages to consider :

# regioneR
# https://www.bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.html

# annotatr
# https://bioconductor.org/packages/release/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html
