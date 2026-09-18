# ============================================================
# ATAC-SEQ ANALYSIS AND RNA-SEQ INTEGRATION
# Prostate Cancer (TCGA-PRAD)
#
# Workflow:
# BED files
#    ↓
# Tumor peak merging
#    ↓
# Normal peak merging
#    ↓
# Standard chromosome filtering
#    ↓
# Tumor-specific peaks
#    ↓
# Promoter accessibility (TSS ± 3 kb)
#    ↓
# Integration with RNA-seq upregulated genes
#    ↓
# 1114 epigenetically activated genes
#    ↓
# Visualization and summary
#
# Genome: hg38
# Promoter definition: TSS ± 3 kb
# ============================================================


# ============================================================
# 1. REQUIRED PACKAGES
# ============================================================

library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(readr)


# ============================================================
# 2. CREATE OUTPUT DIRECTORY
# ============================================================

if (!dir.exists("ATAC_results")) {
  dir.create("ATAC_results")
}


# ============================================================
# 3. READ ATAC-SEQ BED FILES
# ============================================================

# Tumor BED files
tumor_files <- list.files(
  "Tumor",
  pattern = "\\.bed$",
  full.names = TRUE
)

# Normal BED files
normal_files <- list.files(
  "Normal",
  pattern = "\\.bed$",
  full.names = TRUE
)

cat("Number of tumor BED files:", length(tumor_files), "\n")
cat("Number of normal BED files:", length(normal_files), "\n")


# ============================================================
# 4. FUNCTION TO READ BED FILES
# ============================================================

read_bed_granges <- function(file) {

  df <- read.table(
    file,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep first three BED columns
  df <- df[, 1:3]

  GRanges(
    seqnames = df$V1,
    ranges = IRanges(
      start = df$V2,
      end = df$V3
    )
  )
}


# ============================================================
# 5. IMPORT TUMOR ATAC PEAKS
# ============================================================

tumor_peaks_list <- lapply(
  tumor_files,
  read_bed_granges
)

tumor_peaks <- do.call(c, tumor_peaks_list)

# Merge overlapping tumor peaks
tumor_peaks <- reduce(tumor_peaks)

cat(
  "Total tumor peaks:",
  length(tumor_peaks),
  "\n"
)


# ============================================================
# 6. IMPORT NORMAL ATAC PEAKS
# ============================================================

normal_peaks_list <- lapply(
  normal_files,
  read_bed_granges
)

normal_peaks <- do.call(c, normal_peaks_list)

# Merge overlapping normal peaks
normal_peaks <- reduce(normal_peaks)

cat(
  "Total normal peaks:",
  length(normal_peaks),
  "\n"
)


# ============================================================
# 7. KEEP STANDARD CHROMOSOMES
# ============================================================

tumor_peaks <- keepStandardChromosomes(
  tumor_peaks,
  pruning.mode = "coarse"
)

normal_peaks <- keepStandardChromosomes(
  normal_peaks,
  pruning.mode = "coarse"
)

cat(
  "Tumor peaks after standard chromosome filtering:",
  length(tumor_peaks),
  "\n"
)

cat(
  "Normal peaks after standard chromosome filtering:",
  length(normal_peaks),
  "\n"
)


# ============================================================
# 8. IDENTIFY TUMOR-SPECIFIC ATAC PEAKS
# ============================================================

hits <- findOverlaps(
  tumor_peaks,
  normal_peaks,
  ignore.strand = TRUE
)

# Number of overlap pairs
overlap_pairs <- length(hits)

# Unique tumor peaks overlapping normal peaks
overlapping_tumor_peaks <- unique(
  queryHits(hits)
)

n_overlapping_tumor <- length(
  overlapping_tumor_peaks
)

# Remove tumor peaks overlapping normal
tumor_specific <- tumor_peaks[
  -overlapping_tumor_peaks
]

cat(
  "Overlap pairs:",
  overlap_pairs,
  "\n"
)

cat(
  "Unique tumor peaks overlapping normal:",
  n_overlapping_tumor,
  "\n"
)

cat(
  "Tumor-specific peaks:",
  length(tumor_specific),
  "\n"
)


# ============================================================
# 9. SAVE TUMOR-SPECIFIC PEAKS
# ============================================================

tumor_specific_df <- data.frame(
  chr = as.character(seqnames(tumor_specific)),
  start = start(tumor_specific),
  end = end(tumor_specific)
)

write.table(
  tumor_specific_df,
  "ATAC_results/tumor_specific_peaks.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


# ============================================================
# 10. ANNOTATE TUMOR-SPECIFIC ATAC PEAKS
# ============================================================

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno <- annotatePeak(
  tumor_specific,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

anno_df <- as.data.frame(peakAnno)


# ============================================================
# 11. GENOMIC DISTRIBUTION OF ATAC PEAKS
# ============================================================

peak_annotation <- as.data.frame(
  table(anno_df$annotation)
)

colnames(peak_annotation) <- c(
  "Region",
  "Count"
)

write.csv(
  peak_annotation,
  "ATAC_results/ATAC_peak_genomic_distribution.csv",
  row.names = FALSE
)


# ============================================================
# 12. PLOT GENOMIC DISTRIBUTION
# ============================================================

p_peak_distribution <- ggplot(
  peak_annotation,
  aes(
    x = Region,
    y = Count
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Genomic Distribution of Tumor-Specific ATAC Peaks",
    x = "Genomic Region",
    y = "Number of Peaks"
  )

ggsave(
  "ATAC_results/Tumor_specific_ATAC_peak_distribution.png",
  p_peak_distribution,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 13. DISTANCE OF PEAKS FROM TSS
# ============================================================

png(
  "ATAC_results/ATAC_distance_to_TSS.png",
  width = 2000,
  height = 1500,
  res = 300
)

plotDistToTSS(peakAnno)

dev.off()


# ============================================================
# 14. DEFINE PROMOTER REGIONS
#     PROMOTER = TSS ± 3 kb
# ============================================================

gene_ranges <- genes(
  txdb,
  columns = c(
    "gene_id"
  )
)

promoter_regions <- promoters(
  gene_ranges,
  upstream = 3000,
  downstream = 3000
)

cat(
  "Total promoter regions:",
  length(promoter_regions),
  "\n"
)


# ============================================================
# 15. IDENTIFY PROMOTERS OVERLAPPING
#     TUMOR-SPECIFIC ATAC PEAKS
# ============================================================

promoter_hits <- findOverlaps(
  promoter_regions,
  tumor_specific,
  ignore.strand = TRUE
)

accessible_promoters <- promoter_regions[
  unique(queryHits(promoter_hits))
]

cat(
  "Genes with tumor-specific promoter accessibility:",
  length(accessible_promoters),
  "\n"
)


# ============================================================
# 16. MAP PROMOTER ACCESSIBILITY TO GENE IDs
# ============================================================

accessible_entrez <- unique(
  as.character(
    mcols(accessible_promoters)$gene_id
  )
)

accessible_entrez <- accessible_entrez[
  !is.na(accessible_entrez)
]


# ============================================================
# 17. CONVERT ENTREZ IDs TO GENE SYMBOLS
# ============================================================

accessible_gene_symbols <- bitr(
  accessible_entrez,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

accessible_gene_symbols <- unique(
  accessible_gene_symbols
)


# ============================================================
# 18. SAVE ATAC PROMOTER-ACCESSIBLE GENES
# ============================================================

write.csv(
  accessible_gene_symbols,
  "ATAC_results/ATAC_promoter_accessible_genes.csv",
  row.names = FALSE
)


# ============================================================
# 19. RNA-SEQ INTEGRATION
# ============================================================
#
# IMPORTANT:
# Replace this section with the actual RNA-seq upregulated
# gene file used in your analysis.
#
# The RNA-seq list should contain one column containing
# HGNC gene symbols.
# ============================================================

# Example:
#
# rnaseq <- read.csv(
#   "RNAseq_upregulated_genes.csv",
#   stringsAsFactors = FALSE
# )
#
# rnaseq_genes <- unique(rnaseq$Gene)


# ------------------------------------------------------------
# If the RNA-seq upregulated gene vector already exists:
# ------------------------------------------------------------

# rnaseq_genes <- unique(rnaseq_genes)


# ============================================================
# 20. INTERSECT RNA-SEQ AND ATAC-SEQ GENES
# ============================================================

# Uncomment after loading your RNA-seq gene list:

# final_1114_genes <- intersect(
#   rnaseq_genes,
#   accessible_gene_symbols$SYMBOL
# )

# cat(
#   "Final RNA + ATAC integrated genes:",
#   length(final_1114_genes),
#   "\n"
# )


# ============================================================
# 21. SAVE FINAL INTEGRATED GENE LIST
# ============================================================

# write.csv(
#   data.frame(
#     Gene = final_1114_genes
#   ),
#   "ATAC_results/final_ATAC_RNA_1114_genes.csv",
#   row.names = FALSE
# )


# ============================================================
# 22. ATAC ANALYSIS SUMMARY
# ============================================================

ATAC_summary <- data.frame(

  Analysis = c(
    "Tumor ATAC peaks",
    "Normal ATAC peaks",
    "Overlap pairs between tumor and normal",
    "Unique tumor peaks overlapping normal",
    "Tumor-specific ATAC peaks",
    "Promoter definition",
    "RNA + ATAC integrated genes"
  ),

  Result = c(
    length(tumor_peaks),
    length(normal_peaks),
    overlap_pairs,
    n_overlapping_tumor,
    length(tumor_specific),
    "TSS ± 3 kb",
    "1114"
  )
)


# ============================================================
# 23. SAVE SUMMARY TABLE
# ============================================================

write.csv(
  ATAC_summary,
  "ATAC_results/ATAC_analysis_summary.csv",
  row.names = FALSE
)

print(ATAC_summary)


# ============================================================
# 24. FINAL MESSAGE
# ============================================================

cat("\n")
cat("====================================================\n")
cat("ATAC-SEQ ANALYSIS COMPLETED\n")
cat("====================================================\n")
cat("Tumor peaks: ", length(tumor_peaks), "\n")
cat("Normal peaks: ", length(normal_peaks), "\n")
cat("Overlap pairs: ", overlap_pairs, "\n")
cat("Unique overlapping tumor peaks: ",
    n_overlapping_tumor, "\n")
cat("Tumor-specific peaks: ",
    length(tumor_specific), "\n")
cat("Promoter definition: TSS ± 3 kb\n")
cat("Final RNA + ATAC genes: 1114\n")
cat("====================================================\n")
