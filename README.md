# ATAC-seq and RNA-seq Integration – Prostate Cancer (TCGA-PRAD)

Identifies tumor-specific chromatin accessibility regions from ATAC-seq data and integrates them with RNA-seq upregulated genes to find epigenetically activated genes in prostate cancer, in the context of Coronarin D as a therapeutic candidate.

## Workflow
1. Merge tumor and normal ATAC-seq peaks (BED files)
2. Filter to standard chromosomes
3. Identify tumor-specific peaks (not present in normal tissue)
4. Annotate peaks and assess promoter accessibility (TSS ± 3 kb)
5. Integrate with RNA-seq upregulated gene list
6. Output: epigenetically activated gene set + summary visualizations

## Tools & Packages
R, Bioconductor — GenomicRanges, ChIPseeker, clusterProfiler, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, ggplot2, dplyr

## Genome build
hg38
