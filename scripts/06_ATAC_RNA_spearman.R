################################################################################
# 06. Correlation between ATAC-seq and RNA-seq
# Laura Muñoz | 20/07/2026
#
# tests Spearman correlations between chromatin accessibility
# and gene expression in selected OCR-gene pairs.
################################################################################

#-------------------------------------------------------------------------------
# 0. set up
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

# print r version and working directory
message("R version: ", R.version.string)
message("Working directory: ", getwd())

# packages
library(dplyr)
library(tibble)
library(DESeq2)

in_dir_annotation <- "data"
in_dir <- "results"
out_dir <- "results/06_ATAC_RNA_spearman"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
dds_rna <- readRDS(file.path(in_dir, "03_DE_RNA_celltypes", "dds_RNA.rds"))
dds_atac <- readRDS(file.path(in_dir, "04_DE_ATAC_celltypes", "dds_ATAC.rds"))

atac_OCR_annotation <- readRDS(
  file.path(in_dir_annotation, "processed", "atac_OCR_annotation.rds")
  )

common_ocrs_results <- readRDS(
  file.path(in_dir, "05_DE_ATAC_LRT", "common_ocrs_results.rds")
)

# union of significant DEGs from all cell-type transitions
deg_gene_ids <- readRDS(
  file.path(in_dir, "03_DE_RNA_celltypes", "DEG_union.rds")
)

#-------------------------------------------------------------------------------
# 2. VST transformation
#-------------------------------------------------------------------------------
message("RNA design: ", as.character(design(dds_rna)))
message("ATAC design: ", as.character(design(dds_atac)))

vsd_rna <- vst(dds_rna, blind = FALSE)
vsd_atac <- vst(dds_atac, blind = FALSE)

rna_mat <- assay(vsd_rna)
atac_mat <- assay(vsd_atac)

stopifnot(all(deg_gene_ids %in% rownames(rna_mat)))

message(
  "RNA matrix dimensions: ", nrow(rna_mat), 
  " genes x ", ncol(rna_mat), " samples"
  )

message(
  "ATAC matrix dimensions: ", nrow(atac_mat), 
  " OCRs x ", ncol(atac_mat), " samples"
  )

#-------------------------------------------------------------------------------
# 3. match RNA-seq and ATAC-seq samples
#-------------------------------------------------------------------------------
common_samples <- intersect(colnames(rna_mat), colnames(atac_mat))

message("RNA samples: ", ncol(rna_mat))
message("ATAC samples: ", ncol(atac_mat))
message("Common samples: ", length(common_samples))

rna_mat <- rna_mat[, common_samples]
atac_mat <- atac_mat[, common_samples]

stopifnot(identical(colnames(rna_mat), colnames(atac_mat)))

#-------------------------------------------------------------------------------
# 4. filter OCRs (dynamic and with interindividual variability)
#-------------------------------------------------------------------------------
stopifnot(all(common_ocrs_results$OCR %in% rownames(atac_mat)))

message(
  "OCRs that are dynamic and donor-variable: ",
  nrow(common_ocrs_results)
)

atac_mat <- atac_mat[common_ocrs_results$OCR, ]

stopifnot(identical(rownames(atac_mat), common_ocrs_results$OCR))

#-------------------------------------------------------------------------------
# 5. filter OCR-gene pairs 
#-------------------------------------------------------------------------------
# start with every OCR-gene pair in the annotation
all_annotated_pairs <- atac_OCR_annotation %>%
  dplyr::select(ENSEMBL, SYMBOL) %>%
  rownames_to_column("OCR")

message("All annotation rows: ", nrow(all_annotated_pairs))

# remove pairs without an annotated ENSEMBL gene
pairs_with_gene <- all_annotated_pairs %>%
  filter(!is.na(ENSEMBL))

message("Pairs with an annotated gene: ", nrow(pairs_with_gene))

# keep OCRs that are dynamic and donor-variable
pairs_selected_ocrs <- pairs_with_gene %>%
  filter(OCR %in% rownames(atac_mat))

message("Pairs after selected-OCR filter: ", nrow(pairs_selected_ocrs))

# keep genes that are DEGs in at least one cell-type transition
OCR_gene <- pairs_selected_ocrs %>%
  filter(ENSEMBL %in% deg_gene_ids)

message("Pairs after cell-type DEG filter: ", nrow(OCR_gene))
message("Unique DEGs retained: ", n_distinct(OCR_gene$ENSEMBL))

saveRDS(OCR_gene, file.path(out_dir, "OCR_genes_clean.rds"))

#-------------------------------------------------------------------------------
# 6. Spearman correlation for OCR-gene pairs
#-------------------------------------------------------------------------------
spearman_results <- lapply(seq_len(nrow(OCR_gene)), function(i) {
  
  ocr_id <- OCR_gene$OCR[i]
  gene_id <- OCR_gene$ENSEMBL[i]
  gene_symbol <- OCR_gene$SYMBOL[i]
  
  x <- as.numeric(atac_mat[ocr_id, ])
  y <- as.numeric(rna_mat[gene_id, ])
  
  cor_res <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    OCR = ocr_id,
    ENSEMBL = gene_id,
    SYMBOL = gene_symbol,
    rho = unname(cor_res$estimate),
    pvalue = cor_res$p.value
  )
})

spearman_results <- bind_rows(spearman_results) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(padj)

sig_spearman_results <- spearman_results %>%
  filter(padj <= 0.05)

message("Total OCR-gene pairs tested: ", nrow(spearman_results))
message("Significant OCR-gene correlations: ", nrow(sig_spearman_results))

saveRDS(
  spearman_results, 
        file.path(out_dir, "complete_OCR_gene_spearman_results.rds")
  )
saveRDS(
  sig_spearman_results, 
        file.path(out_dir, "significant_OCR_gene_spearman_results.rds")
  )

write.csv(
  spearman_results,
  file.path(out_dir, "complete_OCR_gene_spearman_results.csv"),
  row.names = FALSE
)

write.csv(
  sig_spearman_results,
  file.path(out_dir, "significant_OCR_gene_spearman_results.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 7. extract unique genes
#-------------------------------------------------------------------------------
correlated_genes <- unique(sig_spearman_results$ENSEMBL)

saveRDS(correlated_genes,file.path(out_dir, "unique_correlated_genes.rds"))

#-------------------------------------------------------------------------------
# 8. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
