################################################################################
# 07. Interindividual variability in RNA-seq
# Laura Muñoz | 06/07/2026
#
# evaluates donor-associated expression variability and selects
# variable genes correlated with chromatin accessibility.
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
library(ggplot2) 
library(dplyr)
library(tidyr) 
library(tibble)
library(DESeq2)

in_dir <- "data"
in_dir_genes <- "results/06_ATAC_RNA_spearman"
out_dir <- "results/07_DE_RNA_LRT"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltype_labels <- c(
  HSC = "HSPC", CLP = "CLP", proB = "pro-B", preB = "pre-B",
  ImmatureB = "Immature B", Transitional_B = "Transitional B",
  Naive_CD5pos = "Naive B CD5+", Naive_CD5neg = "Naive B CD5-"
)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
rna_counts <- readRDS(file.path(in_dir, "processed", "rna_counts_clean.rds"))
metadata <- readRDS(file.path(in_dir, "processed", "rna_metadata_clean.rds"))
unique_genes <- readRDS(file.path(in_dir_genes, "unique_correlated_genes.rds"))

#-------------------------------------------------------------------------------
# 2. basic checks
#-------------------------------------------------------------------------------
metadata <- metadata %>%
  mutate(
    Donor = factor(Donor),
    CellType = factor(
      CellType,
      levels = c(
        "HSC", "CLP", "proB", "preB", "ImmatureB",
        "Transitional_B", "Naive_CD5pos", "Naive_CD5neg"
      )
    )
  )


stopifnot(identical(colnames(rna_counts), metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 3. create DESeqDataSet object
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = metadata,
  design = ~ CellType + Donor
)

#-------------------------------------------------------------------------------
# 4. filter low-count genes
#-------------------------------------------------------------------------------
# keep genes with at least 10 counts in at least 5 samples
keep_genes <- rowSums(counts(dds) >= 10) >= 5

message("Genes before filtering: ", nrow(dds))
message("Genes retained after filtering: ", sum(keep_genes))
message("Genes removed by filtering: ", nrow(dds) - sum(keep_genes))

dds <- dds[keep_genes, ]

#-------------------------------------------------------------------------------
# 5. run inter-individual differential expression analysis with LRT
#-------------------------------------------------------------------------------
dds <- DESeq(
  dds,
  test = "LRT",
  reduced = ~ CellType
)

saveRDS(dds, file.path(out_dir, "dds_RNA_LRT.rds"))

#-------------------------------------------------------------------------------
# 6. test correlated genes
#-------------------------------------------------------------------------------
# extract raw LRT p-values, no global multiple-testing correction is used.
res_lrt <- results(
  dds,
  independentFiltering = FALSE,
  pAdjustMethod = "none"
)

# keep only the global LRT columns
res_lrt_tbl <- res_lrt %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  as_tibble() %>%
  dplyr::select(
    ENSEMBL,
    baseMean,
    stat,
    pvalue
  )

message("Total genes tested: ", nrow(res_lrt_tbl))

# restrict the analysis to genes correlated with a selected OCR
correlated_genes_results <- res_lrt_tbl %>%
  filter(ENSEMBL %in% unique_genes) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(padj)

stopifnot(nrow(correlated_genes_results) == length(unique(unique_genes)))

message(
  "Correlated genes included in BH correction: ",
  nrow(correlated_genes_results)
)

#-------------------------------------------------------------------------------
# 7. calculate donor-effect magnitude
#-------------------------------------------------------------------------------
# keep donor coefficients only for the correlated genes
coef_tbl <- coef(dds) %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  as_tibble() %>%
  dplyr::select(
    ENSEMBL,
    starts_with("Donor_")
  ) %>%
  filter(ENSEMBL %in% unique_genes)

# convert donor coefficient columns to long format
coef_tbl_long <- coef_tbl %>%
  pivot_longer(
    cols = starts_with("Donor_"),
    names_to = c("Donor", "reference"),
    names_pattern = "Donor_(.*)_vs_(.*)",
    values_to = "donor_log2FC"
  ) %>%
  dplyr::select(
    ENSEMBL,
    Donor,
    donor_log2FC
  )

# add the reference donor, whose coefficient is zero
reference_donor <- levels(dds$Donor)[1]

reference_donor_tbl <- coef_tbl %>%
  transmute(
    ENSEMBL,
    Donor = reference_donor,
    donor_log2FC = 0
  )

coef_tbl_long <- coef_tbl_long %>%
  bind_rows(reference_donor_tbl) %>%
  arrange(ENSEMBL, Donor)

# save donor coefficients for the correlated genes.
saveRDS(
  coef_tbl_long,
  file.path(out_dir, "donor_coefficients_RNA_LRT.rds")
)

write.csv(
  coef_tbl_long,
  file.path(out_dir, "donor_coefficients_RNA_LRT.csv"),
  row.names = FALSE
)

donor_effect_size_tbl <- coef_tbl_long %>%
  group_by(ENSEMBL) %>%
  summarise(
    # donors with minimum and maximum estimated expression
    donor_min = Donor[which.min(donor_log2FC)],
    donor_max = Donor[which.max(donor_log2FC)],
    
    # minimum and maximum donor coefficients
    min_donor_log2FC = min(
      donor_log2FC,
      na.rm = TRUE
    ),
    
    max_donor_log2FC = max(
      donor_log2FC,
      na.rm = TRUE
    ),
    
    # maximum difference between donors on the log2 scale
    donor_log2FC_range =
      max_donor_log2FC - min_donor_log2FC,
    
    # general dispersion of donor coefficients
    donor_log2FC_sd = sd(
      donor_log2FC,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# add donor-effect magnitude to the correlated-gene results
correlated_genes_results <- correlated_genes_results %>%
  left_join(
    donor_effect_size_tbl,
    by = "ENSEMBL"
  )

#-------------------------------------------------------------------------------
# 8. select genes
#-------------------------------------------------------------------------------
padj_cutoff <- 0.05
min_log2FC_range <- log2(1.5)

selected_correlated_genes <- correlated_genes_results %>%
  filter(
    !is.na(padj),
    padj <= padj_cutoff,
    donor_log2FC_range >= min_log2FC_range
  ) %>%
  arrange(padj)

message(
  "Correlated genes with BH-adjusted p-value <= 0.05: ",
  sum(correlated_genes_results$padj <= padj_cutoff, na.rm = TRUE)
)

message(
  "Correlated genes selected after the effect-size filter: ",
  nrow(selected_correlated_genes)
)

lrt_selection_summary <- tibble(
  n_genes_in_lrt = nrow(res_lrt_tbl),
  n_correlated_genes = nrow(correlated_genes_results),
  n_correlated_genes_padj_at_or_below_0_05 = sum(
    correlated_genes_results$padj <= padj_cutoff,
    na.rm = TRUE
  ),
  n_selected_correlated_genes = nrow(selected_correlated_genes)
)

print(lrt_selection_summary)

#-------------------------------------------------------------------------------
# 9. save results
#-------------------------------------------------------------------------------
# summary of the candidate-gene selection
saveRDS(
  lrt_selection_summary,
  file.path(out_dir, "RNA_LRT_selection_summary.rds")
)

write.csv(
  lrt_selection_summary,
  file.path(out_dir, "RNA_LRT_selection_summary.csv"),
  row.names = FALSE
)

# LRT results for correlated genes
saveRDS(
  correlated_genes_results,
  file.path(
    out_dir,
    "correlated_genes_results_RNA_LRT.rds"
  )
)

write.csv(
  correlated_genes_results,
  file.path(
    out_dir,
    "correlated_genes_results_RNA_LRT.csv"
  ),
  row.names = FALSE
)

# final selected genes
saveRDS(
  selected_correlated_genes,
  file.path(
    out_dir,
    "selected_correlated_genes_RNA_LRT.rds"
  )
)

write.csv(
  selected_correlated_genes,
  file.path(
    out_dir,
    "selected_correlated_genes_RNA_LRT.csv"
  ),
  row.names = FALSE
)

selected_correlated_gene_ids <- selected_correlated_genes$ENSEMBL

saveRDS(
  selected_correlated_gene_ids,
  file.path(
    out_dir,
    "selected_correlated_genes_ids_RNA_LRT.rds"
  )
)

#-------------------------------------------------------------------------------
# 10. PCA with final selected genes
#-------------------------------------------------------------------------------
vsd <- vst(
  dds,
  blind = FALSE
  )
  
vsd_selected <- vsd[
  rownames(vsd) %in% selected_correlated_gene_ids,
]
  
pca_donor <- plotPCA(
  vsd_selected,
  intgroup = "Donor",
  ntop = nrow(vsd_selected)
)

pca_donor <- pca_donor +
  labs(
    title = "PCA of correlated genes with interindividual variability",
    x = pca_donor$labels$x,
    y = pca_donor$labels$y,
    color = "Donor"
  ) +
  theme_minimal()
  
ggsave(
  filename = file.path(
    out_dir,
    "PCA_selected_correlated_genes_by_donor.png"
  ),
  plot = pca_donor,
  width = 7,
  height = 5,
  dpi = 300
)
  
pca_celltype <- plotPCA(
  vsd_selected,
  intgroup = "CellType",
  ntop = nrow(vsd_selected)
)

pca_celltype <- pca_celltype +
  scale_color_discrete(labels = celltype_labels) +
  labs(
    title = "PCA of correlated genes by cell type",
    x = pca_celltype$labels$x,
    y = pca_celltype$labels$y,
    color = "Cell type"
  ) +
  theme_minimal()
  
ggsave(
  filename = file.path(
    out_dir,
    "PCA_selected_correlated_genes_by_celltype.png"
  ),
  plot = pca_celltype,
  width = 7,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 11. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
