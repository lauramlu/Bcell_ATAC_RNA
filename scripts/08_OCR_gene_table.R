################################################################################
# 08. Final OCR-gene selection and annotation
# Laura Muñoz | 23/07/2026
#
# Selects final OCR-gene pairs and combines statistical results
# with genomic annotations, summary tables and figures.
################################################################################

#-------------------------------------------------------------------------------
# 0. set up
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

message("R version: ", R.version.string)
message("Working directory: ", getwd())

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

in_dir_spearman <- "results/06_ATAC_RNA_spearman"
in_dir_rna_lrt <- "results/07_DE_RNA_LRT"
in_dir_annotation <- "data/processed"
in_dir_atac_lrt <- "results/05_DE_ATAC_LRT"
out_dir <- "results/08_final_OCR_gene_table"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
sig_spearman_results <- readRDS(
  file.path(
    in_dir_spearman,
    "significant_OCR_gene_spearman_results.rds"
  )
)

selected_rna_genes <- readRDS(
  file.path(
    in_dir_rna_lrt,
    "selected_correlated_genes_RNA_LRT.rds"
  )
)

atac_annotation <- readRDS(
  file.path(in_dir_annotation, "atac_OCR_annotation.rds")
)

atac_lrt_results <- readRDS(
  file.path(in_dir_atac_lrt, "common_ocrs_results.rds")
)

#-------------------------------------------------------------------------------
# 2. basic checks
#-------------------------------------------------------------------------------
message("Significant Spearman pairs: ", nrow(sig_spearman_results))
message(
  "Unique genes with significant Spearman correlation: ",
  n_distinct(sig_spearman_results$ENSEMBL)
)
message(
  "Correlated genes selected by RNA LRT: ",
  nrow(selected_rna_genes)
)

#-------------------------------------------------------------------------------
# 3. results tables
#-------------------------------------------------------------------------------
# rename statistical columns to indicate the analysis they belong to
sig_spearman_results <- sig_spearman_results %>%
  rename(
    spearman_pvalue = pvalue,
    spearman_padj = padj
  )

selected_rna_genes <- selected_rna_genes %>%
  rename(
    rna_baseMean = baseMean,
    rna_lrt_stat = stat,
    rna_lrt_pvalue = pvalue,
    rna_lrt_padj = padj
  )

#-------------------------------------------------------------------------------
# 4. select OCR-gene pairs
#-------------------------------------------------------------------------------
# keep significant Spearman pairs with gene also selected by the RNA LRT
selected_ocr_gene_pairs <- sig_spearman_results %>%
  inner_join(
    selected_rna_genes,
    by = "ENSEMBL"
  ) %>%
  arrange(rna_lrt_padj, spearman_padj)

selected_ocr_ids <- unique(selected_ocr_gene_pairs$OCR)
selected_gene_ids <- unique(selected_ocr_gene_pairs$ENSEMBL)

message("Selected OCR-gene pairs: ", nrow(selected_ocr_gene_pairs))
message("Unique selected OCRs: ", length(selected_ocr_ids))
message("Unique selected genes: ", length(selected_gene_ids))

#-------------------------------------------------------------------------------
# 5. summarise selected results
#-------------------------------------------------------------------------------
selected_summary <- tibble(
  n_significant_spearman_pairs = nrow(sig_spearman_results),
  n_spearman_correlated_genes = n_distinct(sig_spearman_results$ENSEMBL),
  n_rna_lrt_selected_genes = nrow(selected_rna_genes),
  n_selected_ocr_gene_pairs = nrow(selected_ocr_gene_pairs),
  n_selected_unique_ocrs = length(selected_ocr_ids),
  n_selected_unique_genes = length(selected_gene_ids),
  n_positive_correlations = sum(selected_ocr_gene_pairs$rho > 0, na.rm = TRUE),
  n_negative_correlations = sum(selected_ocr_gene_pairs$rho < 0, na.rm = TRUE)
)

print(selected_summary)

#-------------------------------------------------------------------------------
# 6. prepare RNA and Spearman results
#-------------------------------------------------------------------------------
rna_spearman_results <- selected_ocr_gene_pairs %>%
  transmute(
    OCR,
    ENSEMBL,
    SYMBOL,
    spearman_rho = rho,
    spearman_padj,
    rna_lrt_padj,
    rna_donor_log2FC_range = donor_log2FC_range
  )

#-------------------------------------------------------------------------------
# 7. prepare ATAC results
#-------------------------------------------------------------------------------
atac_lrt_results <- atac_lrt_results %>%
  transmute(
    OCR,
    atac_lrt_padj = padj,
    atac_donor_log2FC_range = donor_log2FC_range
  )

#-------------------------------------------------------------------------------
# 8. prepare genomic annotation
#-------------------------------------------------------------------------------
atac_annotation <- atac_annotation %>%
  as.data.frame() %>%
  rownames_to_column("OCR") %>%
  transmute(
    OCR,
    GENENAME,
    chr = as.character(seqnames),
    start,
    end,
    annotation_simplified,
    distanceToTSS
  )

#-------------------------------------------------------------------------------
# 9. create table
#-------------------------------------------------------------------------------
ocr_gene_table <- rna_spearman_results %>%
  left_join(
    atac_lrt_results,
    by = "OCR"
  ) %>%
  left_join(
    atac_annotation,
    by = "OCR"
  ) %>%
  dplyr::select(
    OCR,
    ENSEMBL,
    SYMBOL,
    GENENAME,
    chr,
    start,
    end,
    annotation_simplified,
    distanceToTSS,
    spearman_rho,
    spearman_padj,
    rna_lrt_padj,
    rna_donor_log2FC_range,
    atac_lrt_padj,
    atac_donor_log2FC_range
  ) %>%
  arrange(
    rna_lrt_padj,
    atac_lrt_padj,
    spearman_padj
  )

stopifnot(
  nrow(ocr_gene_table) == nrow(selected_ocr_gene_pairs)
)

#-------------------------------------------------------------------------------
# 10. summarise OCR-gene table
#-------------------------------------------------------------------------------
ocr_gene_summary <- tibble(
  n_ocr_gene_pairs = nrow(ocr_gene_table),
  n_unique_ocrs = n_distinct(ocr_gene_table$OCR),
  n_unique_genes = n_distinct(ocr_gene_table$ENSEMBL),
  n_positive_correlations = sum(ocr_gene_table$spearman_rho > 0),
  n_negative_correlations = sum(ocr_gene_table$spearman_rho < 0)
)

print(ocr_gene_summary)

#-------------------------------------------------------------------------------
# 11. number of correlated OCRs per gene
#-------------------------------------------------------------------------------
gene_ocr_counts <- ocr_gene_table %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL = dplyr::first(SYMBOL),
    n_OCR = n_distinct(OCR)
    )

ocr_proportions_binary <- gene_ocr_counts %>%
  mutate(
    OCR_group = if_else(
      n_OCR == 1,
      "1 OCR",
      ">1 OCR"
    ),
    OCR_group = factor(
      OCR_group,
      levels = c("1 OCR", ">1 OCR")
    )
  ) %>%
  dplyr::count(
    OCR_group,
    name = "n_genes"
  ) %>%
  mutate(
    proportion = n_genes/sum(n_genes),
    percentage = proportion * 100
  )

print(ocr_proportions_binary)

ocr_proportions_detailed <- gene_ocr_counts %>%
  mutate(
    OCR_group = case_when(
      n_OCR == 1 ~ "1 OCR",
      n_OCR == 2 ~ "2 OCR",
      n_OCR == 3 ~ "3 OCR",
      n_OCR == 4 ~ "4 OCR",
      n_OCR >= 5 ~ "≥5 OCR"
    ),
    OCR_group = factor(
      OCR_group,
      levels = c(
        "1 OCR",
        "2 OCR",
        "3 OCR",
        "4 OCR",
        "≥5 OCR"
      )
    )
  ) %>%
  dplyr::count(
    OCR_group,
    name = "n_genes"
  ) %>%
  mutate(
    proportion = n_genes / sum(n_genes),
    percentage = proportion * 100
  )

print(ocr_proportions_detailed)

ocr_proportions_plot <- ggplot(
  ocr_proportions_detailed,
  aes(
    x = OCR_group,
    y = percentage
  )
) +
  geom_col(
    width = 0.6,
    fill = "#D55C2C"
  ) +
  geom_text(
    aes(
      label = paste0(
        n_genes,
        " genes\n(",
        sprintf("%.1f", percentage),
        "%)"
      )
    ),
    vjust = -0.35,
    size = 3.8
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Number of correlated OCRs per gene",
    x = "Number of correlated OCRs",
    y = "Proportion of genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      hjust = 0.5
    )
  )

print(ocr_proportions_plot)

#-------------------------------------------------------------------------------
# 12. genomic annotation proportions
#-------------------------------------------------------------------------------
ocr_annotation_proportions <- ocr_gene_table %>%
  distinct(
    OCR,
    annotation_simplified
  ) %>%
  dplyr::count(
    annotation_simplified,
    name = "n_OCR"
  ) %>%
  mutate(
    proportion = n_OCR / sum(n_OCR),
    percentage = proportion * 100,
    
    # translate annotation labels
    annotation_label = recode(
      annotation_simplified,
      Intron = "Intron",
      `Distal Intergenic` = "Distal intergenic",
      Exon = "Exon",
      Promoter = "Promoter",
      `3' UTR` = "3' UTR",
      `5' UTR` = "5' UTR",
      `Downstream (<1kb)` = "Downstream (<1 kb)",
      `Downstream (1-2kb)` = "Downstream (1–2 kb)",
      `Downstream (2-3kb)` = "Downstream (2–3 kb)"
    ),
    
    # label displayed next to each bar
    plot_label = paste0(
      n_OCR,
      " OCR (",
      sprintf("%.1f", percentage),
      "%)"
    )
  ) %>%
  arrange(percentage) %>%
  mutate(
    annotation_label = factor(
      annotation_label,
      levels = annotation_label
    )
  )

print(ocr_annotation_proportions)

ocr_annotation_plot <- ggplot(
  ocr_annotation_proportions,
  aes(
    x = annotation_label,
    y = percentage
  )
) +
  geom_col(
    width = 0.7,
    fill = "#D55C2C"
  ) +
  geom_text(
    aes(label = plot_label),
    hjust = -0.1,
    size = 3.7
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Genomic annotation of final OCRs",
    x = NULL,
    y = "Proportion of OCRs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 50,
      b = 10,
      l = 10
    )
  )

print(ocr_annotation_plot)

#-------------------------------------------------------------------------------
# 13. save results
#-------------------------------------------------------------------------------
saveRDS(
  selected_ocr_gene_pairs,
  file.path(out_dir, "selected_OCR_gene_pairs.rds")
)

write.csv(
  selected_ocr_gene_pairs,
  file.path(out_dir, "selected_OCR_gene_pairs.csv"),
  row.names = FALSE
)

saveRDS(
  selected_ocr_ids,
  file.path(out_dir, "selected_OCR_ids.rds")
)

saveRDS(
  selected_gene_ids,
  file.path(out_dir, "selected_gene_ids.rds")
)

saveRDS(
  selected_summary,
  file.path(out_dir, "selected_results_summary.rds")
)

write.csv(
  selected_summary,
  file.path(out_dir, "selected_results_summary.csv"),
  row.names = FALSE
)

saveRDS(
  ocr_gene_table,
  file.path(out_dir, "ocr_gene_table.rds")
)

write.csv(
  ocr_gene_table,
  file.path(out_dir, "ocr_gene_table.csv"),
  row.names = FALSE
)

saveRDS(
  ocr_gene_summary,
  file.path(out_dir, "ocr_gene_table_summary.rds")
)

write.csv(
  ocr_gene_summary,
  file.path(out_dir, "ocr_gene_table_summary.csv"),
  row.names = FALSE
)

write.csv(
  gene_ocr_counts,
  file.path(out_dir, "OCR_count_per_gene.csv"),
  row.names = FALSE
)

write.csv(
  ocr_proportions_binary,
  file.path(out_dir, "OCR_proportions_binary.csv"),
  row.names = FALSE
)

write.csv(
  ocr_proportions_detailed,
  file.path(out_dir, "OCR_proportions_detailed.csv"),
  row.names = FALSE
)

saveRDS(
  gene_ocr_counts,
  file.path(out_dir, "OCR_count_per_gene.rds")
)

saveRDS(
  ocr_proportions_binary,
  file.path(out_dir, "OCR_proportions_binary.rds")
)

saveRDS(
  ocr_proportions_detailed,
  file.path(out_dir, "OCR_proportions_detailed.rds")
)

ggsave(
  filename = file.path(
    out_dir,
    "OCR_number_per_gene.png"
  ),
  plot = ocr_proportions_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  ocr_annotation_proportions,
  file.path(
    out_dir,
    "OCR_genomic_annotation_proportions.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ocr_annotation_proportions,
  file.path(
    out_dir,
    "OCR_genomic_annotation_proportions.rds"
  )
)

ggsave(
  filename = file.path(
    out_dir,
    "OCR_genomic_annotation_proportions.png"
  ),
  plot = ocr_annotation_plot,
  width = 9,
  height = 6,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 14. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
