################################################################################
# 04. Differential accessibility between cell types
# Laura Muñoz | 25/05/2026
#
# identifies differentially accessible regions between consecutive
# b-cell stages, accounting for donor differences.
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
library(pheatmap)

in_dir <- "data"
out_dir <- "results/04_DE_ATAC_celltypes"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltype_labels <- c(
  HSC = "HSPC", CLP = "CLP", proB = "pro-B", preB = "pre-B",
  ImmatureB = "Immature B", Transitional_B = "Transitional B",
  Naive_CD5pos = "Naive B CD5+", Naive_CD5neg = "Naive B CD5-"
)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
atac_counts <- readRDS(file.path(in_dir, "processed", "atac_counts_clean.rds"))
metadata <- readRDS(file.path(in_dir, "processed", "atac_metadata_clean.rds"))

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

colnames(atac_counts) <- metadata$Sample.Name
stopifnot(identical(colnames(atac_counts), metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 3. create DESeqDataSet object
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = atac_counts,
  colData = metadata,
  design = ~ Donor + CellType
)

#-------------------------------------------------------------------------------
# 4. filter low-count OCRs
#-------------------------------------------------------------------------------
# keep OCRs with at least 10 counts in at least 5 samples
keep_ocrs <- rowSums(counts(dds) >= 10) >= 5

message("OCRs before filtering: ", nrow(dds))
message("OCRs retained after filtering: ", sum(keep_ocrs))
message("OCRs removed by filtering: ", nrow(dds) - sum(keep_ocrs))

dds <- dds[keep_ocrs, ]

#-------------------------------------------------------------------------------
# 5. run differential accessibility analysis
#-------------------------------------------------------------------------------
dds <- DESeq(dds)

saveRDS(dds, file.path(out_dir, "dds_ATAC.rds"))

#-------------------------------------------------------------------------------
# 6. extract differential accessibility results
#-------------------------------------------------------------------------------
# define developmental transitions (state1 = reference, state2 = compared)
# positive log2foldchange indicates higher accessibility in state2
# selection criteria: adjusted p-value <= 0.01 and
# absolute log2foldchange >= log2(2)
transitions <- list(
  c("HSC", "CLP"),
  c("CLP", "proB"),
  c("proB", "preB"),
  c("preB", "ImmatureB"),
  c("ImmatureB", "Transitional_B"),
  c("Transitional_B", "Naive_CD5pos"),
  c("Transitional_B", "Naive_CD5neg")
)

res_contrasts <- list()
sig_res <- list()

for (pair in transitions) {
  state1 <- pair[1]
  state2 <- pair[2]
  contrast_name <- paste(state2, "vs", state1, sep = "_")
  
  res <- results(
    dds,
    contrast = c("CellType", state2, state1),
    alpha = 0.01
  )

  res_contrasts[[contrast_name]] <- res
  
  sig_res[[contrast_name]] <- 
    as.data.frame(res_contrasts[[contrast_name]]) %>%
    rownames_to_column(var = "ocr_id") %>%
    filter(
      !is.na(padj),
      padj <= 0.01,
      abs(log2FoldChange) >= log2(2)
      ) %>%
    arrange(padj)
  
  message(
    "Number of OCRs for contrast ", 
    contrast_name,
    ": ",
    nrow(sig_res[[contrast_name]])
    )
}

saveRDS(res_contrasts, file.path(out_dir, "complete_results_DE_ATAC.rds"))

saveRDS(sig_res, file.path(out_dir, "DARs_per_contrast_ATAC.rds"))

#-------------------------------------------------------------------------------
# 7. summarise DARs (differentially accessible regions) results
#-------------------------------------------------------------------------------
dars_summary <- tibble(
  contrast = names(sig_res),
  total_DARs = sapply(sig_res, nrow),
  more_accessible = sapply(sig_res, function(x) sum(x$log2FoldChange > 0)),
  less_accessible = sapply(sig_res, function(x) sum(x$log2FoldChange < 0))
)

saveRDS(dars_summary, file.path(out_dir, "DARs_summary.rds"))

write.csv(
  dars_summary, 
  file.path(out_dir, "DARs_summary_ATAC.csv"), 
  row.names = FALSE
  )

#-------------------------------------------------------------------------------
# 8. export DARs tables
#-------------------------------------------------------------------------------
for (contrast in names(sig_res)){
  write.csv(
    sig_res[[contrast]],
    file.path(out_dir, paste0("DARs_", contrast, ".csv")),
    row.names = FALSE
  )
}

#-------------------------------------------------------------------------------
# 9. PCA with DARs
#-------------------------------------------------------------------------------
# define DARs union
dars_lists <- list()

for (contrast in names(sig_res)){
  dars_lists[[contrast]] <- sig_res[[contrast]]$ocr_id
}

dars_ocrs <- unique(unlist(dars_lists))

# VST transformation
vsd <- vst(dds, blind = FALSE)

# subset VST object to significant DARs
vsd_dars <- vsd[rownames(vsd) %in% dars_ocrs, ]

# PCA colored by cell type
pca_celltype <- plotPCA(
  vsd_dars, 
  intgroup = "CellType",
  ntop = nrow(vsd_dars)
  )

pca_celltype <- pca_celltype +
  scale_color_discrete(labels = celltype_labels) +
  labs(
    title = "PCA of differentially accessible regions (ATAC-seq)",
    x = pca_celltype$labels$x,
    y = pca_celltype$labels$y,
    color = "Cell type"
  ) +
  theme_minimal()
  
ggsave(
  filename = file.path(out_dir, "PCA_DARs_ATAC.png"),
  plot = pca_celltype,
  width = 7,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 10. heatmap of top ATAC DARs
#-------------------------------------------------------------------------------
# select top 3 DARs from each contrast
top_ocrs <- c()

for (contrast in names(sig_res)) {
  ocrs <- head(sig_res[[contrast]]$ocr_id, 3)
  top_ocrs <- c(top_ocrs, ocrs)
}

top_ocrs <- unique(top_ocrs)

message("OCRs selected for heatmap: ", length(top_ocrs))

# save selected OCRs
write.csv(
  data.frame(ocr_id = top_ocrs),
  file.path(out_dir, "top_DARs_heatmap_ATAC.csv"),
  row.names = FALSE
)

# extract VST accessibility values
heatmap_matrix <- assay(vsd)[top_ocrs, ]

# scale by OCR
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# order samples by cell type and donor
celltype_order <- c(
  "HSC", "CLP", "proB", "preB",
  "ImmatureB", "Transitional_B",
  "Naive_CD5pos", "Naive_CD5neg"
)

metadata_heatmap <- metadata %>%
  mutate(CellType = factor(CellType, levels = celltype_order)) %>%
  arrange(CellType, Donor)

heatmap_matrix_scaled <- heatmap_matrix_scaled[
  , metadata_heatmap$Sample.Name
]

# separate cell types in the heatmap
celltype_counts <- table(metadata_heatmap$CellType)
gaps_col <- cumsum(celltype_counts)
gaps_col <- gaps_col[-length(gaps_col)]

# annotation
annotation_col <- metadata_heatmap %>%
  transmute(
    `Cell type` = factor(
      celltype_labels[as.character(CellType)],
      levels = unname(celltype_labels)
    ),
    Donor = factor(Donor)
  )

rownames(annotation_col) <- metadata_heatmap$Sample.Name

# use the same cell type colors in both heatmaps
celltype_colors <- c(
  "HSPC" = "#7CAE00", "CLP" = "#00BFC4", "pro-B" = "#619CFF",
  "pre-B" = "#C5A3FF", "Immature B" = "#F8766D",
  "Transitional B" = "#F564E3", "Naive B CD5+" = "#E69F00",
  "Naive B CD5-" = "#00BA9C"
)

# use the same donor colors in all heatmaps
donor_colors <- c(
  D198 = "#8DD3C7", D199 = "#BEBADA", D200 = "#FB8072",
  D201 = "#80B1D3", D202 = "#FDB462", D204 = "#B3DE69",
  D206 = "#FCCDE5", D210 = "#BC80BD", D211 = "#CCEBC5",
  D212 = "#FFED6F", D213 = "#A65628", D214 = "#1B9E77",
  D215 = "#666666"
)
donor_colors <- donor_colors[sort(unique(as.character(metadata$Donor)))]

annotation_colors <- list(`Cell type` = celltype_colors, Donor = donor_colors)

# heatmap
top_dar_heatmap <- pheatmap(
  heatmap_matrix_scaled,
  breaks = seq(-2.5, 2.5, length.out = 101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = gaps_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = paste(
    "Top differentially accessible regions during B-cell differentiation",
    "Relative accessibility (VST, z-score)",
    sep = "\n"
  ),
  fontsize = 7,
  fontsize_row = 6,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("-2   ", "0   ", "2   ")
)

ggsave(
  filename = file.path(out_dir, "heatmap_top_DARs_ATAC.png"),
  plot = top_dar_heatmap$gtable,
  width = 8,
  height = 10,
  dpi = 300,
  bg = "white"
)

#-------------------------------------------------------------------------------
# 11. export dynamic OCRs union
#-------------------------------------------------------------------------------
dynamic_ocr_ids <- unique(unlist(dars_lists))

dynamic_ocrs_tbl <- tibble(
  OCR = dynamic_ocr_ids
)

message("Total unique dynamic OCRs: ", nrow(dynamic_ocrs_tbl))

saveRDS(
  dynamic_ocr_ids,
  file.path(out_dir, "dynamic_ocrs_ids_ATAC.rds")
)

saveRDS(
  dynamic_ocrs_tbl,
  file.path(out_dir, "dynamic_ocrs_ATAC.rds")
)

write.csv(
  dynamic_ocrs_tbl,
  file.path(out_dir, "dynamic_ocrs_ATAC.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 12. export dynamic OCRs by transition
#-------------------------------------------------------------------------------
dynamic_ocrs_by_transition <- bind_rows(
  lapply(names(sig_res), function(contrast_name) {
    sig_res[[contrast_name]] %>%
      mutate(contrast = contrast_name)
  })
)

saveRDS(
  dynamic_ocrs_by_transition,
  file.path(out_dir, "dynamic_ocrs_by_transition_ATAC.rds")
)

write.csv(
  dynamic_ocrs_by_transition,
  file.path(out_dir, "dynamic_ocrs_by_transition_ATAC.csv"),
  row.names = FALSE
)
#-------------------------------------------------------------------------------
# 13. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
