################################################################################
# 03. Differential expression between cell types
# Laura Muñoz | 18/05/2026
#
# identifies differentially expressed genes between consecutive
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
library(clusterProfiler)
library(org.Hs.eg.db)

in_dir <- "data"
out_dir <- "results/03_DE_RNA_celltypes"

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

#-------------------------------------------------------------------------------
# 2. basic checks
#-------------------------------------------------------------------------------
metadata <- metadata %>%
  mutate(
    CellType = factor(CellType)
  )

stopifnot(identical(colnames(rna_counts), metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 3. create DESeqDataSet object
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = metadata,
  design = ~ Donor + CellType
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
# 5. run differential expression analysis
#-------------------------------------------------------------------------------
dds <- DESeq(dds)

saveRDS(dds, file.path(out_dir, "dds_RNA.rds"))

#-------------------------------------------------------------------------------
# 6. extract differential expression results
#-------------------------------------------------------------------------------
# define developmental transitions (state1 = reference, state2 = compared)
# positive log2foldchange indicates higher expression in state2
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
    alpha = 0.05
  )

  res_contrasts[[contrast_name]] <- res
  
  sig_res[[contrast_name]] <- 
    as.data.frame(res_contrasts[[contrast_name]]) %>%
    rownames_to_column(var = "gene_id") %>%
    filter(
      !is.na(padj),
      padj <= 0.05,
      abs(log2FoldChange) >= log2(1.5)
    ) %>%
    arrange(padj)
  
  message(
    "Number of DEGs for contrast ", 
    contrast_name,
    ": ",
    nrow(sig_res[[contrast_name]])
    )
}

saveRDS(res_contrasts, file.path(out_dir, "complete_results_DE_RNA.rds"))

saveRDS(sig_res, file.path(out_dir, "DEG_per_contrast_RNA.rds"))

#-------------------------------------------------------------------------------
# 7. summarise DEG results
#-------------------------------------------------------------------------------
deg_summary <- tibble(
  contrast = names(sig_res),
  total_DEGs = sapply(sig_res, nrow),
  upregulated = sapply(sig_res, function(x) sum(x$log2FoldChange > 0)),
  downregulated = sapply(sig_res, function(x) sum(x$log2FoldChange < 0))
)

saveRDS(deg_summary, file.path(out_dir, "DEG_summary.rds"))
write.csv(deg_summary, file.path(out_dir, "DEG_summary_RNA.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# 8. export DEG tables
#-------------------------------------------------------------------------------
for (contrast in names(sig_res)){
  write.csv(
    sig_res[[contrast]],
    file.path(out_dir, paste0("DEGs_", contrast, ".csv")),
    row.names = FALSE
  )
}

#-------------------------------------------------------------------------------
# 9. PCA with DEGs
#-------------------------------------------------------------------------------
# define DEG union
deg_lists <- list()

for (contrast in names(sig_res)){
  deg_lists[[contrast]] <- sig_res[[contrast]]$gene_id
}

deg_genes <- unique(unlist(deg_lists))

# save the unique gene identifiers from all contrasts
saveRDS(deg_genes, file.path(out_dir, "DEG_union.rds"))

# VST transformation
vsd <- vst(dds, blind = FALSE)

# subset VST object to significant DEGs
vsd_deg <- vsd[rownames(vsd) %in% deg_genes, ]

# PCA colored by cell type
pca_celltype <- plotPCA(
  vsd_deg, 
  intgroup = "CellType",
  ntop = nrow(vsd_deg)
  )

pca_celltype <- pca_celltype +
  scale_color_discrete(labels = celltype_labels) +
  labs(
    title = "PCA of differentially expressed genes (RNA-seq)",
    x = pca_celltype$labels$x,
    y = pca_celltype$labels$y,
    color = "Cell type"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "PCA_DEGs_RNA.png"),
  plot = pca_celltype,
  width = 7,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 10. heatmap of top RNA-seq DEGs
#-------------------------------------------------------------------------------
# select top 3 DEGs from each contrast
top_genes <- c()

for (contrast in names(sig_res)) {
  genes <- head(sig_res[[contrast]]$gene_id, 3)
  top_genes <- c(top_genes, genes)
}

top_genes <- unique(top_genes)

message("Genes selected for heatmap: ", length(top_genes))

# convert ENSEMBL identifiers to SYMBOL
gene_id_mapping <- bitr(
  top_genes,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

# keep symbols in the same order as the selected genes
gene_symbols <- gene_id_mapping$SYMBOL[
  match(top_genes, gene_id_mapping$ENSEMBL)
]

# save selected genes
write.csv(
  data.frame(gene_id = top_genes, SYMBOL = gene_symbols),
  file.path(out_dir, "top_DEGs_heatmap_RNA.csv"),
  row.names = FALSE
)

# extract VST expression values
heatmap_matrix <- assay(vsd)[top_genes, ]

# scale by gene
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
top_deg_heatmap <- pheatmap(
  heatmap_matrix_scaled,
  breaks = seq(-2.5, 2.5, length.out = 101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = gaps_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = gene_symbols,
  show_colnames = FALSE,
  main = paste(
    "Top differentially expressed genes during B-cell differentiation",
    "Relative expression (VST, z-score)",
    sep = "\n"
  ),
  fontsize = 7,
  fontsize_row = 8,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("-2   ", "0   ", "2   ")
)

ggsave(
  filename = file.path(out_dir, "heatmap_top_DEGs_RNA.png"),
  plot = top_deg_heatmap$gtable,
  width = 8,
  height = 10,
  dpi = 300,
  bg = "white"
)

#-------------------------------------------------------------------------------
# 11. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
