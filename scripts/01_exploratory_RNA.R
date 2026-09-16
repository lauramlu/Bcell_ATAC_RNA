################################################################################
# 01. Exploratory RNA-seq analysis
# Laura Muñoz | 11/05/2026
#
# checks RNA-seq data quality and explores sample variation
# by cell type and donor.
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
out_dir <- "results/01_exploratory_RNA"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (subdir in c("composition", "qc_raw", "qc_deseq2", "pca", "heatmaps")) {
  dir.create(file.path(out_dir, subdir), showWarnings = FALSE, recursive = TRUE)
}

celltype_labels <- c(
  HSC = "HSPC", CLP = "CLP", proB = "pro-B", preB = "pre-B",
  ImmatureB = "Immature B", Transitional_B = "Transitional B",
  Naive_CD5pos = "Naive B CD5+", Naive_CD5neg = "Naive B CD5-"
)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
rna_counts <- readRDS(file.path(in_dir, "raw", "rna_raw_data.rds"))
metadata <- readRDS(file.path(in_dir, "metadata", "rna_metadata.rds"))

#-------------------------------------------------------------------------------
# 2. object type and dimensions
#-------------------------------------------------------------------------------
message("Rna_counts class:", class(rna_counts), "\n")
message("Metadata class:", class(metadata), "\n")

message(
  "RNA_counts dim: ", 
  nrow(rna_counts), 
  " genes x ", 
  ncol(rna_counts), 
  " samples \n"
  )

message(
  "Metadata dim : ", 
  nrow(metadata), 
  " samples x ", 
  ncol(metadata), 
  " features\n"
  )

message("First gene IDs:\n", head(rownames(rna_counts)), "\n")
message("First sample IDs:\n", head(colnames(rna_counts)), "\n")

message("Metadata fields:\n", paste(colnames(metadata), collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# 3. basic checks
#-------------------------------------------------------------------------------
message("Are there NA values in RNA_counts?: ", anyNA(rna_counts), "\n")
message("Are there NA values in metadata?: ", anyNA(metadata), "\n")

message(
  "Expression value range: ", 
  paste(range(rna_counts, na.rm = TRUE), 
        collapse = " .. "), 
  "\n"
  )

message("RNA_counts summary:\n")
print(summary(rna_counts))

message("Example of expression values:\n")
print(rna_counts[1:5, 1:5])

message("Example of metadata:\n")
print(metadata[1:5, ])

#-------------------------------------------------------------------------------
# 4. consistency between rna_counts and metadata
#-------------------------------------------------------------------------------
message("Samples in RNA_counts: ", ncol(rna_counts), "\n",
        "Unique samples: ", length(unique(colnames(rna_counts))), "\n")

message("Samples in metadata: ", nrow(metadata), "\n",
        "Unique samples: ", length(unique(metadata$Sample.Name)), "\n")

message("Any duplicates in metadata$Sample.Name?: ",
        any(duplicated(metadata$Sample.Name)), "\n")

message("Any duplicated genes in RNA_counts rownames?: ",
        any(duplicated(rownames(rna_counts))), "\n")

message("Same sample set (RNA_counts vs metadata): ",
        setequal(colnames(rna_counts), metadata$Alias), "\n")

message("Same sample order: ",
        identical(colnames(rna_counts), metadata$Alias), "\n")

# reorder metadata if needed
if (!identical(colnames(rna_counts), metadata$Alias)) {
  metadata <- metadata[match(colnames(rna_counts), metadata$Alias), ]
  message("Metadata reordered to match rna_counts sample order.\n")
}

stopifnot(identical(colnames(rna_counts), metadata$Alias))

# rename count matrix with readable sample names
colnames(rna_counts) <- metadata$Sample.Name

stopifnot(identical(colnames(rna_counts), rownames(metadata)))

# save processed objects
dir.create(file.path(in_dir, "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rna_counts, file.path(in_dir, "processed", "rna_counts_clean.rds"))
saveRDS(metadata, file.path(in_dir, "processed", "rna_metadata_clean.rds"))

#-------------------------------------------------------------------------------
# 5. dataset composition
#-------------------------------------------------------------------------------
n_patients <- length(unique(metadata$Donor))
message("Number of donors:", n_patients, "\n")

# tables
samples_per_donor <- table(metadata$Donor)
celltype_table <- table(metadata$CellType)
donor_celltype <- table(metadata$Donor, metadata$CellType)

message("\nSamples per donor:\n")
print(samples_per_donor)

message("\nCell type distribution:\n")
print(celltype_table)

message("\nDonor x CellType table:\n")
print(donor_celltype)

# plots
# samples per donor
samples_per_donor_df <- metadata %>%
  dplyr::count(Donor, name = "n_samples")

samples_per_donor_plot <- ggplot(
  samples_per_donor_df,
  aes(x = Donor, y = n_samples)
) +
  geom_col() +
  labs(
    title = "Number of RNA-seq samples per donor",
    x = "Donor",
    y = "Number of samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "composition", "n_samples_per_donor.png"),
  plot = samples_per_donor_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# celltype distribution
celltype_plot <- metadata %>%
  dplyr::count(CellType, name = "n_samples") %>%
  ggplot(aes(x = CellType, y = n_samples)) +
  geom_col() +
  scale_x_discrete(labels = celltype_labels) +
  labs(
    title = "Distribution of RNA-seq samples by cell type",
    x = "Cell type",
    y = "Number of samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "composition", "cell_type_distribution.png"),
  plot = celltype_plot,
  width = 9,
  height = 6,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 6. raw count QC
#-------------------------------------------------------------------------------
# sequencing depth per sample
library_sizes <- colSums(rna_counts)

message("Library size summary:")
print(summary(library_sizes))

# number of genes with at least one read in each sample
detected_genes <- colSums(rna_counts > 0)

message("Detected genes summary:")
print(summary(detected_genes))

# fraction of zero values
zero_fraction <- mean(rna_counts == 0)

message(
  "Fraction of zero values: ",
  round(zero_fraction * 100, 2),
  "%"
)

message("Sample with minimum library size:")
print(library_sizes[which.min(library_sizes)])

message("Sample with minimum detected genes:")
print(detected_genes[which.min(detected_genes)])

qc_raw_df <- tibble(
  Sample.Name = colnames(rna_counts),
  library_size = library_sizes,
  detected_genes = detected_genes
) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType, Donor),
    by = "Sample.Name"
  )

# library size histogram
library_size_hist <- ggplot(qc_raw_df, aes(x = library_size)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of RNA-seq library sizes",
    x = "Total counts per sample",
    y = "Number of samples"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(out_dir, "qc_raw", "library_sizes_histogram.png"),
  plot = library_size_hist,
  width = 10,
  height = 6,
  dpi = 300
)

# detected genes histogram
detected_genes_hist <- ggplot(qc_raw_df, aes(x = detected_genes)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of detected genes per sample",
    x = "Genes with count > 0",
    y = "Number of samples"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(out_dir, "qc_raw", "detected_genes_histogram.png"),
  plot = detected_genes_hist,
  width = 10,
  height = 6,
  dpi = 300
)

# log2 raw counts boxplot
raw_counts_long <- log2(rna_counts + 1) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Sample.Name",
    values_to = "log2_count"
  ) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType),
    by = "Sample.Name"
  )

raw_counts_long <- raw_counts_long %>%
  mutate(
    Sample.Name = factor(Sample.Name, levels = colnames(rna_counts))
  )

raw_counts_boxplot <- ggplot(
  raw_counts_long,
  aes(x = Sample.Name, y = log2_count, fill = CellType)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  scale_fill_discrete(labels = celltype_labels) +
  labs(
    title = "Log2 raw counts across RNA-seq samples",
    x = "Sample",
    y = "log2(count + 1)",
    fill = "Cell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_dir, "qc_raw", "log2_raw_counts_boxplot.png"),
  plot = raw_counts_boxplot,
  width = 18,
  height = 8,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 7. DESeq2 normalization and transformed-data EDA
#-------------------------------------------------------------------------------
# create a DESeq2 object from raw counts and sample metadata for 
# exploratory analysis
dds <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = metadata,
  design = ~ Donor + CellType
)

# filter genes with low counts
keep_genes <- rowSums(counts(dds) >=10) >= 5

message("Genes before filtering: ", nrow(dds))
message("Genes retained after filtering: ", sum(keep_genes))
message("Genes removed by filtering: ", nrow(dds) - sum(keep_genes))

dds <- dds[keep_genes, ]

# estimate DESeq2 size factors to account for differences in sequencing depth
dds <- estimateSizeFactors(dds)

message("Size factor summary:")
print(summary(sizeFactors(dds)))

# apply variance stabilizing transformation.
# blind = true ignores the experimental design during transformation, for
# initial QC and outlier detection.
vsd <- vst(dds, blind = TRUE)

saveRDS(
  dds,
  file.path(out_dir, "qc_deseq2", "dds_filtered_size_factors.rds")
)

saveRDS(
  vsd,
  file.path(out_dir, "qc_deseq2", "vsd_blind_true.rds")
)

#-------------------------------------------------------------------------------
# 8. QC on variance-stabilized expression data
#-------------------------------------------------------------------------------
# prepare VST expression matrix for ggplot
vst_long <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Sample.Name",
    values_to = "vst_expression"
  ) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType),
    by = "Sample.Name"
  ) %>%
  mutate(
    Sample.Name = factor(Sample.Name, levels = colnames(rna_counts))
  )

# boxplot of VST expression values colored by cell type
vst_boxplot <- ggplot(
  vst_long,
  aes(
    x = Sample.Name,
    y = vst_expression,
    fill = CellType
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.2
  ) +
  scale_fill_discrete(labels = celltype_labels) +
  labs(
    title = "VST-transformed expression values",
    x = "Sample",
    y = "VST-transformed expression",
    fill = "Cell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    ),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_dir, "qc_deseq2", "vst_boxplot_by_celltype.png"),
  plot = vst_boxplot,
  width = 18,
  height = 8,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 9. PCA
#-------------------------------------------------------------------------------
# PCA colored by cell type
pca_celltype <- plotPCA(
  vsd, 
  intgroup = "CellType",
  ntop = nrow(vsd)
  )

pca_celltype <- pca_celltype +
  scale_color_discrete(labels = celltype_labels) +
  labs(
    title = "PCA of VST-transformed counts by cell type (RNA-seq)",
    x = pca_celltype$labels$x,
    y = pca_celltype$labels$y,
    color = "Cell type"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "pca", "pca_vst_celltype.png"),
  plot = pca_celltype,
  width = 10,
  height = 6,
  dpi = 300
)

# PCA colored by donor
pca_donor <- plotPCA(
  vsd, 
  intgroup = "Donor",
  ntop = nrow(vsd)
  )

pca_donor <- pca_donor +
  labs(
    title = "PCA of VST-transformed counts by donor (RNA-seq)",
    x = pca_donor$labels$x,
    y = pca_donor$labels$y,
    color = "Donor"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "pca", "pca_vst_donor.png"),
  plot = pca_donor,
  width = 10,
  height = 6,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 10. heatmap
#-------------------------------------------------------------------------------
# euclidean distances on VST transformed data
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

# annotation 
annotation_col <- metadata %>%
  dplyr::select(CellType, Donor)
rownames(annotation_col) <- metadata$Sample.Name

# translate
annotation_col$CellType <- celltype_labels[as.character(annotation_col$CellType)]

# use the same cell type colors in all heatmaps
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

annotation_colors <- list(CellType = celltype_colors, Donor = donor_colors)

sample_distance_heatmap <- pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  main = "Distances between RNA-seq samples (VST)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE,
  fontsize = 7,
  legend = TRUE
)

ggsave(
  filename = file.path(out_dir, "heatmaps", "sample_distance_heatmap_RNA.png"),
  plot = sample_distance_heatmap$gtable,
  width = 8,
  height = 9.5,
  dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

message("\nEDA finished. Outputs saved to: ", out_dir, "\n")
