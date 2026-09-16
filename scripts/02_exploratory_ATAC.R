################################################################################
# 02. Exploratory ATAC-seq analysis
# Laura Muñoz | 12/05/2026
#
# checks ATAC-seq data quality and explores sample variation
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
out_dir <- "results/02_exploratory_ATAC"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (subdir in c("composition", "qc_raw", "qc_deseq2", "pca", "heatmaps")) {
  dir.create(file.path(out_dir, subdir), showWarnings = FALSE, recursive = TRUE)
}
dir.create(file.path(in_dir, "processed"), showWarnings = FALSE, recursive = TRUE)

celltype_labels <- c(
  HSC = "HSPC", CLP = "CLP", proB = "pro-B", preB = "pre-B",
  ImmatureB = "Immature B", Transitional_B = "Transitional B",
  Naive_CD5pos = "Naive B CD5+", Naive_CD5neg = "Naive B CD5-"
)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
atac_counts <- readRDS(file.path(in_dir, "raw", "atac_raw_data.rds"))
metadata <- readRDS(file.path(in_dir, "metadata", "atac_metadata.rds"))

#-------------------------------------------------------------------------------
# 2. object type and dimensions
#-------------------------------------------------------------------------------
message("ATAC_counts class:", class(atac_counts), "\n")
message("Metadata class:", class(metadata), "\n")

message("ATAC_counts dim: ", nrow(atac_counts), " OCRs x ", ncol(atac_counts), " samples \n")
message("Metadata dim : ", nrow(metadata), " samples x ", ncol(metadata), " features\n")

message("First OCR IDs:\n", head(rownames(atac_counts)), "\n")
message("First sample IDs:\n", head(colnames(atac_counts)), "\n")

message("Metadata fields:\n", paste(colnames(metadata), collapse = ", "), "\n")

#-------------------------------------------------------------------------------
# 3. basic checks
#-------------------------------------------------------------------------------
message("Are there NA values in ATAC_counts?: ", anyNA(atac_counts), "\n")
message("Are there NA values in metadata?: ", anyNA(metadata), "\n")

message("ATAC count range: ", paste(range(atac_counts, na.rm = TRUE), collapse = " .. "), "\n")

message("ATAC_counts summary:\n")
print(summary(atac_counts))

message("Example of ATAC count values:\n")
print(atac_counts[1:5, 1:5])

message("Example of metadata:\n")
print(metadata[1:5, ])

#-------------------------------------------------------------------------------
# 4. consistency between ATAC_counts and metadata
#-------------------------------------------------------------------------------
message("Samples in ATAC_counts: ", ncol(atac_counts), "\n",
        "Unique samples: ", length(unique(colnames(atac_counts))), "\n")

message("Samples in metadata: ", nrow(metadata), "\n",
        "Unique samples: ", length(unique(metadata$Sample.Name)), "\n")

message("Any duplicates in metadata$Sample.Name?: ",
        any(duplicated(metadata$Sample.Name)), "\n")

message("Any duplicated OCRs in atac_counts rownames?: ",
        any(duplicated(rownames(atac_counts))), "\n")

message("Same sample set (ATAC_counts vs metadata): ",
        setequal(colnames(atac_counts), rownames(metadata)), "\n")

message("Same sample order: ",
        identical(colnames(atac_counts), rownames(metadata)), "\n")

if (!setequal(colnames(atac_counts), rownames(metadata))) {
  stop("Sample IDs in ATAC_counts do not match metadata rownames.")
}

# reorder metadata if needed
if (!identical(colnames(atac_counts), rownames(metadata))) {
  metadata <- metadata[match(colnames(atac_counts), rownames(metadata)), ]
  message("Metadata reordered to match atac_counts sample order.\n")
}

# keep original ATAC sample names, matching RNA sample naming convention
metadata <- metadata[match(colnames(atac_counts), rownames(metadata)), ]

# use ATAC matrix column names as sample identifiers
metadata$Sample.Name <- colnames(atac_counts)
rownames(metadata) <- colnames(atac_counts)

stopifnot(identical(colnames(atac_counts), metadata$Sample.Name))

# save processed objects
saveRDS(atac_counts, file.path(in_dir, "processed", "atac_counts_clean.rds"))
saveRDS(metadata, file.path(in_dir, "processed", "atac_metadata_clean.rds"))

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
    title = "Number of ATAC-seq samples per donor",
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
    title = "Distribution of ATAC-seq samples by cell type",
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
library_sizes <- colSums(atac_counts)

message("Library size summary:")
print(summary(library_sizes))

# number of OCRs with at least one read in each sample
detected_ocrs <- colSums(atac_counts > 0)

message("Detected OCRs summary:")
print(summary(detected_ocrs))

# fraction of zero values
zero_fraction <- mean(atac_counts == 0)

message(
  "Fraction of zero values: ",
  round(zero_fraction * 100, 2),
  "%"
)

message("Sample with minimum library size:")
print(library_sizes[which.min(library_sizes)])

message("Sample with minimum detected OCRs:")
print(detected_ocrs[which.min(detected_ocrs)])

qc_raw_df <- tibble(
  Sample.Name = colnames(atac_counts),
  library_size = library_sizes,
  detected_ocrs = detected_ocrs
) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType, Donor),
    by = "Sample.Name"
  )

# library size histogram
library_size_hist <- ggplot(qc_raw_df, aes(x = library_size)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of ATAC-seq library sizes",
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

# detected OCRs histogram
detected_ocrs_hist <- ggplot(qc_raw_df, aes(x = detected_ocrs)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of detected OCRs per sample",
    x = "OCRs with count > 0",
    y = "Number of samples"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(out_dir, "qc_raw", "detected_ocrs_histogram.png"),
  plot = detected_ocrs_hist,
  width = 10,
  height = 6,
  dpi = 300
)

# library size by cell type
library_size_celltype_plot <- ggplot(
  qc_raw_df,
  aes(x = CellType, y = library_size)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  scale_x_discrete(labels = celltype_labels) +
  labs(
    title = "ATAC-seq library sizes by cell type",
    x = "Cell type",
    y = "Total counts per sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "qc_raw", "library_sizes_by_celltype.png"),
  plot = library_size_celltype_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# detected OCRs by cell type
detected_ocrs_celltype_plot <- ggplot(
  qc_raw_df,
  aes(x = CellType, y = detected_ocrs)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  scale_x_discrete(labels = celltype_labels) +
  labs(
    title = "Detected OCRs by cell type",
    x = "Cell type",
    y = "OCRs with count > 0"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "qc_raw", "detected_ocrs_by_celltype.png"),
  plot = detected_ocrs_celltype_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# log2 raw counts boxplot
raw_counts_long <- log2(atac_counts + 1) %>%
  as.data.frame() %>%
  rownames_to_column("ocr_id") %>%
  pivot_longer(
    cols = -ocr_id,
    names_to = "Sample.Name",
    values_to = "log2_count"
  ) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType),
    by = "Sample.Name"
  )

raw_counts_long <- raw_counts_long %>%
  mutate(
    Sample.Name = factor(Sample.Name, levels = colnames(atac_counts))
  )


raw_counts_boxplot <- ggplot(
  raw_counts_long,
  aes(x = Sample.Name, y = log2_count, fill = CellType)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  scale_fill_discrete(labels = celltype_labels) +
  labs(
    title = "Log2 raw counts across ATAC-seq samples",
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
  filename = file.path(out_dir, "qc_raw", "log2_raw_counts_ATAC_boxplot.png"),
  plot = raw_counts_boxplot,
  width = 18,
  height = 8,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 7. DESeq2 normalization and transformed-data EDA
#-------------------------------------------------------------------------------
# create a DESeq2 object from raw OCR counts and sample metadata for 
# exploratory analysis
dds <- DESeqDataSetFromMatrix(
  countData = atac_counts,
  colData = metadata,
  design = ~ Donor + CellType
)

# filter OCRs with low counts
keep_ocrs <- rowSums(counts(dds) >=10) >= 5

message("OCRs before filtering: ", nrow(dds))
message("OCRs retained after filtering: ", sum(keep_ocrs))
message("OCRs removed by filtering: ", nrow(dds) - sum(keep_ocrs))

dds <- dds[keep_ocrs, ]

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
  file.path(out_dir, "qc_deseq2", "dds_filtered_size_factors_ATAC.rds")
)

saveRDS(
  vsd,
  file.path(out_dir, "qc_deseq2", "vsd_blind_true_ATAC.rds")
)

#-------------------------------------------------------------------------------
# 8. QC on variance-stabilized ATAC data
#-------------------------------------------------------------------------------
# prepare VST ATAC matrix for ggplot
vst_long <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("OCR_id") %>%
  pivot_longer(
    cols = -OCR_id,
    names_to = "Sample.Name",
    values_to = "vst_accessibility"
  ) %>%
  left_join(
    metadata %>% dplyr::select(Sample.Name, CellType),
    by = "Sample.Name"
  ) %>%
  mutate(
    Sample.Name = factor(Sample.Name, levels = colnames(atac_counts))
  )

# boxplot of VST ATAC values colored by cell type
vst_boxplot <- ggplot(
  vst_long,
  aes(
    x = Sample.Name,
    y = vst_accessibility,
    fill = CellType
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.2
  ) +
  scale_fill_discrete(labels = celltype_labels) +
  labs(
    title = "VST-transformed accessibility values",
    x = "Sample",
    y = "VST-transformed accessibility",
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
  filename = file.path(out_dir, "qc_deseq2", "vst_boxplot_by_celltype_ATAC.png"),
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
    title = "PCA of VST-transformed counts by cell type (ATAC-seq)",
    x = pca_celltype$labels$x,
    y = pca_celltype$labels$y,
    color = "Cell type"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "pca", "pca_vst_celltype_ATAC.png"),
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
    title = "PCA of VST-transformed counts by donor (ATAC-seq)",
    x = pca_donor$labels$x,
    y = pca_donor$labels$y,
    color = "Donor"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "pca", "pca_vst_donor_ATAC.png"),
  plot = pca_donor,
  width = 10,
  height = 6,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 10. heatmap
#-------------------------------------------------------------------------------
# euclidean distances on VST
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
  main = "Distances between ATAC-seq samples (VST)",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE,
  fontsize = 7,
  legend = TRUE
)

ggsave(
  filename = file.path(out_dir, "heatmaps", "sample_distance_heatmap_ATAC.png"),
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
