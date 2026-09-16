################################################################################
# 05. Interindividual variability in ATAC-seq
# Laura Muñoz | 14/07/2026
#
# identifies chromatin regions with donor-associated variability
# and selects those that also change during Bcell differentiation.
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
out_dir <- "results/05_DE_ATAC_LRT"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

stopifnot(identical(colnames(atac_counts), metadata$Sample.Name))

#-------------------------------------------------------------------------------
# 3. create DESeqDataSet object
#-------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = atac_counts,
  colData = metadata,
  design = ~ CellType + Donor
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
# 5. run differential accessibility analysis with LRT
#-------------------------------------------------------------------------------
dds <- DESeq(
  dds,
  test = "LRT",
  reduced = ~ CellType
)

saveRDS(dds, file.path(out_dir, "dds_ATAC_LRT.rds"))

#-------------------------------------------------------------------------------
# 6. extract results
#-------------------------------------------------------------------------------
res_lrt <- results(dds, alpha = 0.01)

res_lrt_tbl <- res_lrt %>%
  as.data.frame() %>%
  rownames_to_column("OCR") %>%
  as_tibble() %>%
  arrange(padj)

message("Total OCRs tested: ", nrow(res_lrt_tbl))

message(
  "OCRs with adjusted p-value <= 0.01: ", 
  sum(res_lrt_tbl$padj <= 0.01, na.rm = TRUE)
  )

#-------------------------------------------------------------------------------
# 7. extract donor coefficients
#-------------------------------------------------------------------------------
# extract model coefficients and keep donor columns
coef_tbl <- coef(dds) %>%
  as.data.frame() %>%
  rownames_to_column("OCR") %>%
  as_tibble() %>%
  dplyr::select(OCR, starts_with("Donor_"))

# move the donor coefficients from columns to rows
coef_tbl_long <- coef_tbl %>%
  pivot_longer(
    cols = starts_with("Donor_"),
    names_to = "comparison",
    values_to = "donor_log2FC"
  ) %>%
  separate(
    col = comparison,
    into = c("variable", "Donor", "versus", "reference"),
    sep = "_"
  ) %>%
  dplyr::select(
    OCR,
    Donor,
    donor_log2FC
  )

# reference donor d200 has coefficient zero
reference_donor <- levels(dds$Donor)[1]

# create one reference-donor row for every OCR
reference_donor_tbl <- coef_tbl %>%
  transmute(
    OCR,
    Donor = reference_donor,
    donor_log2FC = 0
  )

# add the reference donor to the other coefficients
coef_tbl_long <- coef_tbl_long %>%
  bind_rows(reference_donor_tbl) %>%
  arrange(OCR, Donor)

#-------------------------------------------------------------------------------
# 8. calculate donor effect magnitude for each OCR
#-------------------------------------------------------------------------------
# group the 10 donor coefficients belonging to the same OCR
donor_effect_size_tbl <- coef_tbl_long %>%
  group_by(OCR) %>%
  summarise(
    # donors with lowest and highest estimated accessibility
    donor_min = Donor[which.min(donor_log2FC)],
    donor_max = Donor[which.max(donor_log2FC)],
    
    # lowest and highest donor coefficients
    min_donor_log2FC = min(donor_log2FC, na.rm = TRUE),
    max_donor_log2FC = max(donor_log2FC, na.rm = TRUE),
    
    # difference between the highest and lowest donor coefficients
    donor_log2FC_range = max_donor_log2FC - min_donor_log2FC,

    # general dispersion of the 10 donor coefficients
    donor_log2FC_sd = sd(donor_log2FC, na.rm = TRUE),
    .groups = "drop"
  )

# add the donor effect magnitude columns to the LRT results
res_lrt_tbl <- res_lrt_tbl %>%
  left_join(
    donor_effect_size_tbl,
    by = "OCR"
  )

#-------------------------------------------------------------------------------
# 9. extract significant OCRs
#-------------------------------------------------------------------------------
sig_ocrs <- res_lrt_tbl %>%
  filter(!is.na(padj), padj <= 0.01)

message("Significant OCRs associated with Donor effect: ", nrow(sig_ocrs))

#-------------------------------------------------------------------------------
# 10. select OCRs with sufficient donor effect magnitude
#-------------------------------------------------------------------------------
# minimum log2fc range equivalent to a two-fold difference
min_log2FC_range <- log2(2)

# select OCRs with statistical significance and sufficient magnitude
selected_ocrs <- sig_ocrs %>%
  filter(donor_log2FC_range >= min_log2FC_range)

message("Significant OCRs in LRT: ", nrow(sig_ocrs))
message(
  "Significant OCRs with sufficient magnitude: ",
  nrow(selected_ocrs)
)

#-------------------------------------------------------------------------------
# 11. save results
#-------------------------------------------------------------------------------
# save full results
saveRDS(res_lrt_tbl, file.path(
  out_dir, 
  "complete_results_DE_ATAC_individual.rds")
  )

# save significant results
saveRDS(sig_ocrs, file.path(
  out_dir,
  "DARs_ATAC_individual.rds")
  )

# save significant OCR names
sig_ocrs_id <- sig_ocrs$OCR

saveRDS(sig_ocrs_id, file.path(
  out_dir,
  "DARS_ids_ATAC_individual.rds")
)

# save OCRs with significant LRT and sufficient magnitude
saveRDS(selected_ocrs, file.path(
  out_dir,
  "selected_OCRs_ATAC_individual.rds")
)

selected_ocrs_id <- selected_ocrs$OCR

saveRDS(selected_ocrs_id, file.path(
  out_dir,
  "selected_OCRs_ids_ATAC_individual.rds")
)

#-------------------------------------------------------------------------------
# 12. unify dynamic and with interindividual variability OCRs
#-------------------------------------------------------------------------------
dynamic_ocrs_ids <- readRDS(
  "results/04_DE_ATAC_celltypes/dynamic_ocrs_ids_ATAC.rds"
)

dynamic_ocrs_ids <- unique(dynamic_ocrs_ids)

# extract common OCRs IDs (dynamic and with variability)
common_ocrs_ids <- intersect(dynamic_ocrs_ids, selected_ocrs_id)

# final list of results
common_ocrs <- selected_ocrs %>%
  filter(OCR %in% common_ocrs_ids)

message("Dynamic OCRs: ", length(dynamic_ocrs_ids))
message("Selected donor-variable OCRs: ", length(selected_ocrs_id))
message("Common OCRs: ", length(common_ocrs_ids))

# save final results
saveRDS(common_ocrs, file.path(out_dir, "common_ocrs_results.rds"))

#-------------------------------------------------------------------------------
# 13. PCA with selected donor-variable OCRs
#-------------------------------------------------------------------------------
# apply the VST transformation
vsd <- vst(dds, blind = FALSE)

# subset VST object to selected OCRs
vsd_dars <- vsd[rownames(vsd) %in% selected_ocrs_id, ]

# create PCA colored by donor
pca_donor <- plotPCA(
  vsd_dars,
  intgroup = "Donor",
  ntop = nrow(vsd_dars)
)

pca_donor <- pca_donor +
  labs(
    title = "PCA of selected donor-variable OCRs (ATAC-seq)",
    x = pca_donor$labels$x,
    y = pca_donor$labels$y,
    color = "Donor"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "PCA_DARs_ATAC_Individual_effect.png"),
  plot = pca_donor,
  width = 8,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 14. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
