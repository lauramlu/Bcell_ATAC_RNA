################################################################################
# 09. Biological interpretation
# Laura Muñoz | 23/07/2026
#
# explores biological functions and pathways using
# GO and KEGG enrichment analyses.
################################################################################

#-------------------------------------------------------------------------------
# 0. set up
#-------------------------------------------------------------------------------
rm(list = ls())
gc()

message("R version: ", R.version.string)
message("Working directory: ", getwd())

library(ggplot2)
library(dplyr)
library(tibble)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

in_dir <- "results/08_final_OCR_gene_table"
in_dir_GSEA <- "results/07_DE_RNA_LRT"
out_dir <- "results/09_biological_interpretation"

# create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# 1. read data
#-------------------------------------------------------------------------------
ocr_gene_table <- readRDS(
  file.path(
    in_dir,
    "ocr_gene_table.rds"
  )
)

message("Final OCR-gene pairs: ", nrow(ocr_gene_table))
message("Unique final genes: ",n_distinct(ocr_gene_table$ENSEMBL))

dds_rna_lrt <- readRDS(
  file.path(in_dir_GSEA, "dds_RNA_LRT.rds")
)

#-------------------------------------------------------------------------------
# 2. GO biological process ORA
#-------------------------------------------------------------------------------
selected_gene_ids <- unique(ocr_gene_table$ENSEMBL)

message("Genes included in GO analysis: ", length(selected_gene_ids))

ora_bp <- enrichGO(
  gene = selected_gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

ora_bp_results <- as.data.frame(ora_bp) %>%
  as_tibble() %>%
  arrange(p.adjust)

message(
  "Significant GO Biological Process terms: ",
  nrow(ora_bp_results)
)

print(
  ora_bp_results %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      p.adjust,
      Count
    ) %>%
    dplyr::slice_head(n = 10)
)

saveRDS(
  ora_bp,
  file.path(out_dir, "GO_BP_ORA_results.rds")
)

write.csv(
  ora_bp_results,
  file.path(out_dir, "GO_BP_ORA_results.csv"),
  row.names = FALSE
)

selected_ora_terms <- c(
  "antigen receptor-mediated signaling pathway",
  "B cell activation",
  "B cell proliferation",
  "positive regulation of lymphocyte activation",
  "leukocyte chemotaxis",
  "pattern recognition receptor signaling pathway",
  "activation of innate immune response",
  "T cell differentiation"
)

selected_ora_results <- ora_bp_results %>%
  dplyr::filter(Description %in% selected_ora_terms) %>%
  dplyr::mutate(
    Description = factor(
      Description,
      levels = rev(selected_ora_terms)
    )
  )

selected_ora_table <- selected_ora_results %>%
  dplyr::mutate(
    term_order = match(
      as.character(Description),
      selected_ora_terms
    )
  ) %>%
  dplyr::arrange(term_order) %>%
  dplyr::transmute(
    `GO ID` = ID,
    `GO-BP term` = as.character(Description),
    `Gene ratio` = GeneRatio,
    `Fold enrichment` = round(FoldEnrichment, 3),
    Genes = Count,
    FDR = signif(p.adjust, 3)
  )

write.csv(
  selected_ora_table,
  file.path(out_dir, "GO_BP_ORA_selected_terms_table.csv"),
  row.names = FALSE
)

selected_go_lollipop <- ggplot(
  selected_ora_results,
  aes(
    x = FoldEnrichment,
    y = Description
  )
) +
  geom_segment(
    aes(
      x = 1,
      xend = FoldEnrichment,
      yend = Description
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      label = paste0(
        "Genes = ", Count,
        "   FDR = ", formatC(p.adjust, format = "e", digits = 1)
      )
    ),
    nudge_x = 0.12,
    hjust = 0,
    size = 5
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.80))
  ) +
  scale_y_discrete(labels = scales::label_wrap(35)) +
  labs(
    title = "Selected biological processes - ORA",
    x = "Fold enrichment",
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 15, lineheight = 0.95),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    plot.title = element_text(size = 17),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  filename = file.path(
    out_dir,
    "GO_BP_ORA_selected_terms_lollipop.png"
  ),
  plot = selected_go_lollipop,
  width = 10,
  height = 7.5,
  dpi = 300
)


#-------------------------------------------------------------------------------
# 3. convert ENSEMBL identifiers to ENTREZID
#-------------------------------------------------------------------------------
gene_id_mapping <- bitr(
  selected_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

selected_entrez_ids <- unique(gene_id_mapping$ENTREZID)

message(
  "Genes with ENTREZID included in KEGG analysis: ",
  length(selected_entrez_ids)
)

write.csv(
  gene_id_mapping,
  file.path(out_dir, "ENSEMBL_to_ENTREZID.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 4. KEGG pathway over-representation analysis
#-------------------------------------------------------------------------------
ora_kegg <- enrichKEGG(
  gene = selected_entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

ora_kegg_results <- as.data.frame(ora_kegg) %>%
  as_tibble() %>%
  arrange(p.adjust)

message(
  "Significant KEGG pathways: ",
  nrow(ora_kegg_results)
)

print(
  ora_kegg_results %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      p.adjust,
      Count
    ) %>%
    dplyr::slice_head(n = 10)
)

saveRDS(
  ora_kegg,
  file.path(out_dir, "KEGG_ORA_results.rds")
)

write.csv(
  ora_kegg_results,
  file.path(out_dir, "KEGG_ORA_results.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 5. GSEA: create the ranked gene list
#-------------------------------------------------------------------------------
lrt_results <- results(
  dds_rna_lrt,
  independentFiltering = FALSE,
  pAdjustMethod = "none"
)

lrt_results_tbl <- lrt_results %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  as_tibble() %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::mutate(
    # negative LRT statistics have p-value 1 and are treated as zero evidence.
    ranking_stat = pmax(stat, 0)
  ) %>%
  dplyr::arrange(dplyr::desc(ranking_stat))

gene_list_go <- lrt_results_tbl$ranking_stat
names(gene_list_go) <- lrt_results_tbl$ENSEMBL
gene_list_go <- sort(gene_list_go, decreasing = TRUE)

write.csv(
  lrt_results_tbl,
  file.path(out_dir, "RNA_LRT_gene_ranking.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 6. GSEA using gene Ontology
#-------------------------------------------------------------------------------
# scoretype = "pos" tests enrichment towards genes with greater donor
# variability. it does not represent pathway activation or repression.
gsea_bp <- gseGO(
  geneList = gene_list_go,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  seed = TRUE,
  scoreType = "pos"
)

write.csv(
  as.data.frame(gsea_bp),
  file.path(out_dir, "GSEA_GO_BP_results.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 7. convert ENSEMBL identifiers to ENTREZID for GSEA
#-------------------------------------------------------------------------------
gene_id_mapping <- bitr(
  lrt_results_tbl$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  # remove ENSEMBL identifiers associated with more than one ENTREZID.
  dplyr::group_by(ENSEMBL) %>%
  dplyr::filter(dplyr::n() == 1) %>%
  dplyr::ungroup()

kegg_ranking_tbl <- lrt_results_tbl %>%
  dplyr::select(ENSEMBL, ranking_stat) %>%
  dplyr::inner_join(gene_id_mapping, by = "ENSEMBL") %>%
  dplyr::arrange(dplyr::desc(ranking_stat)) %>%
  # KEGG requires unique ENTREZID identifiers.
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

gene_list_kegg <- kegg_ranking_tbl$ranking_stat
names(gene_list_kegg) <- kegg_ranking_tbl$ENTREZID
gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)

write.csv(
  kegg_ranking_tbl,
  file.path(out_dir, "KEGG_gene_ranking.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# 8. GSEA using KEGG
#-------------------------------------------------------------------------------
gsea_kegg <- gseKEGG(
  geneList = gene_list_kegg,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  seed = TRUE,
  scoreType = "pos"
)

write.csv(
  as.data.frame(gsea_kegg),
  file.path(out_dir, "GSEA_KEGG_results.csv"),
  row.names = FALSE
)

selected_gsea_kegg_terms <- c(
  "Antigen processing and presentation",
  "Intestinal immune network for IgA production",
  "Cell adhesion molecule (CAM) interaction"
)

selected_gsea_kegg_results <- as.data.frame(gsea_kegg) %>%
  as_tibble() %>%
  dplyr::filter(Description %in% selected_gsea_kegg_terms) %>%
  dplyr::mutate(
    core_count = lengths(
      strsplit(core_enrichment, "/", fixed = TRUE)
    ),
    term_order = match(Description, selected_gsea_kegg_terms)
  ) %>%
  dplyr::arrange(term_order)

if (nrow(selected_gsea_kegg_results) != length(selected_gsea_kegg_terms)) {
  stop("Not all selected GSEA KEGG pathways were found in the results.")
}

selected_gsea_kegg_table <- selected_gsea_kegg_results %>%
  dplyr::transmute(
    `KEGG ID` = ID,
    `KEGG pathway` = Description,
    NES = round(NES, 3),
    `Gene set size` = setSize,
    `Leading-edge genes` = core_count,
    FDR = signif(p.adjust, 3)
  )

write.csv(
  selected_gsea_kegg_table,
  file.path(out_dir, "GSEA_KEGG_selected_terms_table.csv"),
  row.names = FALSE
)

# save the complete GSEA objects for figures
saveRDS(gsea_bp, file.path(out_dir, "GSEA_GO_BP_results.rds"))
saveRDS(gsea_kegg, file.path(out_dir, "GSEA_KEGG_results.rds"))

#-------------------------------------------------------------------------------
# 9. GSEA figures
#-------------------------------------------------------------------------------
selected_gsea_terms <- c(
  "adaptive immune response",
  "regulation of lymphocyte activation",
  "lymphocyte mediated immunity",
  paste0(
    "immune response-activating cell surface receptor ",
    "signaling pathway"
  ),
  "antigen processing and presentation",
  paste0(
    "adaptive immune response based on somatic recombination of immune ",
    "receptors built from immunoglobulin superfamily domains"
  ),
  "regulation of hemopoiesis",
  "response to virus"
)

selected_gsea_results <- as.data.frame(gsea_bp) %>%
  as_tibble() %>%
  dplyr::filter(Description %in% selected_gsea_terms) %>%
  dplyr::mutate(
    core_count = lengths(
      strsplit(core_enrichment, "/", fixed = TRUE)
    ),
    Description = factor(
      Description,
      levels = rev(selected_gsea_terms)
    )
  )

selected_gsea_table <- selected_gsea_results %>%
  dplyr::mutate(
    term_order = match(
      as.character(Description),
      selected_gsea_terms
    )
  ) %>%
  dplyr::arrange(term_order) %>%
  dplyr::transmute(
    `GO ID` = ID,
    `GO-BP term` = as.character(Description),
    NES = round(NES, 3),
    `Gene set size` = setSize,
    `Leading-edge genes` = core_count,
    FDR = signif(p.adjust, 3)
  )

write.csv(
  selected_gsea_table,
  file.path(out_dir, "GSEA_GO_BP_selected_terms_table.csv"),
  row.names = FALSE
)

selected_gsea_lollipop <- ggplot(
  selected_gsea_results,
  aes(
    x = NES,
    y = Description
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = Description
    )
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      label = paste0(
        "Leading edge = ", core_count,
        "   FDR = ", formatC(p.adjust, format = "e", digits = 1)
      )
    ),
    nudge_x = 0.10,
    hjust = 0,
    size = 5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 1.30))
  ) +
  scale_y_discrete(labels = scales::label_wrap(38)) +
  labs(
    title = "Selected biological processes - GSEA",
    x = "NES",
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 15, lineheight = 0.95),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    plot.title = element_text(size = 17),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(
  filename = file.path(
    out_dir,
    "GSEA_GO_BP_selected_terms_lollipop.png"
  ),
  plot = selected_gsea_lollipop,
  width = 10,
  height = 9,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 10. save session information
#-------------------------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
