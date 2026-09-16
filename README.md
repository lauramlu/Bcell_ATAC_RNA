# ATAC-seq and RNA-seq variability during B-cell differentiation

This repository contains the R scripts used in my master's thesis to study interindividual variability in chromatin accessibility and gene expression during human B-cell differentiation.

The analysis identifies changes between cell types, evaluates differences between donors and integrates ATAC-seq and RNA-seq through correlations between open chromatin regions (OCRs) and genes. The final steps combine the selected OCR-gene pairs with genomic annotations and explore their biological functions using GO and KEGG.

## Project structure

```text
TFM_Bcell_ATAC_RNA/
├── TFM_Bcell_ATAC_RNA.Rproj
├── README.md
├── .gitignore
├── scripts/
├── data/       # local input and processed data
└── results/    # generated tables, figures and R objects
```

The `data/` and `results/` folders are excluded from version control. The scripts create their output folders when needed.

## Requirements

The analysis was run using R 4.5.2. The scripts use the following packages:

- Data handling: `dplyr`, `tidyr` and `tibble`.
- Differential analysis: `DESeq2`.
- Figures: `ggplot2`, `pheatmap`, `enrichplot` and `scales`.
- Functional enrichment and gene annotation: `clusterProfiler` and `org.Hs.eg.db`.

To install the packages:

```r
install.packages(c("BiocManager", "dplyr", "tidyr", "tibble",
                   "ggplot2", "pheatmap", "scales"))

BiocManager::install(c("DESeq2", "clusterProfiler",
                      "org.Hs.eg.db", "enrichplot"))
```

Each script saves `sessionInfo.txt` with the R and package versions used for that run. KEGG analyses require internet access, and enrichment results can vary with annotation database versions.

## Input data

The workflow starts from count matrices, sample metadata and an existing OCR annotation. These input files are not included in the repository and must be available locally before running the analysis:

```text
data/
├── raw/
│   ├── rna_raw_data.rds
│   └── atac_raw_data.rds
├── metadata/
│   ├── rna_metadata.rds
│   └── atac_metadata.rds
└── processed/
    └── atac_OCR_annotation.rds
```

Count matrices contain genes or OCRs in rows and samples in columns. Sample metadata include `Sample.Name`, `Donor` and `CellType`; RNA metadata also include `Alias` to match the original count-matrix column names. The OCR annotation provides gene identifiers and genomic information used in scripts 06 and 08.

Scripts 01 and 02 check sample correspondence and save the processed count matrices and metadata for the following analyses. The workflow does not include read alignment, peak calling or generation of the input OCR annotation.

## Running the analysis

Open `TFM_Bcell_ATAC_RNA.Rproj` in RStudio and run the scripts in the order below. All paths are relative to the project folder, not to `scripts/`.

Use a fresh R session for each script to avoid conflicts between functions from different packages. Each script loads the files it needs from the previous steps. The scripts clear the R environment at the start.

| Script | Analysis |
| --- | --- |
| `01_exploratory_RNA.R` | RNA-seq quality control, filtering, VST, PCA and sample distances |
| `02_exploratory_ATAC.R` | ATAC-seq quality control, filtering, VST, PCA and sample distances |
| `03_differential_analysis_RNA_celltypes.R` | Differential expression between consecutive cell types |
| `04_differential_analysis_ATAC_celltypes.R` | Differential accessibility between consecutive cell types |
| `05_differential_analysis_ATAC_donor.R` | Donor-associated accessibility variability and selection of dynamic OCRs |
| `06_ATAC_RNA_spearman.R` | Spearman correlations between selected OCRs and differentially expressed genes |
| `07_differential_analysis_RNA_donor.R` | Donor-associated expression variability in correlated genes |
| `08_OCR_gene_table.R` | Final OCR-gene selection, annotation and summary figures |
| `09_biological_interpretation.R` | GO and KEGG enrichment analyses and selected GO lollipop plots |

Scripts 03 and 04 account for donor differences when comparing cell types. Scripts 05 and 07 use a likelihood ratio test (LRT) to evaluate donor effects while accounting for cell type. In script 07, Benjamini-Hochberg correction is applied within the set of correlated genes. Thresholds and selection criteria are specified in each script.

## Results

Results are saved in numbered subfolders of `results/`. The final annotated table is written to:

```text
results/08_final_OCR_gene_table/ocr_gene_table.csv
results/08_final_OCR_gene_table/ocr_gene_table.rds
```

Script 09 reads this table and saves enrichment results, selected-term tables and the two GO lollipop plots in `results/09_biological_interpretation/`. Its selected terms correspond to those used in the thesis; the script currently stops if any of the three selected KEGG pathways are absent from the GSEA results.

Rerunning a script overwrites files with the same names. Intermediate R objects are needed by later steps and should remain available locally.

## Author

Laura Muñoz
