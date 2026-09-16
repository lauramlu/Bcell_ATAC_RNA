# ATAC-seq and RNA-seq variability during B-cell differentiation

This repository contains the R scripts used in my master's thesis to study interindividual variability in chromatin accessibility and gene expression during human B-cell differentiation. The analysis combines ATAC-seq and RNA-seq data to identify donor-associated variability and select correlated OCR-gene pairs.

## Data source

This work reanalyses matched RNA-seq and ATAC-seq data from the study by Planell et al. on early human B-cell differentiation (https://doi.org/10.1126/sciadv.adw3110), covering eight cell populations along the B-cell differentiation lineage. Input count matrices, sample metadata and OCR annotations are required to run the scripts and are not included in this repository.
Data: https://osf.io/gswpy/

## Running the analysis

Open `TFM_Bcell_ATAC_RNA.Rproj` in RStudio and run the scripts in `scripts/` in numerical order, from 01 to 09, using a fresh R session for each script. Place the input files in `data/`; the scripts create their output folders under `results/`. Required packages are listed at the beginning of each script.

## Selected results

The [selected_results](selected_results/) folder includes the final figures, enrichment tables and R session information. The complete [OCR-gene table](selected_results/ocr_gene_table.csv) contains the 1,888 selected pairs, with their genomic annotations and statistical results.

## Author

Laura Muñoz
