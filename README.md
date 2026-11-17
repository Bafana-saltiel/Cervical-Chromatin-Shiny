# Cervical cancer chromatin explorer

R Shiny explorer for 61 gene loci across SiHa, Ca Ski, and C-33 A. Select a review panel, gene, and cell line to rotate the local 3D model, compare RNA expression and CRISPR gene-effect scores, and download the plotted values.

## Run locally

The repository bundles the original source, tables, scripts, and assets in `chromatin-source.tar.gz`. The root `app.R` extracts and runs that bundle. To develop the app, unpack the archive into a separate directory.

To run locally, open R in the repository directory and run:

```r
install.packages("shiny")
shiny::runApp()
```

## Deploy to Posit Connect Cloud

Push this directory to a GitHub repository with `app.R` at the repository root. In Posit Connect Cloud, open **Publish**, select **Shiny**, and choose the GitHub repository, branch, and `app.R` entry point. Copy the app URL from the published content page. The account profile lists published content separately.

## Data and methods

- Hi-C contact maps are from [GSE195631](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195631). Individual source URLs are in `data/cell_line_sources.csv`. The bundled coordinates are local GRCh38 25 kb cis-contact embeddings from the cervical SCC review-panel analysis. Classical multidimensional scaling provides an independently inferred geometry for each locus and line. Axes and distances are arbitrary. Rotation and reflection can change between models. Use the rank correlation table for the descriptive comparison of contacts within distance bands.
- [DepMap Public 24Q4](https://doi.org/10.25452/figshare.plus.27993248.v1) provides gene expression in log2(TPM + 1) and CRISPR gene-effect scores. Inputs and checksums are recorded in `data/depmap_sources.csv`. The model joins use ACH identifiers and the features use gene symbols. Blank scores are missing. RNA values cover all 61 loci per line and CRISPR values cover 59 per line.
- From this directory, run `Rscript scripts/depmap.R` to rebuild the compact DepMap tables from the pinned release. The raw matrix cache requires about 900 MB. The running app reads only bundled tables.
- The three cell-line Hi-C maps and DepMap assays are independent observations. Use the app for exploratory comparisons across cell lines. Patient mRNA, protein, treatment and outcome data are outside the viewer's scope.

The public Hi-C contact windows and locus definitions are bundled under `data/`. Run `Rscript scripts/geometry.R` and `Rscript scripts/concordance.R` from the app directory to reproduce the coordinates, diagnostics and contact comparisons.
