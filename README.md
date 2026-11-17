# Cervical cancer chromatin explorer

R Shiny explorer for all DepMap 24Q4 Cervix cell lines and measured gene symbols. Local 3D reconstructions are computed on demand for annotated genes with sufficient Hi-C coverage in SiHa, Ca Ski, C-33 A and HeLa. Upload a normalized contact CSV to reconstruct a locus in another cell line. Select a review panel, gene, and cell line to rotate the local 3D model, compare RNA expression and CRISPR gene-effect scores, and download the plotted values.

## Run locally

Extract every repository archive into a separate folder to inspect the source. The root loader uses the bundled archives directly. From the repository directory in R, run:

```r
install.packages("shiny")
shiny::runApp()
```

## Deploy to Posit Connect Cloud

Push this directory to a GitHub repository with `app.R` at the repository root. In Posit Connect Cloud, open **Publish**, select **Shiny**, and choose the GitHub repository, branch, and `app.R` entry point. Copy the app URL from the published content page. The account profile lists published content separately.

## Data and methods

- All four Hi-C maps are 25 kb GRCh38 KR-normalized cis contacts from [GSE195631](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195631). Accessions are listed in `data/hic_models.csv`. The app stores contact values at offsets up to 1 Mb as compressed chromosome bands. For each selected gene it reads the ±500 kb window around the annotated transcription start site and fits a local geometry by classical multidimensional scaling. Axes and distances are arbitrary.
- Transcription start sites use UCSC's GRCh38 [RefSeq Select and Curated RefSeq tracks](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/). RefSeq Select takes precedence, followed by curated RefSeq with protein-coding transcripts preferred. The chosen transcript and source are shown with the result. Sex chromosomes or sparse regions may lack enough usable contacts.
- For a model with no bundled map, upload a CSV for the selected gene and line. Required columns are `chr,start1,start2,balanced`. Use one row per pair, zero-based GRCh38 25 kb bin starts, and normalized positive contact values from the matching cell line. Uploads are used by the current Shiny session. The concordance table compares the bundled maps.
- [DepMap Public 24Q4](https://doi.org/10.25452/figshare.plus.27993248.v1) provides gene expression in log2(TPM + 1) and CRISPR gene-effect scores. Inputs and checksums are recorded in `data/depmap_sources.csv`. The model joins use ACH identifiers and the features preserve the release gene symbols and Entrez IDs. The model list follows DepMap lineage and model-type metadata. Assay coverage varies by model and gene. Blank scores are missing.
- From this directory, run `Rscript scripts/depmap.R` to rebuild the compact DepMap tables from the pinned release. The raw matrix cache requires about 900 MB. The running app reads only bundled tables.
- The three cell-line Hi-C maps and DepMap assays are independent observations. Use the app for exploratory comparisons across cell lines. Patient mRNA, protein, treatment and outcome data are outside the viewer's scope.

Rebuild the contact bands with `Rscript scripts/hic_bands.R MODEL path/to/map.hic` (requires strawr) and the gene annotation with `Rscript scripts/gene_loci.R ncbiRefSeqSelect.txt.gz ncbiRefSeqCurated.txt.gz`. The smaller original review-panel data and scripts remain available for reproducing the initial figures.

## Patient cohort tab

Upload a CSV containing normalized local Hi-C contacts for one or more coded patients. Use columns `patient_id,chr,start1,start2,balanced` with zero-based GRCh38 25 kb or hg19 40 kb bin starts and one record per contact pair per patient. Use the same normalization method across the cohort. The hosted app uses the file in the current Shiny session and supports uploads up to 5 MB. Exclude identifiers and identifying clinical information.

Select a gene and patient to view the fitted structure and download coordinates. The summary CSV shows coverage and fit diagnostics for all patients at the selected locus.

## Public cohorts

- **TCGA-CESC.** Gene-level quartiles by histology for primary tumours with data in the [UCSC Xena HiSeqV2 expression matrix](https://xenabrowser.net/datapages/?dataset=TCGA.CESC.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net). The build script matches the [Xena phenotype file](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.CESC.sampleMap%2FCESC_clinicalMatrix) by sample ID and excludes metastatic and normal samples. Only group summaries are bundled. Run `Rscript scripts/tcga_expression.R XENA_EXPRESSION.tsv.gz XENA_CLINICAL.tsv` to regenerate.
- **CRA001401.** [GSA study](https://ngdc.cncb.ac.cn/gsa/browse/CRA001401) and [reported TAD and sequencing QC counts](https://doi.org/10.3389/fcell.2021.706375). Study labels use experiment accessions because the archive and article number the two tumours differently. Bundled summaries contain no normalized contact matrix. The Patient cohort tab can reconstruct uploaded hg19 40 kb matrices using the UCSC hg19 RefSeq Select track.
