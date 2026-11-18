library(shiny)

# Connect Cloud can keep application source read-only.
source_dir <- tempfile("chromatin-source-")
dir.create(source_dir)
archives <- c("chromatin-source.tar.gz", list.files(pattern = "^hic-.*[.]tar[.]gz$"))
for (archive in archives) {
  stopifnot(utils::untar(archive, exdir = source_dir) == 0L)
}

# Extend the locus overview in the unpacked source.
reanalysis_file <- file.path(source_dir, "R", "chromatin_reanalysis.R")
reanalysis <- readLines(reanalysis_file, warn = FALSE)
# Add a public, bundled cervical study to the existing reanalysis data sources.
cervical_zip <- "GSE195631_RAD51B_HiC_upload.zip"
cervical_dir <- file.path(source_dir, "data", "GSE195631")
dir.create(cervical_dir, recursive = TRUE, showWarnings = FALSE)
cervical_files <- c("RAD51B_cervical_contacts_25kb.csv.gz", "samples_GRCh38.csv",
                    "RAD51B_GRCh38.csv", "HeLa_H3K27ac_GRCh38.bed")
utils::unzip(cervical_zip, files = cervical_files, exdir = cervical_dir)
stopifnot(all(file.exists(file.path(cervical_dir, cervical_files))))
replace_reanalysis_line <- function(lines, old, replacement) {
  i <- which(trimws(lines) == old)
  stopifnot(length(i) == 1L)
  append(lines[-i], replacement, after = i - 1L)
}
reanalysis <- replace_reanalysis_line(reanalysis,
  'radioButtons("re_source","Data source", c("HeLa-S3 study","Upload study"), inline=TRUE),',
  c('      radioButtons("re_source", "Data source", c("HeLa-S3 study", "Cervical RAD51B (GSE195631)", "Upload study"), inline=TRUE),',
    '      conditionalPanel("input.re_source === \'Cervical RAD51B (GSE195631)\'",',
    '        p("Five baseline cervical cell-line maps. GRCh38 chr14:65.25–70.50 Mb at 25 kb. The shared overlay is HeLa H3K27ac, and compartment scores are unphased.")),'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'ex <- identical(input$re_source,"HeLa-S3 study")',
  c('    ex <- identical(input$re_source, "HeLa-S3 study")',
    '    cervical <- identical(input$re_source, "Cervical RAD51B (GSE195631)")',
    '    cervical_path <- function(file) file.path(reanalysis_root, "data", "GSE195631", file)'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'samples <- read.csv(experiment_example_samples_path)',
  c('        samples <- read.csv(experiment_example_samples_path)',
    '      } else if (cervical) {',
    '        contacts <- read.csv(cervical_path("RAD51B_cervical_contacts_25kb.csv.gz"))',
    '        samples <- read.csv(cervical_path("samples_GRCh38.csv"))'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'g <- if(ex) read.csv(example) else {',
  c('      g <- if (ex) read.csv(example) else if (cervical) {',
    '        read.csv(cervical_path("RAD51B_GRCh38.csv"))',
    '      } else {'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'if (!ex && !is.null(input$re_phase)) {',
  '      if (!ex && !cervical && !is.null(input$re_phase)) {')
reanalysis <- replace_reanalysis_line(reanalysis,
  'if (!ex && !is.null(input$re_features)) {',
  '      if (cervical || (!ex && !is.null(input$re_features))) {')
reanalysis <- replace_reanalysis_line(reanalysis,
  'features <- read.delim(input$re_features$datapath,header=FALSE,comment.char="#",quote="")',
  c('        feature_path <- if (cervical) cervical_path("HeLa_H3K27ac_GRCh38.bed") else input$re_features$datapath',
    '        features <- read.delim(feature_path, header=FALSE, comment.char="#", quote="")'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'c(d,list(genes=g,phase=phase,features=features,example=ex))',
  '      c(d, list(genes=g, phase=phase, features=features, example=ex, cervical=cervical))')
reanalysis <- replace_reanalysis_line(reanalysis,
  'labels <- setNames(s$sample_id,paste(s$condition,s$time_hours,"h ·",s$sample_id))',
  c('    label <- if (stored()$cervical) paste(s$condition, "·", s$sample_id) else paste(s$condition, s$time_hours, "h ·", s$sample_id)',
    '    labels <- setNames(s$sample_id, label)'))
reanalysis <- replace_reanalysis_line(reanalysis,
  'list(ids=ids,labels=setNames(paste(s$condition,s$time_hours,"h ·",ids),ids))',
  c('    label <- if (d$cervical) paste(s$condition, "·", ids) else paste(s$condition, s$time_hours, "h ·", ids)',
    '    list(ids=ids, labels=setNames(label, ids))'))
writeLines(reanalysis, reanalysis_file)
map_row <- grep('h4("Local contact maps"),plotOutput("re_locus_maps",height=400),', reanalysis, fixed = TRUE)
track_output <- grep('  output$re_locus_tracks <-', reanalysis, fixed = TRUE)
if (length(map_row) == 1L && length(track_output) == 1L) {
  reanalysis <- append(reanalysis, c(
    '  output$re_locus_3d <- renderPlot({',
    '    d <- locus()',
    '    draw_locus_3d(d$view, d$annotation$tss0[1], selected()$labels)',
    '  }, alt="Inferred 3D models for both samples and their aligned overlay.")'
  ), after = track_output - 1L)
  reanalysis <- append(reanalysis, c(
    '          h4("Inferred 3D chromatin"),',
    '          plotOutput("re_locus_3d", height=300),',
    '          p("Each bead is a genomic bin. Outlined beads mark the promoter. Teal is the reference and orange is the comparison. The overlay aligns shared bins. Coordinates are inferred from contacts, in arbitrary units. Missing contacts in sufficiently covered windows are assigned the contact floor."),'
  ), after = map_row)
  writeLines(reanalysis, reanalysis_file)
}

# Add the RAD51B evidence tab to the unpacked Shiny UI.
ui_file <- file.path(source_dir, "app.R")
ui_lines <- readLines(ui_file, warn = FALSE)
insert_after <- grep('reanalysis_ui(),', ui_lines, fixed = TRUE)
stopifnot(length(insert_after) == 1L)
ui_lines <- append(ui_lines, c(
  '      tabPanel("cooltools results",',
  '        tags$section(class="card",',
  '          h3("Chromosome 14 cooltools analysis"),',
  '          p(tags$a("Run status and downloadable results", href="https://github.com/Bafana-saltiel/Cervical-Chromatin-Shiny/actions/workflows/cervical-cooltools.yml", target="_blank", rel="noopener")),',
  '          p("Open a successful run, download the HeLa or CaSki artifact, extract it, and select its TSV files below. The analysis covers full GRCh38 chromosome 14. The plots focus on the RAD51B region."),',
  '          fileInput("ct_ins", "Insulation TSV", accept=".tsv"),',
  '          fileInput("ct_eig", "Compartment eigenvectors TSV", accept=".tsv"),',
  '          fileInput("ct_dots", "Dots TSV", accept=".tsv"),',
  '          p("Keep each cell line in its own upload set. Compartment orientation uses that cell line’s H3K27ac signal. Dot significance describes enrichment within a map. Cross-line differences are descriptive."),',
  '          h4("Insulation and boundaries"), plotOutput("ct_ins_plot", height=260), tableOutput("ct_ins_table"),',
  '          h4("Compartment eigenvector"), plotOutput("ct_eig_plot", height=260),',
  '          h4("RAD51B promoter-associated dots"), tableOutput("ct_dot_table"))),',
  '      tabPanel("Cervical RAD51B workflow",',
  '        tags$section(class="card",',
  '          h3("Compare RAD51B promoter contacts in HeLa and CaSki"),',
  '          p("GSE195631 · GRCh38 · chr14:65.25–70.50 Mb · 25 kb bins"),',
  '          p(actionButton("open_rad51b_example", "Load cervical example"), " ", tags$a("Example files on GitHub", href="https://github.com/Bafana-saltiel/Cervical-Chromatin-Shiny/blob/main/GSE195631_RAD51B_HiC_upload.zip", target="_blank", rel="noopener")),',
  '          tags$ol(',
  '            tags$li("Load the bundled data to open the Hi-C reanalysis with the cervical example selected."),',
  '            tags$li("Start with HeLa as the reference, CaSki as the comparison, and RAD51B as the promoter."),',
  '            tags$li("Both views show the bundled HeLa H3K27ac track. Overlap in the CaSki view identifies a location annotated in HeLa.")),',
  '          p("The regional analysis provides contact maps, candidate TAD boundaries, unphased compartment scores, and promoter-centred virtual 4C. Compare RNA expression using the separate CSV. WGBS links are included for follow-up. Methylation has yet to be summarized."),',
  '          p(tags$a("GSE195631 study and source files", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195631", target="_blank", rel="noopener"))),',
  '        tags$section(class="card",',
  '          h3("Genome-wide analysis tools"),',
  '          tags$ul(',
  '            tags$li(tags$a("cooltools CLI", href="https://cooltools.readthedocs.io/en/latest/cli.html", target="_blank", rel="noopener"), " · Run eigs-cis on a full balanced .cool or .mcool map with a matching phasing track for oriented A/B scores. Use insulation for boundaries or dots for focal enrichments."),',
  '            tags$li(tags$a("Juicer HiCCUPS", href="https://github.com/aidenlab/juicer/wiki/HiCCUPS", target="_blank", rel="noopener"), " · Call loops with a full .hic map at a depth-supported resolution. Loop significance applies within a sample. Treatment or cell-line differences require a separate replicate-aware analysis.")),',
  '          p("These tools run outside Shiny on the full maps. The supplied regional CSV supports descriptive reanalysis."))),',
  '      tabPanel("Model datasets",',
  '        tags$section(class="card",',
  '          h3("Cervical cell-line multi-omics"),',
  '          p("Use GSE195631 to compare RAD51B promoter contacts, expression, DNA methylation, and H3K27ac across cervical cell lines. Assays are aligned by cell-line identity. The deposit does not establish pairing by biological replicate across assays."),',
  '          tags$table(class="table table-striped",',
  '            tags$thead(tags$tr(tags$th("Cell line"), tags$th("Hi-C"), tags$th("RNA-seq"), tags$th("WGBS"), tags$th("RAD51B (TPM)"))),',
  '            tags$tbody(',
  '              tags$tr(tags$th("HeLa"), tags$td("GSM6658912"), tags$td("GSM5841997"), tags$td("GSM6658935"), tags$td("2.90")),',
  '              tags$tr(tags$th("CaSki"), tags$td("GSM6658910"), tags$td("GSM5841996"), tags$td("GSM6658934"), tags$td("9.60")),',
  '              tags$tr(tags$th("SiHa"), tags$td("GSM6658915"), tags$td("GSM5841999"), tags$td("GSM6658937"), tags$td("4.92")),',
  '              tags$tr(tags$th("C33A"), tags$td("GSM6658908"), tags$td("GSM5841995"), tags$td("GSM6658933"), tags$td("24.44")),',
  '              tags$tr(tags$th("S12"), tags$td("GSM6658914"), tags$td("GSM5841998"), tags$td("GSM6658936"), tags$td("4.13")))),',
  '          p("Expression values are from the deposited gene-level tables for ENSG00000182185. Use GRCh38 and retain CpG coverage when summarizing methylation at the promoter and distal contact. Cross-line comparisons are descriptive."),',
  '          tags$a("GSE195631 assay records and downloads", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195631")),',
  '        tags$section(class="card",',
  '          h3("HeLa and HeLa-S3 perturbation model"),',
  '          p("Use GSE176066 for conventional Hi-C and RNA-seq in baseline HeLa, and for contact–expression comparisons in HeLa-S3 after TNF-α stimulation. The perturbation contact assay is NicE-C, which profiles interactions involving accessible chromatin."),',
  '          tags$table(class="table table-striped",',
  '            tags$thead(tags$tr(tags$th("Condition"), tags$th("Contact assay and accessions"), tags$th("RNA-seq accessions"))),',
  '            tags$tbody(',
  '              tags$tr(tags$th("HeLa baseline"), tags$td("Hi-C: GSM5354334–GSM5354335"), tags$td("GSM5354338–GSM5354339")),',
  '              tags$tr(tags$th("HeLa-S3 control"), tags$td("NicE-C: GSM5354344–GSM5354345"), tags$td("GSM5354348–GSM5354349")),',
  '              tags$tr(tags$th("HeLa-S3 TNF-α"), tags$td("NicE-C: GSM5354342–GSM5354343"), tags$td("GSM5354346–GSM5354347")))),',
  '          p("Keep TNF-α and control in a within-study comparison. This study uses hg19. Map the RAD51B promoter and distal anchor to hg19 before extracting contacts, and use replicate RNA-seq for differential expression."),',
  '          tags$a("GSE176066 assay records and downloads", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176066")),',
  '        tags$section(class="card",',
  '          h3("Interpretation across studies"),',
  '          p("The GSE198610 ICRF-193 comparison remains a separate HeLa-S3 Hi-C analysis. These model datasets can provide baseline multi-omics and a TNF-α example. RAD51B transcription after ICRF-193 still requires treatment-matched RNA."))),',
  '      tabPanel("RAD51B evidence",',
  '        tags$section(class="card",',
  '          h3("Candidate distal regulatory contact"),',
  '          p("RAD51B promoter bin chr14:67,750,000–68,000,000 and distal bin chr14:68,750,000–69,000,000 (GRCh38, zero-based, end-exclusive)."),',
  '          h4("Enhancer marks in the distal bin"),',
  '          p("Seven HeLa-S3 H3K27ac peak intervals from ENCODE2 overlap the bin after mapping from hg19 to GRCh38. The ENCODE combined registry lists 166 distal enhancer-like elements in the bin, 94 of which intersect those H3K27ac intervals. The registry aggregates cell types. The histone marks come from a separate HeLa-S3 experiment and do not establish activity under either treatment."),',
  '          tags$details(tags$summary("HeLa-S3 H3K27ac overlap intervals (GRCh38)"),',
  '            p("68,785,420–68,787,532; 68,789,920–68,798,246; 68,813,497–68,821,443; 68,824,398–68,824,513; 68,824,596–68,825,054; 68,911,873–68,912,379; 68,919,493–68,990,265 (chr14).")),',
  '          tags$a("ENCODE H3K27ac experiment", href="https://www.encodeproject.org/experiments/ENCSR000AOC/"), " · ",',
  '          tags$a("ENCODE cCRE track", href="https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=encodeCcreCombined")),',
  '        tags$section(class="card",',
  '          h4("Contact enrichment by replicate"),',
  '          p("Unsorted HeLa-S3 mitotic release at 4 h. Observed/expected uses the balanced contact divided by the mean balanced chr14 contact at 1 Mb in the same sample."),',
  '          tags$table(class="table table-striped",',
  '            tags$thead(tags$tr(tags$th("Replicate"), tags$th("DMSO O/E"), tags$th("ICRF-193 O/E"), tags$th("Change"))),',
  '            tags$tbody(',
  '              tags$tr(tags$th("R1"), tags$td("2.819 (GSM5953578)"), tags$td("2.410 (GSM5953577)"), tags$td("−14.5%")),',
  '              tags$tr(tags$th("R2"), tags$td("2.636 (GSM5953572)"), tags$td("2.242 (GSM5953571)"), tags$td("−15.0%")),',
  '              tags$tr(tags$th("R3"), tags$td("2.745 (GSM5953542)"), tags$td("2.534 (GSM5953540)"), tags$td("−7.7%")))),',
  '          p("Enrichment exceeded 1 in all six samples and was lower with ICRF-193 in each pair. This is a descriptive comparison of the selected bin pair."),',
  '          tags$a("GEO GSE198610", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198610")),',
  '        tags$section(class="card",',
  '          h4("RAD51B expression in HeLa"),',
  '          p("DepMap Public 24Q4 HeLa (ACH-001086) baseline RAD51B expression is 1.952 log2(TPM + 1), equivalent to approximately 2.87 TPM. This measurement comes from the HeLa cell-line data bundled with the app."),',
  '          p(tags$a("DepMap Public 24Q4", href="https://doi.org/10.25452/figshare.plus.27993248.v1")),',
  '          h5("GDSC2 drug response in HeLa"),',
  '          p("DepMap ACH-001086 maps to GDSC2 HeLa SIDM00846 (COSMIC 1298134). The table shows topoisomerase II–related compounds with fitted response data for this model in GDSC2 release 8.5."),',
  '          tags$table(class="table table-striped",',
  '            tags$thead(tags$tr(tags$th("Compound (GDSC ID)"), tags$th("Fitted IC50 (µM)"), tags$th("Maximum screened (µM)"), tags$th("AUC"))),',
  '            tags$tbody(',
  '              tags$tr(tags$th("Epirubicin (1511)"), tags$td("1.262"), tags$td("1"), tags$td("0.918")),',
  '              tags$tr(tags$th("Teniposide (1809)"), tags$td("3.113"), tags$td("3"), tags$td("0.905")),',
  '              tags$tr(tags$th("Mitoxantrone (1810)"), tags$td("3.704"), tags$td("3"), tags$td("0.917")))),',
  '          p("IC50 was calculated as exp(LN_IC50). All three estimates exceed the maximum concentration screened. GDSC2 release 8.5 contains no fitted response entry for ICRF-193. These viability assays provide pharmacological context for other compounds. RAD51B transcription after treatment remains unmeasured."),',
  '          p(tags$a("GDSC2 dose-response data", href="https://cellmodelpassports.sanger.ac.uk/downloads/gdsc"), " · ", tags$a("HeLa model identity", href="https://cellmodelpassports.sanger.ac.uk/passports/SIDM00846")),',
  '          h5("Treatment-matched expression") ,',
  '          p("Untested. The GSE198610 deposit contains chromatin-contact assays. We did not identify matched RNA measurements for these six Hi-C samples. The regulatory effect on RAD51B expression remains unknown."))),'
), after = insert_after)
server_hook <- grep('server <- function(input,output,session){', ui_lines, fixed = TRUE)
stopifnot(length(server_hook) == 1L)
ui_lines <- append(ui_lines, c(
  '  ct_read <- function(file, columns) {',
  '    req(file)',
  '    x <- read.delim(file$datapath, check.names=FALSE)',
  '    validate(need(all(columns %in% names(x)), paste("Required columns", paste(columns, collapse=", "))))',
  '    x',
  '  }',
  '  ct_region <- function(x) {',
  '    x[x$chrom %in% c("chr14", "14") & x$start < 70500000 & x$end > 65250000, , drop=FALSE]',
  '  }',
  '  ct_ins <- reactive({',
  '    x <- ct_read(input$ct_ins, c("chrom", "start", "end", "log2_insulation_score_250000", "is_boundary_250000"))',
  '    ct_region(x)',
  '  })',
  '  output$ct_ins_plot <- renderPlot({',
  '    x <- ct_ins()',
  '    validate(need(nrow(x)>1, "Upload chromosome 14 insulation results."))',
  '    plot((x$start+x$end)/2e6, x$log2_insulation_score_250000, type="l", col="#008898", xlab="chr14 (Mb)", ylab="Log2 insulation", main="250 kb insulation window")',
  '    b <- x$is_boundary_250000 %in% c(TRUE, "True", "TRUE")',
  '    abline(v=(x$start[b]+x$end[b])/2e6, col="#d65e3a", lty=3)',
  '    abline(v=67.819778, lty=2)',
  '  })',
  '  output$ct_ins_table <- renderTable({',
  '    x <- ct_ins()',
  '    x[x$is_boundary_250000 %in% c(TRUE, "True", "TRUE"), intersect(c("chrom", "start", "end", "boundary_strength_250000"), names(x)), drop=FALSE]',
  '  })',
  '  output$ct_eig_plot <- renderPlot({',
  '    x <- ct_region(ct_read(input$ct_eig, c("chrom", "start", "end", "E1")))',
  '    validate(need(nrow(x)>1, "Upload chromosome 14 eigenvectors."))',
  '    plot((x$start+x$end)/2e6, x$E1, type="h", col=ifelse(x$E1>=0, "#008898", "#d65e3a"), xlab="chr14 (Mb)", ylab="E1", main="H3K27ac-phased eigenvector at 250 kb")',
  '    abline(h=0, col="grey60")',
  '    abline(v=67.819778, lty=2)',
  '  })',
  '  output$ct_dot_table <- renderTable({',
  '    x <- ct_read(input$ct_dots, c("chrom1", "start1", "end1", "chrom2", "start2", "end2"))',
  '    a <- x$chrom1 %in% c("chr14", "14") & x$start1 <= 67819778 & x$end1 > 67819778',
  '    b <- x$chrom2 %in% c("chr14", "14") & x$start2 <= 67819778 & x$end2 > 67819778',
  '    x <- x[a | b, , drop=FALSE]',
  '    validate(need(nrow(x)>0, "No called dot overlaps the RAD51B promoter bin in this file."))',
  '    x',
  '  })',
  '  observeEvent(input$open_rad51b_example, {',
  '    updateTabsetPanel(session, "analysis_tab", selected = "Hi-C reanalysis")',
  '    updateRadioButtons(session, "re_source", selected = "Cervical RAD51B (GSE195631)")',
  '  })'
), after = server_hook)
writeLines(ui_lines, ui_file)

addResourcePath("chromatin-assets", file.path(source_dir, "www"))
app_environment <- new.env(parent = globalenv())
evalq({
  draw_locus_3d <- function(view, tss, labels) {
    coverage <- lapply(view$maps, function(m) {
      positive <- is.finite(m) & m > 0
      diag(positive) <- FALSE
      rowSums(positive) > 0L
    })
    keep <- Reduce("&", coverage)
    bins <- view$bins[keep]
    anchor <- which(bins <= tss & tss < bins + view$width)
    validate(need(length(bins) >= 8L && length(anchor) == 1L,
                  "The 3D view needs eight shared covered bins, including the promoter. Try a wider flank."))
    fits <- lapply(seq_along(view$maps), function(i) {
      m <- view$maps[[i]][keep, keep, drop = FALSE]
      if (mean(is.finite(m[upper.tri(m)])) < .8) return(NULL)
      m[!is.finite(m)] <- 0
      fit_hic_window(list(bins = seq_along(bins), contacts = m, anchor = anchor,
                          bin_size = view$width), view$ids[i],
                     list(gene = "", chr = view$chr))
    })
    validate(need(all(vapply(fits, Negate(is.null), logical(1))),
                  "The shared bins have insufficient contact coverage for 3D reconstruction."))
    xyz <- lapply(fits, function(f) as.matrix(f$coordinates[c("x", "y", "z")]))
    alignment <- svd(crossprod(xyz[[2]], xyz[[1]]))
    xyz[[2]] <- xyz[[2]] %*% alignment$u %*% t(alignment$v)
    camera <- rbind(c(.8, 0, .6), c(-.3, .866, .4), c(-.52, -.5, .693))
    coords <- lapply(xyz, function(x) x %*% t(camera))
    extent <- max(abs(unlist(coords))) * 1.12
    colors <- c("#087e89", "#c45c38")
    j <- which(diff(bins) == view$width)
    par(mfrow = c(1, 3), mar = c(1, .5, 3, .5))
    for (panel in seq_len(3L)) {
      plot.new()
      plot.window(c(-extent, extent), c(-extent, extent), asp = 1)
      title(main = c(unname(labels), "Aligned overlay")[panel], cex.main = .78)
      ids <- if (panel == 3L) 1:2 else panel
      pieces <- do.call(rbind, lapply(ids, function(i) {
        x <- coords[[i]]
        data.frame(x0 = x[j, 1], y0 = x[j, 2], x1 = x[j+1L, 1], y1 = x[j+1L, 2],
                   depth = (x[j, 3] + x[j+1L, 3]) / 2, color = rep(colors[i], length(j)))
      }))
      pieces <- pieces[order(pieces$depth), ]
      if (nrow(pieces)) with(pieces, segments(x0, y0, x1, y1, col = color, lwd = 2.2, lend = "round"))
      for (i in ids) {
        x <- coords[[i]]
        points(x[, 1], x[, 2], pch = 16, cex = .48, col = colors[i])
        points(x[anchor, 1], x[anchor, 2], pch = 21, cex = 1.25, bg = colors[i], col = "#172532", lwd = 1.3)
      }
    }
  }
}, app_environment)
previous <- options(chromatin.assetPrefix = "chromatin-assets/")
app <- tryCatch(
  source(file.path(source_dir, "app.R"), local = app_environment, chdir = TRUE)$value,
  finally = options(previous)
)
app
