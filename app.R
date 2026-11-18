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
  '      tabPanel("RAD51B evidence",',
  '        tags$section(class="card",',
  '          h3("Candidate distal regulatory contact"),',
  '          p("RAD51B promoter bin chr14:67,750,000–68,000,000 and distal bin chr14:68,750,000–69,000,000 (GRCh38, zero-based, end-exclusive)."),',
  '          h4("Enhancer marks in the distal bin"),',
  '          p("Seven HeLa-S3 H3K27ac peak intervals from ENCODE2 overlap the bin after mapping from hg19 to GRCh38. The ENCODE combined registry lists 166 distal enhancer-like elements in the bin, 94 of which intersect those H3K27ac intervals. The registry aggregates cell types. The histone marks come from a separate HeLa-S3 experiment and do not establish activity under either treatment."),',
  '          tags$details(tags$summary("HeLa-S3 H3K27ac overlap intervals (GRCh38)"),',
  '            p("68,785,420–68,787,532; 68,789,920–68,798,246; 68,813,497–68,821,443; 68,824,398–68,824,513; 68,824,596–68,825,054; 68,911,873–68,912,379; 68,919,493–68,990,265 (chr14).")),',
  '          a("ENCODE H3K27ac experiment", href="https://www.encodeproject.org/experiments/ENCSR000AOC/"), " · ",',
  '          a("ENCODE cCRE track", href="https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=encodeCcreCombined")),',
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
  '          a("GEO GSE198610", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198610")),',
  '        tags$section(class="card",',
  '          h4("RAD51B expression in HeLa"),
  '          p("DepMap Public 24Q4 HeLa (ACH-001086) baseline RAD51B expression is 1.952 log2(TPM + 1), equivalent to approximately 2.87 TPM. This measurement comes from the HeLa cell-line data bundled with the app."),',
  '          p(a("DepMap Public 24Q4", href="https://doi.org/10.25452/figshare.plus.27993248.v1")),',
  '          h5("Treatment-matched expression") ,',
  '          p("Untested. The GSE198610 deposit contains chromatin-contact assays. We did not identify matched RNA measurements for these six Hi-C samples. The regulatory effect on RAD51B expression remains unknown."))),'
), after = insert_after)
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
