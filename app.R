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
