library(shiny)

# Connect Cloud can keep application source read-only.
source_dir <- tempfile("chromatin-source-")
dir.create(source_dir)
archives <- c("chromatin-source.tar.gz", list.files(pattern = "^hic-.*[.]tar[.]gz$"))
for (archive in archives) {
  stopifnot(utils::untar(archive, exdir = source_dir) == 0L)
}
addResourcePath("chromatin-assets", file.path(source_dir, "www"))
app_environment <- new.env(parent = globalenv())
previous <- options(chromatin.assetPrefix = "chromatin-assets/")
app <- tryCatch(
  source(file.path(source_dir, "app.R"), local = app_environment, chdir = TRUE)$value,
  finally = options(previous)
)
app
