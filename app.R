library(shiny)

# Extract the reproducible source bundle into a writable location on the host.
source_dir <- tempfile("chromatin-source-")
dir.create(source_dir)
utils::untar("chromatin-source.tar.gz", exdir = source_dir)

# The app reads the data during startup. Shiny serves the bundled assets
# through an explicit resource path so they remain available across sessions.
addResourcePath("chromatin-assets", file.path(source_dir, "www"))
app_environment <- new.env(parent = globalenv())
previous <- options(chromatin.assetPrefix = "chromatin-assets/")
app <- tryCatch(
  source(file.path(source_dir, "app.R"), local = app_environment, chdir = TRUE)$value,
  finally = options(previous)
)
app
