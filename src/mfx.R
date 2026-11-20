# Paths for R scripts: reads the flat key: value lines of config/paths.local.yaml (or paths.yaml,
# or $MFX_CONFIG), same resolution order as mfx.paths in Python.
mfx_config <- function() {
  root <- Sys.getenv("MFX_ROOT", unset = ".")
  env <- Sys.getenv("MFX_CONFIG", unset = "")
  cands <- c(if (nzchar(env)) env, file.path(root, "config/paths.local.yaml"), file.path(root, "config/paths.yaml"))
  f <- cands[file.exists(cands)][1]
  if (is.na(f)) stop("no paths config found")
  lines <- readLines(f)
  lines <- lines[grepl("^[A-Za-z_]+:", lines)]
  key <- sub(":.*$", "", lines)
  val <- trimws(sub("#.*$", "", sub("^[A-Za-z_]+:", "", lines)))
  setNames(as.list(val), key)
}
mfx_results <- function(name) {
  p <- file.path(mfx_config()$work, "results", name)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
