# Compatibility wrapper. The modular QMD chapters are the source of truth.
builder <- file.path("docs", "build_in_app_guide.R")
if (!file.exists(builder)) {
  stop("Run this script from the VARGtools_Copy project root.", call. = FALSE)
}

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(rscript, builder)
if (!identical(status, 0L)) {
  stop("Guide rebuild failed with status ", status, call. = FALSE)
}
