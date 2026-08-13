if (!requireNamespace("geochem", quietly = TRUE)) {
  stop("The geochem package is unavailable.", call. = FALSE)
}
if (!"tas_diagram" %in% getNamespaceExports("geochem")) {
  stop("geochem does not export tas_diagram().", call. = FALSE)
}

tas_diagram <- getExportedValue("geochem", "tas_diagram")
plot <- ggplot2::ggplot() + tas_diagram(fields = TRUE, labels = TRUE)
invisible(ggplot2::ggplot_build(plot))

cat("tas_overlay_ok\n")
