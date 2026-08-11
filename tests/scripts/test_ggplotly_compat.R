script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Cannot determine the compatibility-test location.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
app_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(app_dir, "R", "utils", "package_compatibility.R"))
varg_assert_package_compatibility()

cat(
  sprintf(
    "R=%s ggplot2=%s plotly=%s\n",
    getRversion(),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("plotly"))
  )
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(plotly)
})

plot <- ggplot(
  data.frame(x = 1:3, y = 1:3),
  aes(x = x, y = y)
) +
  geom_point()

widget <- ggplotly(plot)
stopifnot(inherits(widget, "plotly"))
cat("ggplotly_compat_ok\n")
