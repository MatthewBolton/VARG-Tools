args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: test_processing_plotly_renderers.R <app-directory>", call. = FALSE)
}

app_directory <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
old_directory <- setwd(app_directory)
on.exit(setwd(old_directory), add = TRUE)

source("R/utils/package_compatibility.R")
varg_assert_package_compatibility()
source("global.R")
source("R/utils/tooltip_helpers.R")
source("R/utils/utils_ui.R")
source("R/utils/heavy_job_limiter.R")
source("R/functions/fct_analysis.R")
source("R/modules/mod_processing.R")

set.seed(42)
row_count <- 50L
test_data <- data.frame(
  UID = seq_len(row_count),
  sample_point = sprintf("test_%03d", seq_len(row_count)),
  UMAP_VARG26_2D_1 = stats::rnorm(row_count),
  UMAP_VARG26_2D_2 = stats::rnorm(row_count),
  UMAP_new_2D_1 = stats::rnorm(row_count),
  UMAP_new_2D_2 = stats::rnorm(row_count),
  gmm_cluster = rep(seq_len(5L), length.out = row_count),
  population = rep(paste0("Population ", LETTERS[1:5]), length.out = row_count),
  stringsAsFactors = FALSE
)

global_state <- shiny::reactiveValues(
  data = test_data,
  umap_mode_ran = "new",
  mclust_result = NULL,
  data_stale = FALSE,
  pipeline_config = NULL,
  user_umap_model = NULL,
  source_filename = "synthetic-render-check.csv",
  data_generation = 1L,
  heavy_job = NULL
)

shiny::testServer(
  mod_processing_server,
  args = list(global_rv = global_state),
  {
    session$setInputs(
      popsel_xvar = "UMAP_VARG26_2D_1",
      popsel_yvar = "UMAP_VARG26_2D_2",
      popsel_hover = "sample_point",
      filter_col = "",
      filter_val = "All"
    )
    session$flushReact()

    umap_render <- output$umap_plot_varg26_2d
    population_render <- output$featplot
    stopifnot(!is.null(umap_render), !is.null(population_render))
  }
)

cat("processing_plotly_renderers_ok\n")
