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

stopifnot(
  identical(umap_result_tab_name("pretrained", "1d"), "VARG26 1D"),
  identical(umap_result_tab_name("pretrained", "both"), "VARG26 2D"),
  identical(umap_result_tab_name("new", 1L), "New 1D"),
  identical(umap_result_tab_name("new", c(1L, 2L)), "New 2D"),
  identical(umap_result_tab_name("loaded", 2L), "Loaded 2D"),
  is.null(umap_result_tab_name("unknown", 2L)),
  is.null(umap_result_tab_name("new", 3L))
)

both_umap_report <- processing_report_umap_html(list(
  mode = "new",
  n_components = c(1L, 2L),
  n_neighbors = 15L,
  min_dist = 0.1,
  dens_scale = 0,
  semisupervised = FALSE,
  columns_used = c("pivot_A", "pivot_B")
))

stopifnot(
  length(both_umap_report) == 1L,
  grepl(">1D and 2D<", both_umap_report, fixed = TRUE),
  lengths(regmatches(both_umap_report, gregexpr("4. UMAP Dimensionality Reduction", both_umap_report, fixed = TRUE))) == 1L,
  !grepl("c(\"", both_umap_report, fixed = TRUE),
  identical(processing_report_imputation_method_label("ltsReg"), "Robust regression (ltsReg)"),
  identical(processing_report_imputation_method_label("lm"), "Standard regression (lm)")
)

report_fixture <- data.frame(
  core_id = c(rep("AP97", 3L), rep("AP17", 4L)),
  MnO = c(NA, NA, NA, 0.05, 0.06, 0.07, 0.08),
  TiO2 = c(0, rep(0.8, 6L)),
  stringsAsFactors = FALSE
)
missingness_report <- processing_report_imputation_summary(report_fixture, c("MnO", "TiO2"))
structured_report <- processing_report_structured_missingness(report_fixture, c("MnO", "TiO2"))

stopifnot(
  missingness_report$Missing[missingness_report$Analyte == "MnO"] == 3L,
  missingness_report$Affected[missingness_report$Analyte == "TiO2"] == 1L,
  missingness_report[["Affected (%)"]][missingness_report$Analyte == "MnO"] == 42.9,
  nrow(structured_report) == 1L,
  structured_report[["Grouping column"]] == "core_id",
  structured_report$Group == "AP97",
  structured_report$Analyte == "MnO",
  structured_report$Rows == 3L
)

auto_robust_html <- processing_report_imputation_auto_html(missingness_report, nrow(report_fixture), "ltsReg")
high_missingness_report <- missingness_report
high_missingness_report[["Affected (%)"]][high_missingness_report$Analyte == "MnO"] <- 60
auto_standard_html <- processing_report_imputation_auto_html(high_missingness_report, nrow(report_fixture), "lm")

stopifnot(
  grepl("not above the 50% threshold", auto_robust_html, fixed = TRUE),
  grepl("Robust regression (ltsReg)", auto_robust_html, fixed = TRUE),
  grepl("above the 50% threshold", auto_standard_html, fixed = TRUE),
  grepl("Standard regression (lm)", auto_standard_html, fixed = TRUE)
)

composition_fixture <- data.frame(
  SiO2 = c(60, 62, 70, 72),
  MnO = c(NA, NA, 0.08, 0.10),
  SiO2_imp = c(60, 62, 70, 72),
  MnO_imp = c(0.05, 0.06, 0.08, 0.10),
  gmm_cluster = c(1L, 1L, 2L, 2L),
  stringsAsFactors = FALSE
)
composition_config <- list(preprocessing = list(comp_cols = c("SiO2", "MnO")))
selected_composition <- processing_gmm_selected_compositional_columns(composition_config, NULL)
composition_summary <- processing_gmm_cluster_composition(composition_fixture, selected_composition)
composition_html <- processing_report_gmm_composition_html(
  composition_summary,
  function(data) paste0("<table><tr><td>", paste(names(data), collapse = "|"), "</td></tr></table>")
)

stopifnot(
  identical(selected_composition, c("SiO2", "MnO")),
  nrow(composition_summary) == 2L,
  identical(composition_summary$Cluster, c("1", "2")),
  identical(composition_summary$n, c(2L, 2L)),
  identical(attr(composition_summary, "data_source"), "imputed"),
  identical(composition_summary$SiO2, c("61.00 \u00B1 1.41", "71.00 \u00B1 1.41")),
  identical(composition_summary$MnO, c("0.06 \u00B1 0.01", "0.09 \u00B1 0.01")),
  grepl("Cluster Composition Summary", composition_html, fixed = TRUE),
  grepl("imputed values", composition_html, fixed = TRUE),
  grepl("Cluster|n|SiO2|MnO", composition_html, fixed = TRUE)
)

composition_state <- shiny::reactiveValues(
  data = composition_fixture,
  umap_mode_ran = NULL,
  mclust_result = NULL,
  data_stale = FALSE,
  pipeline_config = composition_config,
  user_umap_model = NULL,
  source_filename = "composition-regression.csv",
  data_generation = 1L,
  heavy_job = NULL
)

shiny::testServer(
  mod_processing_server,
  args = list(global_rv = composition_state),
  {
    reactive_summary <- gmm_cluster_composition_data()
    stopifnot(
      !is.null(reactive_summary),
      nrow(reactive_summary) == 2L,
      identical(names(reactive_summary), c("Cluster", "n", "SiO2", "MnO"))
    )
  }
)

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
