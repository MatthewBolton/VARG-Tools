args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: test_ap97_fixture_preprocessing.R <app-directory>", call. = FALSE)
}

app_directory <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
old_directory <- setwd(app_directory)
on.exit(setwd(old_directory), add = TRUE)
Sys.setenv(VARG_BG_PROCESSES = "false")

source("R/utils/package_compatibility.R")
varg_assert_package_compatibility()
source("global.R")
source("R/utils/tooltip_helpers.R")
source("R/utils/utils_ui.R")
source("R/utils/heavy_job_limiter.R")
source("R/functions/fct_analysis.R")
source("R/modules/mod_processing.R")

fixture_archive <- "Data for Beta Testers.zip"
fixture_entry <- utils::unzip(fixture_archive, list = TRUE)$Name
fixture_entry <- fixture_entry[grepl("AP_Geochem_to test[.]csv$", fixture_entry)]
if (length(fixture_entry) != 1L) {
  stop("Expected exactly one AP_Geochem_to test.csv file in the tester archive.", call. = FALSE)
}
fixture <- utils::read.csv(
  unz(fixture_archive, fixture_entry),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
oxide_columns <- c("SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "Cl")
ap97 <- fixture$core_id == "AP97"
ap17 <- fixture$core_id == "AP17"
source_zero <- fixture$sample_point == "97AP-M_011"
zero_locations <- which(as.matrix(fixture[oxide_columns]) == 0, arr.ind = TRUE)

stopifnot(
  nrow(fixture) == 770L,
  sum(ap97) == 171L,
  sum(ap17) == 599L,
  all(is.na(fixture$MnO[ap97])),
  !any(is.na(fixture$MnO[ap17])),
  !any(is.na(fixture$MgO)),
  nrow(zero_locations) == 1L,
  fixture$sample_point[zero_locations[, "row"]] == "97AP-M_011",
  oxide_columns[zero_locations[, "col"]] == "TiO2",
  identical(fixture$TiO2[source_zero], 0),
  isTRUE(all.equal(fixture$TiO2[fixture$sample_point == "97AP-E-a_002"], 0.8, tolerance = 1e-10)),
  isTRUE(all.equal(fixture$CaO[fixture$sample_point == "97AP-G_001"], 0.9, tolerance = 1e-10)),
  isTRUE(all.equal(fixture$Na2O[fixture$sample_point == "97AP-K_006"], 3.83, tolerance = 1e-10)),
  isTRUE(all.equal(fixture$FeO[fixture$sample_point == "97AP-N-a_003"], 1.39263252, tolerance = 1e-10))
)

missingness_report <- processing_report_imputation_summary(fixture, oxide_columns)
structured_report <- processing_report_structured_missingness(fixture, oxide_columns)
auto_report <- processing_report_imputation_auto_html(missingness_report, nrow(fixture), "ltsReg")

stopifnot(
  missingness_report$Missing[missingness_report$Analyte == "MnO"] == 171L,
  missingness_report[["Affected (%)"]][missingness_report$Analyte == "MnO"] == 22.2,
  any(
    structured_report[["Grouping column"]] == "core_id" &
      structured_report$Group == "AP97" &
      structured_report$Analyte == "MnO" &
      structured_report$Rows == 171L
  ),
  grepl("not above the 50% threshold", auto_report, fixed = TRUE),
  grepl("Robust regression (ltsReg)", auto_report, fixed = TRUE)
)

global_state <- shiny::reactiveValues(
  data = fixture,
  original_cols = names(fixture),
  umap_mode_ran = NULL,
  mclust_result = NULL,
  data_stale = FALSE,
  pipeline_config = NULL,
  user_umap_model = NULL,
  source_filename = basename(fixture_entry),
  data_generation = 1L,
  heavy_job = NULL
)

shiny::testServer(
  mod_processing_server,
  args = list(global_rv = global_state),
  {
    session$setInputs(
      comp_cols = oxide_columns,
      noncomp_cols = c("core_id", "sample_id", "sample_point", "age_mean", "age_sd"),
      do_impute = TRUE,
      impute_method = "auto",
      transform_type = "ilr",
      pivot_var = "SiO2"
    )
    session$flushReact()
    session$setInputs(apply_transforms = 1L)
    session$flushReact()

    stopifnot(
      "MnO_imp" %in% names(global_state$data),
      all(is.finite(global_state$data$MnO_imp[ap97])),
      all(is.finite(global_state$data$MnO_imp[ap17])),
      all(is.na(global_state$data$MnO[ap97])),
      identical(global_state$data$TiO2[source_zero], 0),
      !any(is.na(global_state$data$MgO)),
      sum(is.na(global_state$data$MnO)) == 171L,
      all(as.matrix(global_state$data[paste0(oxide_columns, "_imp")]) > 0),
      all(rowSums(global_state$data[paste0(oxide_columns, "_imp")]) > 99.999),
      all(rowSums(global_state$data[paste0(oxide_columns, "_imp")]) < 100.001),
      length(grep("^pivot_", names(global_state$data))) == 9L,
      all(is.finite(as.matrix(global_state$data[grep("^pivot_", names(global_state$data))])))
    )
  }
)

cat("ap97_fixture_preprocessing_ok missing_MnO=171 zero_TiO2=1\n")
