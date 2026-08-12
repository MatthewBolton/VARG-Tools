args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: test_processing_parallel_stability.R <app-directory>", call. = FALSE)
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

old_environment <- Sys.getenv(
  c(
    "VARG_MAX_CORES", "VARG_HEAVY_CPU_TOKENS", "VARG_HEAVY_SLOTS",
    "VARG_HEAVY_LOCK_DIR", "VARG_HEAVY_LOCK_STALE_MINUTES"
  ),
  unset = NA_character_
)
on.exit({
  for (name in names(old_environment)) {
    value <- old_environment[[name]]
    if (is.na(value)) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(value), name))
  }
}, add = TRUE)

lock_directory <- tempfile("vargtools_parallel_locks_")
on.exit(unlink(lock_directory, recursive = TRUE, force = TRUE), add = TRUE)
Sys.setenv(
  VARG_MAX_CORES = "8",
  VARG_HEAVY_CPU_TOKENS = "12",
  VARG_HEAVY_SLOTS = "3",
  VARG_HEAVY_LOCK_DIR = lock_directory,
  VARG_HEAVY_LOCK_STALE_MINUTES = "75"
)

stopifnot(processing_heavy_worker_cap() == 8L)
blas_environment <- processing_blas_worker_env(4L)
stopifnot(
  identical(unname(blas_environment[["OPENBLAS_NUM_THREADS"]]), "4"),
  identical(unname(blas_environment[["OMP_NUM_THREADS"]]), "4"),
  identical(unname(blas_environment[["MKL_NUM_THREADS"]]), "4")
)

first_lock <- heavy_job_limiter_acquire("parallel test 1", min_tokens = 1L, max_tokens = processing_heavy_worker_cap())
second_lock <- heavy_job_limiter_acquire("parallel test 2", min_tokens = 1L, max_tokens = processing_heavy_worker_cap())
third_lock <- heavy_job_limiter_acquire("parallel test 3", min_tokens = 1L, max_tokens = processing_heavy_worker_cap())
stopifnot(
  isTRUE(first_lock$acquired), first_lock$workers == 8L,
  isTRUE(second_lock$acquired), second_lock$workers == 4L,
  !isTRUE(third_lock$acquired)
)
stopifnot(
  isTRUE(heavy_job_limiter_release(first_lock$lock)),
  isTRUE(heavy_job_limiter_release(second_lock$lock))
)

environment_probe <- callr::r(
  function() Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
  env = processing_blas_worker_env(4L)
)
stopifnot(all(unname(environment_probe) == "4"))

set.seed(2468)
gmm_data <- rbind(
  matrix(stats::rnorm(600, mean = -2, sd = 0.45), ncol = 3),
  matrix(stats::rnorm(600, mean = 0.5, sd = 0.55), ncol = 3),
  matrix(stats::rnorm(600, mean = 3.5, sd = 0.50), ncol = 3)
)

fit_gmm <- function(data, workers) {
  callr::r(
    function(data) {
      library(mclust)
      set.seed(42)
      fit <- mclust::Mclust(
        data,
        G = 1:6,
        initialization = list(noise = TRUE),
        prior = NULL,
        verbose = FALSE
      )
      list(
        G = fit$G,
        modelName = fit$modelName,
        classification = fit$classification,
        bic = unname(fit$bic)
      )
    },
    args = list(data = data),
    env = processing_blas_worker_env(workers)
  )
}

gmm_single <- fit_gmm(gmm_data, 1L)
gmm_parallel <- fit_gmm(gmm_data, 4L)
stopifnot(
  identical(gmm_single$G, gmm_parallel$G),
  identical(gmm_single$modelName, gmm_parallel$modelName),
  identical(gmm_single$classification, gmm_parallel$classification),
  isTRUE(all.equal(gmm_single$bic, gmm_parallel$bic, tolerance = 1e-8))
)

set.seed(1357)
umap_data <- matrix(stats::rnorm(1800), ncol = 6)
umap_single <- uwot::umap(
  umap_data,
  n_neighbors = 15,
  n_components = 2,
  ret_model = TRUE,
  seed = 42,
  n_threads = 1,
  n_sgd_threads = 1,
  verbose = FALSE
)
umap_parallel <- uwot::umap(
  umap_data,
  n_neighbors = 15,
  n_components = 2,
  ret_model = TRUE,
  seed = 42,
  n_threads = 4,
  n_sgd_threads = 1,
  verbose = FALSE
)
stopifnot(identical(umap_single$embedding, umap_parallel$embedding))

projection_data <- umap_data[1:40, , drop = FALSE]
projection_single <- uwot::umap_transform(
  projection_data, umap_single,
  n_threads = 1,
  n_sgd_threads = 1
)
projection_parallel <- uwot::umap_transform(
  projection_data, umap_single,
  n_threads = 4,
  n_sgd_threads = 1
)
stopifnot(identical(projection_single, projection_parallel))

cat("processing_parallel_stability_ok gmm_workers=4 umap_workers=4 allocation=8+4\n")
