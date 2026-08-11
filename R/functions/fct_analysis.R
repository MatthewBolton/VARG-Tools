# R/functions/fct_analysis.R

#' Run GMM Clustering
#'
#' @param data A data frame or matrix of numeric data for clustering.
#' @param g_min Minimum number of clusters (G).
#' @param g_max Maximum number of clusters (G).
#' @param noise_init Logical, whether to include noise initialization.
#' @param use_prior Logical, whether to use a conjugate prior.
#' @return An mclust object.
#' @export
run_gmm_analysis <- function(data, g_min, g_max, noise_init = TRUE, use_prior = FALSE) {
    library(mclust)

    # Validate inputs
    if (is.null(data) || nrow(data) == 0) stop("No data provided for GMM.")
    if (g_min > g_max) stop("g_min cannot be greater than g_max.")

    # Run Mclust
    res <- mclust::Mclust(
        data,
        G = g_min:g_max,
        initialization = if (noise_init) list(noise = TRUE) else NULL,
        prior = if (use_prior) mclust::priorControl() else NULL,
        verbose = FALSE
    )

    return(res)
}

#' Run UMAP Projection (New Model)
#'
#' @param data A data frame or matrix of numeric data.
#' @param n_neighbors Number of neighbors.
#' @param min_dist Minimum distance.
#' @param n_components Number of dimensions (1 or 2).
#' @param metric Distance metric (default "euclidean").
#' @param set_op_mix_ratio Set operation mix ratio (default 1).
#' @param local_connectivity Local connectivity (default 1).
#' @param dens_scale Density scale (default 0).
#' @param y Optional target data for semi-supervised learning.
#' @param target_weight Weight for target data (0 to 1).
#' @return A matrix of UMAP embeddings.
#' @export
run_umap_new <- function(data, n_neighbors = 15, min_dist = 0.1, n_components = 2,
                         metric = "euclidean", set_op_mix_ratio = 1, local_connectivity = 1,
                         dens_scale = 0, y = NULL, target_weight = 0.5) {
    library(uwot)

    # Validate inputs
    if (is.null(data) || nrow(data) == 0) stop("No data provided for UMAP.")

    # Run UMAP
    res <- uwot::umap(
        data,
        n_neighbors = n_neighbors,
        min_dist = min_dist,
        n_components = n_components,
        metric = metric,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        dens_scale = dens_scale,
        y = y,
        target_weight = target_weight,
        ret_model = TRUE # Return model to allow future projection if needed
    )

    return(res)
}

#' Project Data onto Pre-trained UMAP Model
#'
#' @param data A data frame or matrix of numeric data.
#' @param model_path Path to the saved uwot model file.
#' @return A matrix of UMAP embeddings.
#' @export
project_umap_pretrained <- function(data, model_path) {
    library(uwot)

    if (!file.exists(model_path)) stop("Model file not found.")

    model <- uwot::load_uwot(model_path)
    res <- uwot::umap_transform(data, model)

    return(res)
}


#' Normalize Saved UMAP Dimension Metadata
#'
#' @param umap_config A saved UMAP configuration.
#' @return An integer vector containing one or both of `1L` and `2L`.
normalize_umap_dimensions <- function(umap_config) {
    if (is.null(umap_config)) return(integer(0))

    raw <- umap_config$dimensions
    if (is.null(raw)) raw <- umap_config$n_components
    if (is.null(raw)) return(integer(0))

    raw <- tolower(trimws(as.character(unlist(raw, use.names = FALSE))))
    dimensions <- integer(0)
    for (value in raw) {
        if (identical(value, "both")) {
            dimensions <- c(dimensions, 1L, 2L)
        } else if (value %in% c("1", "1d")) {
            dimensions <- c(dimensions, 1L)
        } else if (value %in% c("2", "2d")) {
            dimensions <- c(dimensions, 2L)
        }
    }

    sort(unique(dimensions))
}

#' Validate the Reproducible Saved-UMAP Projection Contract
#'
#' Saved UMAP projection uses already prepared training columns. It deliberately
#' does not re-estimate normalization or imputation from the projection batch.
#'
#' @param pipeline_config The saved pipeline configuration.
#' @return Validation details, required columns, and saved dimensions.
validate_umap_projection_contract <- function(pipeline_config) {
    errors <- character(0)
    warnings <- character(0)

    if (is.null(pipeline_config) || is.null(pipeline_config$umap)) {
        return(list(
            valid = FALSE,
            errors = "Saved model configuration is missing UMAP settings.",
            warnings = warnings,
            columns = character(0),
            dimensions = integer(0)
        ))
    }

    umap_config <- pipeline_config$umap
    contract <- umap_config$projection_contract
    if (is.null(contract)) {
        errors <- c(
            errors,
            paste0(
                "This legacy UMAP bundle does not contain a reproducible projection contract. ",
                "Recreate and save the model with the current VARG-Tools version."
            )
        )
    } else {
        version_raw <- unlist(contract$version, use.names = FALSE)
        version <- if (length(version_raw) == 0) NA_integer_ else {
            suppressWarnings(as.integer(version_raw[[1]]))
        }
        type_raw <- unlist(contract$type, use.names = FALSE)
        type <- if (length(type_raw) == 0) NA_character_ else as.character(type_raw[[1]])
        if (is.na(version) || version != 2L || !identical(type, "direct_columns")) {
            errors <- c(errors, "Unsupported saved UMAP projection contract.")
        }
    }

    columns <- if (is.null(contract$columns)) character(0) else {
        as.character(unlist(contract$columns, use.names = FALSE))
    }
    if (length(columns) == 0 || anyNA(columns) || any(!nzchar(columns))) {
        errors <- c(errors, "Saved UMAP projection contract has no valid training columns.")
    } else if (anyDuplicated(columns)) {
        errors <- c(errors, "Saved UMAP projection contract contains duplicate training columns.")
    }

    configured_columns <- if (is.null(umap_config$columns_used)) character(0) else {
        as.character(unlist(umap_config$columns_used, use.names = FALSE))
    }
    if (length(configured_columns) > 0 && !identical(configured_columns, columns)) {
        errors <- c(errors, "Saved UMAP training-column metadata is internally inconsistent.")
    }

    expected_raw <- unlist(umap_config$n_input_cols, use.names = FALSE)
    expected_ncol <- if (length(expected_raw) == 0) NA_integer_ else {
        suppressWarnings(as.integer(expected_raw[[1]]))
    }
    if (is.na(expected_ncol) || expected_ncol != length(columns)) {
        errors <- c(errors, "Saved UMAP input-column count does not match its training-column metadata.")
    }

    dimensions <- normalize_umap_dimensions(umap_config)
    dimension_fields <- list(
        dimensions = umap_config$dimensions,
        n_components = umap_config$n_components
    )
    supported_dimension_tokens <- c("both", "1", "1d", "2", "2d")
    for (field_name in names(dimension_fields)) {
        raw_field <- dimension_fields[[field_name]]
        if (is.null(raw_field)) next
        tokens <- tolower(trimws(as.character(unlist(raw_field, use.names = FALSE))))
        unsupported <- setdiff(tokens, supported_dimension_tokens)
        if (length(unsupported) > 0) {
            errors <- c(
                errors,
                paste0(
                    "Saved UMAP ", field_name,
                    " contains unsupported output dimensions: ",
                    paste(unique(unsupported), collapse = ", ")
                )
            )
        }
    }
    if (!is.null(umap_config$dimensions) && !is.null(umap_config$n_components)) {
        normalized_dimensions <- normalize_umap_dimensions(list(dimensions = umap_config$dimensions))
        normalized_components <- normalize_umap_dimensions(list(dimensions = umap_config$n_components))
        if (!identical(normalized_dimensions, normalized_components)) {
            errors <- c(errors, "Saved UMAP output-dimension metadata is internally inconsistent.")
        }
    }
    if (length(dimensions) == 0) {
        errors <- c(errors, "Saved UMAP configuration has no valid output dimension metadata.")
    }

    list(
        valid = length(errors) == 0,
        errors = unique(errors),
        warnings = warnings,
        columns = columns,
        dimensions = dimensions
    )
}

#' Prepare Data for Reproducible Saved-UMAP Projection
#'
#' @param new_data A data frame containing the exact prepared training columns.
#' @param pipeline_config The saved pipeline configuration.
#' @return A list with projection matrix, valid row indices, warnings, and error.
prepare_saved_umap_projection <- function(new_data, pipeline_config) {
    tryCatch({
        contract <- validate_umap_projection_contract(pipeline_config)
        if (!contract$valid) {
            return(list(
                data = NULL,
                valid_idx = integer(0),
                warnings = contract$warnings,
                error = paste(contract$errors, collapse = "\n")
            ))
        }

        if (!is.data.frame(new_data)) {
            stop("Projection input must be a data frame.")
        }
        if (nrow(new_data) == 0) {
            stop("Projection input is empty.")
        }

        missing_columns <- setdiff(contract$columns, names(new_data))
        if (length(missing_columns) > 0) {
            stop(paste("Missing required prepared columns:", paste(missing_columns, collapse = ", ")))
        }

        non_numeric <- contract$columns[!vapply(new_data[contract$columns], is.numeric, logical(1))]
        if (length(non_numeric) > 0) {
            stop(paste("Required projection columns must be numeric:", paste(non_numeric, collapse = ", ")))
        }

        projection_matrix <- as.matrix(new_data[, contract$columns, drop = FALSE])
        excluded <- rep(FALSE, nrow(new_data))
        if ("row_excluded" %in% names(new_data)) {
            excluded_values <- suppressWarnings(as.logical(new_data$row_excluded))
            excluded <- !is.na(excluded_values) & excluded_values
        }

        finite_rows <- apply(projection_matrix, 1, function(row) all(is.finite(row)))
        valid <- !excluded & complete.cases(projection_matrix) & finite_rows
        valid_idx <- which(valid)
        if (length(valid_idx) == 0) {
            stop("No non-excluded rows have complete finite values in every required prepared column.")
        }

        warnings <- character(0)
        skipped <- nrow(new_data) - length(valid_idx)
        if (skipped > 0) {
            warnings <- c(warnings, sprintf(
                "%d row(s) will be skipped because they are excluded or lack complete finite prepared values.",
                skipped
            ))
        }
        if (length(valid_idx) < 5) {
            warnings <- c(warnings, "Fewer than five rows are available for projection.")
        }

        list(
            data = projection_matrix[valid_idx, , drop = FALSE],
            valid_idx = valid_idx,
            warnings = warnings,
            error = NULL
        )
    }, error = function(e) {
        list(data = NULL, valid_idx = integer(0), warnings = character(0), error = e$message)
    })
}

#' Validate Data for UMAP Projection
#'
#' @param new_data A data frame to validate.
#' @param pipeline_config The saved pipeline configuration from a UMAP model.
#' @return A list with `valid`, `errors`, and `warnings`.
#' @export
validate_for_projection <- function(new_data, pipeline_config) {
    result <- prepare_saved_umap_projection(new_data, pipeline_config)
    list(
        valid = is.null(result$error),
        errors = if (is.null(result$error)) character(0) else result$error,
        warnings = result$warnings
    )
}

#' Apply a Saved UMAP Projection Contract
#'
#' This compatibility wrapper no longer re-estimates preprocessing from the new
#' batch. It returns the exact prepared columns recorded when the model was fit.
#'
#' @param new_data A data frame to transform.
#' @param pipeline_config The saved pipeline configuration.
#' @return A list with `data`, `valid_idx`, `warnings`, and `error`.
#' @export
apply_saved_pipeline <- function(new_data, pipeline_config) {
    prepare_saved_umap_projection(new_data, pipeline_config)
}

#' Normalize a Saved UMAP Model Collection
#'
#' @param models A named `1d`/`2d` model collection or one legacy model.
#' @param umap_config The matching UMAP configuration.
#' @return A named list of models.
normalize_umap_model_collection <- function(models, umap_config = NULL) {
    if (is.null(models)) return(list())

    model_names <- names(models)
    is_collection <- is.list(models) && length(models) > 0 &&
        !is.null(model_names) && all(nzchar(model_names)) &&
        all(model_names %in% c("1d", "2d"))
    if (is_collection) return(models)

    dimensions <- normalize_umap_dimensions(umap_config)
    if (length(dimensions) != 1) {
        stop("A legacy UMAP model cannot represent multiple saved dimensions.")
    }
    stats::setNames(list(models), paste0(dimensions, "d"))
}

#' Saved UMAP Model Filenames
#'
#' @param umap_config The matching UMAP configuration.
#' @return A named character vector of bundle filenames.
umap_model_filenames <- function(umap_config) {
    dimensions <- normalize_umap_dimensions(umap_config)
    if (length(dimensions) == 0) stop("UMAP configuration has no saved dimensions.")
    stats::setNames(paste0("model_", dimensions, "d.uwot"), paste0(dimensions, "d"))
}

#' Save a UMAP Model Collection
#'
#' @param models A named model collection.
#' @param umap_config The matching UMAP configuration.
#' @param directory Existing destination directory.
#' @return The saved model basenames.
save_umap_model_collection <- function(models, umap_config, directory) {
    collection <- normalize_umap_model_collection(models, umap_config)
    filenames <- umap_model_filenames(umap_config)
    if (!setequal(names(collection), names(filenames))) {
        stop("Saved UMAP models do not match the configured output dimensions.")
    }

    written <- character(0)
    tryCatch({
        for (name in names(filenames)) {
            path <- file.path(directory, filenames[[name]])
            uwot::save_uwot(collection[[name]], path, unload = FALSE)
            written <- c(written, path)
        }
        unname(filenames)
    }, error = function(e) {
        if (length(written) > 0) unlink(written, recursive = TRUE, force = TRUE)
        stop(e)
    })
}

#' Load a UMAP Model Collection
#'
#' @param directory Bundle directory.
#' @param umap_config The matching UMAP configuration.
#' @return A named list of loaded models.
load_umap_model_collection <- function(directory, umap_config) {
    filenames <- umap_model_filenames(umap_config)
    paths <- file.path(directory, filenames)
    missing <- filenames[!file.exists(paths)]
    if (length(missing) > 0) {
        stop(paste("Saved UMAP bundle is missing:", paste(missing, collapse = ", ")))
    }

    loaded <- list()
    tryCatch({
        for (name in names(filenames)) {
            loaded[[name]] <- uwot::load_uwot(file.path(directory, filenames[[name]]))
        }
        loaded
    }, error = function(e) {
        unload_umap_model_collection(loaded)
        stop(e)
    })
}

#' Unload a UMAP Model Collection
#'
#' @param models A model collection.
#' @param umap_config Optional matching configuration for legacy models.
#' @return Invisibly returns `NULL`.
unload_umap_model_collection <- function(models, umap_config = NULL) {
    collection <- tryCatch(
        normalize_umap_model_collection(models, umap_config),
        error = function(e) list()
    )
    for (model in collection) {
        try(uwot::unload_uwot(model), silent = TRUE)
    }
    invisible(NULL)
}

#' Create a Versioned Token for Mutable Analysis Data
#'
#' @param data Current analysis data.
#' @param generation Monotonic state generation.
#' @param context Optional analysis settings that must still match at completion.
#' @return A list containing generation and a serialized data fingerprint.
varg_make_data_token <- function(data, generation, context = NULL) {
    raw_generation <- unlist(generation, use.names = FALSE)
    generation <- if (length(raw_generation) == 0) 0L else {
        suppressWarnings(as.integer(raw_generation[[1]]))
    }
    if (length(generation) != 1 || is.na(generation)) generation <- 0L
    list(
        generation = generation,
        fingerprint = digest::digest(data, algo = "xxhash64", serialize = TRUE),
        context_fingerprint = digest::digest(context, algo = "xxhash64", serialize = TRUE)
    )
}

#' Check Whether a Background Result Still Matches Live Data
#'
#' @param token Token captured when work started.
#' @param data Current analysis data.
#' @param generation Current state generation.
#' @param context Current analysis settings.
#' @return Logical scalar.
varg_data_token_matches <- function(token, data, generation, context = NULL) {
    if (is.null(token) || is.null(token$generation) || is.null(token$fingerprint) ||
        is.null(token$context_fingerprint)) return(FALSE)
    current <- varg_make_data_token(data, generation, context)
    identical(as.integer(token$generation), current$generation) &&
        identical(as.character(token$fingerprint), current$fingerprint) &&
        identical(as.character(token$context_fingerprint), current$context_fingerprint)
}
