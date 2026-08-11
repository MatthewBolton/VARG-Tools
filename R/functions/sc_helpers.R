# sc_helpers.R
# Helper functions for Stratigraphic Correlation Module

# ==============================================================================
# DATA IMPORT AND VALIDATION
# ==============================================================================

#' Load data file (CSV or XLSX)
#' @param path File path
#' @param file_name Original upload name, used when Shiny stores the upload in
#'   a temporary path without an extension
#' @return data.frame or NULL on error
sc_load_file <- function(path, file_name = path) {
    tryCatch(
        {
            file_ext <- tolower(tools::file_ext(file_name))

            data <- if (file_ext == "csv") {
                header <- readLines(path, n = 1, warn = FALSE)
                tab_count <- lengths(regmatches(header, gregexpr("\t", header, fixed = TRUE)))
                comma_count <- lengths(regmatches(header, gregexpr(",", header, fixed = TRUE)))
                if (length(tab_count) > 0 && tab_count[[1]] > comma_count[[1]]) {
                    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
                } else {
                    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
                }
            } else if (file_ext == "xlsx") {
                readxl::read_excel(path)
            } else {
                stop("Unsupported file type. Please upload .csv or .xlsx")
            }

            # Convert to data.frame if tibble
            as.data.frame(data)
        },
        error = function(e) {
            warning(paste("Error loading file:", e$message))
            NULL
        }
    )
}

#' Smart column detection for SC module
#' Attempts to match columns to expected types using pattern matching
#' @param data data.frame with columns to detect
#' @return Named list with site_id, sample_name, z_var, x_var (or NULL for undetected)
sc_detect_columns <- function(data) {
    cols <- names(data)
    cols_lower <- tolower(cols)
    nums <- names(data)[sapply(data, is.numeric)]
    nums_lower <- tolower(nums)
    
    result <- list(site_id = NULL, sample_name = NULL, z_var = NULL, x_var = NULL)
    
    # Patterns for each column type (in priority order)
    site_patterns <- c("^coreid$", "^core_id$", "^core$", "^site$", "^siteid$", "^site_id$", 
                       "^location$", "^locality$", "^section$", "^unit$")
    sample_patterns <- c("^sample$", "^sample_id$", "^sampleid$", "^sample_name$", 
                         "^samplename$", "^id$", "^name$", "^specimen$", "^label$")
    depth_patterns <- c("^depth$", "^depth_m$", "^depth_cm$", "^depth_mm$", "^z$", 
                        "^age$", "^age_ka$", "^age_ma$", "^age_bp$", "^elevation$", 
                        "^height$", "^level$", "^position$", "^strat_position$")
    x_patterns <- c("^x$", "^umap1$", "^umap_1$", "^pc1$", "^pca1$", "^comp1$", 
                    "^component1$", "^axis1$", "^dim1$", "^score1$")
    
    # Function to find first matching column
    find_match <- function(patterns, col_names, col_names_lower) {
        for (pattern in patterns) {
            matches <- grep(pattern, col_names_lower, value = FALSE)
            if (length(matches) > 0) {
                return(col_names[matches[1]])  # Return original case
            }
        }
        NULL
    }
    
    # Detect each column type
    result$site_id <- find_match(site_patterns, cols, cols_lower)
    result$sample_name <- find_match(sample_patterns, cols, cols_lower)
    result$z_var <- find_match(depth_patterns, nums, nums_lower)
    result$x_var <- find_match(x_patterns, nums, nums_lower)
    
    # Fallbacks: use first available column if not detected
    if (is.null(result$site_id) && length(cols) >= 1) {
        result$site_id <- cols[1]
    }
    if (is.null(result$sample_name) && length(cols) >= 2) {
        # Use second column, or first if only one
        result$sample_name <- if (length(cols) >= 2) cols[2] else cols[1]
    }
    if (is.null(result$z_var) && length(nums) >= 1) {
        result$z_var <- nums[1]
    }
    if (is.null(result$x_var) && length(nums) >= 2) {
        result$x_var <- nums[2]
    }
    
    result
}

#' Validate column mappings
#' @param data data.frame
#' @param mappings Named list with site_id, sample_name, x_var, z_var
#' @return List with valid (TRUE/FALSE) and message
sc_validate_columns <- function(data, mappings) {
    required <- c("site_id", "sample_name", "x_var", "z_var")
    missing <- setdiff(required, names(mappings))

    if (length(missing) > 0) {
        return(list(
            valid = FALSE,
            message = paste("Missing mappings:", paste(missing, collapse = ", "))
        ))
    }

    # Check all mapped columns exist
    mapped_cols <- unlist(mappings)
    missing_cols <- setdiff(mapped_cols, names(data))

    if (length(missing_cols) > 0) {
        return(list(
            valid = FALSE,
            message = paste("Columns not found in data:", paste(missing_cols, collapse = ", "))
        ))
    }

    # Check x and z are numeric
    if (!is.numeric(data[[mappings$x_var]])) {
        return(list(
            valid = FALSE,
            message = paste("X variable", mappings$x_var, "must be numeric")
        ))
    }

    if (!is.numeric(data[[mappings$z_var]])) {
        return(list(
            valid = FALSE,
            message = paste("Z variable", mappings$z_var, "must be numeric")
        ))
    }

    list(valid = TRUE, message = "All columns validated successfully")
}

#' Normalize identifier columns used by stratigraphic plots and joins
#' @param data Input data.frame
#' @param mappings Column mappings containing site_id and sample_name
#' @return A copy with mapped identifier columns stored as character vectors
sc_normalize_identifier_columns <- function(data, mappings) {
    identifier_columns <- unique(unname(unlist(mappings[c("site_id", "sample_name")])))
    identifier_columns <- identifier_columns[
        !is.na(identifier_columns) & identifier_columns %in% names(data)
    ]
    for (column in identifier_columns) {
        data[[column]] <- as.character(data[[column]])
    }
    data
}

#' Convert the stratigraphic plot-height mode to a CSS height
#' @param mode One of standard, tall, or double
#' @return CSS height string
sc_plot_height_css <- function(mode = "standard") {
    if (is.null(mode) || length(mode) == 0 || !mode[[1]] %in% c("standard", "tall", "double")) {
        mode <- "standard"
    }
    switch(
        mode[[1]],
        tall = "112.5vh",
        double = "150vh",
        "75vh"
    )
}

#' Calculate a display-only affine depth transformation
#' @param data Original stratigraphic data
#' @param mappings Column mappings
#' @param ref_core Reference core ID
#' @param target_core Target core ID
#' @param mode Either auto or manual
#' @param manual_scale Manual target-depth multiplier
#' @param manual_shift Manual target-depth offset
#' @param ref_direction Whether reference values increase down or up
#' @param target_direction Whether target values increase down or up
#' @return List containing scale, shift, and source depth ranges
sc_affine_depth_parameters <- function(data, mappings, ref_core, target_core,
                                       mode = "auto", manual_scale = 1,
                                       manual_shift = 0,
                                       ref_direction = "down",
                                       target_direction = "down") {
    data <- sc_normalize_identifier_columns(data, mappings)
    mode <- if (identical(mode, "manual")) "manual" else "auto"

    normalize_direction <- function(value, label) {
        value <- as.character(value)[[1]]
        if (!value %in% c("down", "up")) {
            stop(label, " direction must be either 'down' or 'up'.")
        }
        value
    }
    ref_direction <- normalize_direction(ref_direction, "Reference")
    target_direction <- normalize_direction(target_direction, "Target")
    direction_sign <- if (identical(ref_direction, target_direction)) 1 else -1

    finite_range <- function(values, label) {
        values <- suppressWarnings(as.numeric(values))
        values <- values[is.finite(values)]
        if (length(values) == 0) {
            stop(label, " has no finite depth/age values.")
        }
        range(values)
    }

    site_values <- as.character(data[[mappings$site_id]])
    ref_range <- finite_range(
        data[[mappings$z_var]][site_values == as.character(ref_core)],
        "Reference core"
    )
    target_range <- finite_range(
        data[[mappings$z_var]][site_values == as.character(target_core)],
        "Target core"
    )

    if (identical(mode, "manual")) {
        scale_values <- suppressWarnings(as.numeric(manual_scale))
        shift_values <- suppressWarnings(as.numeric(manual_shift))
        scale_magnitude <- if (length(scale_values) > 0) scale_values[[1]] else NA_real_
        shift <- if (length(shift_values) > 0) shift_values[[1]] else NA_real_
        if (!is.finite(scale_magnitude) || scale_magnitude <= 0) {
            stop("Target Z scale magnitude must be a finite value greater than zero.")
        }
        if (!is.finite(shift)) {
            stop("Target Z shift must be finite.")
        }
        scale <- direction_sign * scale_magnitude
    } else {
        ref_span <- diff(ref_range)
        target_span <- diff(target_range)
        if (!is.finite(ref_span) || ref_span <= sqrt(.Machine$double.eps)) {
            stop("Reference core depth/age range is too small for automatic range matching.")
        }
        if (!is.finite(target_span) || target_span <= sqrt(.Machine$double.eps)) {
            stop("Target core depth/age range is too small for automatic range matching.")
        }
        scale <- direction_sign * ref_span / target_span
        ref_top <- if (identical(ref_direction, "down")) ref_range[[1]] else ref_range[[2]]
        target_top <- if (identical(target_direction, "down")) target_range[[1]] else target_range[[2]]
        shift <- ref_top - scale * target_top
    }

    list(
        mode = mode,
        scale = unname(scale),
        shift = unname(shift),
        ref_range = unname(ref_range),
        target_range = unname(target_range),
        ref_direction = ref_direction,
        target_direction = target_direction
    )
}

#' Return the Plotly autorange setting for a stratigraphic direction
#' @param direction Either down or up
#' @return "reversed" when values increase downward, otherwise TRUE
sc_plot_y_autorange <- function(direction = "down") {
    if (identical(direction, "up")) TRUE else "reversed"
}

#' Resolve a Plotly click to canonical stratigraphic coordinates
#'
#' Plot coordinates may be warped or jittered for display. Tie points must
#' always retain the original depth from the uploaded data.
#' @param click Plotly click event
#' @param data Original stratigraphic data
#' @param mappings Column mappings
#' @param ref_core Reference core ID
#' @param target_core Target core ID
#' @param show_points Whether point-cloud traces are displayed
#' @return List with sample_name, role, is_ref, and original z
sc_resolve_plot_click <- function(click, data, mappings, ref_core, target_core,
                                  show_points = FALSE) {
    if (is.null(click) || is.null(click$key) || length(click$key) == 0) {
        stop("Click a plotted sample point, not a line or annotation.")
    }

    sample_name <- as.character(click$key[[1]])
    if (is.na(sample_name) || !nzchar(sample_name)) {
        stop("The selected point does not have a sample identifier.")
    }

    role <- NA_character_
    if (!is.null(click$customdata)) {
        custom_values <- as.character(unlist(click$customdata, use.names = FALSE))
        recognized_roles <- custom_values[custom_values %in% c("ref", "target")]
        if (length(recognized_roles) > 0) {
            role <- recognized_roles[[1]]
        }
    }

    if (is.na(role)) {
        curve_num <- suppressWarnings(as.integer(click$curveNumber[[1]]))
        if (length(curve_num) == 0 || is.na(curve_num)) {
            stop("The selected point could not be assigned to a core.")
        }
        is_ref_curve <- if (isTRUE(show_points)) {
            curve_num %in% c(0L, 2L)
        } else {
            curve_num == 0L
        }
        role <- if (is_ref_curve) "ref" else "target"
    }

    core_id <- if (identical(role, "ref")) ref_core else target_core
    site_values <- as.character(data[[mappings$site_id]])
    sample_values <- as.character(data[[mappings$sample_name]])
    matched <- !is.na(site_values) & !is.na(sample_values) &
        site_values == as.character(core_id) & sample_values == sample_name
    depths <- data[[mappings$z_var]][matched]
    depths <- depths[is.finite(depths)]

    if (length(depths) == 0) {
        stop(paste0(
            "Original depth was not found for sample '", sample_name,
            "' in core '", core_id, "'."
        ))
    }

    list(
        sample_name = sample_name,
        role = role,
        is_ref = identical(role, "ref"),
        z = median(depths)
    )
}

# ==============================================================================
# TIE-POINT MANAGEMENT
# ==============================================================================

#' Create a new tie point
#' @param id Unique ID
#' @return Single-row data.frame
sc_create_tiepoint <- function(id = 1) {
    data.frame(
        id = id,
        ref_sample = NA_character_,
        ref_z = NA_real_,
        target_sample = NA_character_,
        target_z = NA_real_,
        custom_name = NA_character_,
        use_in_warp = TRUE,
        stringsAsFactors = FALSE
    )
}

#' Validate tie points for warping
#' @param tiepoints data.frame with tie point data
#' @param ref_direction Whether reference values increase down or up
#' @param target_direction Whether target values increase down or up
#' @return List with valid (TRUE/FALSE), message, valid_count, warnings
sc_validate_tiepoints <- function(tiepoints, ref_direction = "down",
                                  target_direction = "down") {
    if (is.null(tiepoints) || nrow(tiepoints) == 0) {
        return(list(
            valid = FALSE,
            message = "No tie points defined",
            valid_count = 0,
            warnings = character()
        ))
    }

    # Filter to those marked for use in warp
    valid_ties <- tiepoints[
        !is.na(tiepoints$use_in_warp) & tiepoints$use_in_warp &
            !is.na(tiepoints$ref_sample) & nzchar(trimws(as.character(tiepoints$ref_sample))) &
            !is.na(tiepoints$target_sample) & nzchar(trimws(as.character(tiepoints$target_sample))) &
            is.finite(tiepoints$ref_z) &
            is.finite(tiepoints$target_z),
        , drop = FALSE
    ]

    valid_count <- nrow(valid_ties)
    warnings <- character()
    errors <- character()

    if (valid_count < 2) {
        return(list(
            valid = FALSE,
            message = paste("Need at least 2 valid tie points for warping. Currently have:", valid_count),
            valid_count = valid_count,
            warnings = warnings
        ))
    }

    direction_sign <- function(direction, label) {
        if (!direction %in% c("down", "up")) {
            stop(label, " direction must be either 'down' or 'up'.")
        }
        if (identical(direction, "down")) 1 else -1
    }
    ref_sign <- direction_sign(ref_direction, "Reference")
    target_sign <- direction_sign(target_direction, "Target")

    # Check for crossing tie points in physical top-to-bottom order.
    if (valid_count >= 2) {
        ref_order <- ref_sign * valid_ties$ref_z
        sorted_ties <- valid_ties[order(ref_order), ]
        target_z_ordered <- target_sign * sorted_ties$target_z
        
        # Check for inversions (crossing lines)
        for (i in 2:length(target_z_ordered)) {
            if (target_z_ordered[i] < target_z_ordered[i-1]) {
                # Found a crossing - identify the problematic tie points
                tp1_id <- sorted_ties$id[i-1]
                tp2_id <- sorted_ties$id[i]
                errors <- c(errors, paste0(
                    "Crossing tie points detected: TP", tp1_id, " and TP", tp2_id,
                    " violate stratigraphic order. Tie points must not cross."
                ))
            }
        }
    }

    # Check for duplicate sample assignments (one sample linked to multiple others)
    if (valid_count >= 2) {
        # Check reference samples - each should only appear once
        ref_sample_counts <- table(valid_ties$ref_sample[!is.na(valid_ties$ref_sample) & valid_ties$ref_sample != ""])
        dup_refs <- names(ref_sample_counts[ref_sample_counts > 1])
        if (length(dup_refs) > 0) {
            errors <- c(errors, paste0(
                "Duplicate reference sample(s): ", paste(dup_refs, collapse = ", "),
                ". Each reference sample can only be linked once."
            ))
        }
        
        # Check target samples - each should only appear once
        target_sample_counts <- table(valid_ties$target_sample[!is.na(valid_ties$target_sample) & valid_ties$target_sample != ""])
        dup_targets <- names(target_sample_counts[target_sample_counts > 1])
        if (length(dup_targets) > 0) {
            errors <- c(errors, paste0(
                "Duplicate target sample(s): ", paste(dup_targets, collapse = ", "),
                ". Each target sample can only be linked once."
            ))
        }
    }

    # If we have errors, validation fails
    if (length(errors) > 0) {
        return(list(
            valid = FALSE,
            message = paste(errors, collapse = " | "),
            valid_count = valid_count,
            warnings = warnings
        ))
    }

    list(
        valid = TRUE,
        message = paste("Ready to warp with", valid_count, "tie points"),
        valid_count = valid_count,
        warnings = warnings
    )
}

# ==============================================================================
# WARP CALCULATIONS
# ==============================================================================

#' Calculate warp function from tie points
#' @param tiepoints data.frame with ref_z and target_z columns
#' @param method "auto", "linear", "monotonic", or legacy "gam"
#' @param ref_direction Optional reference coordinate direction ("down" or
#'   "up") used to check the expected mapping direction
#' @param target_direction Optional target coordinate direction ("down" or
#'   "up") used to check the expected mapping direction
#' @param extrapolation How to handle target values outside the tie-point
#'   range: "linear" continues the selected model (the nearest tied interval
#'   for monotonic piecewise warping or the fitted slope for least-squares
#'   warping), "error" stops, and "na" returns NA
#' @return List with warp_func (function), model_type, model_object, and diagnostics
sc_calculate_warp <- function(tiepoints, method = "auto",
                              ref_direction = NULL,
                              target_direction = NULL,
                              extrapolation = c("linear", "error", "na")) {
    # Filter to valid tie points
    valid_ties <- tiepoints[
        !is.na(tiepoints$use_in_warp) & tiepoints$use_in_warp &
            !is.na(tiepoints$ref_sample) & nzchar(trimws(as.character(tiepoints$ref_sample))) &
            !is.na(tiepoints$target_sample) & nzchar(trimws(as.character(tiepoints$target_sample))) &
            is.finite(tiepoints$ref_z) &
            is.finite(tiepoints$target_z),
        , drop = FALSE
    ]

    if (nrow(valid_ties) < 2) {
        stop("Need at least 2 valid tie points for warping")
    }

    method <- match.arg(method, c("auto", "linear", "monotonic", "gam"))
    extrapolation <- match.arg(extrapolation)

    normalize_direction <- function(direction, label) {
        if (is.null(direction)) {
            return(NULL)
        }
        value <- as.character(direction)[[1]]
        if (is.na(value) || !value %in% c("down", "up")) {
            stop(label, " direction must be either 'down' or 'up'.")
        }
        value
    }

    ref_direction <- normalize_direction(ref_direction, "Reference")
    target_direction <- normalize_direction(target_direction, "Target")

    # Automatic and legacy GAM requests use an exact-knot monotonic mapping.
    # Least-squares alignment remains available only when explicitly requested.
    n_points <- nrow(valid_ties)
    use_method <- if (method == "auto") {
        "monotonic"
    } else if (method == "gam") {
        "monotonic"
    } else {
        method
    }

    # Fit model
    if (use_method == "linear") {
        if (!all(is.finite(valid_ties$target_z)) ||
            !all(is.finite(valid_ties$ref_z))) {
            stop("Linear warping requires finite tie-point coordinates.")
        }

        ordered_ties <- valid_ties[order(valid_ties$target_z), , drop = FALSE]
        target_knots <- ordered_ties$target_z
        ref_knots <- ordered_ties$ref_z
        ref_tolerance <- 100 * .Machine$double.eps *
            max(1, max(abs(ref_knots)))

        duplicate_targets <- unique(target_knots[duplicated(target_knots)])
        if (length(duplicate_targets) > 0) {
            for (target_value in duplicate_targets) {
                duplicate_refs <- ref_knots[target_knots == target_value]
                if (max(duplicate_refs) - min(duplicate_refs) > ref_tolerance) {
                    stop("Duplicate target tie coordinates have conflicting reference values.")
                }
            }
        }

        ref_differences <- diff(ref_knots)
        nonzero_differences <- ref_differences[
            abs(ref_differences) > ref_tolerance
        ]
        expected_sign <- if (!is.null(ref_direction) &&
            !is.null(target_direction)) {
            if (identical(ref_direction, target_direction)) 1 else -1
        } else if (length(nonzero_differences) == 0) {
            1
        } else if (all(nonzero_differences > 0)) {
            1
        } else if (all(nonzero_differences < 0)) {
            -1
        } else {
            0
        }
        if (expected_sign == 0 ||
            any(expected_sign * ref_differences < -ref_tolerance)) {
            stop("Tie points must preserve a monotonic target-to-reference order.")
        }

        target_range <- range(target_knots)
        fit <- lm(ref_z ~ target_z, data = valid_ties)
        warp_func <- function(z) {
            if (!is.numeric(z)) {
                stop("Warp coordinates must be numeric.")
            }

            z_names <- names(z)
            z <- as.numeric(z)
            result <- rep(NA_real_, length(z))
            finite <- is.finite(z)
            if (any(finite)) {
                finite_indices <- which(finite)
                finite_values <- z[finite_indices]
                outside <- finite_values < target_range[[1]] |
                    finite_values > target_range[[2]]
                if (extrapolation == "error" && any(outside)) {
                    stop("Warp coordinate is outside the tie-point target range; extrapolation is disabled.")
                }
                keep <- if (extrapolation == "na") !outside else rep(TRUE, length(outside))
                if (any(keep)) {
                    result[finite_indices[keep]] <- unname(predict(
                        fit,
                        newdata = data.frame(target_z = finite_values[keep])
                    ))
                }
            }
            if (!is.null(z_names)) names(result) <- z_names
            result
        }
        attr(warp_func, "target_range") <- unname(target_range)
        attr(warp_func, "extrapolation") <- extrapolation
        model_type <- "Linear (LM)"
        model_obj <- fit
    } else {
        if (!all(is.finite(valid_ties$target_z)) ||
            !all(is.finite(valid_ties$ref_z))) {
            stop("Monotonic warping requires finite tie-point coordinates.")
        }

        ordered_ties <- valid_ties[order(valid_ties$target_z), , drop = FALSE]
        target_knots <- ordered_ties$target_z
        ref_knots <- ordered_ties$ref_z
        ref_tolerance <- 100 * .Machine$double.eps *
            max(1, max(abs(ref_knots)))

        # A function cannot map one target coordinate to two different
        # reference coordinates. Identical pairs are harmless and collapse to
        # one knot; conflicting pairs are rejected rather than averaged.
        duplicate_targets <- unique(target_knots[duplicated(target_knots)])
        if (length(duplicate_targets) > 0) {
            for (target_value in duplicate_targets) {
                duplicate_refs <- ref_knots[target_knots == target_value]
                if (max(duplicate_refs) - min(duplicate_refs) > ref_tolerance) {
                    stop("Duplicate target tie coordinates have conflicting reference values.")
                }
            }
            keep <- !duplicated(target_knots)
            target_knots <- target_knots[keep]
            ref_knots <- ref_knots[keep]
        }

        if (length(target_knots) < 2) {
            stop("Monotonic warping requires at least two distinct target coordinates.")
        }

        ref_differences <- diff(ref_knots)
        nonzero_differences <- ref_differences[
            abs(ref_differences) > ref_tolerance
        ]
        expected_sign <- if (!is.null(ref_direction) &&
            !is.null(target_direction)) {
            if (identical(ref_direction, target_direction)) 1 else -1
        } else if (length(nonzero_differences) == 0) {
            1
        } else if (all(nonzero_differences > 0)) {
            1
        } else if (all(nonzero_differences < 0)) {
            -1
        } else {
            0
        }

        if (expected_sign == 0 ||
            any(expected_sign * ref_differences < -ref_tolerance)) {
            stop("Tie points must preserve a monotonic target-to-reference order.")
        }

        target_steps <- diff(target_knots)
        slopes <- diff(ref_knots) / target_steps
        target_range <- range(target_knots)

        map_finite_values <- function(values) {
            segment <- findInterval(values, target_knots, all.inside = FALSE)
            segment <- pmin(pmax(segment, 1L), length(target_knots) - 1L)
            ref_knots[segment] + slopes[segment] *
                (values - target_knots[segment])
        }

        warp_func <- function(z) {
            if (!is.numeric(z)) {
                stop("Warp coordinates must be numeric.")
            }

            z_names <- names(z)
            z <- as.numeric(z)
            result <- rep(NA_real_, length(z))
            finite <- is.finite(z)

            if (any(finite)) {
                finite_indices <- which(finite)
                finite_values <- z[finite_indices]
                outside <- finite_values < target_range[[1]] |
                    finite_values > target_range[[2]]

                if (extrapolation == "error" && any(outside)) {
                    stop("Warp coordinate is outside the tie-point target range; extrapolation is disabled.")
                }

                keep <- if (extrapolation == "na") !outside else rep(TRUE, length(outside))
                if (any(keep)) {
                    result[finite_indices[keep]] <- map_finite_values(
                        finite_values[keep]
                    )
                }
            }

            if (!is.null(z_names)) {
                names(result) <- z_names
            }
            result
        }

        attr(warp_func, "target_range") <- unname(target_range)
        attr(warp_func, "extrapolation") <- extrapolation
        model_type <- "Monotonic piecewise-linear"
        model_obj <- list(
            target_knots = unname(target_knots),
            ref_knots = unname(ref_knots),
            slopes = unname(slopes),
            target_range = unname(target_range),
            extrapolation = extrapolation,
            ref_direction = ref_direction,
            target_direction = target_direction
        )
        class(model_obj) <- c("stratigraphic_monotonic_warp", "list")
    }

    # Calculate diagnostics
    valid_ties$pred_z <- warp_func(valid_ties$target_z)
    valid_ties$residual <- valid_ties$ref_z - valid_ties$pred_z

    rmse <- sqrt(mean(valid_ties$residual^2))
    mae <- mean(abs(valid_ties$residual))

    list(
        warp_func = warp_func,
        model_type = model_type,
        model_object = model_obj,
        diagnostics = list(
            rmse = rmse,
            mae = mae,
            n_points = n_points,
            target_range = if (exists("target_range")) {
                unname(target_range)
            } else {
                unname(range(valid_ties$target_z))
            },
            extrapolation = extrapolation,
            residuals = valid_ties[, c(
                "id", "ref_sample", "target_sample",
                "ref_z", "target_z", "pred_z", "residual"
            )]
        )
    )
}

#' Apply warp to target data
#' @param data data.frame with target core data
#' @param z_col Name of depth/age column
#' @param warp_func Warping function
#' @return data.frame with added z_warped column
sc_apply_warp <- function(data, z_col, warp_func) {
    data$z_warped <- warp_func(data[[z_col]])
    data
}

# ==============================================================================
# OUTPUT OBJECT CONSTRUCTION
# ==============================================================================

#' Create structured output object for chronology module
#' @param tiepoints data.frame of tie points
#' @param warp_result Result from sc_calculate_warp
#' @param metadata List with cores, variables, timestamp
#' @return S3 object of class "stratigraphic_warp"
sc_create_output <- function(tiepoints, warp_result, metadata) {
    obj <- list(
        tiepoints = tiepoints[
            !is.na(tiepoints$use_in_warp) & tiepoints$use_in_warp &
                !is.na(tiepoints$ref_sample) & nzchar(trimws(as.character(tiepoints$ref_sample))) &
                !is.na(tiepoints$target_sample) & nzchar(trimws(as.character(tiepoints$target_sample))) &
                is.finite(tiepoints$ref_z) &
                is.finite(tiepoints$target_z),
            , drop = FALSE
        ],
        warp_function = warp_result$warp_func,
        model_type = warp_result$model_type,
        model_object = warp_result$model_object,
        diagnostics = warp_result$diagnostics,
        metadata = c(metadata, list(created = Sys.time()))
    )

    class(obj) <- c("stratigraphic_warp", "list")
    obj
}

#' Print method for stratigraphic_warp
#' @param x stratigraphic_warp object
#' @param ... Additional arguments
print.stratigraphic_warp <- function(x, ...) {
    cat("Stratigraphic Warp Object\n")
    cat("========================\n\n")
    cat("Reference Core:", x$metadata$reference_core, "\n")
    cat("Target Core:", x$metadata$target_core, "\n")
    cat("Model Type:", x$model_type, "\n")
    cat("Tie Points:", nrow(x$tiepoints), "\n")
    cat("RMSE:", round(x$diagnostics$rmse, 4), "\n")
    cat("MAE:", round(x$diagnostics$mae, 4), "\n")
    cat("Created:", format(x$metadata$created, "%Y-%m-%d %H:%M:%S"), "\n")
}

#' Convert to JSON for clipboard
#' @param obj stratigraphic_warp object
#' @return JSON string
sc_to_json <- function(obj) {
    # Create simplified version for JSON export
    export_obj <- list(
        reference_core = obj$metadata$reference_core,
        target_core = obj$metadata$target_core,
        model_type = obj$model_type,
        tiepoints = obj$tiepoints[, c("ref_sample", "ref_z", "target_sample", "target_z")],
        diagnostics = list(
            rmse = obj$diagnostics$rmse,
            mae = obj$diagnostics$mae,
            n_points = obj$diagnostics$n_points
        ),
        created = format(obj$metadata$created, "%Y-%m-%d %H:%M:%S")
    )

    jsonlite::toJSON(export_obj, pretty = TRUE, auto_unbox = TRUE)
}

# ==============================================================================
# PLOT GENERATION HELPERS
# ==============================================================================

#' Generate initial alignment plot
#' @param data data.frame with all core data
#' @param mappings Column mappings
#' @param ref_core Reference core ID
#' @param target_core Target core ID
#' @param tiepoints data.frame of tie points
#' @param selected_ids Vector of selected tie point IDs (for highlighting)
#' @param show_points Logical, whether to show raw point cloud behind median lines
#' @param plot_source Character, source ID for plotly click events
#' @param highlight_sample Character, sample name to highlight with a larger point marker (for click-to-add mode)
#' @param target_offset Numeric, horizontal offset to apply to target data for visual separation
#' @param line_width Numeric width of the reference and target profile lines
#' @param point_size Numeric size of the reference and target profile markers
#' @param reference_direction Whether reference values increase down or up
#' @return plotly object
sc_plot_initial <- function(data, mappings, ref_core, target_core,
                            tiepoints = NULL, selected_ids = NULL,
                            show_points = FALSE, plot_source = "sc_initial",
                            highlight_sample = NULL,
                            target_offset = 0,
                            line_width = 2,
                            point_size = 8,
                            jitter_x = 0,
                            jitter_z = 0,
                            reference_direction = "down") {
    library(plotly)
    library(dplyr)

    data <- sc_normalize_identifier_columns(data, mappings)

    # Filter data
    ref_data <- data %>%
        filter(.data[[mappings$site_id]] == ref_core) %>%
        arrange(.data[[mappings$z_var]])

    target_data <- data %>%
        filter(.data[[mappings$site_id]] == target_core) %>%
        arrange(.data[[mappings$z_var]])

    # Apply target offset to target data x values for display
    target_data_display <- target_data
    if (target_offset != 0) {
        target_data_display <- target_data %>%
            mutate(!!mappings$x_var := .data[[mappings$x_var]] + target_offset)
    }
    
    # Apply jitter if specified (for display only)
    set.seed(42)  # Reproducible jitter
    if (jitter_x > 0 || jitter_z > 0) {
        ref_data <- ref_data %>%
            mutate(
                !!mappings$x_var := .data[[mappings$x_var]] + runif(n(), -jitter_x, jitter_x),
                !!mappings$z_var := .data[[mappings$z_var]] + runif(n(), -jitter_z, jitter_z)
            )
        target_data_display <- target_data_display %>%
            mutate(
                !!mappings$x_var := .data[[mappings$x_var]] + runif(n(), -jitter_x, jitter_x),
                !!mappings$z_var := .data[[mappings$z_var]] + runif(n(), -jitter_z, jitter_z)
            )
    }

    # Create summary data (median per sample) - z_med is ORIGINAL coordinates for tie point matching
    ref_summary <- ref_data %>%
        group_by(.data[[mappings$sample_name]]) %>%
        summarise(
            x_med = median(.data[[mappings$x_var]], na.rm = TRUE),
            z_med = median(.data[[mappings$z_var]], na.rm = TRUE),  # Original coordinates
            .groups = "drop"
        ) %>%
        arrange(z_med)

    target_summary <- target_data_display %>%
        group_by(.data[[mappings$sample_name]]) %>%
        summarise(
            x_med = median(.data[[mappings$x_var]], na.rm = TRUE),  # Includes offset
            z_med = median(.data[[mappings$z_var]], na.rm = TRUE),  # Original coordinates
            .groups = "drop"
        ) %>%
        arrange(z_med)

    # Use z_med directly for display (no scaling)
    ref_summary_display <- ref_summary %>%
        mutate(z_display = z_med)
    target_summary_display <- target_summary %>%
        mutate(z_display = z_med)

    # Create plot with source for click events
    p <- plot_ly(source = plot_source)

    # Add point cloud layers first (behind everything) if enabled
    if (show_points) {
        # Reference point cloud - NO scaling
        ref_data_display <- ref_data %>%
            mutate(.z_display = .data[[mappings$z_var]])
        # Target point cloud
        target_data_display_pts <- target_data_display %>%
            mutate(.z_display = .data[[mappings$z_var]])
        
        # Reference point cloud
        p <- p %>% add_trace(
            data = ref_data_display,
            x = ~ .data[[mappings$x_var]], 
            y = ~ .z_display,
            type = "scatter", mode = "markers",
            name = paste("Ref points:", ref_core),
            marker = list(size = max(1, point_size * 0.625), color = "#1f77b4", opacity = 0.25, symbol = "circle"),
            text = ~ .data[[mappings$sample_name]],
            key = ~ .data[[mappings$sample_name]],
            customdata = "ref",
            hovertemplate = paste0(
                "<b>%{text}</b><br>",
                mappings$x_var, ": %{x:.2f}<br>",
                mappings$z_var, ": %{y:.2f}<extra>Reference</extra>"
            ),
            showlegend = FALSE
        )
        
        # Target point cloud
        p <- p %>% add_trace(
            data = target_data_display_pts,
            x = ~ .data[[mappings$x_var]], 
            y = ~ .z_display,
            type = "scatter", mode = "markers",
            name = paste("Target points:", target_core),
            marker = list(size = max(1, point_size * 0.625), color = "#ff7f0e", opacity = 0.25, symbol = "triangle-up"),
            text = ~ .data[[mappings$sample_name]],
            key = ~ .data[[mappings$sample_name]],
            customdata = "target",
            hovertemplate = paste0(
                "<b>%{text}</b><br>",
                mappings$x_var, ": %{x:.2f}<br>",
                mappings$z_var, ": %{y:.2f}<extra>Target</extra>"
            ),
            showlegend = FALSE
        )
    }

    # Add reference profile (NO scaling)
    p <- p %>% add_trace(
        data = ref_summary_display,
        x = ~x_med, y = ~z_display,
        type = "scatter", mode = "lines+markers",
        name = paste("Reference:", ref_core),
        line = list(color = "#1f77b4", width = line_width),
        marker = list(size = point_size, color = "#1f77b4", symbol = "circle"),
        text = ~ .data[[mappings$sample_name]],
        key = ~ .data[[mappings$sample_name]],
        customdata = "ref",
        hovertemplate = paste0(
            "<b>%{text}</b><br>",
            mappings$x_var, ": %{x:.2f}<br>",
            mappings$z_var, ": %{y:.2f}<extra>Reference</extra>"
        )
    )

    # Add target profile with any display-only X offset applied
    p <- p %>% add_trace(
        data = target_summary_display,
        x = ~x_med, y = ~z_display,
        type = "scatter", mode = "lines+markers",
        name = paste("Target:", target_core),
        line = list(color = "#ff7f0e", width = line_width),
        marker = list(size = point_size, color = "#ff7f0e", symbol = "triangle-up"),
        text = ~ .data[[mappings$sample_name]],
        key = ~ .data[[mappings$sample_name]],
        customdata = "target",
        hovertemplate = paste0(
            "<b>%{text}</b><br>",
            mappings$x_var, ": %{x:.2f}<br>",
            mappings$z_var, ": %{y:.2f}<extra>Target</extra>"
        )
    )

    # Add tie point lines and labels
    tie_point_annotations <- list()
    if (!is.null(tiepoints) && nrow(tiepoints) > 0) {
        valid_ties <- tiepoints %>% filter(!is.na(ref_z) & !is.na(target_z))

        if (nrow(valid_ties) > 0) {
            # Join with summary data to get x coordinates (using original z_med for matching)
            join_by_ref <- setNames(mappings$sample_name, "ref_sample")
            join_by_target <- setNames(mappings$sample_name, "target_sample")

            ties_with_coords <- valid_ties %>%
                left_join(
                    ref_summary %>%
                        rename(ref_x = x_med),
                    by = join_by_ref
                ) %>%
                left_join(
                    target_summary %>%
                        rename(target_x = x_med),
                    by = join_by_target
                ) %>%
                filter(!is.na(ref_x) & !is.na(target_x))

            if (nrow(ties_with_coords) > 0) {
                for (i in seq_len(nrow(ties_with_coords))) {
                    tie <- ties_with_coords[i, ]

                    # Determine if this tie point is selected
                    is_selected <- !is.null(selected_ids) && tie$id %in% selected_ids

                    # Style based on selection and use_in_warp
                    if (is_selected) {
                        line_color <- "#00ff00"
                        line_width <- 3
                        line_alpha <- 1.0
                    } else if (tie$use_in_warp) {
                        line_color <- "#ff0000"
                        line_width <- 1.5
                        line_alpha <- 0.6
                    } else {
                        line_color <- "#999999"
                        line_width <- 1
                        line_alpha <- 0.3
                    }

                    p <- p %>% add_segments(
                        x = tie$ref_x, xend = tie$target_x,
                        y = tie$ref_z, yend = tie$target_z,
                        line = list(color = line_color, width = line_width, dash = "dash"),
                        opacity = line_alpha,
                        showlegend = FALSE,
                        hoverinfo = "text",
                        text = paste0(
                            "Tie Point #", tie$id, "<br>",
                            "Ref: ", tie$ref_sample, " (", round(tie$ref_z, 2), ")<br>",
                            "Target: ", tie$target_sample, " (", round(tie$target_z, 2), ")"
                        )
                    )
                    
                    # Add label annotation at midpoint of tie line
                    mid_x <- (tie$ref_x + tie$target_x) / 2
                    mid_y <- (tie$ref_z + tie$target_z) / 2
                    tie_point_annotations <- c(tie_point_annotations, list(
                        list(
                            x = mid_x, y = mid_y,
                            text = paste0("TP", tie$id),
                            showarrow = FALSE,
                            font = list(size = 10, color = if (tie$use_in_warp) "#cc0000" else "#666666"),
                            bgcolor = "rgba(255,255,255,0.7)",
                            borderpad = 2
                        )
                    ))
                }
            }
        }
    }

    # Build annotations list for highlight
    annotations_list <- tie_point_annotations
    
    # Add enlarged point marker for pending reference sample (reference is NOT scaled)
    if (!is.null(highlight_sample)) {
        highlight_point <- ref_summary_display %>% 
            filter(.data[[mappings$sample_name]] == highlight_sample)
        if (nrow(highlight_point) > 0) {
            # Add a larger green point marker for the selected sample
            p <- p %>% add_trace(
                x = highlight_point$x_med[1],
                y = highlight_point$z_display[1],
                type = "scatter", mode = "markers",
                name = "Selected Reference",
                marker = list(size = max(point_size + 4, point_size * 2.25), color = "#00cc00", symbol = "circle",
                              line = list(color = "#006600", width = 2)),
                text = paste0("Selected: ", highlight_sample),
                hovertemplate = paste0("<b>", highlight_sample, "</b><br>Click target sample to create tie point<extra></extra>"),
                showlegend = FALSE
            )
            annotations_list <- c(annotations_list, list(
                list(
                    x = highlight_point$x_med[1],
                    y = highlight_point$z_display[1],
                    text = paste0("Selected: ", highlight_sample),
                    showarrow = TRUE,
                    arrowhead = 2,
                    arrowcolor = "#00cc00",
                    ax = 40, ay = -30,
                    font = list(color = "#00cc00", size = 11)
                )
            ))
        }
    }

    # Layout - optimized for depth/age data
    p <- p %>% layout(
        title = list(text = "Initial Core Alignment", y = 0.98),
        xaxis = list(title = mappings$x_var),
        yaxis = list(
            title = mappings$z_var, 
            autorange = sc_plot_y_autorange(reference_direction),
            automargin = TRUE
        ),
        hovermode = "closest",
        legend = list(orientation = "h", y = -0.1, x = 0.5, xanchor = "center"),
        margin = list(t = 40, b = 60, l = 60, r = 20),
        autosize = TRUE,
        annotations = annotations_list
    )
    
    # Register click event
    p <- p %>% event_register("plotly_click")
    
    p
}

#' Generate a display-only affine depth preview
#' @param data Original stratigraphic data
#' @param mappings Column mappings
#' @param ref_core Reference core ID
#' @param target_core Target core ID
#' @param affine_parameters Output from sc_affine_depth_parameters
#' @inheritParams sc_plot_initial
#' @return plotly object whose sample keys resolve against original data
sc_plot_affine_preview <- function(data, mappings, ref_core, target_core,
                                   affine_parameters, tiepoints = NULL,
                                   selected_ids = NULL, show_points = FALSE,
                                   plot_source = "sc_affine",
                                   highlight_sample = NULL,
                                   target_offset = 0, line_width = 2,
                                   point_size = 8, jitter_x = 0,
                                   jitter_z = 0,
                                   reference_direction = "down") {
    preview_data <- data
    preview_z_name <- paste("Preview", mappings$z_var)
    preview_data[[preview_z_name]] <- preview_data[[mappings$z_var]]

    target_rows <- as.character(preview_data[[mappings$site_id]]) == as.character(target_core)
    preview_data[[preview_z_name]][target_rows] <-
        preview_data[[mappings$z_var]][target_rows] * affine_parameters$scale +
        affine_parameters$shift

    preview_mappings <- mappings
    preview_mappings$z_var <- preview_z_name

    preview_ties <- tiepoints
    if (!is.null(preview_ties) && nrow(preview_ties) > 0) {
        preview_ties$target_z <-
            preview_ties$target_z * affine_parameters$scale + affine_parameters$shift
    }

    p <- sc_plot_initial(
        preview_data, preview_mappings, ref_core, target_core,
        tiepoints = preview_ties, selected_ids = selected_ids,
        show_points = show_points, plot_source = plot_source,
        highlight_sample = highlight_sample, target_offset = target_offset,
        line_width = line_width, point_size = point_size,
        jitter_x = jitter_x, jitter_z = jitter_z,
        reference_direction = reference_direction
    )

    plotly::layout(
        p,
        title = list(text = "Affine Coordinate Preview (display only)", y = 0.98)
    )
}

#' Generate warped alignment plot
#' @param data data.frame with all core data
#' @param warped_data data.frame with warped target data
#' @param mappings Column mappings
#' @param ref_core Reference core ID
#' @param target_core Target core ID
#' @param tiepoints data.frame of tie points
#' @param selected_ids Vector of selected tie point IDs
#' @param show_points Logical, whether to show raw point cloud behind median lines
#' @param plot_source Character, source ID for plotly click events
#' @param highlight_sample Character, sample name to highlight with a larger point marker (for click-to-add mode)
#' @param target_offset Numeric, horizontal offset to apply to target data for visual separation
#' @param line_width Numeric width of the reference and target profile lines
#' @param point_size Numeric size of the reference and target profile markers
#' @param reference_direction Whether reference values increase down or up
#' @return plotly object
sc_plot_warped <- function(data, warped_data, mappings, ref_core, target_core,
                           tiepoints = NULL, selected_ids = NULL,
                           show_points = FALSE, plot_source = "sc_warped",
                           highlight_sample = NULL,
                           target_offset = 0,
                           line_width = 2,
                           point_size = 8,
                           jitter_x = 0,
                           jitter_z = 0,
                           reference_direction = "down") {
    library(plotly)
    library(dplyr)

    data <- sc_normalize_identifier_columns(data, mappings)
    warped_data <- sc_normalize_identifier_columns(warped_data, mappings)

    if (target_offset != 0) {
        warped_data <- warped_data %>%
            mutate(!!mappings$x_var := .data[[mappings$x_var]] + target_offset)
    }

    # Filter reference data
    ref_data <- data %>%
        filter(.data[[mappings$site_id]] == ref_core) %>%
        arrange(.data[[mappings$z_var]])
    
    # Apply jitter if specified (for display only)
    set.seed(42)  # Reproducible jitter
    if (jitter_x > 0 || jitter_z > 0) {
        ref_data <- ref_data %>%
            mutate(
                !!mappings$x_var := .data[[mappings$x_var]] + runif(n(), -jitter_x, jitter_x),
                !!mappings$z_var := .data[[mappings$z_var]] + runif(n(), -jitter_z, jitter_z)
            )
        warped_data <- warped_data %>%
            mutate(
                !!mappings$x_var := .data[[mappings$x_var]] + runif(n(), -jitter_x, jitter_x),
                z_warped = z_warped + runif(n(), -jitter_z, jitter_z)
            )
    }

    # Create summary data - sorted by depth for proper line connections
    ref_summary <- ref_data %>%
        group_by(.data[[mappings$sample_name]]) %>%
        summarise(
            x_med = median(.data[[mappings$x_var]], na.rm = TRUE),
            z_med = median(.data[[mappings$z_var]], na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(z_med)

    target_summary <- warped_data %>%
        group_by(.data[[mappings$sample_name]]) %>%
        summarise(
            x_med = median(.data[[mappings$x_var]], na.rm = TRUE),
            z_med = median(z_warped, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(z_med)

    # Create plot with source for click events
    p <- plot_ly(source = plot_source)

    # Add point cloud layers first (behind everything) if enabled
    if (show_points) {
        # Reference point cloud
        p <- p %>% add_trace(
            data = ref_data,
            x = ~ .data[[mappings$x_var]], 
            y = ~ .data[[mappings$z_var]],
            type = "scatter", mode = "markers",
            name = paste("Ref points:", ref_core),
            marker = list(size = max(1, point_size * 0.625), color = "#1f77b4", opacity = 0.25, symbol = "circle"),
            text = ~ .data[[mappings$sample_name]],
            key = ~ .data[[mappings$sample_name]],
            customdata = "ref",
            hovertemplate = paste0(
                "<b>%{text}</b><br>",
                mappings$x_var, ": %{x:.2f}<br>",
                mappings$z_var, ": %{y:.2f}<extra>Reference</extra>"
            ),
            showlegend = FALSE
        )
        
        # Warped target point cloud
        p <- p %>% add_trace(
            data = warped_data,
            x = ~ .data[[mappings$x_var]], 
            y = ~z_warped,
            type = "scatter", mode = "markers",
            name = paste("Target points:", target_core),
            marker = list(size = max(1, point_size * 0.625), color = "#ff7f0e", opacity = 0.25, symbol = "triangle-up"),
            text = ~ .data[[mappings$sample_name]],
            key = ~ .data[[mappings$sample_name]],
            customdata = "target",
            hovertemplate = paste0(
                "<b>%{text}</b><br>",
                mappings$x_var, ": %{x:.2f}<br>",
                "Warped ", mappings$z_var, ": %{y:.2f}<extra>Target</extra>"
            ),
            showlegend = FALSE
        )
    }

    # Add reference profile
    p <- p %>% add_trace(
        data = ref_summary,
        x = ~x_med, y = ~z_med,
        type = "scatter", mode = "lines+markers",
        name = paste("Reference:", ref_core),
        line = list(color = "#1f77b4", width = line_width),
        marker = list(size = point_size, color = "#1f77b4", symbol = "circle"),
        text = ~ .data[[mappings$sample_name]],
        key = ~ .data[[mappings$sample_name]],
        customdata = "ref",
        hovertemplate = paste0(
            "<b>%{text}</b><br>",
            mappings$x_var, ": %{x:.2f}<br>",
            "Warped ", mappings$z_var, ": %{y:.2f}<extra>Reference</extra>"
        )
    )

    # Add warped target profile
    p <- p %>% add_trace(
        data = target_summary,
        x = ~x_med, y = ~z_med,
        type = "scatter", mode = "lines+markers",
        name = paste("Warped Target:", target_core),
        line = list(color = "#ff7f0e", width = line_width),
        marker = list(size = point_size, color = "#ff7f0e", symbol = "triangle-up"),
        text = ~ .data[[mappings$sample_name]],
        key = ~ .data[[mappings$sample_name]],
        customdata = "target",
        hovertemplate = paste0(
            "<b>%{text}</b><br>",
            mappings$x_var, ": %{x:.2f}<br>",
            "Warped ", mappings$z_var, ": %{y:.2f}<extra>Target</extra>"
        )
    )

    # Add tie point lines and labels (similar to initial plot)
    tie_point_annotations <- list()
    if (!is.null(tiepoints) && nrow(tiepoints) > 0) {
        valid_ties <- tiepoints %>% filter(!is.na(ref_z) & !is.na(target_z))

        if (nrow(valid_ties) > 0) {
            # Create join vector properly
            join_by_ref <- setNames(mappings$sample_name, "ref_sample")
            join_by_target <- setNames(mappings$sample_name, "target_sample")

            ties_with_coords <- valid_ties %>%
                left_join(
                    ref_summary %>%
                        rename(ref_x = x_med, ref_z_check = z_med),
                    by = join_by_ref
                ) %>%
                left_join(
                    target_summary %>%
                        rename(target_x = x_med, target_z_warped = z_med),
                    by = join_by_target
                ) %>%
                # Filter out ties with missing coordinates (failed joins)
                filter(!is.na(ref_x) & !is.na(target_x) & !is.na(target_z_warped))

            if (nrow(ties_with_coords) > 0) {
                for (i in seq_len(nrow(ties_with_coords))) {
                    tie <- ties_with_coords[i, ]
                    is_selected <- !is.null(selected_ids) && tie$id %in% selected_ids

                    if (is_selected) {
                        line_color <- "#00ff00"
                        line_width <- 3
                        line_alpha <- 1.0
                    } else if (tie$use_in_warp) {
                        line_color <- "#ff0000"
                        line_width <- 1.5
                        line_alpha <- 0.6
                    } else {
                        line_color <- "#999999"
                        line_width <- 1
                        line_alpha <- 0.3
                    }

                    p <- p %>% add_segments(
                        x = tie$ref_x, xend = tie$target_x,
                        y = tie$ref_z, yend = tie$target_z_warped,
                        line = list(color = line_color, width = line_width, dash = "dash"),
                        opacity = line_alpha,
                        showlegend = FALSE,
                        hoverinfo = "text",
                        text = paste0(
                            "Tie Point #", tie$id, "<br>",
                            "Ref: ", tie$ref_sample, " (", round(tie$ref_z, 2), ")<br>",
                            "Target: ", tie$target_sample, " (warped: ", round(tie$target_z_warped, 2), ")"
                        )
                    )
                    
                    # Add label annotation at midpoint of tie line
                    mid_x <- (tie$ref_x + tie$target_x) / 2
                    mid_y <- (tie$ref_z + tie$target_z_warped) / 2
                    tie_point_annotations <- c(tie_point_annotations, list(
                        list(
                            x = mid_x, y = mid_y,
                            text = paste0("TP", tie$id),
                            showarrow = FALSE,
                            font = list(size = 10, color = if (tie$use_in_warp) "#cc0000" else "#666666"),
                            bgcolor = "rgba(255,255,255,0.7)",
                            borderpad = 2
                        )
                    ))
                }
            }
        }
    }

    # Build annotations list for highlight
    annotations_list <- tie_point_annotations
    
    # Add enlarged point marker for pending reference sample
    if (!is.null(highlight_sample)) {
        highlight_point <- ref_summary %>% 
            filter(.data[[mappings$sample_name]] == highlight_sample)
        if (nrow(highlight_point) > 0) {
            # Add a larger green point marker for the selected sample
            p <- p %>% add_trace(
                x = highlight_point$x_med[1],
                y = highlight_point$z_med[1],
                type = "scatter", mode = "markers",
                name = "Selected Reference",
                marker = list(size = max(point_size + 4, point_size * 2.25), color = "#00cc00", symbol = "circle",
                              line = list(color = "#006600", width = 2)),
                text = paste0("Selected: ", highlight_sample),
                hovertemplate = paste0("<b>", highlight_sample, "</b><br>Click target sample to create tie point<extra></extra>"),
                showlegend = FALSE
            )
            annotations_list <- c(annotations_list, list(
                list(
                    x = highlight_point$x_med[1],
                    y = highlight_point$z_med[1],
                    text = paste0("Selected: ", highlight_sample),
                    showarrow = TRUE,
                    arrowhead = 2,
                    arrowcolor = "#00cc00",
                    ax = 40, ay = -30,
                    font = list(color = "#00cc00", size = 11)
                )
            ))
        }
    }

    # Layout - optimized for depth/age data
    p <- p %>% layout(
        title = list(text = "Warped Core Alignment", y = 0.98),
        xaxis = list(title = mappings$x_var),
        yaxis = list(
            title = paste("Warped", mappings$z_var), 
            autorange = sc_plot_y_autorange(reference_direction),
            automargin = TRUE
        ),
        hovermode = "closest",
        legend = list(orientation = "h", y = -0.1, x = 0.5, xanchor = "center"),
        margin = list(t = 40, b = 60, l = 60, r = 20),
        autosize = TRUE,
        annotations = annotations_list
    )
    
    # Register click event
    p <- p %>% event_register("plotly_click")
    
    p
}

#' Generate warp fit plot
#' @param warp_result Result from sc_calculate_warp
#' @return plotly object
sc_plot_fit <- function(warp_result) {
    library(plotly)

    residuals <- warp_result$diagnostics$residuals

    # Generate smooth prediction line
    z_range <- range(residuals$target_z)
    z_seq <- seq(z_range[1], z_range[2], length.out = 200)
    pred_df <- data.frame(
        target_z = z_seq,
        pred_z = warp_result$warp_func(z_seq)
    )

    # Create plot
    p <- plot_ly()

    # Add prediction line
    p <- p %>% add_trace(
        data = pred_df,
        x = ~target_z, y = ~pred_z,
        type = "scatter", mode = "lines",
        name = "Warp Function",
        line = list(color = "#1f77b4", width = 2)
    )

    # Add tie points
    p <- p %>% add_trace(
        data = residuals,
        x = ~target_z, y = ~ref_z,
        type = "scatter", mode = "markers",
        name = "Tie Points",
        marker = list(size = 10, color = "#ff7f0e"),
        text = ~ paste0(
            "Tie Point #", id, "<br>",
            "Target: ", target_sample, "<br>",
            "Target Z: ", round(target_z, 2), "<br>",
            "Reference Z: ", round(ref_z, 2), "<br>",
            "Predicted Z: ", round(pred_z, 2), "<br>",
            "Residual: ", round(residual, 4)
        ),
        hovertemplate = "%{text}<extra></extra>"
    )

    # Layout
    p %>% layout(
        title = paste("Warp Fit -", warp_result$model_type),
        xaxis = list(title = "Target Depth/Age"),
        yaxis = list(title = "Reference Depth/Age"),
        hovermode = "closest",
        legend = list(orientation = "h", y = -0.15)
    )
}
