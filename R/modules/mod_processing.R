# mod_processing.R

processing_heavy_worker_cap <- function() {
  limiter_cap <- tryCatch(
    as.integer(heavy_job_limiter_config()$cpu_tokens),
    error = function(e) NA_integer_
  )
  configured_cap <- suppressWarnings(as.integer(Sys.getenv("VARG_MAX_CORES", unset = "0")))
  if (length(configured_cap) != 1L || is.na(configured_cap) || configured_cap < 1L) {
    available <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (length(available) != 1L || is.na(available) || available < 2L) available <- 2L
    configured_cap <- max(1L, available - 1L)
  }
  if (length(limiter_cap) != 1L || is.na(limiter_cap) || limiter_cap < 1L) {
    limiter_cap <- configured_cap
  }
  max(1L, min(limiter_cap, configured_cap))
}

processing_blas_worker_env <- function(workers) {
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) workers <- 1L
  value <- as.character(workers)
  c(
    OPENBLAS_NUM_THREADS = value,
    OMP_NUM_THREADS = value,
    OMP_THREAD_LIMIT = value,
    MKL_NUM_THREADS = value,
    BLIS_NUM_THREADS = value,
    VECLIB_MAXIMUM_THREADS = value
  )
}

processing_report_imputation_method_label <- function(method) {
  method <- if (length(method) > 0L && !is.null(method)) as.character(method[[1L]]) else "unknown"
  switch(
    method,
    ltsReg = "Robust regression (ltsReg)",
    lm = "Standard regression (lm)",
    auto = "Auto",
    unknown = "Unknown",
    method
  )
}

processing_report_imputation_auto_html <- function(missingness, row_count, resolved_method) {
  if (nrow(missingness) == 0L) return("")

  maximum_row <- missingness[which.max(missingness[["Affected (%)"]]), , drop = FALSE]
  maximum_fraction <- maximum_row[["Affected (%)"]]
  threshold_result <- if (maximum_fraction > 50) {
    "was above the 50% threshold, so Auto selected standard regression"
  } else {
    "was not above the 50% threshold, so Auto retained robust regression"
  }

  paste0(
    "<p class='meta'>Auto evaluated missingness across the full ", row_count,
    "-row dataset. The maximum affected fraction was ", maximum_fraction, "% in ",
    htmltools::htmlEscape(maximum_row$Analyte), "; this ", threshold_result,
    ". Method used: ", processing_report_imputation_method_label(resolved_method), ".</p>"
  )
}

processing_report_dimensions_label <- function(dimensions) {
  if (is.null(dimensions) || length(dimensions) == 0L) return("Unknown")

  if (is.character(dimensions)) {
    normalized <- tolower(dimensions)
    if (any(normalized == "both")) return("1D and 2D")
    dimensions <- sub("d$", "", normalized)
  }

  dimensions <- suppressWarnings(as.integer(dimensions))
  dimensions <- unique(dimensions[dimensions %in% c(1L, 2L)])
  if (length(dimensions) == 0L) return("Unknown")
  paste(paste0(dimensions, "D"), collapse = " and ")
}

processing_report_umap_html <- function(umap_config) {
  if (is.null(umap_config)) {
    return("<h2>4. UMAP Dimensionality Reduction</h2><p class='muted'>UMAP was not run in this session.</p>")
  }

  esc <- function(value) htmltools::htmlEscape(as.character(value))
  mode <- if (length(umap_config$mode) > 0L) as.character(umap_config$mode[[1L]]) else "unknown"
  html <- paste0(
    "<h2>4. UMAP Dimensionality Reduction</h2><table><tbody>",
    "<tr><td><strong>Mode</strong></td><td>", esc(mode), "</td></tr>"
  )

  if (identical(mode, "pretrained")) {
    html <- paste0(
      html,
      "<tr><td><strong>VARG26 dimensions</strong></td><td>",
      esc(processing_report_dimensions_label(umap_config$VARG26_dims)), "</td></tr>",
      "<tr><td><strong>VARG26 oxides</strong></td><td>",
      esc(paste(umap_config$VARG26_oxides, collapse = ", ")), "</td></tr>"
    )
  } else {
    dimensions <- umap_config$n_components
    if (is.null(dimensions)) dimensions <- umap_config$dimensions
    html <- paste0(
      html,
      "<tr><td><strong>Components</strong></td><td>",
      esc(processing_report_dimensions_label(dimensions)), "</td></tr>",
      "<tr><td><strong>n_neighbors</strong></td><td>", esc(umap_config$n_neighbors), "</td></tr>",
      "<tr><td><strong>min_dist</strong></td><td>", esc(umap_config$min_dist), "</td></tr>",
      "<tr><td><strong>dens_scale</strong></td><td>", esc(umap_config$dens_scale), "</td></tr>",
      "<tr><td><strong>Semi-supervised</strong></td><td>",
      if (isTRUE(umap_config$semisupervised)) "Yes" else "No", "</td></tr>",
      "<tr><td><strong>Columns used</strong></td><td>",
      esc(paste(umap_config$columns_used, collapse = ", ")), "</td></tr>",
      "</tbody></table>"
    )
  }

  if (identical(mode, "pretrained")) html <- paste0(html, "</tbody></table>")
  html
}

processing_report_imputation_summary <- function(data, compositional_columns) {
  compositional_columns <- intersect(compositional_columns, names(data))
  if (length(compositional_columns) == 0L || nrow(data) == 0L) return(data.frame())

  rows <- lapply(compositional_columns, function(column) {
    values <- suppressWarnings(as.numeric(as.character(data[[column]])))
    missing <- is.na(values) | !is.finite(values)
    non_positive <- !missing & values <= 1e-10
    affected <- missing | non_positive
    data.frame(
      Analyte = column,
      Missing = sum(missing),
      `Zero/non-positive` = sum(non_positive),
      Total = length(values),
      Affected = sum(affected),
      `Affected (%)` = round(mean(affected) * 100, 1),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, rows)
  summary[summary$Affected > 0L, , drop = FALSE]
}

processing_report_structured_missingness <- function(data, compositional_columns) {
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(data)))
  candidate_keys <- c("coreid", "siteid", "recordid", "core", "site")
  group_index <- match(candidate_keys, normalized_names, nomatch = 0L)
  group_index <- group_index[group_index > 0L]
  if (length(group_index) == 0L || nrow(data) == 0L) return(data.frame())

  group_column <- names(data)[group_index[[1L]]]
  groups <- as.character(data[[group_column]])
  valid_groups <- !is.na(groups) & nzchar(trimws(groups))
  compositional_columns <- intersect(compositional_columns, names(data))
  findings <- list()

  for (column in compositional_columns) {
    values <- suppressWarnings(as.numeric(as.character(data[[column]])))
    affected <- is.na(values) | !is.finite(values) | values <= 1e-10
    if (!any(!affected)) next

    for (group in unique(groups[valid_groups])) {
      in_group <- valid_groups & groups == group
      if (sum(in_group) >= 2L && all(affected[in_group])) {
        findings[[length(findings) + 1L]] <- data.frame(
          `Grouping column` = group_column,
          Group = group,
          Analyte = column,
          Rows = sum(in_group),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(findings) == 0L) return(data.frame())
  do.call(rbind, findings)
}

processing_gmm_selected_compositional_columns <- function(pipeline_config, ui_columns = NULL) {
  configured <- pipeline_config$preprocessing$comp_cols
  if (!is.null(configured) && length(configured) > 0L) {
    return(as.character(configured))
  }
  if (is.null(ui_columns)) character(0) else as.character(ui_columns)
}

processing_gmm_cluster_composition <- function(data, compositional_columns) {
  if (is.null(data) || !"gmm_cluster" %in% names(data)) return(NULL)

  valid_clusters <- !is.na(data$gmm_cluster)
  if (!any(valid_clusters)) return(NULL)

  compositional_columns <- unique(as.character(compositional_columns))
  compositional_columns <- compositional_columns[nzchar(compositional_columns)]
  if (length(compositional_columns) == 0L) return(NULL)

  df <- data[valid_clusters, , drop = FALSE]
  cluster_labels <- as.character(df$gmm_cluster)
  cluster_levels <- unique(cluster_labels)
  numeric_levels <- suppressWarnings(as.numeric(cluster_levels))
  if (all(!is.na(numeric_levels))) {
    cluster_levels <- as.character(sort(unique(numeric_levels)))
  } else {
    cluster_levels <- sort(unique(cluster_levels))
  }

  source_columns <- vapply(
    compositional_columns,
    function(column) {
      imputed <- paste0(column, "_imp")
      if (imputed %in% names(df)) imputed else if (column %in% names(df)) column else NA_character_
    },
    character(1)
  )
  available <- !is.na(source_columns)
  source_columns <- source_columns[available]
  display_names <- compositional_columns[available]
  if (length(source_columns) == 0L) return(NULL)

  uses_imputed <- grepl("_imp$", source_columns)
  data_source <- if (all(uses_imputed)) {
    "imputed"
  } else if (any(uses_imputed)) {
    "mixed"
  } else {
    "raw"
  }

  format_mean_sd <- function(values) {
    values <- suppressWarnings(as.numeric(as.character(values)))
    if (all(is.na(values))) return("NA \u00B1 NA")
    mean_value <- mean(values, na.rm = TRUE)
    sd_value <- stats::sd(values, na.rm = TRUE)
    if (!is.finite(mean_value)) mean_value <- NA_real_
    if (!is.finite(sd_value)) sd_value <- NA_real_
    paste0(
      formatC(mean_value, format = "f", digits = 2),
      " \u00B1 ",
      formatC(sd_value, format = "f", digits = 2)
    )
  }

  summary_df <- data.frame(
    Cluster = cluster_levels,
    n = vapply(cluster_levels, function(cluster) sum(cluster_labels == cluster), integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (index in seq_along(source_columns)) {
    summary_df[[display_names[index]]] <- vapply(
      cluster_levels,
      function(cluster) format_mean_sd(df[[source_columns[index]]][cluster_labels == cluster]),
      character(1)
    )
  }

  attr(summary_df, "data_source") <- data_source
  summary_df
}

processing_gmm_composition_source_label <- function(data_source) {
  switch(
    data_source,
    imputed = "imputed values (after compositional imputation, before log-ratio transformation)",
    mixed = "available imputed values, with raw values used where no imputed column was available",
    "raw uploaded values (imputation was not applied)"
  )
}

processing_report_gmm_composition_html <- function(composition_data, df_to_html) {
  if (is.null(composition_data) || nrow(composition_data) == 0L) return("")
  source_label <- processing_gmm_composition_source_label(attr(composition_data, "data_source") %||% "raw")
  paste0(
    "<h3>Cluster Composition Summary (Mean &plusmn; SD)</h3>",
    "<p class='meta'>Computed from ", htmltools::htmlEscape(source_label),
    " for the compositional columns used in preprocessing.</p>",
    df_to_html(composition_data)
  )
}

mod_processing_ui <- function(id) {
  ns <- NS(id)

  tagList(
    # CSS adjustments (Required for Optimized Tab 3)
    # Moved to www/custom.css
    bslib::navset_card_underline(
      id = ns("proc_wizard"),
      title = "Processing Workflow",

      # --- Step 1: Data Import ---
      bslib::nav_panel(
        title = "1. Import Data",
        icon = icon("file-import"),
        module_banner(
          goal = "Load your geochemical dataset and filter out unwanted rows.",
          inputs = "CSV or Excel file with oxide analyses and metadata.",
          outputs = "Clean data table ready for preprocessing.",
          why = "Every analysis starts with clean, validated data. This step catches blank rows, lets you exclude unwanted samples, and sets up your dataset for the transformation pipeline."
        ),
        uiOutput(ns("step1_status_ui")),
        layout_columns(
          col_widths = c(4, 8),
          varg_card(
            title = tagList(icon("upload"), " Upload"),
            # Instructional Help Box
            help_box(
              title = "Data Requirements",
              content = "Upload your geochemical dataset (CSV or Excel). Your data should include:
              <ul>
                <li><b>Compositional Columns:</b> Major oxides (e.g., SiO2, Al2O3, CaO, FeO, etc.) that sum to ~100%. <b>Required for VARG26 projection:</b> SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O. <b>Recommended but optional:</b> P2O5, Cl, F.</li>
                <li><b>Non-Compositional Columns:</b> Metadata like Sample ID, Depth, Age, coordinates. These don't need statistical transformation but can be used in Advanced Mode for custom analyses.</li>
                <li><b>Sample Identifiers:</b> A unique ID for each row is recommended to track results. If not present, one will be generated.</li>
              </ul>
              <small><i>Tip: Use the template below to ensure correct formatting. Check which analytes your region typically reports. For example, Icelandic basalts may benefit from including P2O5 even though it's optional for Alaskan glass.</i></small>"
            ),
            fileInput(ns("file_upload"), "Upload CSV or Excel", accept = c(".csv", ".xlsx", ".xls")),
            uiOutput(ns("sheet_ui")),
            actionButton(ns("load_data"), "Load Data", class = "btn-primary w-100", icon = icon("arrow-right")),
            
            # Row filtering options (shown after data is loaded)
            uiOutput(ns("row_filter_ui")),
            
            hr(),
            div(
              class = "text-muted small",
              icon("info-circle"), " Need a template?",
              downloadButton(ns("download_template"), "Download Template", class = "btn-outline-secondary btn-sm mt-2 w-100")
            )
          ),
          div(
            varg_card(
              title = tagList(icon("table"), " Column Summary"),
              uiOutput(ns("data_summary_ui"))
            ),
            varg_card(
              title = tagList(icon("eye"), " Data Preview"),
              DTOutput(ns("preview_data_step1"))
            )
          )
        )
      ),

      # --- Step 2: Preprocessing ---
      bslib::nav_panel(
        title = "2. Preprocessing",
        icon = icon("filter"),
        module_banner(
          goal = "Transform compositional data into a statistically valid form for downstream analysis.",
          inputs = "Loaded data from Step 1 with compositional columns identified.",
          outputs = "Imputed, normalized, and transformed columns (ILR pivot or CLR).",
          why = "Raw oxide data has a closure problem (values must sum to 100%), which distorts standard statistics. These transformations open data into proper Euclidean space so GMM clustering and UMAP projections give meaningful results."
        ),
        uiOutput(ns("step2_status_ui")),
        layout_columns(
          col_widths = c(4, 8),
          varg_card(
            title = tagList(icon("cogs"), " Configuration"),
            div(
              class = "scrollable-config",
            help_box(
              title = "What does preprocessing do?",
              content = "Before clustering or UMAP, compositional data needs special treatment. Here's why:
              <ul>
                <li><b>The problem:</b> Geochemical oxides sum to ~100% (they are 'closed'). This creates artificial correlations that mislead standard statistics. We need to 'open' the data before analysis.</li>
                <li><b>Compositional Columns:</b> Select your oxide analytes (e.g., SiO2, Al2O3, CaO). These can be transformed to break closure using ILR (recommended) or CLR.</li>
                <li><b>Non-Compositional Columns:</b> Metadata like Sample ID, Depth, or Age. These are carried through without transformation but can be used in Advanced Mode for clustering/UMAP if desired (e.g., mixing temporal and compositional data).</li>
                <li><b>Imputation:</b> Replaces zeros and missing values using regression-based compositional methods implemented through robCompositions::impCoda() (Hron et al. 2010). Essential when some analytes weren't used in all samples, or if you want to derive ratios (or log-ratios) from your data.</li>
                <li><b>Log-Ratio Transformation (ILR or CLR):</b> ILR (pivot coordinates) is preferred for GMM and UMAP. CLR keeps one transformed value per oxide and is useful for exploratory analysis and visualization, but can produce singular covariance for GMM.</li>
              </ul>
              <small><i>If you need more detail, see the User Guide (accessible from the header).</i></small>"
            ),
            uiOutput(ns("compcol_ui")),
            uiOutput(ns("noncomp_ui")),
            checkboxInput(ns("do_impute"), div(style = "display: inline;", "Impute Zeros/Missing", help_icon("<strong>Fills missing and zero values before compositional transforms.</strong><details><summary>Learn more</summary>Uses regression-based compositional imputation through <code>robCompositions::impCoda()</code> (Hron et al. 2010). Auto mode selects robust or standard regression from column missingness and can fall back if a fit fails.<br><br><b>Why use it?</b> Many datasets include analytes that were not measured, or zeros caused by detection-limit artefacts.<br><br>Imputation prevents these values from breaking log-ratio transformations <em>and</em> derived ratio variables (for example, SiO\u2082/Al\u2082O\u2083).<br><br><b>Default:</b> ON. Leave this checked unless you are certain compositional columns contain no zeros or missing values.<br><br><b>If unchecked:</b> Zeros can cause log-ratio transforms to fail or return -Inf values.</details>")), value = TRUE),
            conditionalPanel(
              condition = paste0("input['", ns("do_impute"), "']"),
              selectInput(ns("impute_method"),
                div(style = "display: inline;", "Imputation method",
                  help_icon("<strong>Controls how missing oxide values are estimated.</strong><details><summary>Learn more</summary><b>Auto (recommended):</b> Uses Least Trimmed Squares (LTS) high-breakdown robust regression for best compatibility with the VARG26 reference database. Automatically falls back to standard ordinary least-squares regression if your data has columns with very high missingness (>50%), which can cause robust methods to produce unreliable estimates.<br><br><b>Robust:</b> Least Trimmed Squares high-breakdown regression. It is resistant to outliers because it discards the most extreme residuals before fitting. This matches the method used to build the VARG26 reference database. It can produce unreliable estimates when a large proportion of a column is missing.<br><br><b>Standard:</b> Ordinary least-squares linear regression. This is always stable, but slightly less resistant to outliers.<br><br><b>Most users should leave this on Auto.</b></details>")),
                choices = c("Auto (recommended)" = "auto", "Robust" = "ltsReg", "Standard" = "lm"),
                selected = "auto")
            ),
            selectInput(ns("transform_type"), div(style = "display: inline;", "Log-ratio transformation", help_icon("<strong>Select ILR (recommended), CLR, or no log-ratio transform.</strong><details><summary>Learn more</summary><b>ILR (Pivot Coordinates):</b> Isometric log-ratio transform using pivot coordinates. Recommended for compositional GMM clustering and UMAP because covariance is well-behaved in Euclidean space.<br><br><b>CLR (Centered Log-Ratio):</b> Divides each component by the row geometric mean, then takes the log. CLR keeps the same number of columns as the input and each coordinate maps directly to one oxide, so it is often easier to interpret.<br><br><b>Important:</b> CLR coordinates sum to zero, so the covariance matrix is singular and GMM may fail for some datasets.<br><br><b>None:</b> Skip log-ratio transformation (use only for non-compositional analyses or teaching/demo workflows).<br><br><b>Recommendation:</b> Use ILR (Pivot) for GMM and UMAP. CLR is acceptable for exploratory analysis and visualization.</details>")), choices = c("ILR (Pivot Coordinates)" = "ilr", "CLR (Centered Log-Ratio)" = "clr", "None" = "none"), selected = "ilr"),
            uiOutput(ns("pivot_var_ui")),
            actionButton(ns("apply_transforms"), "Apply Transforms", class = "btn-success w-100", icon = icon("check")),
            shinyjs::hidden(actionButton(ns("cancel_preprocess"), "Cancel", class = "btn-danger w-100 mt-2", icon = icon("stop")))
            ) # end scrollable-config
          ),
          div(
            varg_card(
              title = tagList(icon("columns"), " Generated Columns"),
              uiOutput(ns("column_summary_ui"))
            ),
            varg_card(
              title = tagList(icon("eye"), " Data Preview"),
              DTOutput(ns("preview_data_step2"))
            )
          )
        )
      ),

      # --- Step 3: Clustering (GMM) ---
      bslib::nav_panel(
        title = "3. Clustering (GMM)",
        icon = icon("layer-group"),
        module_banner(
          goal = "Automatically identify distinct geochemical populations using Gaussian Mixture Modeling.",
          inputs = "Transformed data (preferably ILR pivot coordinates; CLR supported for exploratory use).",
          outputs = "Cluster assignments (gmm_cluster column) and BIC model selection plot.",
          why = "GMM finds natural groupings in multivariate geochemical space that may be difficult to see in simple bivariate plots. This helps identify distinct tephra populations, detect mixed samples or detrital grains, and flag outliers that don't belong to known populations."
        ),
        uiOutput(ns("step3_status_ui")),
        layout_columns(
          col_widths = c(4, 8),
          varg_card(
            title = tagList(icon("sliders-h"), " GMM Settings"),
            div(
              class = "scrollable-config",
            # Instructional Help Box
            help_box(
              title = "Gaussian Mixture Model (GMM) Clustering",
              content = "Automatically identify distinct geochemical populations in your data.
              <ul>
                <li><b>What it does:</b> Fits multiple Gaussian (bell-curve) distributions to your data to find natural groupings. Each group represents a potential population with its own geochemical signature (mean composition, variance).</li>
                <li><b>Choosing how many clusters:</b> You set a <em>range</em> (Min to Max clusters). The algorithm tests every number in that range and uses BIC (Bayesian Information Criterion) to determine which number of clusters best describes your data.</li>
                <li><b>Interpreting BIC:</b> In this implementation (mclust), <b>higher BIC values are better</b>. The plot shows how model quality changes with the number of clusters. Look for the peak; that's the optimal number of groups.</li>
                <li><b>What if BIC looks noisy?</b> If you see 'jaggy' lines that don't extend across the full range, some models failed to converge. Try enabling <b>Stabilize Estimation</b> (conjugate prior), or reduce the Max clusters. This is normal for complex or small datasets.</li>
              </ul>
              <small><i>Reference: Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, The R Journal, 8/1, 289-317.</i></small>"
            ),
            conditionalPanel(
              condition = "input.gmm_advanced_mode == false", ns = ns,
              uiOutput(ns("gmm_data_source_ui"))
            ),

            # Advanced Column Selection
            checkboxInput(ns("gmm_advanced_mode"), div("Advanced Mode (Custom Columns)", help_icon("<strong>Lets you manually choose exactly which columns enter GMM.</strong><details><summary>Learn more</summary>Overrides the default data source and lets you hand-pick clustering columns.<br><br><b>Example use cases:</b> combine compositional data with temporal information (for example, depth), cluster a targeted subset of analytes, or exclude noisy variables.<br><br><b>Default:</b> OFF. Most users should leave this off and use the data source selector above.<br><br><b>Effect when ON:</b> A multi-select list appears so you choose exactly which columns enter the GMM.<br><br><b>Note:</b> Only selected variables are used. For example, to combine depth with normalized oxides, you must select all needed columns manually.</details>")), value = FALSE),
            conditionalPanel(
              condition = "input.gmm_advanced_mode == true", ns = ns,
              selectizeInput(ns("gmm_custom_cols"), "Select Columns:", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button", "drag_drop")))
            ),
            hr(),
            # Standard Controls
            numericInput(ns("gmin"), label = div("Min Clusters (G):", help_icon("<strong>Sets the minimum number of clusters tested by GMM.</strong><details><summary>Learn more</summary>Defines the minimum number of populations to consider.<br><br><b>Default:</b> 2 (at least two groups).<br><br>Setting this higher skips simpler models, which is useful if you already know there are at least N groups and want to save computation time.<br><br><b>When to change:</b> If you are confident your data contains many populations, setting Min to 5+ can speed analysis by skipping trivial solutions.</details>")), value = 2, min = 1, step = 1),
            numericInput(ns("gmax"), label = div("Max Clusters (G):", help_icon("<strong>Sets the maximum number of clusters tested by GMM.</strong><details><summary>Learn more</summary>Defines the maximum number of populations to test.<br><br>The algorithm tests every value from Min to Max and selects the best via BIC.<br><br><b>Default:</b> 15.<br><br><b>When to increase:</b> If the BIC curve is still climbing at the right edge, your data may contain more populations than tested; increase Max and re-run.<br><br><b>When to decrease:</b> If computation is slow or the BIC plot is very noisy, a lower Max can help.<br><br><b>Tip:</b> Start with 15–20 for most datasets.</details>")), value = 15, min = 1, step = 1),
            checkboxInput(ns("noise_init"), div("Detect Outliers", help_icon("<strong>Adds a noise cluster to catch points that do not fit well.</strong><details><summary>Learn more</summary>Enables <b>mclust noise initialization</b> by adding a noise component (Cluster 0).<br><br>This catches outlier points that do not fit neatly into any population.<br><br><b>Default:</b> ON.<br><br><b>Why?</b> Real-world geochemical data often includes scattered or mixed-source analyses that do not belong to clean clusters.<br><br><b>When to turn OFF:</b> If the model fails to converge (for example, incomplete lines in the BIC plot), or if you want every point forced into its nearest cluster instead of flagged as noise.</details>")), value = TRUE),
            checkboxInput(ns("use_prior"), div("Stabilize Estimation", help_icon("<strong>Adds prior regularization to stabilize difficult GMM fits.</strong><details><summary>Learn more</summary>Adds a <b>conjugate prior</b> regularization term to stabilize GMM fitting.<br><br><b>Default:</b> OFF.<br><br><b>When to turn ON:</b> If the BIC plot is jagged and lines do not extend across the full cluster range, or if you get convergence errors.<br><br>The prior smooths estimation and can rescue unstable models, especially for small datasets or when many model types fail.<br><br><b>Effect:</b> Results may be slightly less precise but more reliable.</details>")), value = FALSE),
            checkboxInput(ns("include_uncertainty"), div("Include Uncertainty/Probabilities", help_icon("<strong>Adds per-cluster probabilities and an overall uncertainty score.</strong><details><summary>Learn more</summary>Adds extra output columns: one probability column per cluster (for example, gmm_prob_1, gmm_prob_2, ...) plus an overall uncertainty score.<br><br><b>Default:</b> OFF.<br><br><b>When to turn ON:</b> If you want to identify ambiguous samples between populations, or quantify clustering confidence for each point.<br><br><b>Trade-off:</b> If G is large, many columns are added and exports may become cluttered.</details>")), value = FALSE),
            div(
              class = "d-flex gap-2",
              actionButton(ns("run_gmm"), "Run GMM", class = "btn-primary flex-grow-1", icon = icon("play")),
              shinyjs::hidden(
                actionButton(ns("cancel_gmm"), "Cancel", class = "btn-danger", icon = icon("stop"))
              )
            ),

            # Progress Bar for GMM
            shinyjs::hidden(
              div(
                id = ns("gmm_progress_container"), class = "mt-3",
                tags$div(
                  class = "progress",
                  tags$div(
                    class = "progress-bar progress-bar-striped progress-bar-animated bg-info",
                    role = "progressbar", style = "width: 100%", "Clustering in progress..."
                  )
                )
              )
            )
            ) # end scrollable-config
          ),
          varg_card(
            title = tagList(icon("chart-bar"), " BIC Plot"),
            shinycssloaders::withSpinner(plotOutput(ns("bic_plot"), height = "400px"), type = 4, color = "#2C3E50"),
            textOutput(ns("bic_summary"))
          )
        ),
        # Cluster Interpretation Report (appears after GMM completes)
        conditionalPanel(
          condition = sprintf("output['%s']", ns("gmm_has_result")),
          varg_card(
            title = tagList(icon("microscope"), " Cluster Interpretation Report"),
            p(class = "text-muted", style = "font-size: 13px; margin-bottom: 12px;",
              "Mean composition and standard deviation per cluster. Use this to understand what each cluster represents geochemically.",
              "Clusters with similar means across all variables may warrant merging; clusters with very few samples may be outlier groups."
            ),
            DTOutput(ns("gmm_cluster_report")),
            tags$h5("Cluster Composition Summary", class = "mt-3 mb-2"),
            DTOutput(ns("gmm_cluster_composition_summary")),
            tags$p(
              class = "text-muted small mt-2 mb-0",
              "Means and SDs computed on untransformed oxide values. Compositional data are non-Gaussian; these are approximate summaries."
            ),
            hr(),
            div(
              class = "d-flex justify-content-end",
              downloadButton(ns("download_cluster_report"), "Download Report", class = "btn-outline-primary btn-sm")
            ),
            tags$p(
              class = "text-muted small mt-2 mb-0",
              "Tip: These cluster assignments can be used as a starting point in the Population Definition tab (Step 5). Use 'Initialize from Column' to load them, then refine manually."
            )
          )
        )
      ),

      # --- Step 4: Projection (UMAP) ---
      bslib::nav_panel(
        title = "4. Projection (UMAP)",
        icon = icon("project-diagram"),
        module_banner(
          goal = "Reduce high-dimensional geochemical data to a visual 2D or 1D map.",
          inputs = "Transformed data. Optionally, labels from GMM for semi-supervised mode.",
          outputs = "UMAP coordinate columns for visualization and population definition.",
          why = "With 10+ oxide dimensions, you can't visualize all relationships at once. UMAP compresses this into a 2D map where similar samples cluster together, making it easy to spot populations, compare unknowns against the VARG26 reference database, and define groups for tephra correlation."
        ),
        uiOutput(ns("step4_status_ui")),
        layout_columns(
          col_widths = c(4, 8),

          # Left column: Controls
          div(
            class = "scrollable-config",
            varg_card(
              title = tagList(icon("map"), " VARG26 Projection (Pre-trained)"),
              # Instructional Help Box
              help_box(
                title = "What is UMAP and why use it?",
                content = "UMAP (Uniform Manifold Approximation and Projection) takes your high-dimensional geochemical data (10+ oxide columns) and creates a 2D or 1D map where similar compositions plot close together. <b>The X and Y axes of a UMAP plot do not have individual chemical meanings</b>; they represent compressed dimensions of the full compositional space.
                <ul>
                  <li><b>VARG26 Projection:</b> Projects your unknowns into a <em>fixed, reproducible</em> coordinate space built from the VARG26 regional reference database. <b>Use this when</b> you want to compare your samples against known tephra populations. The coordinates are consistent across sessions and datasets. A tiny internal noise term is added before projection to avoid exact self-matches for reference-identical rows.</li>
                  <li><b>New UMAP:</b> Creates a custom projection optimized for the internal structure of <em>your data only</em>. <b>Use this when</b> you want to explore relationships and clusters within your own dataset without the reference context.</li>
                  <li><b>You can use both!</b> Create a VARG26 projection for regional comparison AND a New UMAP to explore your dataset's structure. Both sets of coordinates will be available for visualization.</li>
                </ul>"
              ),
              selectInput(ns("VARG26_dims"), "Dimensions:", choices = c("2D" = "2d", "1D" = "1d", "Both" = "both")),
              actionButton(ns("run_umap_pre"), "Project onto VARG26", class = "btn-info w-100", icon = icon("map")),
              shinyjs::hidden(actionButton(ns("cancel_varg26"), "Cancel", class = "btn-danger w-100 mt-2", icon = icon("stop")))
            ),
            varg_card(
              title = tagList(icon("plus-circle"), " Create New UMAP"),
              selectInput(ns("umap_dims"), "Dimensions:", choices = c("2D" = "2d", "1D" = "1d", "Both" = "both"), selected = "2d"),
              conditionalPanel(
                condition = "input.umap_advanced_mode == false", ns = ns,
                uiOutput(ns("umap_data_source_ui"))
              ),
              checkboxInput(ns("use_semisupervised"), div("Semi-supervised?", help_icon("<strong>Uses known labels to guide how UMAP groups samples.</strong><details><summary>Learn more</summary>Semi-supervised UMAP uses labels you already know (for example, source volcano names or GMM clusters) to <em>guide</em> the projection and pull same-label samples closer together.<br><br><b>Default:</b> OFF.<br><br><b>When to turn ON:</b> When at least some samples are labeled and you want the map to emphasize known groupings.<br><br><b>Example:</b> The VARG26 UMAP used volcanic region, volcanic source (for example, volcano), and tephra where known.<br><br>Blank labels are allowed, which is why this is <em>semi</em>-supervised.<br><br><b>Tip:</b> You can select multiple target variables at once; this is especially useful when combining labeled reference data with unknowns.<br><br><b>Caution:</b> Over-reliance on labels can create artificially clean separations; use Target Weight to balance this.</details>")), value = FALSE),
              conditionalPanel(
                condition = "input.use_semisupervised == true", ns = ns,
                uiOutput(ns("semisupervised_ui")),
                sliderInput(ns("umap_target_weight"), label = div("Target Weight:", help_icon("<strong>Controls how strongly labels steer semi-supervised UMAP.</strong><details><summary>Learn more</summary>Sets how much label information influences the projection.<br><br><b>0</b> = labels ignored (equivalent to unsupervised).<br><b>1</b> = labels dominate the projection.<br><br><b>Default:</b> 0.5 (balanced blend).<br><br><b>When to adjust:</b> Increase toward 0.7–0.8 if known groups are not separating enough; decrease toward 0.2–0.3 if the map looks artificially structured or circular.</details>")), 0, 1, 0.5, step = 0.05)
              ),

              # Advanced Column Selection
              checkboxInput(ns("umap_advanced_mode"), div("Advanced Mode (Custom Columns)", help_icon("<strong>Lets you manually choose exactly which columns UMAP uses.</strong><details><summary>Learn more</summary>Hand-picks which columns enter the UMAP projection.<br><br><b>Default:</b> OFF.<br><br><b>When to turn ON:</b> To mix temporal data (for example, depth or age) with compositional data, focus on specific analytes, or exclude variables.<br><br><b>Note:</b> Only selected variables are used. For example, if combining depth with normalized oxides, select all needed columns manually.<br><br>UMAP requires at least 2 input dimensions; selecting only one column causes an error.</details>")), value = FALSE),
              conditionalPanel(
                condition = "input.umap_advanced_mode == true", ns = ns,
                selectizeInput(ns("umap_custom_cols"), "Select Columns:", choices = NULL, multiple = TRUE, options = list(plugins = list("remove_button", "drag_drop")))
              ),
              hr(),
              tags$details(
                class = "mt-1 mb-2",
                tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("sliders-h"), " Advanced Settings"),
                div(
                  class = "mt-2 p-2 border rounded bg-light",
                  numericInput(ns("umap_n_neighbors"), label = div("Neighbors:", help_icon("<strong>Sets local versus global structure emphasis in UMAP.</strong><details><summary>Learn more</summary>Controls how many nearby points UMAP considers when building the map.<br><br><b>Default:</b> 15, a good balance for most geochemical datasets.<br><br><b>Low values (5–10):</b> Emphasize fine local clusters and tightly grouped sub-populations, but can look noisy or fragmented.<br><br><b>High values (30–50+):</b> Emphasize broad relationships among major groups, creating smoother but less detailed maps.<br><br><b>When to adjust:</b> Lower if small populations are swallowed by larger ones; increase if the map looks too scattered.<br><br><b>Tip:</b> In semi-supervised mode, increasing neighbors is usually advisable.</details>")), 15, min = 2, max = 100),
                  numericInput(ns("umap_min_dist"), label = div("Min Distance:", help_icon("<strong>Controls how tightly points are packed in UMAP space.</strong><details><summary>Learn more</summary>Sets point packing in the projected map.<br><br><b>Default:</b> 0.1, good for distinguishing populations.<br><br><b>Low values (0.0–0.05):</b> Create very tight clumps and stronger separation.<br><br><b>Higher values (0.3–0.5):</b> Spread points out more, preserving broader topology but with less sharp cluster boundaries.<br><br><b>When to adjust:</b> Lower if clusters overlap; increase if you want to show more gradational relationships between groups.</details>")), 0.1, min = 0, max = 0.99, step = 0.01),
                  numericInput(ns("umap_dens_scale"), label = div("Density Scale:", help_icon("<strong>Adjusts LEOPOLD density preservation in the UMAP map.</strong><details><summary>Learn more</summary>Fine-tunes cluster density using LEOPOLD (Lightweight Estimate of Preservation of Local Density).<br><br>Instead of a uniform scale, attractive forces vary by each point's local density, keeping high-density clusters compact and low-density clusters more spread out.<br><br><b>Default:</b> 0 (standard UMAP, no density correction).<br><br><b>Range 0–1:</b> 0.2–0.5 is a good starting point to reveal density differences; 1.0 applies the full density range.<br><br><b>When to adjust:</b> Only if clusters still appear overlaid despite appropriate neighbor and min-distance settings.</details>")), 0, min = 0, max = 1, step = 0.1),
                  numericInput(ns("umap_seed"), label = div("Random Seed:", help_icon("<strong>Sets the random seed for reproducible UMAP outputs.</strong><details><summary>Learn more</summary>UMAP uses pseudorandom number generation: deterministic, but outputs can look random.<br><br>The seed sets the start point in that sequence, so the same seed always reproduces the same result.<br><br><b>Default:</b> 42 (the value itself is arbitrary).<br><br><b>When to change:</b> For sensitivity analysis, try multiple seeds (1, 2, 3, ...). Robust populations should appear regardless of seed.<br><br>If a group appears only for certain seeds, it may not be a reliable population.</details>")), 42, min = 1)
                )
              ),
              actionButton(ns("run_umap_new"), "Create UMAP", class = "btn-primary w-100", icon = icon("play"))
            ),
            varg_card(
              title = tagList(icon("save"), " Save/Load UMAP Model"),
              p(class = "text-muted small", "Save your trained UMAP model to project new data later."),
              uiOutput(ns("umap_model_status_ui")),
              downloadButton(ns("save_umap_model"), "Save UMAP Model", class = "btn-success w-100 mb-2"),
              hr(),
              fileInput(ns("load_umap_model"), "Load UMAP Model", accept = ".varg_umap"),
              actionButton(ns("project_loaded_model"), "Project Using Loaded Model", class = "btn-warning w-100", icon = icon("upload"))
            )
          ),

          # Right column: Plots
          varg_card(
            title = tagList(icon("chart-area"), " UMAP Visualization"),
            shinycssloaders::withSpinner(uiOutput(ns("umap_plots_ui")), type = 4, color = "#2C3E50")
          )
        )
      ),

      # --- Step 5: Define Populations ---
      bslib::nav_panel(
        title = "5. Define Populations",
        icon = icon("users"),
        module_banner(
          goal = "Assign meaningful population names to your samples using interactive selection.",
          inputs = "Imported data, typically preprocessed, clustered, and UMAPed. You can also jump here directly from Step 1 if your data is ready.",
          outputs = "A 'population' column with named identities for each sample.",
          why = "While GMM gives you numbered clusters, this step is where you assign interpretive names (e.g., 'Katmai 1912', 'Hayes F', or simply 'Pop. A', 'Pop. B', 'Feldspar contam.', or 'outlier'). This is especially helpful if you notice repeated populations between samples; there may be a detrital signal. Here, you can easily flag shards that are likely non-primary."
        ),
        uiOutput(ns("step5_status_ui")),
        layout_columns(
          col_widths = c(3, 9),
          varg_card(
            title = tagList(icon("mouse-pointer"), " Selection Tools"),
            div(
              class = "scrollable-config",
            # Instructional Help Box
            help_box(
              title = "Population Assignment",
              content = "Assign final population identities to your samples based on visual clustering.
              <ul>
                <li><b>How it works:</b> Use the <em>lasso</em> or <em>box select</em> tools on the plot to highlight a group of points. Then type a population name and click 'Assign Selection'.</li>
                <li><b>Visualizing GMM Clusters:</b> When clusters are selected on an axis, they are shown as jittered points. This makes it easy to see how individual samples within a GMM cluster relate to other variables or existing labels.</li>
                <li><b>Refinement:</b> You can re-select and re-label any points at any time. To revert points back to unassigned, select them and type 'unassigned'.</li>
                <li><b>Filter:</b> Use the filter dropdowns below the label controls to focus on a specific sample or group, making it easier to work with subsets of your data.</li>
              </ul>
              <small><i>Tip: For the best view, plot your data in UMAP 2D space. Switch to oxide pairs (e.g., SiO2 vs CaO) to validate assignments in compositional space.</i></small>"
            ),
            uiOutput(ns("popsel_xvar_ui")),
            uiOutput(ns("popsel_yvar_ui")),
            uiOutput(ns("popsel_hover_label_ui")),
            
            # Derived Variable Creator
            tags$details(
              class = "mt-2 mb-2",
              tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("calculator"), " Create Ratio or Sum Variable"),
              div(
                class = "mt-2 p-2 border rounded bg-light",
                p(class = "small text-muted mb-2", "Create a new variable from two existing ones (e.g., SiO\u2082/Al\u2082O\u2083 or Na\u2082O+K\u2082O). It will appear in the X/Y dropdowns above."),
                uiOutput(ns("derived_var_a_ui")),
                selectInput(ns("derived_op"), "Operation:", choices = c("A / B (ratio)" = "ratio", "A + B (sum)" = "sum"), selected = "ratio", width = "100%"),
                uiOutput(ns("derived_var_b_ui")),
                actionButton(ns("create_derived_var"), "Create Variable", class = "btn-outline-primary btn-sm w-100 mt-1", icon = icon("plus")),
                uiOutput(ns("derived_var_status"))
              )
            ),
            hr(),
            # Initialize populations from an existing column
            tags$details(
              class = "mb-2",
              tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("magic"), " Initialize from Column"),
              div(
                class = "mt-2 p-2 border rounded bg-light",
                p(class = "small text-muted mb-2", "Bulk-assign population labels from an existing column (e.g., GMM cluster, sample ID, volcanic source). Rows with NA or empty values become 'Unassigned'."),
                uiOutput(ns("init_pop_col_ui")),
                actionButton(ns("init_pop_from_col"), "Initialize Populations", class = "btn-outline-info btn-sm w-100 mt-1", icon = icon("magic")),
                uiOutput(ns("init_pop_status"))
              )
            ),
            selectizeInput(
              ns("feature_label"),
              label = div("Assign Label:", help_icon("<strong>Assigns selected points to a new or existing population label.</strong><details><summary>Learn more</summary>Type a new population name or choose an existing one from the dropdown, then click <b>Assign Selection</b>.<br><br>Use the <b>Revert</b> option below to return points to <b>Unassigned</b>.<br><br>You can also use <b>Initialize from Column</b> to bulk-assign labels from an existing variable.</details>")),
              choices = NULL,
              selected = NULL,
              options = list(
                create = TRUE,
                dropdownParent = "body",
                placeholder = "Start typing a population name...",
                allowEmptyOption = TRUE,
                plugins = list('clear_button'),
                onInitialize = I('function() { this.clear(); }')
              )
            ),
            div(
              class = "d-flex gap-2",
              actionButton(ns("assign_feature"), "Assign Selection", class = "btn-warning flex-grow-1", icon = icon("tag")),
              shinyjs::disabled(
                actionButton(ns("undo_assign_feature"), "Undo Last Assignment", class = "btn-outline-secondary", icon = icon("undo"))
              )
            ),
            
            # Revert entire population to Unassigned
            tags$details(
              class = "mt-2 mb-2",
              tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("undo"), " Revert a population to Unassigned"),
              div(
                class = "mt-2 p-2 border rounded bg-light",
                uiOutput(ns("revert_pop_ui")),
                actionButton(ns("revert_population"), "Revert All Points", class = "btn-outline-danger btn-sm w-100 mt-1", icon = icon("undo"))
              )
            ),
            # Rename a population
            tags$details(
              class = "mb-2",
              tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("pen"), " Rename a population"),
              div(
                class = "mt-2 p-2 border rounded bg-light",
                p(class = "small text-muted mb-2", "Rename all points assigned to a population. Useful for replacing GMM cluster numbers with interpretive names (e.g., '2' → 'Redoubt')."),
                uiOutput(ns("rename_pop_ui")),
                textInput(ns("rename_pop_new"), "New name:", placeholder = "e.g., Redoubt"),
                actionButton(ns("rename_population"), "Rename", class = "btn-outline-primary btn-sm w-100 mt-1", icon = icon("pen"))
              )
            ),
            hr(),
            uiOutput(ns("filter_col_ui")),
            uiOutput(ns("filter_val_ui"))
            ) # end scrollable-config
          ),
          varg_card(
            title = tagList(icon("braille"), " Interactive Plot"),
            tags$div(
              class = "alert alert-info py-1 px-2 mb-2 small",
              style = "border-left: 3px solid #0d6efd;",
              icon("mouse-pointer"), " ",
              "Draw around points with the ",
              tags$b("Lasso Select"),
              " tool (active by default), then type a name and click ",
              tags$b("Assign Selection"),
              " on the left."
            ),
            plotlyOutput(ns("featplot"), height = "600px")
          )
        )
      ),

      # --- Step 6: Export ---
      bslib::nav_panel(
        title = "6. Export",
        icon = icon("file-export"),
        module_banner(
          goal = "Download your processed data with all generated columns.",
          inputs = "Your working dataset with transformations, clusters, UMAP, and population labels.",
          outputs = "CSV or XLSX file with original data plus all new columns.",
          why = "Export your results for use in other software (Excel, R, Python), for publication, or to re-import into the Visualization or Chronology modules."
        ),
        uiOutput(ns("step6_status_ui")),
        uiOutput(ns("export_preview_ui")),
        varg_card(
          title = tagList(icon("download"), " Download Results"),
          tags$p(class = "text-muted small mb-3", 
            icon("info-circle"), 
            " Generated columns (normalized, imputed, UMAP, clusters, populations) are added to the ", 
            tags$strong("right"), 
            " of your original data. UID column is added at the left for row matching."
          ),
          div(
            class = "d-flex gap-2",
            downloadButton(ns("download_csv"), "Download CSV", class = "btn-success"),
            downloadButton(ns("download_xlsx"), "Download XLSX", class = "btn-primary"),
            downloadButton(ns("download_report"), "Download Processing Report", class = "btn-outline-secondary", icon = icon("file-alt"))
          )
        )
      )
    )
  )
}

umap_result_tab_name <- function(mode, dimensions) {
  mode_prefix <- switch(
    tolower(as.character(mode)[1]),
    pretrained = "VARG26",
    new = "New",
    loaded = "Loaded",
    NULL
  )
  if (is.null(mode_prefix)) return(NULL)

  dimension_values <- if (length(dimensions) == 1L && is.character(dimensions)) {
    switch(
      tolower(dimensions),
      "1d" = 1L,
      "2d" = 2L,
      "both" = c(1L, 2L),
      integer(0)
    )
  } else {
    suppressWarnings(as.integer(dimensions))
  }
  dimension_values <- dimension_values[dimension_values %in% c(1L, 2L)]
  if (length(dimension_values) == 0L) return(NULL)

  paste(mode_prefix, if (2L %in% dimension_values) "2D" else "1D")
}

mod_processing_server <- function(id, global_rv = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    focus_umap_result_tab <- function(mode, dimensions) {
      selected_tab <- umap_result_tab_name(mode, dimensions)
      if (is.null(selected_tab)) return(invisible(NULL))

      session$onFlushed(
        function() updateTabsetPanel(session, "umap_tab", selected = selected_tab),
        once = TRUE
      )
      invisible(selected_tab)
    }

    # --- Reactive Values ---
    # Use global_rv if provided, otherwise create local
    rv <- if (!is.null(global_rv)) {
      global_rv
    } else {
      reactiveValues(
        data = NULL,
        umap_mode_ran = "new",
        mclust_result = NULL,
        data_stale = FALSE,
        pipeline_config = NULL,
        user_umap_model = NULL,
        source_filename = NULL,
        data_generation = 0L
      )
    }

    # Local reactive for GMM running state (UI only)
    gmm_running <- reactiveVal(FALSE)
    gmm_bg_process <- reactiveVal(NULL)  # callr background process handle
    gmm_bg_context <- reactiveValues(
      valid = NULL, src = NULL, notif_id = NULL, include_uncertainty = FALSE,
      noise_init = FALSE, use_prior = FALSE, gmin = NULL, gmax = NULL, token = NULL
    )
    output$gmm_is_running <- reactive({
      gmm_running()
    })
    outputOptions(output, "gmm_is_running", suspendWhenHidden = FALSE)

    # Local reactive for preprocessing background process
    preprocess_bg_process <- reactiveVal(NULL)
    preprocess_bg_context <- reactiveValues(notif_id = NULL, token = NULL)

    # Local reactive for VARG26 projection background process
    varg26_bg_process <- reactiveVal(NULL)
    varg26_bg_context <- reactiveValues(notif_id = NULL, dims_choice = NULL, token = NULL)

    # Warn if user requests very high G (may be slow)
    observeEvent(list(input$gmin, input$gmax), {
      gvals <- c(input$gmin, input$gmax)
      if (any(!is.null(gvals) & gvals > 30, na.rm = TRUE)) {
        showNotification("Warning: G values above 30 may take a while to compute.", type = "warning", duration = 5)
      }
    })

    # Initialize slots if using global_rv (they may not exist yet)
    if (!is.null(global_rv)) {
      if (is.null(isolate(rv$pipeline_config))) rv$pipeline_config <- NULL
      if (is.null(isolate(rv$user_umap_model))) rv$user_umap_model <- NULL
      if (is.null(isolate(rv$source_filename))) rv$source_filename <- NULL
      if (is.null(isolate(rv$heavy_job))) rv$heavy_job <- NULL
      if (is.null(isolate(rv$data_generation))) rv$data_generation <- 0L
    }

    heavy_job_lock <- NULL

    begin_heavy_job <- function(label, min_workers = 1L, max_workers = min_workers) {
      active <- isolate(rv$heavy_job)
      if (!is.null(active) && nzchar(active)) {
        showNotification(
          paste0("Another analysis is already running: ", active, ". Cancel it or wait for it to finish before starting a new heavy job."),
          type = "warning",
          duration = 8
        )
        return(FALSE)
      }

      global_slot <- heavy_job_limiter_acquire(
        label,
        min_tokens = min_workers,
        max_tokens = max_workers
      )
      if (!isTRUE(global_slot$acquired)) {
        showNotification(global_slot$message, type = "warning", duration = 10)
        return(FALSE)
      }

      heavy_job_lock <<- global_slot$lock
      rv$heavy_job <- label
      TRUE
    }

    finish_heavy_job <- function() {
      rv$heavy_job <- NULL
      if (!is.null(heavy_job_lock)) {
        heavy_job_limiter_release(heavy_job_lock)
        heavy_job_lock <<- NULL
      }
    }

    pop_styles <- reactiveVal(list(colors = setNames("grey50", "Unassigned"), shapes = setNames("circle", "Unassigned")))
    last_feature_label <- reactiveVal(NULL) # Track last-selected label to preserve dropdown selection
    undo_history <- reactiveVal(list())
    max_undo_levels <- 10L

    push_undo_snapshot <- function(uids, old_labels) {
      history <- undo_history()
      history[[length(history) + 1]] <- list(uids = uids, old_labels = old_labels)
      if (length(history) > max_undo_levels) {
        history <- tail(history, max_undo_levels)
      }
      undo_history(history)
    }

    # Output for conditional panel
    output$data_is_stale <- reactive({
      return(rv$data_stale)
    })
    outputOptions(output, "data_is_stale", suspendWhenHidden = FALSE)

    # Store filter state separately to prevent reset
    filter_state <- reactiveValues(col = "", val = "All", initialized = FALSE)
    
    # Store popsel x/y variable state to prevent reset on data changes
    popsel_state <- reactiveValues(xvar = NULL, yvar = NULL, hover_label = NULL, initialized = FALSE)
    
    # Store preprocessing column selections to prevent reset after Apply Transforms
    preproc_state <- reactiveValues(comp_cols = NULL, noncomp_cols = NULL, transform_type = "ilr", pivot_var = NULL, initialized = FALSE)

    # Store GMM settings to prevent reset after running
    gmm_ui_state <- reactiveValues(data_source = NULL, initialized = FALSE)

    # Store UMAP settings to prevent reset after running
    umap_ui_state <- reactiveValues(data_source = NULL, dims = NULL, semisupervised = NULL, 
                                     y_vars = NULL, target_weight = NULL, n_neighbors = NULL,
                                     min_dist = NULL, dens_scale = NULL, seed = NULL, initialized = FALSE)

    # Track if we're restoring from a loaded project (to avoid triggering stale state)
    restoring_state <- reactiveVal(FALSE)

    # Constants
    ADD_NOISE_LEVEL <- 1e-5

    # --- Step Status Banners ---
    # Helper to create status banner
    status_banner <- function(complete, message, warning = NULL) {
      if (!is.null(warning)) {
        div(
          class = "alert alert-warning py-2 mb-3",
          icon("exclamation-triangle"), " ", warning
        )
      } else if (complete) {
        div(
          class = "alert alert-success py-2 mb-3",
          icon("check-circle"), " ", message
        )
      } else {
        div(
          class = "alert alert-info py-2 mb-3",
          icon("info-circle"), " ", message
        )
      }
    }

    # Step 1: Data Import
    output$step1_status_ui <- renderUI({
      if (!is.null(rv$data) && nrow(rv$data) > 0) {
        status_banner(TRUE, paste("Data loaded:", nrow(rv$data), "rows,", ncol(rv$data), "columns"))
      } else {
        status_banner(FALSE, "Upload a CSV/Excel file or load a saved project to begin.")
      }
    })

    # Step 2: Preprocessing
    output$step2_status_ui <- renderUI({
      if (is.null(rv$data)) {
        return(status_banner(FALSE, "Please import data in Step 1 first.", warning = "No data loaded"))
      }

      has_transforms <- !is.null(rv$original_cols) && length(setdiff(names(rv$data), rv$original_cols)) > 0

      if (rv$data_stale) {
        status_banner(FALSE, NULL, warning = "Settings changed. Click 'Apply Transforms' to update.")
      } else if (has_transforms) {
        pivot_cols <- grep("^pivot_", names(rv$data), value = TRUE)
        clr_cols <- grep("^clr_", names(rv$data), value = TRUE)
        status_banner(TRUE, paste("Transforms applied.", length(pivot_cols), "pivot coordinates and", length(clr_cols), "CLR coordinates generated."))
      } else {
        status_banner(FALSE, "Select compositional columns and apply transforms.")
      }
    })

    # Step 3: Clustering (GMM)
    output$step3_status_ui <- renderUI({
      if (is.null(rv$data)) {
        return(status_banner(FALSE, NULL, warning = "No data loaded. Start at Step 1."))
      }

      has_transforms <- !is.null(rv$original_cols) && length(setdiff(names(rv$data), rv$original_cols)) > 0
      if (!has_transforms) {
        return(status_banner(FALSE, NULL, warning = "Apply preprocessing in Step 2 first."))
      }

      if (!is.null(rv$mclust_result)) {
        n_clusters <- rv$mclust_result$G
        status_banner(TRUE, paste("GMM complete:", n_clusters, "clusters identified."))
      } else {
        status_banner(FALSE, "Configure GMM settings and run clustering.")
      }
    })

    # Step 4: UMAP Projection
    output$step4_status_ui <- renderUI({
      if (is.null(rv$data)) {
        return(status_banner(FALSE, NULL, warning = "No data loaded. Start at Step 1."))
      }

      has_umap <- any(grepl("^(V1|V2|UMAP)", names(rv$data)))
      has_transforms <- !is.null(rv$original_cols) && length(setdiff(names(rv$data), rv$original_cols)) > 0

      if (has_umap) {
        umap_mode <- if (!is.null(rv$umap_mode_ran)) rv$umap_mode_ran else "unknown"
        status_banner(TRUE, paste("UMAP coordinates generated (", umap_mode, " mode).", sep = ""))
      } else if (!has_transforms) {
        status_banner(FALSE, "Apply preprocessing in Step 2 first. For VARG26, you also need the 9 major oxides.")
      } else {
        status_banner(FALSE, "Choose VARG26 pretrained projection or create a new UMAP.")
      }
    })

    # Step 5: Define Populations
    output$step5_status_ui <- renderUI({
      if (is.null(rv$data)) {
        return(status_banner(FALSE, NULL, warning = "No data loaded. Start at Step 1."))
      }

      has_pop <- "population" %in% names(rv$data)
      has_umap <- any(grepl("^(V1|V2|UMAP)", names(rv$data)))

      if (has_pop) {
        n_pops <- length(unique(na.omit(rv$data$population)))
        status_banner(TRUE, paste(n_pops, "populations defined."))
      } else if (!has_umap) {
        status_banner(FALSE, "Use lasso selection to define population groups. Tip: Run UMAP in Step 4 first for easier grouping.")
      } else {
        status_banner(FALSE, "Use lasso selection to define population groups.")
      }
    })

    # Step 6: Export
    export_ready <- reactiveVal(FALSE)
    # Reset export state when data changes
    observeEvent(rv$data, {
      export_ready(FALSE)
    }, ignoreNULL = FALSE, priority = 10)
    output$step6_status_ui <- renderUI({
      if (is.null(rv$data)) {
        return(status_banner(FALSE, NULL, warning = "No data to export. Complete earlier steps first."))
      }
      if (!export_ready()) {
        return(status_banner(FALSE, "Preparing data for export\u2026"))
      }
      status_banner(TRUE, "Data ready for export.")
    })
    # Pre-render export outputs even when tab is hidden so they're ready on switch
    outputOptions(output, "step6_status_ui", suspendWhenHidden = FALSE)

    output$export_preview_ui <- renderUI({
      req(rv$data)

      export_df <- tryCatch(prepare_export_data(rv$data), error = function(e) NULL)
      if (is.null(export_df)) {
        export_ready(FALSE)
        return(
          div(
            class = "alert alert-warning py-2 mb-3",
            icon("exclamation-triangle"),
            " Export preview unavailable right now."
          )
        )
      }

      export_ready(TRUE)

      original_cols <- if (!is.null(rv$original_cols)) intersect(rv$original_cols, names(export_df)) else character(0)
      norm_cols <- grep("_norm$", names(export_df), value = TRUE)
      imp_cols <- grep("_imp$", names(export_df), value = TRUE)
      pivot_cols <- grep("^pivot_", names(export_df), value = TRUE)
      clr_cols <- grep("^clr_", names(export_df), value = TRUE)
      umap_cols <- grep("^UMAP_", names(export_df), value = TRUE)
      gmm_cols <- grep("^gmm_", names(export_df), value = TRUE)
      population_cols <- intersect("population", names(export_df))

      column_group_ui <- function(label, cols) {
        tags$details(
          class = "mb-1",
          tags$summary(
            class = "small",
            style = "cursor: pointer;",
            paste0(label, " (", length(cols), ")")
          ),
          if (length(cols) > 0) {
            tags$ul(
              class = "small text-muted mb-1",
              style = "padding-left: 1.25rem;",
              lapply(cols, tags$li)
            )
          } else {
            tags$div(class = "small text-muted", "None")
          }
        )
      }

      varg_card(
        title = tagList(icon("table"), " Export Preview"),
        tags$p(
          class = "small mb-2",
          sprintf("%s rows x %s columns ready for export", nrow(export_df), ncol(export_df))
        ),
        DTOutput(ns("export_preview_table")),
        tags$details(
          class = "mt-2",
          tags$summary(class = "small text-muted", style = "cursor: pointer;", icon("list"), " Column groups"),
          div(
            class = "mt-2",
            column_group_ui("Original", original_cols),
            column_group_ui("_norm", norm_cols),
            column_group_ui("_imp", imp_cols),
            column_group_ui("pivot_", pivot_cols),
            column_group_ui("clr_", clr_cols),
            column_group_ui("UMAP_", umap_cols),
            column_group_ui("gmm_", gmm_cols),
            column_group_ui("population", population_cols)
          )
        )
      )
    })
    outputOptions(output, "export_preview_ui", suspendWhenHidden = FALSE)

    output$export_preview_table <- renderDT({
      req(rv$data)
      export_df <- prepare_export_data(rv$data)
      datatable(
        head(export_df, 5),
        options = list(
          dom = "t",
          scrollX = TRUE,
          paging = FALSE
        ),
        rownames = FALSE
      )
    })

    # --- Restore UI State from Loaded Project ---
    reset_processing_ui_state <- function() {
      preproc_state$comp_cols <- NULL
      preproc_state$noncomp_cols <- NULL
      preproc_state$transform_type <- "ilr"
      preproc_state$pivot_var <- NULL
      preproc_state$initialized <- FALSE
      gmm_ui_state$data_source <- NULL
      gmm_ui_state$initialized <- FALSE
      umap_ui_state$data_source <- NULL
      umap_ui_state$dims <- NULL
      umap_ui_state$semisupervised <- NULL
      umap_ui_state$y_vars <- NULL
      umap_ui_state$target_weight <- NULL
      umap_ui_state$n_neighbors <- NULL
      umap_ui_state$min_dist <- NULL
      umap_ui_state$dens_scale <- NULL
      umap_ui_state$seed <- NULL
      umap_ui_state$initialized <- FALSE

      updateCheckboxInput(session, "do_impute", value = TRUE)
      updateSelectInput(session, "impute_method", selected = "auto")
      updateSelectInput(session, "transform_type", selected = "ilr")
      updateCheckboxInput(session, "gmm_advanced_mode", value = FALSE)
      updateSelectizeInput(session, "gmm_custom_cols", selected = character(0))
      updateNumericInput(session, "gmin", value = 2)
      updateNumericInput(session, "gmax", value = 15)
      updateCheckboxInput(session, "noise_init", value = TRUE)
      updateCheckboxInput(session, "use_prior", value = FALSE)
      updateCheckboxInput(session, "include_uncertainty", value = FALSE)
      updateSelectInput(session, "VARG26_dims", selected = "2d")
      updateSelectInput(session, "umap_dims", selected = "2d")
      updateCheckboxInput(session, "umap_advanced_mode", value = FALSE)
      updateSelectizeInput(session, "umap_custom_cols", selected = character(0))
      updateCheckboxInput(session, "use_semisupervised", value = FALSE)
      updateNumericInput(session, "umap_n_neighbors", value = 15)
      updateNumericInput(session, "umap_min_dist", value = 0.1)
      updateNumericInput(session, "umap_dens_scale", value = 0)
      updateNumericInput(session, "umap_seed", value = 42)
      updateSliderInput(session, "umap_target_weight", value = 0.5)
      invisible(NULL)
    }

    # Watch for restore_trigger changes (from project load) and restore inputs
    observeEvent(rv$restore_trigger,
      {
        req(rv$restore_trigger > 0) # Skip initial value
        config <- rv$pipeline_config
        if (is.null(rv$data)) return()

        restoring_state(TRUE)
        on.exit(
          session$onFlushed(function() restoring_state(FALSE), once = TRUE),
          add = TRUE
        )
        reset_processing_ui_state()

        if (is.null(config)) {
          rv$data_stale <- FALSE
          showNotification(
            "Loaded data-only project; processing controls were reset to current defaults.",
            type = "warning",
            duration = 6
          )
          return()
        }

        prep <- config$preprocessing
        if (!is.null(prep)) {
          # Restore compositional columns
          if (!is.null(prep$comp_cols)) {
            valid_cols <- intersect(prep$comp_cols, names(rv$data))
            if (length(valid_cols) > 0) {
              updateCheckboxGroupInput(session, "comp_cols", selected = valid_cols)
            }
          }

          # Restore non-compositional columns
          if (!is.null(prep$noncomp_cols)) {
            valid_cols <- intersect(prep$noncomp_cols, names(rv$data))
            if (length(valid_cols) > 0) {
              updateCheckboxGroupInput(session, "noncomp_cols", selected = valid_cols)
            }
          }

          if (!is.null(prep$do_impute)) {
            updateCheckboxInput(session, "do_impute", value = prep$do_impute)
          }
          restored_impute_method <- prep$impute_method %||% prep$imputation_method %||% "auto"
          if (!(restored_impute_method %in% c("auto", "ltsReg", "lm"))) {
            restored_impute_method <- "auto"
          }
          updateSelectInput(session, "impute_method", selected = restored_impute_method)

          restored_transform_type <- prep$transform_type
          if (is.null(restored_transform_type)) {
            if (!is.null(prep$do_pivot)) {
              restored_transform_type <- if (isTRUE(prep$do_pivot)) "ilr" else "none"
            } else {
              restored_transform_type <- "ilr"
            }
          }
          if (!(restored_transform_type %in% c("ilr", "clr", "none"))) {
            restored_transform_type <- "ilr"
          }
          updateSelectInput(session, "transform_type", selected = restored_transform_type)

          if (!is.null(prep$pivot_var) && identical(restored_transform_type, "ilr")) {
            shinyjs::delay(100, {
              updateSelectInput(session, "pivot_var", selected = prep$pivot_var)
            })
          }
        }

        gmm <- config$gmm
        if (!is.null(gmm)) {
          gmm_ui_state$data_source <- gmm$data_source
          gmm_ui_state$initialized <- TRUE
          updateCheckboxInput(session, "gmm_advanced_mode", value = isTRUE(gmm$advanced_mode))
          if (!is.null(gmm$columns_used)) {
            updateSelectizeInput(
              session,
              "gmm_custom_cols",
              selected = intersect(gmm$columns_used, names(rv$data))
            )
          }
          if (!is.null(gmm$gmin)) updateNumericInput(session, "gmin", value = gmm$gmin)
          if (!is.null(gmm$gmax)) updateNumericInput(session, "gmax", value = gmm$gmax)
          if (!is.null(gmm$noise_init)) updateCheckboxInput(session, "noise_init", value = isTRUE(gmm$noise_init))
          if (!is.null(gmm$use_prior)) updateCheckboxInput(session, "use_prior", value = isTRUE(gmm$use_prior))
          if (!is.null(gmm$include_uncertainty)) {
            updateCheckboxInput(session, "include_uncertainty", value = isTRUE(gmm$include_uncertainty))
          }
        }

        umap <- config$umap
        if (!is.null(umap)) {
          if (identical(umap$mode, "pretrained")) {
            if (!is.null(umap$VARG26_dims)) {
              updateSelectInput(session, "VARG26_dims", selected = umap$VARG26_dims)
            }
          } else if (identical(umap$mode, "new")) {
            dimensions <- normalize_umap_dimensions(umap)
            dims_choice <- if (identical(dimensions, c(1L, 2L))) {
              "both"
            } else if (identical(dimensions, 1L)) {
              "1d"
            } else {
              "2d"
            }
            updateSelectInput(session, "umap_dims", selected = dims_choice)
            updateCheckboxInput(session, "umap_advanced_mode", value = isTRUE(umap$advanced_mode))
            if (!is.null(umap$columns_used)) {
              updateSelectizeInput(
                session,
                "umap_custom_cols",
                selected = intersect(umap$columns_used, names(rv$data))
              )
            }
            umap_ui_state$data_source <- umap$data_source
            umap_ui_state$y_vars <- umap$y_vars
            umap_ui_state$initialized <- TRUE
            updateCheckboxInput(session, "use_semisupervised", value = isTRUE(umap$semisupervised))
            if (!is.null(umap$n_neighbors)) updateNumericInput(session, "umap_n_neighbors", value = umap$n_neighbors)
            if (!is.null(umap$min_dist)) updateNumericInput(session, "umap_min_dist", value = umap$min_dist)
            if (!is.null(umap$dens_scale)) updateNumericInput(session, "umap_dens_scale", value = umap$dens_scale)
            if (!is.null(umap$seed)) updateNumericInput(session, "umap_seed", value = umap$seed)
            if (!is.null(umap$target_weight)) {
              updateSliderInput(session, "umap_target_weight", value = umap$target_weight)
            }
          }
        }

        # Mark data as not stale since we just loaded a complete project
        rv$data_stale <- FALSE

        showNotification("Pipeline settings restored from project.", type = "message", duration = 3)
      },
      ignoreInit = TRUE,
      ignoreNULL = TRUE
    )

    # --- Session Cleanup ---
    # Unload any user UMAP models when session ends to free resources
    session$onSessionEnded(function() {
      for (process_getter in list(preprocess_bg_process, gmm_bg_process, varg26_bg_process)) {
        proc <- tryCatch(isolate(process_getter()), error = function(e) NULL)
        if (!is.null(proc)) {
          try({
            if (proc$is_alive()) {
              proc$kill()
              proc$wait(timeout = 2000)
            }
          }, silent = TRUE)
        }
      }
      try(finish_heavy_job(), silent = TRUE)
      tryCatch(
        {
          if (!is.null(isolate(rv$user_umap_model))) {
            unload_umap_model_collection(
              isolate(rv$user_umap_model),
              if (is.null(isolate(rv$pipeline_config))) NULL else isolate(rv$pipeline_config$umap)
            )
          }
        },
        error = function(e) {
          # Ignore cleanup errors silently
        }
      )
    })

    # --- Helpers ---
    show_loading <- function(msg) {
      showModal(modalDialog(title = "Processing", msg, easyClose = FALSE, footer = NULL))
    }
    hide_loading <- function() {
      removeModal()
    }
    clear_user_umap_model <- function() {
      model <- isolate(rv$user_umap_model)
      if (!is.null(model)) {
        config <- isolate(rv$pipeline_config)
        unload_umap_model_collection(model, if (is.null(config)) NULL else config$umap)
      }
      rv$user_umap_model <- NULL
      invisible(NULL)
    }
    bump_data_generation <- function() {
      current <- isolate(rv$data_generation)
      if (is.null(current) || length(current) != 1 || is.na(current)) current <- 0L
      rv$data_generation <- as.integer(current) + 1L
      invisible(rv$data_generation)
    }
    capture_data_token <- function(context = NULL) {
      varg_make_data_token(isolate(rv$data), isolate(rv$data_generation), context)
    }
    data_token_is_current <- function(token, context = NULL) {
      varg_data_token_matches(token, isolate(rv$data), isolate(rv$data_generation), context)
    }
    current_preprocess_context <- function() {
      transform_type <- input$transform_type %||% "ilr"
      do_impute <- isTRUE(input$do_impute)
      list(
        comp_cols = input$comp_cols,
        noncomp_cols = input$noncomp_cols,
        do_impute = do_impute,
        impute_method = if (do_impute) (input$impute_method %||% "auto") else "auto",
        transform_type = transform_type,
        pivot_var = if (identical(transform_type, "ilr")) input$pivot_var else NULL
      )
    }
    current_gmm_context <- function() {
      list(
        advanced_mode = isTRUE(input$gmm_advanced_mode),
        data_source = input$gmm_data_source,
        custom_cols = input$gmm_custom_cols,
        include_uncertainty = isTRUE(input$include_uncertainty),
        gmin = input$gmin,
        gmax = input$gmax,
        noise_init = isTRUE(input$noise_init),
        use_prior = isTRUE(input$use_prior)
      )
    }
    current_varg26_context <- function() {
      list(dimensions = input$VARG26_dims)
    }
    cancel_processing_background_jobs <- function() {
      process_slots <- list(
        list(get = preprocess_bg_process, set = preprocess_bg_process),
        list(get = gmm_bg_process, set = gmm_bg_process),
        list(get = varg26_bg_process, set = varg26_bg_process)
      )
      for (slot in process_slots) {
        proc <- tryCatch(isolate(slot$get()), error = function(e) NULL)
        if (!is.null(proc)) {
          try({
            if (proc$is_alive()) {
              proc$kill()
              proc$wait(timeout = 2000)
            }
          }, silent = TRUE)
          slot$set(NULL)
        }
      }
      gmm_running(FALSE)
      shinyjs::hide("cancel_preprocess")
      shinyjs::hide("cancel_gmm")
      shinyjs::hide("cancel_varg26")
      shinyjs::hide("gmm_progress_container")
      shinyjs::enable("apply_transforms")
      shinyjs::enable("run_gmm")
      shinyjs::enable("run_umap_pre")
      for (notification_id in list(
        isolate(preprocess_bg_context$notif_id),
        isolate(gmm_bg_context$notif_id),
        isolate(varg26_bg_context$notif_id)
      )) {
        if (!is.null(notification_id)) try(removeNotification(notification_id), silent = TRUE)
      }
      finish_heavy_job()
      invisible(NULL)
    }
    if (!is.null(global_rv)) rv$cancel_processing_jobs <- cancel_processing_background_jobs

    # --- Download Handlers ---
    output$download_template <- downloadHandler(
      filename = function() {
        "VARG_Tools_ProcessViz_Template.xlsx"
      },
      content = function(file) {
        sheets <- generate_template("geochem")
        writexl::write_xlsx(sheets, path = file)
      }
    )

    # --- UI Outputs ---
    output$sheet_ui <- renderUI({
      req(input$file_upload)
      file_ext <- tolower(tools::file_ext(input$file_upload$name))
      if (file_ext %in% c("xlsx", "xls")) {
        sheets <- tryCatch(readxl::excel_sheets(input$file_upload$datapath), error = function(e) NULL)
        if (!is.null(sheets)) {
          selectInput(ns("sheet_select"), label = tags$span("Select worksheet:", help_icon("<strong>Select the worksheet to import; only that sheet is read.</strong>")), choices = sheets)
        }
      }
    })
    
    # Row filtering UI - shown after data is loaded
    output$row_filter_ui <- renderUI({
      req(rv$data)
      
      # Detect blank rows (all analyte-type numeric columns are NA)
      # Exclude common non-analyte numeric columns from the check
      non_analyte_cols <- c("UID", "Depth", "depth", "Age", "age", "Latitude", "Longitude",
                            "latitude", "longitude", "lat", "lon", "Elevation", "elevation",
                            "row_number", "X", "x", "Y", "y", "Z", "z")
      numeric_cols <- names(rv$data)[sapply(rv$data, is.numeric) & !names(rv$data) %in% c("UID", non_analyte_cols)]
      if (length(numeric_cols) > 0) {
        blank_rows <- apply(rv$data[, numeric_cols, drop = FALSE], 1, function(x) all(is.na(x)))
        n_blank <- sum(blank_rows)
      } else {
        n_blank <- 0
      }
      
      # Get text columns for filtering
      text_cols <- names(rv$data)[sapply(rv$data, function(x) is.character(x) || is.factor(x))]
      
      div(
        class = "mt-3 p-2 border rounded bg-light",
        tags$h6(icon("filter"), " Row Filtering", class = "text-primary mb-2"),
        
        # Blank row detection
        if (n_blank > 0) {
          div(
            class = "alert alert-warning py-2 mb-2",
            style = "font-size: 0.85rem;",
            icon("exclamation-triangle"),
            sprintf(" %d blank row%s detected (all numeric values are NA).", n_blank, if (n_blank > 1) "s" else ""),
            div(
              class = "mt-2",
              checkboxInput(ns("exclude_blank_rows"), 
                            label = "Exclude blank rows from processing", 
                            value = TRUE)
            )
          )
        } else {
          div(class = "text-success small mb-2", icon("check"), " No blank rows detected.")
        },
        
        # Text pattern filter
        if (length(text_cols) > 0) {
          tagList(
            div(
              class = "small text-muted mb-1",
              "Filter out summary rows (e.g., Mean, StDev):"
            ),
            div(
              style = "max-width: 100%;",
              selectInput(ns("filter_text_col"), 
                          label = "Column to check:",
                          choices = c("None" = "", text_cols),
                          selected = "",
                          width = "100%")
            ),
            div(
              style = "max-width: 100%;",
              textInput(ns("filter_text_pattern"), 
                        label = "Pattern to exclude:",
                        placeholder = "e.g., Mean|StDev|Average",
                        width = "100%")
            ),
            div(
              class = "text-muted small",
              style = "font-size: 0.75rem;",
              icon("info-circle"), " Use | to separate patterns (no spaces around |)"
            ),
            uiOutput(ns("filter_preview_ui"))
          )
        },
        
        # Apply filter button
        actionButton(ns("apply_row_filter"), "Apply Row Filter", 
                     class = "btn-outline-secondary btn-sm w-100 mt-2", 
                     icon = icon("filter"))
      )
    })
    
    # Preview matching rows for text filter
    output$filter_preview_ui <- renderUI({
      req(input$filter_text_col, input$filter_text_col != "")
      req(input$filter_text_pattern, nchar(trimws(input$filter_text_pattern)) > 0)
      req(rv$data)
      
      tryCatch({
        pattern <- trimws(input$filter_text_pattern)
        col_data <- as.character(rv$data[[input$filter_text_col]])
        # Split on | and use fixed matching to avoid regex errors
        # Case-insensitive via tolower (fixed=TRUE ignores ignore.case arg)
        parts <- strsplit(pattern, "\\|")[[1]]
        parts <- trimws(parts)
        parts <- parts[nchar(parts) > 0]
        if (length(parts) == 0) return(NULL)
        col_lower <- tolower(col_data)
        matches <- Reduce(`|`, lapply(parts, function(p) {
          grepl(tolower(p), col_lower, fixed = TRUE)
        }))
        n_matches <- sum(matches)
        
        if (n_matches > 0) {
          # Show sample of matching values
          matching_vals <- unique(col_data[matches])
          sample_vals <- head(matching_vals, 3)
          
          div(
            class = "alert alert-info py-2 mt-2",
            style = "font-size: 0.8rem;",
            sprintf("%d row%s match%s: ", n_matches, 
                    if (n_matches > 1) "s" else "", 
                    if (n_matches == 1) "es" else ""),
            tags$em(paste(sample_vals, collapse = ", ")),
            if (length(matching_vals) > 3) paste0(" + ", length(matching_vals) - 3, " more...")
          )
        }
      }, error = function(e) NULL)
    })

    output$compcol_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric)]
      nums <- nums[!grepl("UID|^pivot_|^clr_|^UMAP_|_imp$|_norm$|_raw$|gmm_cluster", nums)]
      
      # Use stored state if initialized and valid, otherwise use defaults
      if (preproc_state$initialized && !is.null(preproc_state$comp_cols)) {
        # Keep only selections that are still valid
        selected <- intersect(preproc_state$comp_cols, nums)
      } else {
        # Default to VARG26 oxides
        selected <- intersect(nums, VARG26_OXIDES)
      }
      
      checkboxGroupInput(ns("comp_cols"), "Compositional columns:", choices = nums, selected = selected)
    })
    
    # Track comp_cols changes
    observeEvent(input$comp_cols, {
      preproc_state$comp_cols <- input$comp_cols
      preproc_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    output$noncomp_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric)]
      nums <- nums[!grepl("UID|^pivot_|^clr_|^UMAP_|_imp$|_norm$|_raw$|gmm_cluster", nums)]
      # Exclude already selected compositional columns
      if (!is.null(input$comp_cols)) {
        nums <- setdiff(nums, input$comp_cols)
      }
      if (length(nums) > 0) {
        # Use stored state if initialized and valid
        if (preproc_state$initialized && !is.null(preproc_state$noncomp_cols)) {
          selected <- intersect(preproc_state$noncomp_cols, nums)
        } else {
          selected <- NULL
        }
        checkboxGroupInput(ns("noncomp_cols"), "Non-compositional columns:", choices = nums, selected = selected)
      } else {
        p(class = "info-box-text", style = "font-style: italic; color: #999;", "No additional numeric columns available.")
      }
    })
    
    # Track noncomp_cols changes
    observeEvent(input$noncomp_cols, {
      preproc_state$noncomp_cols <- input$noncomp_cols
      preproc_state$initialized <- TRUE
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    output$pivot_var_ui <- renderUI({
      req(input$transform_type, input$comp_cols)
      if (!identical(input$transform_type, "ilr")) {
        return(NULL)
      }
      
      # Use stored state if valid, otherwise auto-guess SiO2
      if (preproc_state$initialized && !is.null(preproc_state$pivot_var) && 
          preproc_state$pivot_var %in% input$comp_cols) {
        selected <- preproc_state$pivot_var
      } else {
        guess <- grep("SiO2", input$comp_cols, ignore.case = TRUE, value = TRUE)
        selected <- if (length(guess) > 0) guess[1] else input$comp_cols[1]
      }
      
      div(
        selectInput(ns("pivot_var"), label = div("Select pivot variable:", help_icon("<strong>Chooses which oxide anchors the first ILR pivot coordinate.</strong><details><summary>Learn more</summary>The pivot variable is the oxide placed in the <b>numerator</b> of the first ILR log-ratio. The denominator is the geometric mean of all remaining parts.<br><br>Mathematically: z\u2081 = constant \u00d7 ln(pivot / geometric_mean(others)).<br><br><b>For glass data, use SiO\u2082</b>; this is the convention used by the VARG26 reference model.<br><br><b>Why it matters:</b> Pivot choice defines the transformed coordinate system. The first coordinate captures the pivot oxide's relative abundance versus everything else.<br><br>To compare datasets, or compare with VARG26, use the same pivot consistently.<br><br>If you are only exploring your own data, any abundant analyte can work, but SiO\u2082 is the recommended default.</details>")), choices = input$comp_cols, selected = selected)
      )
    })
    
    # Track pivot_var changes
    observeEvent(input$pivot_var, {
      preproc_state$pivot_var <- input$pivot_var
      preproc_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    output$gmm_data_source_ui <- renderUI({
      req(rv$data)
      choices <- c()
      if (any(grepl("^pivot_", names(rv$data)))) choices <- c(choices, "pivot")
      if (any(grepl("^clr_", names(rv$data)))) choices <- c(choices, "clr")
      if (any(grepl("_imp$", names(rv$data)))) choices <- c(choices, "imp")
      if (any(grepl("_norm$", names(rv$data)))) choices <- c(choices, "norm")
      if (any(grepl("_raw$", names(rv$data)))) choices <- c(choices, "raw")
      req(length(choices) > 0)
      # Preserve user's selection across data updates
      selected <- if (gmm_ui_state$initialized && !is.null(gmm_ui_state$data_source) && gmm_ui_state$data_source %in% choices) gmm_ui_state$data_source else choices[1]
      selectInput(ns("gmm_data_source"), label = div("Data for clustering:", help_icon("<strong>Chooses which data version GMM uses for clustering.</strong><details><summary>Learn more</summary>Pick the data representation used by GMM:<br><br><b>pivot</b> (recommended): ILR pivot coordinates; this is the preferred representation for compositional GMM.<br><b>clr</b>: CLR-transformed coordinates. Easier to interpret, but covariance is singular (columns sum to zero), so GMM may fail for some datasets.<br><b>imp</b>: Imputed compositional data before log-ratio transform; may work, but is less reliable for oxides.<br><b>norm</b>: Normalized data.<br><b>raw</b>: Non-compositional data only; useful only when clustering non-oxide variables.<br><br><b>For standard geochemical clustering, use pivot.</b></details>")), choices = choices, selected = selected)
    })

    # Track GMM data source changes
    observeEvent(input$gmm_data_source, {
      gmm_ui_state$data_source <- input$gmm_data_source
      gmm_ui_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    output$umap_data_source_ui <- renderUI({
      req(rv$data)
      choices <- c()
      if (any(grepl("^pivot_", names(rv$data)))) choices <- c(choices, "pivot")
      if (any(grepl("^clr_", names(rv$data)))) choices <- c(choices, "clr")
      if (any(grepl("_imp$", names(rv$data)))) choices <- c(choices, "imp")
      if (any(grepl("_norm$", names(rv$data)))) choices <- c(choices, "norm")
      if (any(grepl("_raw$", names(rv$data)))) choices <- c(choices, "raw")
      if (length(choices) == 0) choices <- c("raw" = "raw")
      # Preserve user's selection across data updates
      selected <- if (umap_ui_state$initialized && !is.null(umap_ui_state$data_source) && umap_ui_state$data_source %in% choices) umap_ui_state$data_source else choices[1]
      selectInput(ns("umap_data_source"), label = div("Data for UMAP:", help_icon("<strong>Chooses which data version UMAP uses for projection.</strong><details><summary>Learn more</summary>Select the data representation for UMAP:<br><br><b>pivot</b> (recommended): ILR-transformed coordinates; best for compositional distance geometry.<br><b>clr</b>: CLR-transformed coordinates; interpretable and works well for exploratory analysis and visualization, including UMAP.<br><b>imp</b>: Imputed data before log-ratio transform; distances are less reliable for oxides summing to 100%.<br><b>norm</b>: Normalized data.<br><b>raw</b>: Untransformed data.<br><br><b>Default:</b> pivot.<br><br>For consistency with the VARG26 reference model, always use pivot. Switch to clr for exploratory interpretation or to raw/norm only for non-compositional variables.<br><br><b>Note:</b> When using VARG26 UMAP, the correct variables are selected for you automatically, including re-imputation and ILR, even if not done manually beforehand.</details>")), choices = choices, selected = selected)
    })

    # Track UMAP data source changes
    observeEvent(input$umap_data_source, {
      umap_ui_state$data_source <- input$umap_data_source
      umap_ui_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    # --- Advanced Column Selectors ---
    observe({
      req(rv$data)
      # Get all numeric columns, excluding metadata
      nums <- names(rv$data)[sapply(rv$data, is.numeric)]
      nums <- nums[!grepl("^UID$|^gmm_cluster$|^UMAP_", nums)]

      # Update GMM custom columns
      # Preserve selection if possible
      current_gmm <- isolate(input$gmm_custom_cols)
      updateSelectInput(session, "gmm_custom_cols", choices = nums, selected = current_gmm)

      # Update UMAP custom columns
      current_umap <- isolate(input$umap_custom_cols)
      updateSelectInput(session, "umap_custom_cols", choices = nums, selected = current_umap)
    })

    output$semisupervised_ui <- renderUI({
      req(rv$data)
      # Exclude UID and transform columns, keep potential categorical labels or numeric targets
      potential_cols <- names(rv$data)[!grepl("^UID$|_orig$|_norm$|_imp$|^pivot_|^clr_|^UMAP_|^gmm_cluster$", names(rv$data))]
      # Preserve user's selection
      selected <- if (umap_ui_state$initialized && !is.null(umap_ui_state$y_vars)) {
        intersect(umap_ui_state$y_vars, potential_cols)
      } else {
        NULL
      }
      selectInput(ns("y_vars"), "Select label (Y) columns:", choices = potential_cols, selected = selected, multiple = TRUE)
    })

    # Track semi-supervised Y column selection
    observeEvent(input$y_vars, {
      umap_ui_state$y_vars <- input$y_vars
      umap_ui_state$initialized <- TRUE
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    output$popsel_xvar_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric) & !names(rv$data) %in% "UID"]
      # Use stored state if initialized and still valid, otherwise use default
      if (popsel_state$initialized && !is.null(popsel_state$xvar) && popsel_state$xvar %in% nums) {
        selected <- popsel_state$xvar
      } else {
        # Prefer 2D UMAP columns (both axes should be 2D)
        umap_2d <- grep("^UMAP_.*2D_1$", nums, value = TRUE)
        selected <- if (length(umap_2d) > 0) umap_2d[1] else if (length(nums) > 0) nums[1] else NULL
      }
      selectInput(ns("popsel_xvar"), label = div("X variable:", help_icon("<strong>Chooses the X-axis variable for population selection plots.</strong><details><summary>Learn more</summary>Select the variable shown on the X-axis.<br><br><b>Tip:</b> Start with UMAP 2D coordinates for selection, then switch to oxide pairs (for example, SiO\u2082 vs CaO) to validate assignments in compositional space.<br><br>When <b>gmm_cluster</b> is selected, points are jittered to reveal overlap.</details>")), choices = nums, selected = selected)
    })

    output$popsel_yvar_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric) & !names(rv$data) %in% "UID"]
      # Use stored state if initialized and still valid, otherwise use default
      if (popsel_state$initialized && !is.null(popsel_state$yvar) && popsel_state$yvar %in% nums) {
        selected <- popsel_state$yvar
      } else {
        # Prefer 2D UMAP columns (second axis)
        umap_2d <- grep("^UMAP_.*2D_2$", nums, value = TRUE)
        selected <- if (length(umap_2d) > 0) umap_2d[1] else if (length(nums) > 1) nums[2] else if (length(nums) > 0) nums[1] else NULL
      }
      selectInput(ns("popsel_yvar"), label = div("Y variable:", help_icon("<strong>Chooses the Y-axis variable for population selection plots.</strong><details><summary>Learn more</summary>Select the variable shown on the Y-axis.<br><br>Use the same variable space as X (for example, both UMAP or both oxides) for meaningful plots.<br><br>Mixing spaces (for example, UMAP on X and an oxide on Y) can still be informative, but distances are not interpretable.</details>")), choices = nums, selected = selected)
    })

    # Hover label column selector for Population Definition plot
    output$popsel_hover_label_ui <- renderUI({
      req(rv$data)
      all_cols <- names(rv$data)
      # Exclude UID (always shown) and internal columns
      label_cols <- setdiff(all_cols, c("UID"))
      # Use stored state if initialized and still valid, otherwise use default
      if (popsel_state$initialized && !is.null(popsel_state$hover_label) && popsel_state$hover_label %in% c("", label_cols)) {
        selected <- popsel_state$hover_label
      } else {
        # Default to sample_point or sample_id if available
        default <- intersect(c("sample_point", "sample_id", "Samplepop"), label_cols)
        selected <- if (length(default) > 0) default[1] else NULL
      }
      selectInput(ns("popsel_hover_label"),
        label = div("Hover Label:", help_icon("<strong>Chooses which column is shown in plot hover tooltips.</strong><details><summary>Learn more</summary>Select a column to display when hovering over points.<br><br>This helps identify which samples belong to a cluster by showing metadata such as sample names, volcanic sources, or site ID directly in the tooltip.<br><br><b>Tip:</b> Sample ID or volcanic source columns usually provide the most useful hover labels.</details>")),
        choices = c("(none)" = "", label_cols),
        selected = selected
      )
    })

    # --- Initialize Populations from Column ---
    output$init_pop_col_ui <- renderUI({
      req(rv$data)
      # Show text/character columns + gmm_cluster (integer but categorical)
      all_cols <- names(rv$data)
      # Include character columns, factor columns, and gmm_cluster
      candidate_cols <- all_cols[sapply(rv$data, function(x) is.character(x) || is.factor(x))]
      # Also include gmm_cluster if it exists (stored as integer)
      if ("gmm_cluster" %in% all_cols && !"gmm_cluster" %in% candidate_cols) {
        candidate_cols <- c("gmm_cluster", candidate_cols)
      }
      # Exclude internal columns
      candidate_cols <- setdiff(candidate_cols, c("UID", "population", ".hover_text"))
      req(length(candidate_cols) > 0)
      selectInput(ns("init_pop_col"), "Column to use:", choices = candidate_cols, selected = candidate_cols[1])
    })

    observeEvent(input$init_pop_from_col, {
      req(rv$data, input$init_pop_col)
      col <- input$init_pop_col
      req(col %in% names(rv$data))

      # Copy column values as population labels
      values <- as.character(rv$data[[col]])
      # Replace NA, empty, or whitespace-only values with "Unassigned"
      values[is.na(values) | trimws(values) == ""] <- "Unassigned"
      rv$data$population <- values
      bump_data_generation()

      # Rebuild styles for all the new population names
      rebuilt_styles <- rebuild_pop_styles(rv$data$population)
      pop_styles(rebuilt_styles)

      # Update the selectize dropdown with the new labels
      existing_labels <- unique(rv$data$population)
      updateSelectizeInput(session, "feature_label", choices = existing_labels, selected = "", server = TRUE)

      n_pops <- length(setdiff(unique(values), "Unassigned"))
      showNotification(
        sprintf("Initialized %d populations from '%s'.", n_pops, col),
        type = "message"
      )

      output$init_pop_status <- renderUI({
        tags$p(class = "small text-success mt-1", icon("check"),
          sprintf("%d populations created from '%s'", n_pops, col))
      })
    })

    # Track popsel variable changes and store them
    observeEvent(input$popsel_xvar, {
      popsel_state$xvar <- input$popsel_xvar
      popsel_state$initialized <- TRUE
    }, ignoreInit = TRUE)
    
    observeEvent(input$popsel_yvar, {
      popsel_state$yvar <- input$popsel_yvar
      popsel_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    observeEvent(input$popsel_hover_label, {
      popsel_state$hover_label <- input$popsel_hover_label
      popsel_state$initialized <- TRUE
    }, ignoreInit = TRUE)

    # --- Derived Variable Creator ---
    output$derived_var_a_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric) & !names(rv$data) %in% "UID"]
      selectInput(ns("derived_var_a"), "Variable A:", choices = nums, selected = nums[1], width = "100%")
    })
    
    output$derived_var_b_ui <- renderUI({
      req(rv$data)
      nums <- names(rv$data)[sapply(rv$data, is.numeric) & !names(rv$data) %in% "UID"]
      selectInput(ns("derived_var_b"), "Variable B:", choices = nums, selected = if (length(nums) > 1) nums[2] else nums[1], width = "100%")
    })
    
    observeEvent(input$create_derived_var, {
      req(rv$data, input$derived_var_a, input$derived_var_b, input$derived_op)
      
      var_a <- input$derived_var_a
      var_b <- input$derived_var_b
      op <- input$derived_op
      
      if (var_a == var_b) {
        output$derived_var_status <- renderUI({
          div(class = "text-danger small mt-1", icon("exclamation-circle"), "Variables A and B must be different.")
        })
        return()
      }
      
      if (!var_a %in% names(rv$data) || !var_b %in% names(rv$data)) {
        output$derived_var_status <- renderUI({
          div(class = "text-danger small mt-1", icon("exclamation-circle"), "Selected variable(s) not found in data.")
        })
        return()
      }
      
      # Build column name and compute
      if (op == "ratio") {
        col_name <- paste0(var_a, "_div_", var_b)
        values <- rv$data[[var_a]] / rv$data[[var_b]]
        # Handle division by zero
        values[is.infinite(values)] <- NA
      } else {
        col_name <- paste0(var_a, "_plus_", var_b)
        values <- rv$data[[var_a]] + rv$data[[var_b]]
      }
      
      # Check if column already exists
      if (col_name %in% names(rv$data)) {
        output$derived_var_status <- renderUI({
          div(class = "text-warning small mt-1", icon("info-circle"), paste0("'", col_name, "' already exists. Values updated."))
        })
      } else {
        output$derived_var_status <- renderUI({
          div(class = "text-success small mt-1", icon("check-circle"), paste0("Created '", col_name, "'. Select it from X or Y above."))
        })
      }
      
      rv$data[[col_name]] <- values
      bump_data_generation()
      showNotification(paste0("Derived variable '", col_name, "' is now available."), type = "message", duration = 4)
    })

    # Filter UI for sample/group filtering
    output$filter_col_ui <- renderUI({
      req(rv$data)
      char_cols <- names(rv$data)[sapply(rv$data, function(x) is.character(x) || is.factor(x))]
      char_cols <- char_cols[!char_cols %in% c("population")]
      if (length(char_cols) == 0) {
        return(NULL)
      }

      # Use stored state if initialized, otherwise default to ""
      selected <- if (filter_state$initialized) filter_state$col else ""
      selectizeInput(ns("filter_col"), label = div("Filter column:", help_icon("<strong>Filters the plot to one category so assignments are easier.</strong><details><summary>Learn more</summary>Optionally filters the plot to show only points that match a selected category.<br><br>Click the \u00d7 to clear the filter and show all data again.<br><br><b>When to use:</b> Focus on one sample or site at a time to assign populations without visual clutter from other groups.</details>")), choices = c("None" = "", char_cols), selected = selected, options = list(plugins = list('clear_button')))
    })

    output$filter_val_ui <- renderUI({
      req(input$filter_col)
      if (input$filter_col == "") {
        return(NULL)
      }
      req(rv$data)
      vals <- sort(unique(as.character(rv$data[[input$filter_col]])))

      # Use stored state if available
      selected <- if (filter_state$initialized && filter_state$val %in% c("All", vals)) filter_state$val else "All"
      selectizeInput(ns("filter_val"), "Show only:", choices = c("All", vals), selected = selected, options = list(plugins = list('clear_button')))
    })

    # Track filter changes and store them
    observeEvent(input$filter_col,
      {
        filter_state$col <- input$filter_col
        filter_state$initialized <- TRUE
      },
      ignoreInit = TRUE
    )

    observeEvent(input$filter_val,
      {
        filter_state$val <- input$filter_val
      },
      ignoreInit = TRUE
    )

    # --- Stale State Observers ---
    observeEvent(c(input$comp_cols, input$noncomp_cols, input$do_impute, input$impute_method, input$transform_type, input$pivot_var),
      {
        # Don't mark as stale if we're restoring from a loaded project
        if (!restoring_state()) {
          rv$data_stale <- TRUE
        }
      },
      ignoreInit = TRUE
    )

    # --- Data Loading ---
    observeEvent(input$load_data, {
      req(input$file_upload)
      show_loading("Reading file...")
      on.exit(hide_loading())

      tryCatch(
        {
          file_ext <- tolower(tools::file_ext(input$file_upload$name))
          df <- if (file_ext == "xlsx") {
            readxl::read_xlsx(input$file_upload$datapath, sheet = input$sheet_select)
          } else if (file_ext == "xls") {
            readxl::read_xls(input$file_upload$datapath, sheet = input$sheet_select)
          } else {
            read.csv(input$file_upload$datapath, stringsAsFactors = FALSE, check.names = FALSE)
          }

          # Normalize column names (e.g. FeOT -> FeO)
          df <- normalize_geochem_data(df)

          if (!"UID" %in% names(df)) df <- df %>% mutate(UID = row_number(), .before = 1)
          
          # Normalize population column name (handle both "Population" and "population")
          # Check for capitalized version first and rename if found
          if ("Population" %in% names(df) && !"population" %in% names(df)) {
            names(df)[names(df) == "Population"] <- "population"
          }
          if (!"population" %in% names(df)) df$population <- "Unassigned"

          cancel_processing_background_jobs()
          clear_user_umap_model()
          rv$mclust_result <- NULL
          rv$umap_mode_ran <- NULL
          rv$data_stale <- FALSE
          rv$pipeline_config <- NULL
          rv$original_cols <- names(df)
          rv$source_filename <- input$file_upload$name
          rv$data <- df
          bump_data_generation()
          undo_history(list())
          showNotification("Data loaded.", type = "message")
        },
        error = function(e) {
          showNotification(
            paste0(
              "Could not load the uploaded file. Please upload a CSV, XLSX, or XLS file with column headers. ",
              "Technical detail: ", e$message
            ),
            type = "error"
          )
        }
      )
    })
    
    # --- Row Filtering ---
    observeEvent(input$apply_row_filter, {
      req(rv$data)
      
      df <- rv$data
      original_nrow <- nrow(df)
      rows_to_exclude <- rep(FALSE, nrow(df))
      
      # 1. Exclude blank rows if checkbox is checked
      if (isTRUE(input$exclude_blank_rows)) {
        non_analyte_cols <- c("UID", "Depth", "depth", "Age", "age", "Latitude", "Longitude",
                              "latitude", "longitude", "lat", "lon", "Elevation", "elevation",
                              "row_number", "X", "x", "Y", "y", "Z", "z")
        numeric_cols <- names(df)[sapply(df, is.numeric) & !names(df) %in% c("UID", non_analyte_cols)]
        if (length(numeric_cols) > 0) {
          blank_rows <- apply(df[, numeric_cols, drop = FALSE], 1, function(x) all(is.na(x)))
          rows_to_exclude <- rows_to_exclude | blank_rows
        }
      }
      
      # 2. Exclude rows matching text pattern
      if (!is.null(input$filter_text_col) && input$filter_text_col != "" &&
          !is.null(input$filter_text_pattern) && nchar(trimws(input$filter_text_pattern)) > 0) {
        tryCatch({
          pattern <- trimws(input$filter_text_pattern)
          col_data <- as.character(df[[input$filter_text_col]])
          # Split on | and use fixed matching to avoid regex errors
          # Case-insensitive via tolower (fixed=TRUE ignores ignore.case arg)
          parts <- strsplit(pattern, "\\|")[[1]]
          parts <- trimws(parts)
          parts <- parts[nchar(parts) > 0]
          if (length(parts) > 0) {
            col_lower <- tolower(col_data)
            pattern_matches <- Reduce(`|`, lapply(parts, function(p) {
              grepl(tolower(p), col_lower, fixed = TRUE)
            }))
            rows_to_exclude <- rows_to_exclude | pattern_matches
          }
        }, error = function(e) {
          showNotification(paste("Invalid filter pattern:", e$message), type = "warning")
        })
      }
      
      n_excluded <- sum(rows_to_exclude)
      
      if (n_excluded > 0) {
        # Mark excluded rows - set a flag column for exclusion tracking
        # The excluded rows will have NA in derived columns but keep original data
        df$row_excluded <- rows_to_exclude
        
        rv$data <- df
        bump_data_generation()
        showNotification(
          sprintf("Marked %d row%s for exclusion. These will be skipped during processing but preserved in output.", 
                  n_excluded, if (n_excluded > 1) "s" else ""),
          type = "message",
          duration = 5
        )
      } else {
        # Remove exclusion flag if no rows to exclude
        df$row_excluded <- FALSE
        rv$data <- df
        bump_data_generation()
        showNotification("No rows to filter.", type = "message")
      }
    })

    # Data summary for Step 1 - shows column types and counts
    output$data_summary_ui <- renderUI({
      req(rv$data)

      # Create column type summary
      col_info <- data.frame(
        Column = names(rv$data),
        Type = sapply(rv$data, function(x) {
          if (is.numeric(x)) {
            "Numeric"
          } else if (is.factor(x)) {
            "Factor"
          } else if (is.character(x)) {
            "Text"
          } else if (inherits(x, "Date")) {
            "Date"
          } else {
            class(x)[1]
          }
        }),
        `Non-NA` = sapply(rv$data, function(x) sum(!is.na(x))),
        Example = sapply(rv$data, function(x) {
          vals <- na.omit(x)[1:3]
          if (length(vals) == 0) {
            return("(all NA)")
          }
          paste(head(vals, 2), collapse = ", ")
        }),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      n_numeric <- sum(col_info$Type == "Numeric")
      n_text <- sum(col_info$Type %in% c("Text", "Factor"))

      div(
        div(
          class = "alert alert-success py-2 mb-2",
          icon("check-circle"), " ",
          strong(nrow(rv$data)), " rows × ", strong(ncol(rv$data)), " columns loaded"
        ),
        div(
          class = "small text-muted mb-2",
          icon("hashtag"), " ", n_numeric, " numeric | ",
          icon("font"), " ", n_text, " text/factor"
        ),
        div(
          style = "max-height: 200px; overflow-y: auto;",
          tags$table(
            class = "table table-sm table-striped",
            tags$thead(
              tags$tr(
                tags$th("Column", style = "width: 35%"),
                tags$th("Type", style = "width: 20%"),
                tags$th("Non-NA", style = "width: 15%"),
                tags$th("Example", style = "width: 30%")
              )
            ),
            tags$tbody(
              lapply(seq_len(nrow(col_info)), function(i) {
                tags$tr(
                  tags$td(col_info$Column[i], style = "font-family: monospace; font-size: 0.85em;"),
                  tags$td(
                    span(
                      class = if (col_info$Type[i] == "Numeric") "badge bg-primary" else "badge bg-secondary",
                      col_info$Type[i]
                    )
                  ),
                  tags$td(col_info$`Non-NA`[i]),
                  tags$td(col_info$Example[i], style = "font-size: 0.85em; color: #666;")
                )
              })
            )
          )
        )
      )
    })

    # Data preview for Step 1
    output$preview_data_step1 <- renderDT({
      req(rv$data)
      datatable(rv$data, options = list(scrollX = TRUE, scrollY = "250px", pageLength = 5), rownames = FALSE)
    })

    # Data preview for Step 2
    output$preview_data_step2 <- renderDT({
      req(rv$data)
      datatable(rv$data, options = list(scrollX = TRUE, scrollY = "250px", pageLength = 5), rownames = FALSE)
    })

    # Column summary - shows what columns were generated
    output$column_summary_ui <- renderUI({
      req(rv$data)

      # Track original columns when first loaded
      if (is.null(rv$original_cols)) {
        return(p("Load data and apply transforms to see generated columns.", class = "text-muted", style = "font-style: italic;"))
      }

      all_cols <- names(rv$data)
      new_cols <- setdiff(all_cols, rv$original_cols)

      if (length(new_cols) == 0) {
        return(div(
          icon("info-circle", class = "text-info"),
          span(" No new columns yet. Click 'Apply Transforms' to process data.", style = "margin-left: 5px;")
        ))
      }

      # Group columns by type
      norm_cols <- grep("_norm$", new_cols, value = TRUE)
      imp_cols <- grep("_imp$", new_cols, value = TRUE)
      raw_cols <- grep("_raw$", new_cols, value = TRUE)
      user_pivot_cols <- grep("^pivot_(?!VARG26_)", new_cols, value = TRUE, perl = TRUE)
      user_clr_cols <- grep("^clr_", new_cols, value = TRUE)
      varg26_pivot_cols <- grep("^pivot_VARG26_", new_cols, value = TRUE)

      tagList(
        div(
          style = "margin-bottom: 10px;",
          icon("check-circle", class = "text-success"),
          strong(paste(" Generated", length(new_cols), "new columns"), style = "margin-left: 5px;")
        ),
        tags$hr(style = "margin: 10px 0;"),
        if (length(norm_cols) > 0) {
          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("balance-scale"), paste0(" Normalized (", length(norm_cols), ")"), style = "font-size: 0.9rem;"),
            tags$div(
              style = "margin-left: 20px; font-size: 0.85rem; color: #666;",
              paste(sub("_norm$", "", norm_cols), collapse = ", ")
            )
          )
        },
        if (length(imp_cols) > 0) {
          # Build per-column imputation breakdown
          imputed_details <- tryCatch({
            detail_lines <- character(0)
            total_filled <- 0
            for (ic in imp_cols) {
              base <- sub("_imp$", "", ic)
              norm_cn <- paste0(base, "_norm")
              if (norm_cn %in% names(rv$data) && ic %in% names(rv$data)) {
                # Count values that were NA/zero in norm but got filled in imp
                n_filled <- sum(
                  !is.na(rv$data[[ic]]) &
                  (is.na(rv$data[[norm_cn]]) | rv$data[[norm_cn]] <= 1e-10),
                  na.rm = TRUE
                )
                if (n_filled > 0) {
                  detail_lines <- c(detail_lines, paste0(base, ": ", n_filled, " values"))
                  total_filled <- total_filled + n_filled
                }
              }
            }
            if (length(detail_lines) > 0) {
              list(
                summary = paste0(length(detail_lines), " column(s) needed imputation (", total_filled, " total values filled)"),
                details = paste(detail_lines, collapse = " | ")
              )
            } else {
              list(summary = "No missing values detected. All columns were complete", details = NULL)
            }
          }, error = function(e) {
            list(summary = paste(sub("_imp$", "", imp_cols), collapse = ", "), details = NULL)
          })

          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("fill-drip"), paste0(" Imputed (", length(imp_cols), ")"), style = "font-size: 0.9rem;"),
            tags$div(
              style = "margin-left: 20px; font-size: 0.85rem; color: #666;",
              imputed_details$summary
            ),
            if (!is.null(imputed_details$details)) {
              tags$div(
                style = "margin-left: 20px; font-size: 0.8rem; color: #999; font-style: italic;",
                imputed_details$details
              )
            }
          )
        },
        if (length(user_pivot_cols) > 0) {
          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("project-diagram"), paste0(" Pivot Coordinates (", length(user_pivot_cols), ")"), style = "font-size: 0.9rem;"),
            tags$div(
              style = "margin-left: 20px; font-size: 0.85rem; color: #666;",
              paste0(length(user_pivot_cols), " ILR coordinates from ", length(norm_cols), " compositional variables")
            )
          )
        },
        if (length(user_clr_cols) > 0) {
          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("project-diagram"), paste0(" CLR Coordinates (", length(user_clr_cols), ")"), style = "font-size: 0.9rem;"),
            tags$div(
              style = "margin-left: 20px; font-size: 0.85rem; color: #666;",
              paste0(length(user_clr_cols), " CLR coordinates from ", length(norm_cols), " compositional variables")
            )
          )
        },
        if (length(varg26_pivot_cols) > 0) {
          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("project-diagram"), paste0(" VARG26 Pivot Coordinates (", length(varg26_pivot_cols), ")"), style = "font-size: 0.9rem;")
          )
        },
        if (length(raw_cols) > 0) {
          div(
            style = "margin-bottom: 8px;",
            tags$strong(icon("table"), paste0(" Raw Non-compositional (", length(raw_cols), ")"), style = "font-size: 0.9rem;"),
            tags$div(
              style = "margin-left: 20px; font-size: 0.85rem; color: #666;",
              paste(sub("_raw$", "", raw_cols), collapse = ", ")
            )
          )
        }
      )
    })

    # --- Processing (background with cancel support) ---
    observeEvent(input$apply_transforms, {
      req(rv$data)
      if (is.null(input$comp_cols) && is.null(input$noncomp_cols)) {
        showNotification("Please select at least one compositional or non-compositional column.", type = "warning")
        return()
      }
      if (!begin_heavy_job("preprocessing")) return()

      # Capture all parameters before launching background
      df <- rv$data
      comp_cols <- input$comp_cols
      noncomp_cols <- input$noncomp_cols
      do_impute <- input$do_impute
      impute_method <- if (isTRUE(do_impute)) (input$impute_method %||% "auto") else "auto"
      transform_type <- input$transform_type %||% "ilr"
      do_pivot <- identical(transform_type, "ilr")
      do_clr <- identical(transform_type, "clr")
      pivot_var <- if (do_pivot) input$pivot_var else NULL

      # Remove existing transform columns before passing to background
      existing_transform_cols <- grep("_norm$|_imp$|_raw$|^pivot_(?!VARG26_)|^clr_|^gmm_cluster$|^UMAP_", names(df), value = TRUE, perl = TRUE)
      if (length(existing_transform_cols) > 0) {
        df <- df[, !names(df) %in% existing_transform_cols, drop = FALSE]
      }

      # Update original_cols to exclude any transform columns from imported data
      # so the readout correctly shows ALL generated columns as "new"
      if (!is.null(rv$original_cols)) {
        rv$original_cols <- rv$original_cols[!grepl("_norm$|_imp$|_raw$|^pivot_(?!VARG26_)|^clr_|^gmm_cluster$|^UMAP_", rv$original_cols, perl = TRUE)]
      }

      # Clear downstream reactive values immediately
      rv$mclust_result <- NULL
      rv$umap_mode_ran <- NULL
      clear_user_umap_model()

      # Pre-validate: check for all-zero/blank compositional columns before processing
      if (!is.null(comp_cols) && length(comp_cols) > 0 && (do_impute || do_pivot || do_clr)) {
        rows_check <- if ("row_excluded" %in% names(df)) !df$row_excluded else rep(TRUE, nrow(df))
        zero_cols <- character(0)
        for (col in comp_cols) {
          vals <- df[[col]][rows_check]
          if (all(is.na(vals) | vals <= 1e-10)) {
            zero_cols <- c(zero_cols, col)
          }
        }
        if (length(zero_cols) > 0) {
          showNotification(
            HTML(paste0(
              "<strong>Cannot proceed:</strong> The following compositional column(s) have no non-zero values: <em>",
              paste(zero_cols, collapse = ", "),
              "</em>.<br>Please remove them from the Compositional Columns selector before applying transforms."
            )),
            type = "error", duration = 15
          )
          finish_heavy_job()
          return()
        }
      }

      # Define the preprocessing function (shared by both paths)
        preprocess_func <- function(df, comp_cols, noncomp_cols, do_impute, transform_type, pivot_var, impute_method = "auto") {
          imputation_method <- if (isTRUE(do_impute)) impute_method else NULL
          imputation_fallback <- FALSE
          imputation_auto_reason <- NULL
          imputation_floor_msg <- NULL
          transform_skip_msg <- NULL

          # Process compositional columns
          if (!is.null(comp_cols) && length(comp_cols) > 0) {
            rows_to_process <- if ("row_excluded" %in% names(df)) !df$row_excluded else rep(TRUE, nrow(df))

            # 1. Normalize
            df_norm <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(comp_cols)))
            names(df_norm) <- paste0(comp_cols, "_norm")
            if (sum(rows_to_process) > 0) {
              df_subset <- df[rows_to_process, comp_cols, drop = FALSE]
              rs <- rowSums(df_subset, na.rm = TRUE)
              df_subset_norm <- df_subset / rs * 100
              df_norm[rows_to_process, ] <- df_subset_norm
            }
            df <- dplyr::bind_cols(df, df_norm)

            # 2. Impute (on RAW data, not normalized — lets impCoda handle closure properly)
            if (do_impute) {
              prep <- df[, comp_cols, drop = FALSE]
              if (sum(rows_to_process) > 0) {
                # Only work with non-excluded rows
                prep_sub <- prep[rows_to_process, , drop = FALSE]
              } else {
                prep_sub <- prep
              }
              prep_sub[prep_sub <= 1e-10] <- NA

              # Exclude columns that are entirely NA (e.g. all-zero columns like P2O5)
              col_all_na <- apply(prep_sub, 2, function(x) all(is.na(x)))
              skipped_cols <- names(prep_sub)[col_all_na]
              impute_cols <- names(prep_sub)[!col_all_na]

              valid_for_impute <- rows_to_process & !apply(prep[, impute_cols, drop = FALSE], 1, function(x) all(is.na(x)))

              if (length(impute_cols) >= 2 && sum(valid_for_impute) >= 2) {
                prep_valid <- prep_sub[valid_for_impute[rows_to_process], impute_cols, drop = FALSE]
                set.seed(42)

                # --- Assess column missingness ---
                col_miss_frac <- colMeans(is.na(prep_valid))
                max_miss <- max(col_miss_frac)
                max_miss_col <- names(which.max(col_miss_frac))
                high_missingness <- max_miss > 0.50

                # --- Resolve method ---
                # "auto" switches intelligently; explicit user choices are respected
                resolved_method <- impute_method
                if (impute_method == "auto") {
                  if (high_missingness) {
                    resolved_method <- "lm"
                    imputation_fallback <- TRUE
                    imputation_auto_reason <- paste0(
                      "Auto selected standard regression because ",
                      max_miss_col, " has ", round(max_miss * 100, 1),
                      "% missing values. Robust regression can produce ",
                      "unreliable estimates when >50% of a column is missing."
                    )
                    message("Auto imputation: ", imputation_auto_reason)
                  } else {
                    resolved_method <- "ltsReg"
                    message("Auto imputation: using robust regression (max column missingness = ",
                            round(max_miss * 100, 1), "% in ", max_miss_col, ")")
                  }
                } else if (impute_method == "ltsReg" && high_missingness) {
                  # User explicitly chose Robust — warn but respect their choice
                  imputation_auto_reason <- paste0(
                    "<u>Warning:</u> ", max_miss_col, " has ", round(max_miss * 100, 1),
                    "% missing values. Robust regression may produce unreliable ",
                    "estimates when >50% of a column is missing. ",
                    "<em>Consider switching to Auto or Standard if results look unexpected.</em>"
                  )
                  imputation_fallback <- TRUE
                  message("Imputation warning: ", imputation_auto_reason)
                }

                # Run imputation with tryCatch fallback
                used_method <- resolved_method
                res <- tryCatch(
                  robCompositions::impCoda(prep_valid, method = resolved_method, maxit = 5)$xImp,
                  error = function(e) {
                    fallback <- if (resolved_method == "ltsReg") "lm" else "ltsReg"
                    message("impCoda ", resolved_method, " failed (", e$message,
                            "), falling back to method='", fallback, "'")
                    used_method <<- fallback
                    imputation_fallback <<- TRUE
                    if (is.null(imputation_auto_reason)) {
                      imputation_auto_reason <<- paste0(
                        "Robust regression failed: ", e$message,
                        ". Switched to standard regression."
                      )
                    }
                    robCompositions::impCoda(prep_valid, method = fallback, maxit = 5)$xImp
                  }
                )

                # --- Post-check: validate imputed values (ltsReg only) ---
                # Only relevant when robust regression was used — checks for
                # biased imputation caused by high-missingness KNN contamination.
                # When lm is used (either by Auto fallback or user choice), skip.
                #
                # Compare within `res` only (impCoda may rescale internally),
                # using non-missing rows as the reference baseline on the same scale.
                if (used_method == "ltsReg") {
                  post_check_warnings <- character(0)
                  for (cn in impute_cols) {
                    orig_vals <- prep_valid[[cn]]
                    imp_vals <- res[[cn]]
                    non_missing <- !is.na(orig_vals) & orig_vals > 1e-10
                    was_missing <- is.na(orig_vals) | orig_vals <= 1e-10

                    if (sum(non_missing) >= 5 && sum(was_missing) >= 1) {
                      # Use res values for BOTH medians so they're on the same scale
                      med_ref <- median(imp_vals[non_missing])
                      med_imp <- median(imp_vals[was_missing])

                      if (med_ref > 0 && (med_imp / med_ref > 10 || med_imp / med_ref < 0.1)) {
                        warn_msg <- paste0(cn, " imputed median deviates >10x from ",
                                           "measured values (ratio: ",
                                           round(med_imp / med_ref, 2), ")")
                        message("Post-check: ", warn_msg)

                        if (impute_method == "auto") {
                          # Auto mode: re-run with standard regression
                          used_method <- "lm"
                          imputation_fallback <- TRUE
                          imputation_auto_reason <- paste0(
                            "Post-check detected unreliable robust imputation for ", cn,
                            ". Re-ran with standard regression."
                          )
                          set.seed(42)
                          res <- robCompositions::impCoda(prep_valid, method = "lm", maxit = 5)$xImp
                          break
                        } else {
                          # User explicitly chose Robust — warn but don't override
                          post_check_warnings <- c(post_check_warnings, warn_msg)
                        }
                      }
                    }
                  }

                  # Attach post-check warnings for explicit Robust users
                  if (length(post_check_warnings) > 0) {
                    pc_msg <- paste0(
                      "<u>Post-check warning:</u> ", paste(post_check_warnings, collapse = "; "),
                      ". <em>Consider switching to Auto or Standard.</em>"
                    )
                    imputation_fallback <- TRUE
                    imputation_auto_reason <- if (is.null(imputation_auto_reason)) {
                      pc_msg
                    } else {
                      paste0(imputation_auto_reason, "<br><br>", pc_msg)
                    }
                  }
                }

                imputation_method <- used_method
                # Build _imp columns: imputed columns get results, skipped columns get NA
                res_df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(comp_cols)))
                names(res_df) <- paste0(comp_cols, "_imp")
                # Map imputed results back to their corresponding _imp columns
                for (cn in impute_cols) {
                  imp_col <- paste0(cn, "_imp")
                  if (imp_col %in% names(res_df) && cn %in% names(res)) {
                    res_df[[imp_col]][valid_for_impute] <- res[[cn]]
                  }
                }
                # Normalize imputed values to sum to 100 (same scale as _norm columns)
                imp_numeric_cols <- paste0(comp_cols, "_imp")
                for (i in which(valid_for_impute)) {
                  row_vals <- as.numeric(res_df[i, imp_numeric_cols])
                  rs <- sum(row_vals, na.rm = TRUE)
                  if (!is.na(rs) && rs > 0) {
                    res_df[i, imp_numeric_cols] <- row_vals / rs * 100
                  }
                }
                # Floor near-zero imputed values to prevent log-ratio failures
                floored_cols <- character(0)
                floored_count <- 0L
                for (cn in imp_numeric_cols) {
                  vals <- res_df[[cn]]
                  hit <- !is.na(vals) & vals > 0 & vals < 1e-4
                  n_hit <- sum(hit, na.rm = TRUE)
                  if (n_hit > 0) {
                    floored_cols <- c(floored_cols, sub("_imp$", "", cn))
                    floored_count <- floored_count + n_hit
                    res_df[[cn]] <- ifelse(hit, 1e-4, vals)
                  }
                }
                if (floored_count > 0) {
                  imputation_floor_msg <- paste0(
                    floored_count, " near-zero imputed value",
                    if (floored_count > 1) "s" else "",
                    " floored to 0.0001 wt% in ",
                    paste(floored_cols, collapse = ", "),
                    " to prevent log-ratio failures."
                  )
                  message("Imputation floor: ", imputation_floor_msg)
                }
                df <- dplyr::bind_cols(df, res_df)

                if (length(skipped_cols) > 0) {
                  message("Imputation skipped for all-zero/NA columns: ", paste(skipped_cols, collapse = ", "))
                }
              } else {
                # Not enough columns or rows for imputation — create empty _imp columns
                res_df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(comp_cols)))
                names(res_df) <- paste0(comp_cols, "_imp")
                df <- dplyr::bind_cols(df, res_df)
              }
            }

            # 3. Log-ratio transform
            if (identical(transform_type, "ilr")) {
              suffix <- if (do_impute) "_imp" else "_norm"
              p_cols <- paste0(comp_cols, suffix)
              prep <- df[, p_cols, drop = FALSE]
              valid <- rows_to_process & complete.cases(prep) & apply(prep, 1, function(r) all(r > 1e-10))
              if (sum(valid) >= 2) {
                idx <- match(pivot_var, comp_cols)
                res <- robCompositions::pivotCoord(as.matrix(prep[valid, ]), pivot = idx)
                res_df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(res)))
                names(res_df) <- paste0("pivot_", colnames(res))
                res_df[valid, ] <- res
                df <- dplyr::bind_cols(df, res_df)
              }
              # Track rows that couldn't be transformed
              n_skipped <- sum(rows_to_process) - sum(valid)
              if (n_skipped > 0 && !do_impute) {
                # Find which columns had zeros/NAs causing the skip
                problem_cols <- character(0)
                for (pc in p_cols) {
                  vals <- df[[pc]][rows_to_process & !valid]
                  if (any(is.na(vals) | vals <= 1e-10, na.rm = TRUE)) {
                    problem_cols <- c(problem_cols, sub("_norm$", "", pc))
                  }
                }
                transform_skip_msg <- paste0(
                  n_skipped, " row", if (n_skipped > 1) "s" else "",
                  " did not receive ILR pivot coordinates because ",
                  if (length(problem_cols) > 0)
                    paste0(paste(problem_cols, collapse = ", "), " contained")
                  else "compositional columns contained",
                  " zeros or missing values. ",
                  "Consider enabling imputation to fill these gaps."
                )
              }
            } else if (identical(transform_type, "clr")) {
              suffix <- if (do_impute) "_imp" else "_norm"
              p_cols <- paste0(comp_cols, suffix)
              prep <- df[, p_cols, drop = FALSE]
              valid <- rows_to_process & complete.cases(prep) & apply(prep, 1, function(r) all(r > 1e-10))
              if (sum(valid) >= 1) {
                prep_mat <- as.matrix(prep[valid, , drop = FALSE])
                log_mat <- log(prep_mat)
                clr_mat <- log_mat - rowMeans(log_mat)
                clr_df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(clr_mat)))
                names(clr_df) <- paste0("clr_", comp_cols)
                clr_df[valid, ] <- clr_mat
                df <- dplyr::bind_cols(df, clr_df)
              }
              # Track rows that couldn't be transformed
              n_skipped <- sum(rows_to_process) - sum(valid)
              if (n_skipped > 0 && !do_impute) {
                problem_cols <- character(0)
                for (pc in p_cols) {
                  vals <- df[[pc]][rows_to_process & !valid]
                  if (any(is.na(vals) | vals <= 1e-10, na.rm = TRUE)) {
                    problem_cols <- c(problem_cols, sub("_norm$", "", pc))
                  }
                }
                transform_skip_msg <- paste0(
                  n_skipped, " row", if (n_skipped > 1) "s" else "",
                  " did not receive CLR coordinates because ",
                  if (length(problem_cols) > 0)
                    paste0(paste(problem_cols, collapse = ", "), " contained")
                  else "compositional columns contained",
                  " zeros or missing values. ",
                  "Consider enabling imputation to fill these gaps."
                )
              }
            }
          }

          # Process non-compositional columns
          if (!is.null(noncomp_cols) && length(noncomp_cols) > 0) {
            df_raw <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(noncomp_cols)))
            names(df_raw) <- paste0(noncomp_cols, "_raw")
            rows_to_process <- if ("row_excluded" %in% names(df)) !df$row_excluded else rep(TRUE, nrow(df))
            if (sum(rows_to_process) > 0) {
              df_raw[rows_to_process, ] <- df[rows_to_process, noncomp_cols, drop = FALSE]
            }
            df <- dplyr::bind_cols(df, df_raw)
          }

          return(list(
            data = df,
            imputation_method = imputation_method,
            imputation_fallback = imputation_fallback,
            imputation_auto_reason = imputation_auto_reason,
            imputation_floor_msg = imputation_floor_msg,
            transform_skip_msg = transform_skip_msg
          ))
      }
      environment(preprocess_func) <- globalenv()  # Prevent callr from serializing reactive env

      if (USE_BG_PROCESSES) {
        # --- Background path: cancellable via callr ---
        shinyjs::show("cancel_preprocess")
        shinyjs::disable("apply_transforms")
        preprocess_bg_context$notif_id <- showNotification(
          "Preprocessing started in background. You can cancel anytime.",
          type = "message", duration = NULL, closeButton = FALSE
        )
        preprocess_bg_context$comp_cols <- comp_cols
        preprocess_bg_context$noncomp_cols <- noncomp_cols
        preprocess_bg_context$do_impute <- do_impute
        preprocess_bg_context$impute_method <- impute_method
        preprocess_bg_context$transform_type <- transform_type
        preprocess_bg_context$do_pivot <- do_pivot
        preprocess_bg_context$pivot_var <- pivot_var
        preprocess_bg_context$token <- capture_data_token(current_preprocess_context())

        proc <- tryCatch(
          callr::r_bg(
            func = preprocess_func,
            args = list(df = df, comp_cols = comp_cols, noncomp_cols = noncomp_cols,
                        do_impute = do_impute, transform_type = transform_type, pivot_var = pivot_var, impute_method = impute_method),
            supervise = TRUE
          ),
          error = function(e) {
            shinyjs::hide("cancel_preprocess")
            shinyjs::enable("apply_transforms")
            finish_heavy_job()
            removeNotification(preprocess_bg_context$notif_id)
            showNotification(paste("Could not start preprocessing:", e$message), type = "error", duration = 10)
            NULL
          }
        )
        if (is.null(proc)) return()
        preprocess_bg_process(proc)
      } else {
        # --- Synchronous path: no extra process, lower memory ---
        shinyjs::disable("apply_transforms")
        withProgress(message = "Preprocessing data...", value = 0.5, {
          result <- tryCatch(
            preprocess_func(df, comp_cols, noncomp_cols, do_impute, transform_type, pivot_var, impute_method),
            error = function(e) {
              showNotification(
                paste0(
                  "Could not apply preprocessing transforms. Check that selected columns are numeric and contain enough non-missing values. ",
                  "Technical detail: ", e$message
                ),
                type = "error",
                duration = 15
              )
              return(NULL)
            }
          )
          if (!is.null(result)) {
            rv$data <- result$data
            rv$data_stale <- FALSE
            rv$pipeline_config <- list(
              preprocessing = list(
                comp_cols = comp_cols,
                noncomp_cols = noncomp_cols,
                do_impute = do_impute,
                impute_method = impute_method,
                transform_type = transform_type,
                do_pivot = do_pivot,
                pivot_var = pivot_var,
                imputation_method = result$imputation_method
              ),
              timestamp = Sys.time()
            )
            bump_data_generation()
            if (isTRUE(result$imputation_fallback)) {
              fallback_msg <- result$imputation_auto_reason %||%
                "Robust imputation was not possible for this dataset. Standard imputation was used instead."
              showNotification(
                HTML(paste0("<strong>Imputation method note:</strong> ", fallback_msg)),
                type = "warning",
                duration = NULL
              )
            }
            if (!is.null(result$imputation_floor_msg)) {
              showNotification(
                HTML(paste0("<strong>Note:</strong> ", result$imputation_floor_msg)),
                type = "message",
                duration = NULL
              )
            }
            if (!is.null(result$transform_skip_msg)) {
              showNotification(
                HTML(paste0("<strong>Warning:</strong> ", result$transform_skip_msg)),
                type = "warning",
                duration = NULL
              )
            }
            showNotification("Transforms applied.", type = "message")
          }
        })
        shinyjs::enable("apply_transforms")
        finish_heavy_job()
      }
    })

    # Cancel preprocessing
    observeEvent(input$cancel_preprocess, {
      proc <- preprocess_bg_process()
      if (!is.null(proc) && proc$is_alive()) {
        proc$kill()
        try(proc$wait(timeout = 2000), silent = TRUE)
        preprocess_bg_process(NULL)
        shinyjs::hide("cancel_preprocess")
        shinyjs::enable("apply_transforms")
        finish_heavy_job()
        removeNotification(preprocess_bg_context$notif_id)
        showNotification("Preprocessing cancelled.", type = "warning", duration = 5)
      }
    })

    # Poll for preprocessing completion
    observe({
      proc <- preprocess_bg_process()
      req(proc)
      if (proc$is_alive()) {
        invalidateLater(1000)
        return()
      }

      preprocess_bg_process(NULL)
      shinyjs::hide("cancel_preprocess")
      shinyjs::enable("apply_transforms")
      finish_heavy_job()
      removeNotification(preprocess_bg_context$notif_id)

      result <- tryCatch(proc$get_result(), error = function(e) {
        # Extract the actual error from callr subprocess
        err_msg <- e$message
        if (!is.null(e$parent)) {
          err_msg <- paste0(err_msg, " | Cause: ", e$parent$message)
        }
        message("Preprocessing callr error: ", err_msg)  # Log to console
        showNotification(
          paste0(
            "Preprocessing stopped due to a data or transformation issue. Please review selected columns and try again. ",
            "Technical detail: ", err_msg
          ),
          type = "error",
          duration = 30
        )
        return(NULL)
      })

      if (!is.null(result)) {
        if (!data_token_is_current(preprocess_bg_context$token, current_preprocess_context())) {
          showNotification(
            "Preprocessing finished, but its input data or settings changed. The stale result was discarded.",
            type = "warning",
            duration = 8
          )
          return()
        }
        rv$data <- result$data
        rv$data_stale <- FALSE

        # Capture pipeline configuration
        rv$pipeline_config <- list(
          preprocessing = list(
            comp_cols = preprocess_bg_context$comp_cols,
            noncomp_cols = preprocess_bg_context$noncomp_cols,
            do_impute = preprocess_bg_context$do_impute,
            impute_method = preprocess_bg_context$impute_method,
            transform_type = preprocess_bg_context$transform_type %||% if (isTRUE(preprocess_bg_context$do_pivot)) "ilr" else "none",
            do_pivot = preprocess_bg_context$do_pivot,
            pivot_var = preprocess_bg_context$pivot_var,
            imputation_method = result$imputation_method
          ),
          timestamp = Sys.time()
        )
        bump_data_generation()
        if (isTRUE(result$imputation_fallback)) {
          fallback_msg <- result$imputation_auto_reason %||%
            "Robust imputation was not possible for this dataset. Standard imputation was used instead."
          showNotification(
            HTML(paste0("<strong>Imputation method note:</strong> ", fallback_msg)),
            type = "warning",
            duration = NULL
          )
        }
        if (!is.null(result$imputation_floor_msg)) {
          showNotification(
            HTML(paste0("<strong>Note:</strong> ", result$imputation_floor_msg)),
            type = "message",
            duration = NULL
          )
        }
        if (!is.null(result$transform_skip_msg)) {
          showNotification(
            HTML(paste0("<strong>Warning:</strong> ", result$transform_skip_msg)),
            type = "warning",
            duration = NULL
          )
        }

        showNotification("Transforms applied.", type = "message")
      }
    })

    # --- GMM Clustering ---
    observeEvent(input$run_gmm, {
      req(rv$data)

      # Validate G range
      if (!is.null(input$gmin) && !is.null(input$gmax) && input$gmax < input$gmin) {
        showNotification("Max Clusters (G) must be greater than or equal to Min Clusters (G).", type = "error")
        return()
      }

      # 1. Prepare data and parameters in the main session
      cols <- NULL
      src <- "custom" # Default for advanced mode

      if (isTRUE(input$gmm_advanced_mode)) {
        cols <- input$gmm_custom_cols
        if (is.null(cols) || length(cols) == 0) {
          showNotification("Please select at least one column for Advanced GMM.", type = "error")
          return()
        }
      } else {
        req(input$gmm_data_source)
        src <- input$gmm_data_source
        pat <- if (src == "pivot") "^pivot_" else if (src == "clr") "^clr_" else paste0("_", src, "$")
        cols <- grep(pat, names(rv$data), value = TRUE)

        # FIX: Exclude VARG26 pivot columns if selecting generic pivot columns
        if (src == "pivot") {
          cols <- cols[!grepl("^pivot_VARG26_", cols)]
        }
      }

      if (length(cols) == 0) {
        showNotification(
          "Could not start clustering because no columns matched your current GMM data-source selection. Technical detail: No matching columns found for clustering.",
          type = "error"
        )
        return()
      }

      df_sub <- rv$data[, cols, drop = FALSE]
      valid <- complete.cases(df_sub) & apply(df_sub, 1, function(r) all(is.finite(r)))
      
      # Also exclude rows marked for exclusion
      if ("row_excluded" %in% names(rv$data)) {
        valid <- valid & !rv$data$row_excluded
      }

      if (sum(valid) <= length(cols)) {
        showNotification(
          "Could not run clustering because there are too few complete numeric rows. Check missing values and selected columns, then try again. Technical detail: Not enough valid data points for clustering.",
          type = "error"
        )
        return()
      }

      # Extract parameters to pass to future
      data_for_gmm <- df_sub[valid, ]

      # Pre-check for potential issues
      tryCatch(
        {
          # Check for zero variance columns
          col_vars <- apply(data_for_gmm, 2, var, na.rm = TRUE)
          if (any(col_vars < 1e-10, na.rm = TRUE)) {
            zero_var_cols <- names(col_vars)[col_vars < 1e-10]
            showNotification(paste0("Warning: Near-zero variance detected in columns: ", paste(zero_var_cols, collapse = ", "), ". This may cause GMM to fail. Consider using 'pivot' (ILR) transformation."), type = "warning", duration = 10)
          }

          # Check covariance matrix condition for non-pivot data
          if (src != "pivot" && ncol(data_for_gmm) > 1) {
            cov_mat <- tryCatch(cov(data_for_gmm), error = function(e) NULL)
            if (!is.null(cov_mat)) {
              # Check if matrix is nearly singular
              det_val <- tryCatch(det(cov_mat), error = function(e) 0)
              if (abs(det_val) < 1e-10) {
                showNotification(paste0("Warning: Data covariance matrix is nearly singular (determinant = ", format(det_val, scientific = TRUE), "). GMM may fail with 'norm' or 'imp' transformations. Strongly recommend using 'pivot' (ILR) transformation for compositional data."), type = "warning", duration = 15)
              }
            }
          }
        },
        error = function(e) {
          # If pre-check fails, just warn but continue
          showNotification(paste0("Warning during data validation: ", e$message, ". Proceeding with GMM anyway."), type = "warning", duration = 5)
        }
      )

      gmin <- input$gmin
      gmax <- input$gmax
      noise_init <- input$noise_init
      use_prior <- input$use_prior

      # Notify user and set running state
      gmm_worker_cap <- if (isTRUE(USE_BG_PROCESSES)) processing_heavy_worker_cap() else 1L
      if (!begin_heavy_job(
        "GMM clustering",
        min_workers = 1L,
        max_workers = gmm_worker_cap
      )) return()
      gmm_workers <- heavy_job_limiter_workers(heavy_job_lock, default = 1L)
      gmm_running(TRUE)
      shinyjs::show("gmm_progress_container")
      shinyjs::disable("run_gmm")
      rv$mclust_result <- NULL # Clear previous result to trigger spinner

      # Define the GMM function (shared by both paths)
      gmm_func <- function(data, g_min, g_max, noise_init, use_prior) {
        library(mclust)
        mclust::Mclust(
          data,
          G = g_min:g_max,
          initialization = if (noise_init) list(noise = TRUE) else NULL,
          prior = if (use_prior) mclust::priorControl() else NULL,
          verbose = FALSE
        )
      }
      environment(gmm_func) <- globalenv()  # Prevent callr from serializing reactive env

      if (USE_BG_PROCESSES) {
        # --- Background path: cancellable via callr ---
        shinyjs::show("cancel_gmm")
        notif_id <- showNotification(
          paste0(
            "GMM clustering started in background with ", gmm_workers,
            if (gmm_workers == 1L) " CPU worker. " else " CPU workers. ",
            "You can cancel anytime. You'll be notified when finished."
          ),
          type = "message", duration = NULL, closeButton = FALSE
        )

        # Store context for the poller
        gmm_bg_context$valid <- valid
        gmm_bg_context$src <- src
        gmm_bg_context$advanced_mode <- isTRUE(input$gmm_advanced_mode)
        gmm_bg_context$columns_used <- cols
        gmm_bg_context$notif_id <- notif_id
        gmm_bg_context$include_uncertainty <- input$include_uncertainty
        gmm_bg_context$noise_init <- noise_init
        gmm_bg_context$use_prior <- use_prior
        gmm_bg_context$gmin <- gmin
        gmm_bg_context$gmax <- gmax
        gmm_bg_context$workers <- gmm_workers
        gmm_bg_context$token <- capture_data_token(current_gmm_context())

        proc <- tryCatch(
          callr::r_bg(
            func = gmm_func,
            args = list(
              data = data_for_gmm,
              g_min = gmin,
              g_max = gmax,
              noise_init = noise_init,
              use_prior = use_prior
            ),
            env = processing_blas_worker_env(gmm_workers),
            supervise = TRUE
          ),
          error = function(e) {
            gmm_running(FALSE)
            shinyjs::hide("gmm_progress_container")
            shinyjs::hide("cancel_gmm")
            shinyjs::enable("run_gmm")
            finish_heavy_job()
            removeNotification(gmm_bg_context$notif_id)
            showNotification(paste("Could not start GMM clustering:", e$message), type = "error", duration = 10)
            NULL
          }
        )
        if (is.null(proc)) return()
        gmm_bg_process(proc)
      } else {
        # --- Synchronous path: no extra process, lower memory ---
        include_uncertainty <- input$include_uncertainty

        res <- withProgress(message = "Running GMM clustering...", value = 0.5, {
          tryCatch(
            gmm_func(data_for_gmm, gmin, gmax, noise_init, use_prior),
            error = function(e) {
              showNotification(
                paste0(
                  "Clustering failed. Try adjusting the cluster range (Min/Max G), and for compositional data consider using pivot columns. ",
                  "Technical detail: ", e$message
                ),
                type = "error",
                duration = 15
              )
              NULL
            }
          )
        })

        gmm_running(FALSE)
        shinyjs::hide("gmm_progress_container")
        shinyjs::enable("run_gmm")
        finish_heavy_job()

        if (is.null(res) || is.null(res$classification) || length(res$classification) == 0) {
          showNotification(
            "Clustering ran but did not return usable groups. Try a different G range or switch to pivot (ILR) transformed columns. Technical detail: GMM clustering returned no results.",
            type = "error",
            duration = 10
          )
          return()
        }
        if (nrow(rv$data) != length(valid) || length(res$classification) != sum(valid)) {
          showNotification("Clustering returned a result with an unexpected row count; no state was changed.", type = "error")
          return()
        }
        if (isTRUE(include_uncertainty) &&
            (length(res$uncertainty) != sum(valid) || (!is.null(res$z) && nrow(res$z) != sum(valid)))) {
          showNotification("Clustering uncertainty output has an unexpected row count; no state was changed.", type = "error")
          return()
        }

        rv$mclust_result <- res

        # Capture GMM config
        if (is.null(rv$pipeline_config)) {
          rv$pipeline_config <- list(preprocessing = NULL, timestamp = Sys.time())
        }
        rv$pipeline_config$gmm <- list(
          data_source = src,
          advanced_mode = isTRUE(input$gmm_advanced_mode),
          columns_used = cols,
          include_uncertainty = isTRUE(include_uncertainty),
          gmin = gmin,
          gmax = gmax,
          noise_init = noise_init,
          use_prior = use_prior,
          best_model = res$modelName,
          best_G = res$G
        )

        # Update data with clusters
        if ("gmm_cluster" %in% names(rv$data)) rv$data$gmm_cluster <- NULL
        rv$data$gmm_cluster <- NA_integer_

        if (nrow(rv$data) == length(valid)) {
          if (length(res$classification) == sum(valid)) {
            rv$data$gmm_cluster[valid] <- as.integer(res$classification)

            # Add uncertainty and probabilities if requested
            if (include_uncertainty) {
              if ("gmm_uncertainty" %in% names(rv$data)) rv$data$gmm_uncertainty <- NULL
              rv$data$gmm_uncertainty <- as.numeric(NA)
              rv$data$gmm_uncertainty[valid] <- res$uncertainty

              old_prob_cols <- grep("^gmm_prob_", names(rv$data), value = TRUE)
              if (length(old_prob_cols) > 0) {
                rv$data <- rv$data[, !names(rv$data) %in% old_prob_cols]
              }

              if (!is.null(res$z)) {
                n_components <- ncol(res$z)
                for (k in 1:n_components) {
                  if (noise_init && k == n_components) {
                    col_name <- "gmm_prob_outlier"
                  } else {
                    col_name <- paste0("gmm_prob_", k)
                  }
                  rv$data[[col_name]] <- as.numeric(NA)
                  rv$data[[col_name]][valid] <- res$z[, k]
                }
              }

              if (noise_init && !is.null(res$z)) {
                showNotification(
                  "Note: 'gmm_prob_outlier' column represents the probability of belonging to the outlier/noise component.",
                  type = "message", duration = 8
                )
              }
            }

            bump_data_generation()
            success_msg <- paste0("GMM clustering complete! Found ", res$G, " clusters using ", res$modelName, " model.")
            if (noise_init) {
              success_msg <- paste0(success_msg, " (Outlier detection enabled)")
            }
            showNotification(success_msg, type = "message", duration = 5)
          } else {
            showNotification(paste0("Warning: Classification length mismatch. Expected ", sum(valid), " but got ", length(res$classification), ". Model available but not applied to data."), type = "warning", duration = 10)
          }
        } else {
          showNotification("Warning: Data changed during processing. Clusters not applied to table, but model is available.", type = "warning")
        }
      }
    })

    # Cancel GMM
    observeEvent(input$cancel_gmm, {
      proc <- gmm_bg_process()
      if (!is.null(proc) && proc$is_alive()) {
        proc$kill()
        try(proc$wait(timeout = 2000), silent = TRUE)
        gmm_bg_process(NULL)
        gmm_running(FALSE)
        shinyjs::hide("gmm_progress_container")
        shinyjs::hide("cancel_gmm")
        shinyjs::enable("run_gmm")
        finish_heavy_job()
        removeNotification(gmm_bg_context$notif_id)
        showNotification("GMM clustering cancelled.", type = "warning", duration = 5)
      }
    })

    # Poll for GMM completion
    observe({
      proc <- gmm_bg_process()
      req(proc)
      if (proc$is_alive()) {
        invalidateLater(1000)
        return()
      }

      # Process finished — retrieve result
      gmm_bg_process(NULL)
      gmm_running(FALSE)
      shinyjs::hide("gmm_progress_container")
      shinyjs::hide("cancel_gmm")
      shinyjs::enable("run_gmm")
      finish_heavy_job()
      removeNotification(gmm_bg_context$notif_id)

      res <- tryCatch(proc$get_result(), error = function(e) {
        showNotification(
          paste0(
            "Clustering failed in the background. Try changing Min/Max G or using pivot columns, then run again. ",
            "Technical detail: ", e$message
          ),
          type = "error",
          duration = 15
        )
        return(NULL)
      })

      if (is.null(res) || is.null(res$classification) || length(res$classification) == 0) {
        showNotification(
          "Background clustering returned no usable groups. Try a different cluster range or pivot (ILR) transformed inputs. Technical detail: GMM clustering returned no results.",
          type = "error",
          duration = 10
        )
        return()
      }

      if (!data_token_is_current(gmm_bg_context$token, current_gmm_context())) {
        showNotification(
          "Clustering finished, but its input data or settings changed. The stale result was discarded.",
          type = "warning",
          duration = 8
        )
        return()
      }
      valid <- gmm_bg_context$valid
      if (nrow(rv$data) != length(valid) || length(res$classification) != sum(valid)) {
        showNotification("Background clustering returned an unexpected row count; no state was changed.", type = "error")
        return()
      }
      if (isTRUE(gmm_bg_context$include_uncertainty) &&
          (length(res$uncertainty) != sum(valid) || (!is.null(res$z) && nrow(res$z) != sum(valid)))) {
        showNotification("Background clustering uncertainty output has an unexpected row count; no state was changed.", type = "error")
        return()
      }

      rv$mclust_result <- res
      src <- gmm_bg_context$src

      # Capture GMM config
      if (is.null(rv$pipeline_config)) {
        rv$pipeline_config <- list(preprocessing = NULL, timestamp = Sys.time())
      }
      rv$pipeline_config$gmm <- list(
        data_source = src,
        advanced_mode = isTRUE(gmm_bg_context$advanced_mode),
        columns_used = gmm_bg_context$columns_used,
        include_uncertainty = isTRUE(gmm_bg_context$include_uncertainty),
        gmin = gmm_bg_context$gmin,
        gmax = gmm_bg_context$gmax,
        noise_init = gmm_bg_context$noise_init,
        use_prior = gmm_bg_context$use_prior,
        best_model = res$modelName,
        best_G = res$G
      )

      # Update data with clusters
      if ("gmm_cluster" %in% names(rv$data)) rv$data$gmm_cluster <- NULL
      rv$data$gmm_cluster <- NA_integer_

      if (nrow(rv$data) == length(valid)) {
        if (length(res$classification) == sum(valid)) {
          rv$data$gmm_cluster[valid] <- as.integer(res$classification)
          
          # Add uncertainty and probabilities if requested
          if (gmm_bg_context$include_uncertainty) {
            # Add uncertainty column
            if ("gmm_uncertainty" %in% names(rv$data)) rv$data$gmm_uncertainty <- NULL
            rv$data$gmm_uncertainty <- as.numeric(NA)
            rv$data$gmm_uncertainty[valid] <- res$uncertainty
            
            # Add probability columns for each cluster
            old_prob_cols <- grep("^gmm_prob_", names(rv$data), value = TRUE)
            if (length(old_prob_cols) > 0) {
              rv$data <- rv$data[, !names(rv$data) %in% old_prob_cols]
            }
            
            if (!is.null(res$z)) {
              n_components <- ncol(res$z)
              noise_used <- gmm_bg_context$noise_init
              
              for (k in 1:n_components) {
                if (noise_used && k == n_components) {
                  col_name <- "gmm_prob_outlier"
                } else {
                  col_name <- paste0("gmm_prob_", k)
                }
                rv$data[[col_name]] <- as.numeric(NA)
                rv$data[[col_name]][valid] <- res$z[, k]
              }
            }
            
            if (gmm_bg_context$noise_init && !is.null(res$z)) {
              showNotification(
                paste0("Note: 'gmm_prob_outlier' column represents the probability of belonging to the outlier/noise component."),
                type = "message", duration = 8
              )
            }
          }
          
          bump_data_generation()
          success_msg <- paste0("GMM clustering complete! Found ", res$G, " clusters using ", res$modelName, " model.")
          if (gmm_bg_context$noise_init) {
            success_msg <- paste0(success_msg, " (Outlier detection enabled)")
          }
          showNotification(success_msg, type = "message", duration = 5)
        } else {
          showNotification(paste0("Warning: Classification length mismatch. Expected ", sum(valid), " but got ", length(res$classification), ". Model available but not applied to data."), type = "warning", duration = 10)
        }
      } else {
        showNotification("Warning: Data changed during processing. Clusters not applied to table, but model is available.", type = "warning")
      }
    })

    output$bic_plot <- renderPlot({
      req(rv$mclust_result)
      plot(rv$mclust_result, what = "BIC")
    })

    output$bic_summary <- renderText({
      req(rv$mclust_result)
      paste(
        "Best model:", rv$mclust_result$modelName,
        "with G =", rv$mclust_result$G,
        "clusters and BIC =", round(rv$mclust_result$bic, 2)
      )
    })

    # --- Cluster Interpretation Report ---
    output$gmm_has_result <- reactive({
      !is.null(rv$mclust_result) && !is.null(rv$mclust_result$parameters)
    })
    outputOptions(output, "gmm_has_result", suspendWhenHidden = FALSE)

    # Build the cluster report data from model parameters (not recomputed from data)
    gmm_report_data <- reactive({
      req(rv$mclust_result, rv$mclust_result$parameters)
      res <- rv$mclust_result
      params <- res$parameters
      
      # Get the variable names used in GMM
      var_names <- colnames(res$data)
      if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(res$data)))
      
      # Number of clusters
      G <- res$G
      cluster_ids <- seq_len(G)
      
      # Extract model means: matrix [variables x clusters]
      model_means <- params$mean
      if (is.null(dim(model_means))) {
        # 1D case (single variable)
        model_means <- matrix(model_means, nrow = 1)
      }
      
      # Extract model SDs from variance (diagonal of covariance matrices)
      # params$variance$sigma is [p x p x G] array
      model_sds <- matrix(NA, nrow = nrow(model_means), ncol = G)
      if (!is.null(params$variance$sigma)) {
        sigma <- params$variance$sigma
        if (length(dim(sigma)) == 3) {
          for (k in seq_len(G)) {
            model_sds[, k] <- sqrt(diag(sigma[, , k]))
          }
        } else if (is.matrix(sigma)) {
          # Same covariance for all clusters (e.g., EEI model)
          sd_vals <- sqrt(diag(sigma))
          for (k in seq_len(G)) model_sds[, k] <- sd_vals
        } else {
          # Scalar variance (1D)
          for (k in seq_len(G)) model_sds[, k] <- sqrt(sigma)
        }
      }
      
      # Sample counts per cluster - ensure it matches G
      n_per_cluster <- as.integer(table(factor(res$classification, levels = cluster_ids)))
      
      # Mixing proportions - ensure it matches G
      mix_props <- if (!is.null(params$pro)) {
        p <- as.numeric(params$pro)
        if (length(p) != G) {
          # Pad or truncate if there's a mismatch (shouldn't happen with mclust)
          if (length(p) < G) p <- c(p, rep(NA, G - length(p))) else p <- p[1:G]
        }
        round(p * 100, 1)
      } else {
        rep(NA, G)
      }
      
      # Build display data frame: one row per cluster
      display_df <- data.frame(
        Cluster = as.character(cluster_ids),
        n = n_per_cluster,
        `Proportion (%)` = mix_props,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      for (i in seq_along(var_names)) {
        display_df[[var_names[i]]] <- paste0(
          formatC(model_means[i, ], format = "f", digits = 2),
          " \u00B1 ",
          formatC(model_sds[i, ], format = "f", digits = 2)
        )
      }
      
      display_df
    })

    # Build cluster composition summary from untransformed oxide columns in rv$data
    gmm_cluster_composition_data <- reactive({
      req(rv$data)
      selected_comp <- processing_gmm_selected_compositional_columns(
        rv$pipeline_config,
        preproc_state$comp_cols
      )
      processing_gmm_cluster_composition(rv$data, selected_comp)
    })

    output$gmm_cluster_report <- renderDT({
      req(gmm_report_data())
      df <- gmm_report_data()
      
      datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "t",
          ordering = FALSE,
          columnDefs = list(
            list(className = "dt-center", targets = "_all")
          )
        ),
        caption = htmltools::tags$caption(
          style = "caption-side: bottom; font-size: 12px; color: #666;",
          "Values shown as Mean \u00B1 SD from the fitted model parameters (not recomputed from data). Proportion is the model-estimated mixing weight."
        )
      ) %>%
        formatStyle("Cluster", fontWeight = "bold") %>%
        formatStyle("n", fontWeight = "bold") %>%
        formatStyle("Proportion (%)", fontStyle = "italic")
    })

    output$gmm_cluster_composition_summary <- renderDT({
      req(gmm_cluster_composition_data())
      df <- gmm_cluster_composition_data()
      data_source <- attr(df, "data_source") %||% "raw"
      source_label <- processing_gmm_composition_source_label(data_source)
      caption_text <- paste0(
        "Values shown as Mean \u00B1 1 SD computed from the ", source_label,
        " for the compositional columns used in the transformation."
      )

      datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "t",
          ordering = FALSE,
          columnDefs = list(
            list(className = "dt-center", targets = "_all")
          )
        ),
        caption = htmltools::tags$caption(
          style = "caption-side: bottom; font-size: 12px; color: #666;",
          caption_text
        )
      ) %>%
        formatStyle("Cluster", fontWeight = "bold") %>%
        formatStyle("n", fontWeight = "bold")
    })

    output$download_cluster_report <- downloadHandler(
      filename = function() {
        paste0("GMM_Cluster_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv")
      },
      content = function(file) {
        req(rv$mclust_result, rv$mclust_result$parameters)
        res <- rv$mclust_result
        params <- res$parameters
        var_names <- colnames(res$data)
        if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(res$data)))
        
        G <- res$G
        model_means <- params$mean
        if (is.null(dim(model_means))) model_means <- matrix(model_means, nrow = 1)
        
        # Extract SDs from covariance
        model_sds <- matrix(NA, nrow = nrow(model_means), ncol = G)
        if (!is.null(params$variance$sigma)) {
          sigma <- params$variance$sigma
          if (length(dim(sigma)) == 3) {
            for (k in seq_len(G)) model_sds[, k] <- sqrt(diag(sigma[, , k]))
          } else if (is.matrix(sigma)) {
            sd_vals <- sqrt(diag(sigma))
            for (k in seq_len(G)) model_sds[, k] <- sd_vals
          } else {
            for (k in seq_len(G)) model_sds[, k] <- sqrt(sigma)
          }
        }
        
        n_per_cluster <- table(factor(res$classification, levels = seq_len(G)))
        mix_props <- if (!is.null(params$pro)) round(params$pro * 100, 2) else rep(NA, G)
        
        # Long format: one row per cluster-variable combination
        report_rows <- lapply(seq_len(G), function(k) {
          data.frame(
            Cluster = k,
            N_Samples = as.integer(n_per_cluster[k]),
            Proportion_Pct = mix_props[k],
            Variable = var_names,
            Model_Mean = round(model_means[, k], 4),
            Model_SD = round(model_sds[, k], 4),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        })
        
        full_report <- do.call(rbind, report_rows)
        
        # Add model metadata as header rows
        n_cols <- ncol(full_report)
        meta <- data.frame(
          Cluster = c("# Model", "# G", "# BIC", "---"),
          N_Samples = c(res$modelName, res$G, round(res$bic, 2), "---"),
          matrix("", nrow = 4, ncol = n_cols - 2),
          stringsAsFactors = FALSE
        )
        names(meta) <- names(full_report)
        
        write.csv(rbind(meta, full_report), file, row.names = FALSE, na = "")

        comp_report <- tryCatch(gmm_cluster_composition_data(), error = function(e) NULL)
        if (!is.null(comp_report) && nrow(comp_report) > 0) {
          cat("\n\"# Cluster Composition Summary\"\n", file = file, append = TRUE)
          cat(
            "\"# Mean \u00B1 1 SD of the compositional columns used in the transformation (imputed values where available; otherwise raw).\"\n",
            file = file,
            append = TRUE
          )
          utils::write.table(
            comp_report,
            file = file,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            na = "",
            append = TRUE,
            qmethod = "double"
          )
        }
      }
    )

    # --- VARG26 Projection ---
    ensure_VARG26_inputs <- function(data_df) {
      missing <- setdiff(VARG26_OXIDES, names(data_df))
      if (length(missing) > 0) {
        stop(paste0(
          "Missing required VARG26 oxides: ", paste(missing, collapse = ", "),
          ". Go back to Step 1 (Import Data) and confirm these columns exist in your uploaded file, ",
          "or download the data template for the expected format."
        ))
      }

      # Check if user has already imputed all VARG26 oxides
      imp_cols <- paste0(VARG26_OXIDES, "_imp")
      all_imp_exist <- all(imp_cols %in% names(data_df))

      if (all_imp_exist) {
        # User has pre-imputed data - use it directly (just re-normalize)
        showNotification("Using existing imputed VARG26 oxides (faster)", type = "message", duration = 3)
        prep <- data_df[, imp_cols, drop = FALSE]
        colnames(prep) <- VARG26_OXIDES # Remove _imp suffix for processing

        # Re-normalize to 100% (compositional closure)
        row_sums <- rowSums(prep, na.rm = TRUE)
        imputed <- prep / row_sums * 100

        # Replace any remaining issues
        imputed[!is.finite(as.matrix(imputed))] <- NA
      } else {
        # No pre-imputed data - perform imputation from raw oxides
        showNotification("Imputing VARG26 oxides from raw data", type = "message", duration = 3)
        prep <- data_df[, VARG26_OXIDES, drop = FALSE]
        prep[prep <= 1e-10] <- NA

        # Check if there's enough complete data to seed the imputation
        complete_count <- sum(complete.cases(prep))
        if (complete_count < 2) {
          stop("Not enough complete rows to perform imputation. Need at least 2 complete rows.")
        }

        # Perform imputation on all rows
        set.seed(42)  # Ensure reproducible imputation
        imputed <- tryCatch(
          robCompositions::impCoda(prep, method = "ltsReg", maxit = 5)$xImp,
          error = function(e) {
            message("impCoda ltsReg failed (", e$message, "), falling back to method='lm'")
            robCompositions::impCoda(prep, method = "lm", maxit = 5)$xImp
          }
        )
      }

      # Now validate rows after imputation/normalization
      valid_idx <- which(complete.cases(imputed) &
        apply(imputed, 1, function(r) all(is.finite(r) & r > 0)))

      if (length(valid_idx) < 2) {
        stop("Not enough valid rows after imputation.")
      }

      # Apply pivot transformation to valid rows only
      pivot_res <- robCompositions::pivotCoord(as.matrix(imputed[valid_idx, ]), pivot = 1)
      set.seed(42)  # Ensure reproducible noise addition
      # A tiny amount of noise (1e-5) is added before VARG26 projection to ensure that
      # samples already present in the reference database produce consistent nearest-neighbor
      # results rather than exact self-matches. This does not meaningfully affect the projection.
      noisy <- sdcMicro::addNoise(as.data.frame(pivot_res), noise = ADD_NOISE_LEVEL)$xm

      res_df <- as.data.frame(matrix(NA, nrow = nrow(data_df), ncol = ncol(noisy)))
      colnames(res_df) <- paste0("pivot_VARG26_", colnames(noisy))
      res_df[valid_idx, ] <- noisy

      return(bind_cols(data_df, res_df))
    }

    observeEvent(input$run_umap_pre, {
      req(rv$data, input$VARG26_dims)

      if (!check_models()) {
        showNotification(
          "VARG26 projection cannot start because required model files were not found. Check your VARG-Tools installation and model paths. Technical detail: Model files missing.",
          type = "error"
        )
        return()
      }

      dims_choice <- input$VARG26_dims

      # Prepare data in main session
      df <- rv$data
      df <- df[, !grepl("^pivot_VARG26_", names(df))]
      if (dims_choice == "1d") {
        df <- df[, !grepl("^UMAP_VARG26_1D$", names(df))]
      } else if (dims_choice == "2d") {
        df <- df[, !grepl("^UMAP_VARG26_2D_", names(df))]
      } else {
        df <- df[, !grepl("^UMAP_VARG26_", names(df))]
      }

      # Check required oxides
      varg26_oxides <- VARG26_OXIDES
      missing <- setdiff(varg26_oxides, names(df))
      if (length(missing) > 0) {
        showNotification(
          paste0(
            "Missing required VARG26 oxides: ", paste(missing, collapse = ", "),
            ". Go back to Step 1 and confirm these columns exist in the uploaded file, ",
            "or download the data template for the expected format."
          ),
          type = "error",
          duration = 15
        )
        return()
      }

      # Resolve absolute model paths
      model_path_1d <- normalizePath(MODEL_PATH_1D, mustWork = FALSE)
      model_path_2d <- normalizePath(MODEL_PATH_2D, mustWork = FALSE)
      add_noise <- ADD_NOISE_LEVEL
      imp_cols <- paste0(varg26_oxides, "_imp")
      has_imp <- all(imp_cols %in% names(df))

      # Define the VARG26 projection function (shared by both paths)
      varg26_func <- function(df, dims_choice, varg26_oxides, model_path_1d, model_path_2d, add_noise, has_imp, workers = 1L) {
          # --- Inline ensure_VARG26_inputs ---
          imp_cols <- paste0(varg26_oxides, "_imp")
          if (has_imp) {
            prep <- df[, imp_cols, drop = FALSE]
            colnames(prep) <- varg26_oxides
            row_sums <- rowSums(prep, na.rm = TRUE)
            imputed <- prep / row_sums * 100
            imputed[!is.finite(as.matrix(imputed))] <- NA
          } else {
            prep <- df[, varg26_oxides, drop = FALSE]
            prep[prep <= 1e-10] <- NA
            if (sum(complete.cases(prep)) < 2) stop("Not enough complete rows for imputation.")
            set.seed(42)
            imputed <- tryCatch(
              robCompositions::impCoda(prep, method = "ltsReg", maxit = 5)$xImp,
              error = function(e) {
                message("impCoda ltsReg failed (", e$message, "), falling back to method='lm'")
                robCompositions::impCoda(prep, method = "lm", maxit = 5)$xImp
              }
            )
          }

          valid_idx <- which(complete.cases(imputed) &
            apply(imputed, 1, function(r) all(is.finite(r) & r > 0)))
          if (length(valid_idx) < 2) stop("Not enough valid rows after imputation.")

          pivot_res <- robCompositions::pivotCoord(as.matrix(imputed[valid_idx, ]), pivot = 1)
          set.seed(42)
          # A tiny amount of noise (1e-5) is added before VARG26 projection to ensure that
          # samples already present in the reference database produce consistent nearest-neighbor
          # results rather than exact self-matches. This does not meaningfully affect the projection.
          noisy <- sdcMicro::addNoise(as.data.frame(pivot_res), noise = add_noise)$xm

          res_df <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(noisy)))
          colnames(res_df) <- paste0("pivot_VARG26_", colnames(noisy))
          res_df[valid_idx, ] <- noisy
          df <- dplyr::bind_cols(df, res_df)

          # --- Projection ---
          p_cols <- grep("^pivot_VARG26_", names(df), value = TRUE)
          if (length(p_cols) != 8) stop(paste("Expected 8 pivot columns, got", length(p_cols)))
          matrix_in <- as.matrix(df[, p_cols])
          valid <- complete.cases(matrix_in)
          if ("row_excluded" %in% names(df)) valid <- valid & !df$row_excluded
          if (sum(valid) == 0) stop("No valid data for projection.")

          if (dims_choice %in% c("1d", "both")) {
            model_1d <- uwot::load_uwot(model_path_1d)
            emb1 <- uwot::umap_transform(
              matrix_in[valid, ], model_1d,
              n_threads = workers,
              n_sgd_threads = 1
            )
            df$UMAP_VARG26_1D <- NA
            df$UMAP_VARG26_1D[valid] <- emb1[, 1]
          }
          if (dims_choice %in% c("2d", "both")) {
            model_2d <- uwot::load_uwot(model_path_2d)
            emb2 <- uwot::umap_transform(
              matrix_in[valid, ], model_2d,
              n_threads = workers,
              n_sgd_threads = 1
            )
            df$UMAP_VARG26_2D_1 <- NA
            df$UMAP_VARG26_2D_2 <- NA
            df$UMAP_VARG26_2D_1[valid] <- emb2[, 1]
            df$UMAP_VARG26_2D_2[valid] <- emb2[, 2]
          }

          return(df)
      }
      environment(varg26_func) <- globalenv()  # Prevent callr from serializing reactive env

      if (!begin_heavy_job(
        "VARG26 projection",
        min_workers = 1L,
        max_workers = processing_heavy_worker_cap()
      )) return()
      varg26_workers <- heavy_job_limiter_workers(heavy_job_lock, default = 1L)
      if (USE_BG_PROCESSES) {
        # --- Background path: cancellable via callr ---
        shinyjs::show("cancel_varg26")
        shinyjs::disable("run_umap_pre")
        varg26_bg_context$notif_id <- showNotification(
          paste0(
            "VARG26 projection started in background with ", varg26_workers,
            if (varg26_workers == 1L) " CPU worker. " else " CPU workers. ",
            "You can cancel anytime."
          ),
          type = "message", duration = NULL, closeButton = FALSE
        )
        varg26_bg_context$dims_choice <- dims_choice
        varg26_bg_context$token <- capture_data_token(current_varg26_context())

        proc <- tryCatch(
          callr::r_bg(
            func = varg26_func,
            args = list(df = df, dims_choice = dims_choice, varg26_oxides = varg26_oxides,
                        model_path_1d = model_path_1d, model_path_2d = model_path_2d,
                        add_noise = add_noise, has_imp = has_imp,
                        workers = varg26_workers),
            supervise = TRUE
          ),
          error = function(e) {
            shinyjs::hide("cancel_varg26")
            shinyjs::enable("run_umap_pre")
            finish_heavy_job()
            removeNotification(varg26_bg_context$notif_id)
            showNotification(paste("Could not start VARG26 projection:", e$message), type = "error", duration = 10)
            NULL
          }
        )
        if (is.null(proc)) return()
        varg26_bg_process(proc)
      } else {
        # --- Synchronous path: no extra process, lower memory ---
        shinyjs::disable("run_umap_pre")
        result <- withProgress(message = "Running VARG26 projection...", value = 0.5, {
          tryCatch(
            varg26_func(df, dims_choice, varg26_oxides, model_path_1d, model_path_2d, add_noise, has_imp, varg26_workers),
            error = function(e) {
              showNotification(
                paste0(
                  "VARG26 projection failed. Verify required oxide columns and check for enough complete rows before retrying. ",
                  "Technical detail: ", e$message
                ),
                type = "error",
                duration = 15
              )
              NULL
            }
          )
        })

        shinyjs::enable("run_umap_pre")
        finish_heavy_job()

        if (!is.null(result)) {
          rv$data <- result
          rv$umap_mode_ran <- "pretrained"
          clear_user_umap_model()

          if (!is.null(rv$pipeline_config)) {
            rv$pipeline_config$umap <- list(
              mode = "pretrained",
              VARG26_dims = dims_choice,
              VARG26_oxides = VARG26_OXIDES
            )
          }
          bump_data_generation()
          focus_umap_result_tab("pretrained", dims_choice)

          showNotification("VARG26 projection successful!", type = "message")
        }
      }
    })

    # Cancel VARG26 projection
    observeEvent(input$cancel_varg26, {
      proc <- varg26_bg_process()
      if (!is.null(proc) && proc$is_alive()) {
        proc$kill()
        try(proc$wait(timeout = 2000), silent = TRUE)
        varg26_bg_process(NULL)
        shinyjs::hide("cancel_varg26")
        shinyjs::enable("run_umap_pre")
        finish_heavy_job()
        removeNotification(varg26_bg_context$notif_id)
        showNotification("VARG26 projection cancelled.", type = "warning", duration = 5)
      }
    })

    # Poll for VARG26 projection completion
    observe({
      proc <- varg26_bg_process()
      req(proc)
      if (proc$is_alive()) {
        invalidateLater(1000)
        return()
      }

      varg26_bg_process(NULL)
      shinyjs::hide("cancel_varg26")
      shinyjs::enable("run_umap_pre")
      finish_heavy_job()
      removeNotification(varg26_bg_context$notif_id)

      result <- tryCatch(proc$get_result(), error = function(e) {
        showNotification(
          paste0(
            "VARG26 projection failed in the background. Check required oxides and missing values, then try again. ",
            "Technical detail: ", e$message
          ),
          type = "error",
          duration = 15
        )
        return(NULL)
      })

      if (!is.null(result)) {
        if (!data_token_is_current(varg26_bg_context$token, current_varg26_context())) {
          showNotification(
            "VARG26 projection finished, but its input data or settings changed. The stale result was discarded.",
            type = "warning",
            duration = 8
          )
          return()
        }
        rv$data <- result
        rv$umap_mode_ran <- "pretrained"
        clear_user_umap_model()

        if (!is.null(rv$pipeline_config)) {
          rv$pipeline_config$umap <- list(
            mode = "pretrained",
            VARG26_dims = varg26_bg_context$dims_choice,
            VARG26_oxides = VARG26_OXIDES
          )
        }
        bump_data_generation()
        focus_umap_result_tab("pretrained", varg26_bg_context$dims_choice)

        showNotification("VARG26 projection successful!", type = "message")
      }
    })

    # --- New UMAP ---
    observeEvent(input$run_umap_new, {
      req(rv$data, input$umap_dims)
      if (!begin_heavy_job(
        "new UMAP",
        min_workers = 1L,
        max_workers = processing_heavy_worker_cap()
      )) return()
      on.exit(finish_heavy_job(), add = TRUE)
      umap_workers <- heavy_job_limiter_workers(heavy_job_lock, default = 1L)

      # Prepare data in main session first
      tryCatch(
        {
          cols <- NULL
          src <- "custom"

          if (isTRUE(input$umap_advanced_mode)) {
            cols <- input$umap_custom_cols
            if (is.null(cols) || length(cols) == 0) stop("Please select at least one column for Advanced UMAP.")
          } else {
            req(input$umap_data_source)
            src <- input$umap_data_source
            pat <- if (src == "pivot") "^pivot_" else if (src == "clr") "^clr_" else paste0("_", src, "$")
            cols <- grep(pat, names(rv$data), value = TRUE)

            # FIX: Exclude VARG26 pivot columns if selecting generic pivot columns
            if (src == "pivot") {
              cols <- cols[!grepl("^pivot_VARG26_", cols)]
            }
          }

          if (length(cols) == 0) stop("No matching columns found.")

          mat <- as.matrix(rv$data[, cols])
          valid <- complete.cases(mat) & apply(mat, 1, function(r) all(is.finite(r)))
          
          # Also exclude rows marked for exclusion
          if ("row_excluded" %in% names(rv$data)) {
            valid <- valid & !rv$data$row_excluded
          }
          
          if (sum(valid) < input$umap_n_neighbors) stop("Not enough data points.")

          # Prepare Target Data (Y) for Semi-Supervised UMAP
          Y_df <- NULL
          if (input$use_semisupervised && !is.null(input$y_vars)) {
            Y_df <- as.data.frame(rv$data[valid, input$y_vars, drop = FALSE])

            for (col in names(Y_df)) {
              if (is.numeric(Y_df[[col]])) {
                if (any(is.na(Y_df[[col]]))) {
                  stop(paste0("Numeric target variable '", col, "' contains missing values."))
                }
              } else {
                Y_df[[col]][Y_df[[col]] == ""] <- NA
                Y_df[[col]] <- as.factor(Y_df[[col]])
              }
            }

            if (ncol(Y_df) == 1) {
              Y_df <- Y_df[[1]]
            }
          }

          # Extract parameters
          dims_choice <- input$umap_dims
          n_neighbors <- input$umap_n_neighbors
          min_dist <- input$umap_min_dist
          dens_scale <- input$umap_dens_scale
          umap_seed <- input$umap_seed
          target_weight <- if (!is.null(Y_df)) input$umap_target_weight else 0.5
          data_for_umap <- mat[valid, ]

          # Use withProgress for synchronous UMAP
          withProgress(
            message = paste0(
              "Running UMAP with ", umap_workers,
              if (umap_workers == 1L) " CPU worker..." else " CPU workers..."
            ),
            value = 0,
            {
            incProgress(0.1, detail = "Initializing...")

            # Update data — remove only the dimension(s) being replaced (preserve other dims, VARG26, loaded)
            df_processed <- rv$data
            if (dims_choice == "1d") {
              df_processed <- df_processed[, !grepl("^UMAP_new_1D$", names(df_processed))]
            } else if (dims_choice == "2d") {
              df_processed <- df_processed[, !grepl("^UMAP_new_2D_", names(df_processed))]
            } else {
              df_processed <- df_processed[, !grepl("^UMAP_new_", names(df_processed))]
            }

            umap_models <- list()

            # Run 1D if requested
            if (dims_choice %in% c("1d", "both")) {
              incProgress(0.2, detail = "Running 1D UMAP...")
              res_1d <- uwot::umap(
                data_for_umap,
                n_neighbors = n_neighbors,
                min_dist = min_dist,
                n_components = 1,
                metric = "euclidean",
                dens_scale = dens_scale,
                y = Y_df,
                target_weight = target_weight,
                ret_model = TRUE,
                seed = umap_seed,
                n_threads = umap_workers,
                n_sgd_threads = 1
              )
              df_processed$UMAP_new_1D <- NA
              df_processed$UMAP_new_1D[valid] <- res_1d$embedding[, 1]
              umap_models[["1d"]] <- res_1d
            }

            # Run 2D if requested
            if (dims_choice %in% c("2d", "both")) {
              incProgress(0.5, detail = "Running 2D UMAP...")
              res_2d <- uwot::umap(
                data_for_umap,
                n_neighbors = n_neighbors,
                min_dist = min_dist,
                n_components = 2,
                metric = "euclidean",
                dens_scale = dens_scale,
                y = Y_df,
                target_weight = target_weight,
                ret_model = TRUE,
                seed = umap_seed,
                n_threads = umap_workers,
                n_sgd_threads = 1
              )
              df_processed$UMAP_new_2D_1 <- NA
              df_processed$UMAP_new_2D_2 <- NA
              df_processed$UMAP_new_2D_1[valid] <- res_2d$embedding[, 1]
              df_processed$UMAP_new_2D_2[valid] <- res_2d$embedding[, 2]
              umap_models[["2d"]] <- res_2d
            }

            incProgress(0.9, detail = "Finalizing...")

            rv$data <- df_processed
            rv$umap_mode_ran <- "new"
            clear_user_umap_model()
            rv$user_umap_model <- umap_models

            # Capture UMAP config
            if (is.null(rv$pipeline_config)) {
              rv$pipeline_config <- list(preprocessing = NULL, timestamp = Sys.time())
            }
            saved_dimensions <- if (identical(dims_choice, "both")) c(1L, 2L) else {
              if (identical(dims_choice, "1d")) 1L else 2L
            }
            rv$pipeline_config$umap <- list(
              mode = "new",
              data_source = src,
              advanced_mode = isTRUE(input$umap_advanced_mode),
              n_neighbors = n_neighbors,
              min_dist = min_dist,
              dens_scale = dens_scale,
              seed = umap_seed,
              dimensions = saved_dimensions,
              n_components = saved_dimensions,
              n_input_cols = length(cols),
              columns_used = cols,
              semisupervised = input$use_semisupervised,
              y_vars = if (isTRUE(input$use_semisupervised)) input$y_vars else character(0),
              target_weight = target_weight,
              projection_contract = list(
                version = 2L,
                type = "direct_columns",
                columns = cols,
                require_complete = TRUE
              )
            )
            bump_data_generation()
            focus_umap_result_tab("new", saved_dimensions)

            showNotification("UMAP created successfully!", type = "message")
            }
          )
        },
        error = function(e) {
          showNotification(
            paste0(
              "UMAP could not be created from the selected inputs. Check selected columns and missing values, then retry. ",
              "Technical detail: ", e$message
            ),
            type = "error",
            duration = 10
          )
        }
      )
    })

    # --- UMAP Plots ---
    output$umap_plots_ui <- renderUI({
      req(rv$data)
      df <- rv$data
      
      tabs <- list()
      
      # VARG26 pretrained
      if (all(c("UMAP_VARG26_2D_1", "UMAP_VARG26_2D_2") %in% names(df))) {
        tabs <- c(tabs, list(tabPanel("VARG26 2D", plotlyOutput(ns("umap_plot_varg26_2d"), height = "500px"))))
      }
      if ("UMAP_VARG26_1D" %in% names(df)) {
        tabs <- c(tabs, list(tabPanel("VARG26 1D", plotlyOutput(ns("umap_plot_varg26_1d"), height = "500px"))))
      }
      
      # User-trained new UMAP
      if (all(c("UMAP_new_2D_1", "UMAP_new_2D_2") %in% names(df))) {
        tabs <- c(tabs, list(tabPanel("New 2D", plotlyOutput(ns("umap_plot_new_2d"), height = "500px"))))
      }
      if ("UMAP_new_1D" %in% names(df)) {
        tabs <- c(tabs, list(tabPanel("New 1D", plotlyOutput(ns("umap_plot_new_1d"), height = "500px"))))
      }
      
      # Loaded from file
      if (all(c("UMAP_loaded_2D_1", "UMAP_loaded_2D_2") %in% names(df))) {
        tabs <- c(tabs, list(tabPanel("Loaded 2D", plotlyOutput(ns("umap_plot_loaded_2d"), height = "500px"))))
      }
      if ("UMAP_loaded_1D" %in% names(df)) {
        tabs <- c(tabs, list(tabPanel("Loaded 1D", plotlyOutput(ns("umap_plot_loaded_1d"), height = "500px"))))
      }
      
      if (length(tabs) == 0) {
        p("No UMAP results yet. Run UMAP to see visualization.")
      } else if (length(tabs) == 1) {
        tagList(
          tabs[[1]],
          tags$p(
            style = "font-style: italic; font-size: 0.85em;",
            class = "text-muted mt-2 mb-0",
            "Note: UMAP axes represent compressed multivariate relationships and do not correspond to individual oxide concentrations."
          )
        )
      } else {
        tagList(
          tags$style(HTML("
            .inner-tabs .nav-tabs { background-color: #e9ecef !important; border-bottom: 1px solid #dee2e6 !important; }
            .inner-tabs .nav-tabs .nav-link { color: #495057 !important; background-color: transparent !important; border: none !important; }
            .inner-tabs .nav-tabs .nav-link:hover { background-color: rgba(0,0,0,0.05) !important; color: #212529 !important; }
            .inner-tabs .nav-tabs .nav-link.active { background-color: #ffffff !important; color: #2c3e50 !important; font-weight: 700 !important; border-top: 3px solid #1565c0 !important; }
            .inner-tabs .nav-tabs > li > a { color: #495057 !important; background-color: transparent !important; }
            .inner-tabs .nav-tabs > li.active > a { background-color: #ffffff !important; color: #2c3e50 !important; font-weight: 700 !important; border-top: 3px solid #1565c0 !important; }
          ")),
          div(class = "inner-tabs",
            do.call(tabsetPanel, c(list(id = ns("umap_tab")), tabs))
          ),
          tags$p(
            style = "font-style: italic; font-size: 0.85em;",
            class = "text-muted mt-2 mb-0",
            "Note: UMAP axes represent compressed multivariate relationships and do not correspond to individual oxide concentrations."
          )
        )
      }
    })

    # Helper: build a 2D UMAP plotly
    build_umap_2d <- function(df, x_col, y_col) {
      if (!(x_col %in% names(df)) || !(y_col %in% names(df))) return(NULL)
      color_var <- if ("gmm_cluster" %in% names(df)) "gmm_cluster" else "population"
      df[[color_var]] <- as.factor(df[[color_var]])
      p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                           color = .data[[color_var]],
                           text = paste("UID:", UID))) +
        geom_point(alpha = 0.7) + theme_minimal() +
        labs(title = paste0("UMAP 2D"), color = color_var)
      ggplotly(p, tooltip = "text")
    }
    
    # Helper: build a 1D UMAP plotly
    build_umap_1d <- function(df, x_col) {
      if (!(x_col %in% names(df))) return(NULL)
      color_var <- if ("gmm_cluster" %in% names(df)) "gmm_cluster" else "population"
      df[[color_var]] <- as.factor(df[[color_var]])
      p <- ggplot(df, aes(x = .data[[x_col]], y = 0,
                           color = .data[[color_var]],
                           text = paste("UID:", UID))) +
        geom_point(alpha = 0.7, position = position_jitter(width = 0, height = 0.1)) +
        theme_minimal() +
        labs(title = "UMAP 1D", x = x_col, y = NULL, color = color_var) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotly(p, tooltip = "text")
    }

    # Render outputs for each possible UMAP combination
    output$umap_plot_varg26_2d <- renderPlotly({ req(rv$data); build_umap_2d(rv$data, "UMAP_VARG26_2D_1", "UMAP_VARG26_2D_2") })
    output$umap_plot_varg26_1d <- renderPlotly({ req(rv$data); build_umap_1d(rv$data, "UMAP_VARG26_1D") })
    output$umap_plot_new_2d    <- renderPlotly({ req(rv$data); build_umap_2d(rv$data, "UMAP_new_2D_1", "UMAP_new_2D_2") })
    output$umap_plot_new_1d    <- renderPlotly({ req(rv$data); build_umap_1d(rv$data, "UMAP_new_1D") })
    output$umap_plot_loaded_2d <- renderPlotly({ req(rv$data); build_umap_2d(rv$data, "UMAP_loaded_2D_1", "UMAP_loaded_2D_2") })
    output$umap_plot_loaded_1d <- renderPlotly({ req(rv$data); build_umap_1d(rv$data, "UMAP_loaded_1D") })

    # --- UMAP Model Save/Load ---
    loaded_external_model <- reactiveValues(models = NULL, config = NULL, temp_dir = NULL)

    cleanup_loaded_external_model <- function(clear_state = TRUE) {
      models <- isolate(loaded_external_model$models)
      config <- isolate(loaded_external_model$config)
      temp_dir <- isolate(loaded_external_model$temp_dir)
      if (!is.null(models)) {
        unload_umap_model_collection(models, if (is.null(config)) NULL else config$umap)
      }
      if (!is.null(temp_dir) && dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = TRUE)
      }
      if (isTRUE(clear_state)) {
        loaded_external_model$models <- NULL
        loaded_external_model$config <- NULL
        loaded_external_model$temp_dir <- NULL
      }
      invisible(NULL)
    }

    session$onSessionEnded(function() {
      cleanup_loaded_external_model(clear_state = FALSE)
    })

    umap_dimension_label <- function(umap_config) {
      dimensions <- normalize_umap_dimensions(umap_config)
      if (length(dimensions) == 0) return("Unknown-dimension")
      paste0(dimensions, "D", collapse = " + ")
    }

    # UI: Show model status
    output$umap_model_status_ui <- renderUI({
      mode <- rv$umap_mode_ran

      if (!is.null(loaded_external_model$models)) {
        config <- loaded_external_model$config$umap
        n_cols <- if (is.null(config$columns_used)) "?" else length(unlist(config$columns_used))
        div(
          class = "alert alert-info py-2 mb-2",
          icon("download"),
          strong(" Loaded Model Ready: "),
          paste0(umap_dimension_label(config), " UMAP; requires ", n_cols, " prepared columns")
        )
      } else if (!is.null(rv$user_umap_model) && identical(mode, "new")) {
        config <- rv$pipeline_config$umap
        n_cols <- if (is.null(config$columns_used)) "?" else length(config$columns_used)
        div(
          class = "alert alert-success py-2 mb-2",
          icon("check-circle"),
          strong(" Model Ready: "),
          paste0(umap_dimension_label(config), " UMAP trained on ", n_cols, " columns")
        )
      } else if (identical(mode, "pretrained")) {
        div(
          class = "alert alert-warning py-2 mb-2",
          icon("database"),
          strong(" Pretrained: "),
          "Using VARG26 reference UMAP"
        )
      } else {
        div(
          class = "alert alert-secondary py-2 mb-2",
          icon("info-circle"),
          " No user UMAP model available. Create one first."
        )
      }
    })

    # Download handler for saving UMAP model
    output$save_umap_model <- downloadHandler(
      filename = function() {
        paste0("UMAP_Model_", format(Sys.time(), "%Y%m%d_%H%M"), ".varg_umap")
      },
      content = function(file) {
        if (is.null(rv$user_umap_model)) {
          showNotification("No UMAP model to save. Create a UMAP first.", type = "error")
          return(NULL)
        }
        if (is.null(rv$pipeline_config) || is.null(rv$pipeline_config$umap)) {
          showNotification("UMAP configuration not available. Please re-run UMAP.", type = "error")
          return(NULL)
        }

        contract <- validate_umap_projection_contract(rv$pipeline_config)
        if (!contract$valid) {
          showNotification(paste(contract$errors, collapse = "\n"), type = "error", duration = 12)
          return(NULL)
        }

        showModal(modalDialog("Saving UMAP model...", footer = NULL))
        on.exit(removeModal())

        temp_dir <- tempfile("umap_bundle_")
        dir.create(temp_dir)
        on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

        model_files <- save_umap_model_collection(
          rv$user_umap_model,
          rv$pipeline_config$umap,
          temp_dir
        )

        config_file <- file.path(temp_dir, "config.json")
        config_data <- list(
          pipeline_config = rv$pipeline_config,
          source_filename = rv$source_filename,
          metadata = list(
            format = "varg_umap",
            schema_version = 2L,
            created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            app_version = APP_VERSION,
            uwot_version = as.character(utils::packageVersion("uwot"))
          )
        )
        jsonlite::write_json(config_data, config_file, pretty = TRUE, auto_unbox = TRUE)

        old_wd <- setwd(temp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        zip::zip(
          zipfile = file,
          files = c(model_files, "config.json"),
          mode = "cherry-pick"
        )

        showNotification("UMAP model saved successfully!", type = "message")
      }
    )

    # Load UMAP model from file. Stage and validate everything before replacing state.
    observeEvent(input$load_umap_model, {
      req(input$load_umap_model)
      temp_dir <- NULL
      staged_models <- NULL

      tryCatch({
        showModal(modalDialog("Loading UMAP model...", footer = NULL))
        file_path <- input$load_umap_model$datapath
        temp_dir <- tempfile("umap_load_")
        dir.create(temp_dir)
        zip::unzip(file_path, exdir = temp_dir)

        config_file <- file.path(temp_dir, "config.json")
        if (!file.exists(config_file)) stop("Invalid UMAP model file: config.json not found.")
        config_data <- jsonlite::read_json(config_file)

        metadata <- config_data$metadata
        format_name <- if (is.null(metadata$format)) NA_character_ else as.character(metadata$format[[1]])
        schema_version <- if (is.null(metadata$schema_version)) NA_integer_ else {
          suppressWarnings(as.integer(metadata$schema_version[[1]]))
        }
        if (!identical(format_name, "varg_umap") || is.na(schema_version) || schema_version != 2L) {
          stop(paste0(
            "This legacy UMAP bundle cannot be projected reproducibly. ",
            "Recreate and save the model with the current VARG-Tools version."
          ))
        }

        raw_config <- config_data$pipeline_config
        contract <- validate_umap_projection_contract(raw_config)
        if (!contract$valid) stop(paste(contract$errors, collapse = "\n"))

        raw_config$umap$columns_used <- contract$columns
        raw_config$umap$dimensions <- contract$dimensions
        raw_config$umap$n_components <- contract$dimensions
        raw_config$umap$n_input_cols <- length(contract$columns)
        raw_config$umap$projection_contract <- list(
          version = 2L,
          type = "direct_columns",
          columns = contract$columns,
          require_complete = TRUE
        )
        if (!is.null(raw_config$umap$y_vars)) {
          raw_config$umap$y_vars <- as.character(unlist(raw_config$umap$y_vars, use.names = FALSE))
        }

        staged_models <- load_umap_model_collection(temp_dir, raw_config$umap)

        cleanup_loaded_external_model()
        loaded_external_model$models <- staged_models
        loaded_external_model$config <- raw_config
        loaded_external_model$temp_dir <- temp_dir
        staged_models <- NULL
        temp_dir <- NULL

        showNotification(
          "UMAP model loaded. Prepare the required columns, then select Project Using Loaded Model.",
          type = "message",
          duration = 7
        )
      }, error = function(e) {
        if (!is.null(staged_models)) unload_umap_model_collection(staged_models)
        if (!is.null(temp_dir) && dir.exists(temp_dir)) {
          unlink(temp_dir, recursive = TRUE, force = TRUE)
        }
        showNotification(paste("Error loading UMAP model:", e$message), type = "error", duration = 12)
      }, finally = {
        removeModal()
      })
    })

    # Project current data using every model contained in the loaded bundle.
    observeEvent(input$project_loaded_model, {
      req(rv$data, loaded_external_model$models, loaded_external_model$config)
      show_loading("Validating and projecting data...")
      on.exit(hide_loading())

      tryCatch({
        config <- loaded_external_model$config
        projection <- prepare_saved_umap_projection(rv$data, config)
        if (!is.null(projection$error)) stop(projection$error)
        if (length(projection$warnings) > 0) {
          showNotification(paste(projection$warnings, collapse = "\n"), type = "warning", duration = 8)
        }

        models <- normalize_umap_model_collection(loaded_external_model$models, config$umap)
        dimensions <- normalize_umap_dimensions(config$umap)
        expected_names <- paste0(dimensions, "d")
        if (!setequal(names(models), expected_names)) {
          stop("Loaded UMAP models do not match the saved dimension metadata.")
        }

        embeddings <- list()
        for (name in expected_names) {
          dimension <- as.integer(sub("d$", "", name))
          result <- as.matrix(uwot::umap_transform(projection$data, models[[name]]))
          if (nrow(result) != length(projection$valid_idx) || ncol(result) != dimension) {
            stop(paste0("Loaded ", dimension, "D UMAP model returned an unexpected result shape."))
          }
          embeddings[[name]] <- result
        }

        projected_data <- rv$data[, !grepl("^UMAP_loaded_", names(rv$data)), drop = FALSE]
        if ("1d" %in% names(embeddings)) {
          projected_data$UMAP_loaded_1D <- NA_real_
          projected_data$UMAP_loaded_1D[projection$valid_idx] <- embeddings[["1d"]][, 1]
        }
        if ("2d" %in% names(embeddings)) {
          projected_data$UMAP_loaded_2D_1 <- NA_real_
          projected_data$UMAP_loaded_2D_2 <- NA_real_
          projected_data$UMAP_loaded_2D_1[projection$valid_idx] <- embeddings[["2d"]][, 1]
          projected_data$UMAP_loaded_2D_2[projection$valid_idx] <- embeddings[["2d"]][, 2]
        }

        rv$data <- projected_data
        rv$umap_mode_ran <- "loaded"
        bump_data_generation()
        focus_umap_result_tab("loaded", dimensions)
        showNotification(
          sprintf(
            "Projected %d of %d rows using the loaded %s UMAP model.",
            length(projection$valid_idx), nrow(rv$data), umap_dimension_label(config$umap)
          ),
          type = "message"
        )
      }, error = function(e) {
        showNotification(
          paste0(
            "Could not project this dataset. The exact prepared training columns are required. ",
            "Technical detail: ", e$message
          ),
          type = "error",
          duration = 15
        )
      })
    })

    # --- Population Definition ---
    # Helper function to rebuild pop_styles from existing population values
    rebuild_pop_styles <- function(populations) {
      # Get unique population names (excluding NA)
      unique_pops <- unique(na.omit(populations))
      unique_pops <- unique_pops[unique_pops != ""]
      
      # Start with base styles (Unassigned always present)
      new_colors <- list("Unassigned" = "grey50")
      new_shapes <- list("Unassigned" = "circle")
      
      # Add styles for each existing population (excluding Unassigned)
      other_pops <- setdiff(unique_pops, "Unassigned")
      # Sort for deterministic color assignment (same order always gets same colors)
      other_pops <- sort(other_pops)
      for (i in seq_along(other_pops)) {
        pop_name <- other_pops[i]
        # Use consistent 1-based indexing to avoid always starting at color 2
        color_idx <- ((i - 1) %% length(P50)) + 1
        shape_idx <- ((i - 1) %% length(VALID_SHAPES)) + 1
        new_colors[[pop_name]] <- P50[color_idx]
        new_shapes[[pop_name]] <- VALID_SHAPES[shape_idx]
      }
      
      list(colors = new_colors, shapes = new_shapes)
    }
    
    observeEvent(rv$data,
      {
        if ("population" %in% names(rv$data)) {
          rv$data$population <- as.character(rv$data$population)
          rv$data$population[is.na(rv$data$population)] <- "Unassigned"
          
          # Rebuild pop_styles from existing population values
          # This ensures loaded data with pre-existing populations displays correctly
          rebuilt_styles <- rebuild_pop_styles(rv$data$population)
          pop_styles(rebuilt_styles)
          
          # Update the label dropdown with existing population labels (include 'Unassigned' as an option)
          existing_labels <- unique(rv$data$population)
          # Only pre-select a meaningful (non-Unassigned) label the user was recently working with
          selected_label <- ""
          if (!is.null(last_feature_label()) && last_feature_label() != "Unassigned" && last_feature_label() %in% existing_labels) {
            selected_label <- last_feature_label()
          }

          updateSelectizeInput(session, "feature_label", choices = existing_labels, selected = selected_label, server = TRUE)
        } else {
          rv$data$population <- "Unassigned"
          # Reset to default styles for fresh data
          pop_styles(list(colors = setNames("grey50", "Unassigned"), shapes = setNames("circle", "Unassigned")))
          # Reset dropdown to empty choices
          updateSelectizeInput(session, "feature_label", choices = NULL, server = TRUE)
        }
      },
      priority = 10
    )

    plot_data_filtered <- reactive({
      req(rv$data, input$popsel_xvar, input$popsel_yvar)
      df <- rv$data
      
      # Validate that the selected columns still exist in the data
      # (they may have been removed if UMAP was re-run with different settings)
      if (!input$popsel_xvar %in% names(df)) {
        showNotification(paste0("Column '", input$popsel_xvar, "' no longer exists in the data. Please select a different X variable."), type = "warning", duration = 5)
        return(NULL)
      }
      if (!input$popsel_yvar %in% names(df)) {
        showNotification(paste0("Column '", input$popsel_yvar, "' no longer exists in the data. Please select a different Y variable."), type = "warning", duration = 5)
        return(NULL)
      }
      
      # Apply finite filter
      df <- df %>% filter(is.finite(.data[[input$popsel_xvar]]) & is.finite(.data[[input$popsel_yvar]]))
      # Apply group filter if selected
      if (!is.null(input$filter_col) && input$filter_col != "" && !is.null(input$filter_val) && input$filter_val != "All") {
        df <- df %>% filter(.data[[input$filter_col]] == input$filter_val)
      }
      df
    })

    observeEvent(input$assign_feature, {
      req(rv$data)
      sel_data <- event_data("plotly_selected", source = "popselect")
      req(sel_data)
      selected_UIDs <- unique(sel_data$key)
      req(length(selected_UIDs) > 0)

      selected_idx <- which(rv$data$UID %in% selected_UIDs)
      req(length(selected_idx) > 0)

      old_labels <- as.character(rv$data$population[selected_idx])
      old_labels[is.na(old_labels) | trimws(old_labels) == ""] <- "Unassigned"
      push_undo_snapshot(
        uids = rv$data$UID[selected_idx],
        old_labels = old_labels
      )

      feature_name <- trimws(input$feature_label)
      if (feature_name == "") feature_name <- "Unassigned"

      # Remember last chosen label so dropdown stays on it after updates
      last_feature_label(feature_name)

      rv$data$population[selected_idx] <- feature_name
      bump_data_generation()

      current_styles <- pop_styles()
      if (!feature_name %in% names(current_styles$colors)) {
        existing_pops <- setdiff(names(current_styles$colors), "Unassigned")
        next_idx <- length(existing_pops) + 1
        current_styles$colors[[feature_name]] <- P50[next_idx %% length(P50) + 1]
        current_styles$shapes[[feature_name]] <- VALID_SHAPES[next_idx %% length(VALID_SHAPES) + 1]
        pop_styles(current_styles)
      }

      showNotification(paste(length(selected_UIDs), "points assigned to:", feature_name), type = "message")

      # Update choices for the selectize input (keep all labels including 'Unassigned')
      updateSelectizeInput(session, "feature_label", choices = unique(rv$data$population), selected = feature_name, server = TRUE)

      # Note: We do NOT clear the filter selection here, so it persists
    })

    observeEvent(input$undo_assign_feature, {
      req(rv$data)
      history <- undo_history()
      req(length(history) > 0)

      snapshot <- history[[length(history)]]
      history <- history[-length(history)]
      undo_history(history)

      restore_idx <- match(snapshot$uids, rv$data$UID)
      valid <- !is.na(restore_idx)
      if (!any(valid)) {
        showNotification("Undo snapshot does not match current data.", type = "warning")
        return()
      }

      rv$data$population[restore_idx[valid]] <- snapshot$old_labels[valid]
      bump_data_generation()

      rebuilt_styles <- rebuild_pop_styles(rv$data$population)
      pop_styles(rebuilt_styles)

      selected_label <- ""
      if (!is.null(last_feature_label()) && last_feature_label() %in% unique(rv$data$population)) {
        selected_label <- last_feature_label()
      }
      updateSelectizeInput(session, "feature_label", choices = unique(rv$data$population), selected = selected_label, server = TRUE)

      showNotification(sprintf("Undid last assignment for %d point(s).", sum(valid)), type = "message")
    })

    observe({
      shinyjs::toggleState("undo_assign_feature", condition = length(undo_history()) > 0)
    })

    # Revert population dropdown - show all non-Unassigned populations
    output$revert_pop_ui <- renderUI({
      req(rv$data)
      pops <- setdiff(unique(rv$data$population), "Unassigned")
      if (length(pops) == 0) return(tags$p(class = "small text-muted", "No populations defined yet."))
      selectInput(ns("revert_pop_select"), "Population to revert:", choices = pops)
    })

    # Revert all points in a population back to Unassigned
    observeEvent(input$revert_population, {
      req(rv$data, input$revert_pop_select)
      pop_name <- input$revert_pop_select
      n_reverted <- sum(rv$data$population == pop_name, na.rm = TRUE)
      
      if (n_reverted > 0) {
        rv$data$population[rv$data$population == pop_name] <- "Unassigned"
        bump_data_generation()
        showNotification(
          sprintf("Reverted %d points from '%s' back to Unassigned.", n_reverted, pop_name),
          type = "message"
        )
        # Update the selectize choices
        updateSelectizeInput(session, "feature_label", choices = unique(rv$data$population), selected = "", server = TRUE)
      } else {
        showNotification("No points found for that population.", type = "warning")
      }
    })

    # Rename population dropdown - show all non-Unassigned populations
    output$rename_pop_ui <- renderUI({
      req(rv$data)
      pops <- setdiff(unique(rv$data$population), "Unassigned")
      if (length(pops) == 0) return(tags$p(class = "small text-muted", "No populations defined yet."))
      selectInput(ns("rename_pop_select"), "Population to rename:", choices = pops)
    })

    # Rename all points in a population to a new name
    observeEvent(input$rename_population, {
      req(rv$data, input$rename_pop_select)
      old_name <- input$rename_pop_select
      new_name <- trimws(input$rename_pop_new)
      
      # Validate new name
      if (nchar(new_name) == 0) {
        showNotification("Please enter a new population name.", type = "warning")
        return()
      }
      if (tolower(new_name) == "unassigned") {
        showNotification("Use the 'Revert' option to set points back to Unassigned.", type = "warning")
        return()
      }
      if (new_name == old_name) {
        showNotification("New name is the same as the current name.", type = "warning")
        return()
      }
      # If target name already exists, merge into it
      existing_pops <- unique(rv$data$population)
      merging <- new_name %in% existing_pops
      
      n_renamed <- sum(rv$data$population == old_name, na.rm = TRUE)
      if (n_renamed > 0) {
        rv$data$population[rv$data$population == old_name] <- new_name
        bump_data_generation()
        # Rebuild styles so the new name gets proper color/shape
        rebuilt_styles <- rebuild_pop_styles(rv$data$population)
        pop_styles(rebuilt_styles)
        # Update the selectize choices
        updateSelectizeInput(session, "feature_label", choices = unique(rv$data$population), selected = new_name, server = TRUE)
        # Clear the text input
        updateTextInput(session, "rename_pop_new", value = "")
        
        if (merging) {
          showNotification(
            sprintf("Merged %d points from '%s' into existing population '%s'.", n_renamed, old_name, new_name),
            type = "message"
          )
        } else {
          showNotification(
            sprintf("Renamed %d points from '%s' to '%s'.", n_renamed, old_name, new_name),
            type = "message"
          )
        }
      } else {
        showNotification("No points found for that population.", type = "warning")
      }
    })

    output$featplot <- renderPlotly({
      df_plot <- plot_data_filtered()
      req(nrow(df_plot) > 0)
      x_var <- req(input$popsel_xvar)
      y_var <- req(input$popsel_yvar)
      
      # Get current pop_styles for consistent ordering
      current_styles <- pop_styles()
      
      # Get all unique populations in the data
      all_pops <- unique(df_plot$population)
      
      # Ensure all populations have colors/shapes (add missing ones on-the-fly)
      for (pop in all_pops) {
        if (!pop %in% names(current_styles$colors)) {
          existing_count <- length(current_styles$colors) - 1  # exclude Unassigned
          color_idx <- (existing_count %% length(P50)) + 1
          shape_idx <- (existing_count %% length(VALID_SHAPES)) + 1
          current_styles$colors[[pop]] <- P50[color_idx]
          current_styles$shapes[[pop]] <- VALID_SHAPES[shape_idx]
        }
      }
      
      # Update pop_styles if we added any new ones
      if (length(current_styles$colors) != length(pop_styles()$colors)) {
        pop_styles(current_styles)
      }
      
      # Create factor with levels matching pop_styles order for consistent coloring
      # Unassigned first, then others in sorted order
      level_order <- c("Unassigned", sort(setdiff(names(current_styles$colors), "Unassigned")))
      # Only keep levels that exist in the data
      level_order <- intersect(level_order, all_pops)
      df_plot$population <- factor(df_plot$population, levels = level_order)

      # Build hover text with optional label column
      hover_col <- input$popsel_hover_label
      if (!is.null(hover_col) && nzchar(hover_col) && hover_col %in% names(df_plot)) {
        df_plot$.hover_text <- paste0(hover_col, ": ", df_plot[[hover_col]], "<br>Pop: ", df_plot$population)
      } else {
        df_plot$.hover_text <- paste0("UID: ", df_plot$UID, "<br>Pop: ", df_plot$population)
      }

      # If gmm_cluster is on an axis, make it a factor for cleaner discrete plotting
      is_discrete_x <- x_var == "gmm_cluster"
      is_discrete_y <- y_var == "gmm_cluster"
      if (is_discrete_x) df_plot[[x_var]] <- factor(df_plot[[x_var]])
      if (is_discrete_y) df_plot[[y_var]] <- factor(df_plot[[y_var]])

      p <- ggplot(df_plot, aes(
        x = .data[[x_var]], y = .data[[y_var]],
        key = UID, color = population, shape = population,
        text = .hover_text
      ))
      
      if (is_discrete_x || is_discrete_y) {
        # Use jitter to prevent overlap when gmm_cluster is on an axis (simulates beeswarm)
        p <- p + geom_jitter(alpha = 0.7, size = 2.5, 
                             width = if(is_discrete_x) 0.25 else 0, 
                             height = if(is_discrete_y) 0.25 else 0)
      } else {
        p <- p + geom_point(alpha = 0.7, size = 2.5)
      }

      p <- p +
        scale_color_manual(values = unlist(current_styles$colors), name = "Population", drop = FALSE) +
        scale_shape_manual(values = unlist(current_styles$shapes), name = "Population", drop = FALSE) +
        theme_minimal() +
        labs(x = x_var, y = y_var)

      ggplotly(p, source = "popselect", tooltip = "text") %>%
        plotly::layout(dragmode = "lasso", hoverlabel = list(align = "left"))
    })

    # --- Downloads ---
    # Helper to prepare data for export (convert text columns to proper types)
    prepare_export_data <- function(df) {
      # Convert gmm_cluster from text to numeric for proper sorting in external tools
      if ("gmm_cluster" %in% names(df)) {
        df$gmm_cluster <- suppressWarnings(as.integer(df$gmm_cluster))
      }

      # User-facing rename for exported files (internal name remains row_excluded)
      if ("row_excluded" %in% names(df)) {
        names(df)[names(df) == "row_excluded"] <- "Excluded"
      }

      # Reorder key column groups for export readability
      original_cols <- if (!is.null(rv$original_cols)) intersect(rv$original_cols, names(df)) else character(0)
      norm_cols <- grep("_norm$", names(df), value = TRUE)
      imp_cols <- grep("_imp$", names(df), value = TRUE)
      pivot_cols <- grep("^pivot_", names(df), value = TRUE)
      clr_cols <- grep("^clr_", names(df), value = TRUE)
      umap_cols <- grep("^UMAP_", names(df), value = TRUE)
      gmm_col <- intersect("gmm_cluster", names(df))
      pop_col <- intersect("population", names(df))
      excluded_col <- intersect("Excluded", names(df))

      ordered_primary <- unique(c(
        original_cols,
        norm_cols,
        imp_cols,
        pivot_cols,
        clr_cols,
        umap_cols,
        gmm_col,
        pop_col
      ))

      remaining_cols <- setdiff(names(df), c(ordered_primary, excluded_col))
      final_order <- c(ordered_primary, remaining_cols, excluded_col)
      final_order <- intersect(final_order, names(df))

      df[, final_order, drop = FALSE]
    }

    create_column_dictionary <- function(export_df) {
      original_cols <- if (!is.null(rv$original_cols)) intersect(rv$original_cols, names(export_df)) else character(0)

      column_step <- function(col_name) {
        if (col_name %in% original_cols) return("Upload (original data)")
        if (grepl("_norm$", col_name)) return("Preprocessing: normalization")
        if (grepl("_imp$", col_name)) return("Preprocessing: imputation")
        if (grepl("^pivot_", col_name)) return("Preprocessing: ILR transformation")
        if (grepl("^clr_", col_name)) return("Preprocessing: CLR transformation")
        if (grepl("^UMAP_", col_name)) return("Projection: UMAP")
        if (identical(col_name, "gmm_cluster")) return("Clustering: GMM")
        if (identical(col_name, "population")) return("Population definition")
        if (identical(col_name, "Excluded")) return("Row filtering")
        return("Other / carried through")
      }

      column_description <- function(col_name) {
        if (col_name %in% original_cols) return("Original column from uploaded dataset.")
        if (grepl("_norm$", col_name)) return("Normalized compositional value (row closed to 100%).")
        if (grepl("_imp$", col_name)) return("Imputed compositional value after handling zeros/missing entries.")
        if (grepl("^pivot_", col_name)) return("Pivot coordinate (ILR-transformed compositional feature).")
        if (grepl("^clr_", col_name)) return("Centered log-ratio (CLR) transformed compositional feature.")
        if (grepl("^UMAP_", col_name)) return("UMAP projection coordinate.")
        if (identical(col_name, "gmm_cluster")) return("Gaussian Mixture Model cluster assignment.")
        if (identical(col_name, "population")) return("User-defined population label.")
        if (identical(col_name, "Excluded")) return("TRUE for rows flagged to be skipped during processing.")
        return("Additional column preserved or generated by workflow.")
      }

      data.frame(
        `Column Name` = names(export_df),
        Description = vapply(names(export_df), column_description, character(1)),
        `Created In Step` = vapply(names(export_df), column_step, character(1)),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }

    output$download_csv <- downloadHandler(
      filename = function() {
        paste0("VARG-Tools_Processed_Data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(rv$data)
        tryCatch(
          {
            write.csv(prepare_export_data(rv$data), file, row.names = FALSE, na = "")
          },
          error = function(e) {
            showNotification(
              paste0(
                "Export failed while creating the CSV file. Check that data is loaded and the destination is writable. ",
                "Technical detail: ", e$message
              ),
              type = "error",
              duration = 12
            )
            stop(e)
          }
        )
      }
    )

    output$download_xlsx <- downloadHandler(
      filename = function() {
        paste0("VARG-Tools_Processed_Data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        req(rv$data)
        tryCatch(
          {
            export_df <- prepare_export_data(rv$data)
            dict_df <- create_column_dictionary(export_df)

            wb <- openxlsx::createWorkbook()
            openxlsx::addWorksheet(wb, "Data")
            openxlsx::writeData(wb, "Data", export_df)
            openxlsx::addWorksheet(wb, "Column Dictionary")
            openxlsx::writeData(wb, "Column Dictionary", dict_df)
            openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
          },
          error = function(e) {
            showNotification(
              paste0(
                "Export failed while creating the XLSX file. Check that data is loaded and the destination is writable. ",
                "Technical detail: ", e$message
              ),
              type = "error",
              duration = 12
            )
            stop(e)
          }
        )
      }
    )

    # --- Processing Report (self-contained HTML) ---
    output$download_report <- downloadHandler(
      filename = function() {
        paste0("VARG-Tools_Processing_Report_", Sys.Date(), ".html")
      },
      content = function(file) {
        req(rv$data)

        # Helper: encode a plot expression as base64 PNG for embedding
        plot_to_base64 <- function(plot_expr, width = 700, height = 450) {
          tmp <- tempfile(fileext = ".png")
          on.exit(unlink(tmp), add = TRUE)
          png(tmp, width = width, height = height, res = 120)
          tryCatch({ eval(plot_expr) }, error = function(e) NULL)
          dev.off()
          if (file.exists(tmp) && file.info(tmp)$size > 0) {
            raw <- readBin(tmp, "raw", file.info(tmp)$size)
            paste0("data:image/png;base64,", jsonlite::base64_enc(raw))
          } else {
            NULL
          }
        }

        # Helper: data frame -> HTML table
        df_to_html <- function(df, caption = NULL) {
          if (is.null(df) || nrow(df) == 0) return("")
          header <- paste0("<th>", htmltools::htmlEscape(names(df)), "</th>", collapse = "")
          rows <- apply(df, 1, function(r) {
            paste0("<tr>", paste0("<td>", htmltools::htmlEscape(as.character(r)), "</td>", collapse = ""), "</tr>")
          })
          cap <- if (!is.null(caption)) paste0("<caption style='caption-side:top;font-weight:bold;margin-bottom:6px;'>", htmltools::htmlEscape(caption), "</caption>") else ""
          paste0("<table>", cap, "<thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
        }

        # Collect sections
        sections <- list()

        # --- Header ---
        source_name <- if (!is.null(rv$source_filename)) rv$source_filename else "Unknown"
        sections <- c(sections, list(paste0(
          "<h1>VARG-Tools Processing Report</h1>",
          "<p class='meta'>Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
          " &bull; Source file: <strong>", htmltools::htmlEscape(source_name), "</strong></p>"
        )))

        # --- Dataset Overview ---
        n_row <- nrow(rv$data)
        n_col <- ncol(rv$data)
        col_names <- names(rv$data)
        generated_cols <- intersect(col_names, c("UID", "gmm_cluster", "gmm_uncertainty", "population",
                                                  grep("^gmm_prob_|^UMAP_|^norm_|^imp_", col_names, value = TRUE)))
        original_cols <- setdiff(col_names, generated_cols)

        sections <- c(sections, list(paste0(
          "<h2>1. Dataset Overview</h2>",
          "<table><tbody>",
          "<tr><td><strong>Rows</strong></td><td>", n_row, "</td></tr>",
          "<tr><td><strong>Total columns</strong></td><td>", n_col, "</td></tr>",
          "<tr><td><strong>Original columns</strong></td><td>", length(original_cols), "</td></tr>",
          "<tr><td><strong>Generated columns</strong></td><td>", length(generated_cols), "</td></tr>",
          "</tbody></table>"
        )))

        # --- Preprocessing ---
        pre <- if (!is.null(rv$pipeline_config)) rv$pipeline_config$preprocessing else NULL
        if (!is.null(pre)) {
          comp_txt <- if (length(pre$comp_cols) > 0) htmltools::htmlEscape(paste(pre$comp_cols, collapse = ", ")) else "<em>None</em>"
          noncomp_txt <- if (length(pre$noncomp_cols) > 0) htmltools::htmlEscape(paste(pre$noncomp_cols, collapse = ", ")) else "<em>None</em>"
          preprocessing_html <- paste0(
            "<h2>2. Preprocessing</h2>",
            "<table><tbody>",
            "<tr><td><strong>Compositional columns</strong></td><td>", comp_txt, "</td></tr>",
            "<tr><td><strong>Non-compositional columns</strong></td><td>", noncomp_txt, "</td></tr>",
            "<tr><td><strong>Imputation</strong></td><td>",
            if (isTRUE(pre$do_impute)) {
              processing_report_imputation_method_label(pre$imputation_method)
            } else {
              "No"
            },
            "</td></tr>",
            "<tr><td><strong>Transform</strong></td><td>",
            if (identical(pre$transform_type, "ilr")) {
              paste0("ILR (pivot", if (!is.null(pre$pivot_var)) paste0(" by ", htmltools::htmlEscape(pre$pivot_var)) else "", ")")
            } else if (identical(pre$transform_type, "clr")) {
              "CLR"
            } else if (isTRUE(pre$do_pivot)) {
              paste0("ILR (pivot", if (!is.null(pre$pivot_var)) paste0(" by ", htmltools::htmlEscape(pre$pivot_var)) else "", ")")
            } else {
              "None"
            },
            "</td></tr>",
            "</tbody></table>"
          )

          if (isTRUE(pre$do_impute)) {
            missingness <- processing_report_imputation_summary(rv$data, pre$comp_cols)
            if (nrow(missingness) > 0L) {
              preprocessing_html <- paste0(
                preprocessing_html,
                "<h3>Values requiring imputation</h3>",
                df_to_html(missingness)
              )

              if (identical(pre$impute_method, "auto")) {
                preprocessing_html <- paste0(
                  preprocessing_html,
                  processing_report_imputation_auto_html(
                    missingness,
                    nrow(rv$data),
                    pre$imputation_method
                  )
                )
              }

              structured_missingness <- processing_report_structured_missingness(rv$data, pre$comp_cols)
              if (nrow(structured_missingness) > 0L) {
                preprocessing_html <- paste0(
                  preprocessing_html,
                  "<div class='caution'><strong>Structured missingness:</strong> The following core/site groups have no measured values for an analyte. Their imputed values are estimated from relationships learned from other rows; neither robust nor standard regression can validate group-specific behavior without measured values.</div>",
                  df_to_html(structured_missingness)
                )
              }
            }
          }

          sections <- c(sections, list(preprocessing_html))
        } else {
          sections <- c(sections, list("<h2>2. Preprocessing</h2><p class='muted'>No preprocessing was performed.</p>"))
        }

        # --- GMM Clustering ---
        gmm_cfg <- if (!is.null(rv$pipeline_config)) rv$pipeline_config$gmm else NULL
        if (!is.null(gmm_cfg) && !is.null(rv$mclust_result)) {
          gmm_html <- paste0(
            "<h2>3. GMM Clustering</h2>",
            "<table><tbody>",
            "<tr><td><strong>Data source</strong></td><td>", htmltools::htmlEscape(as.character(gmm_cfg$data_source)), "</td></tr>",
            "<tr><td><strong>G range tested</strong></td><td>", gmm_cfg$gmin, " &ndash; ", gmm_cfg$gmax, "</td></tr>",
            "<tr><td><strong>Noise initialisation</strong></td><td>", if (isTRUE(gmm_cfg$noise_init)) "Yes" else "No", "</td></tr>",
            "<tr><td><strong>Prior regularisation</strong></td><td>", if (isTRUE(gmm_cfg$use_prior)) "Yes" else "No", "</td></tr>",
            "<tr><td><strong>Best model</strong></td><td>", htmltools::htmlEscape(as.character(gmm_cfg$best_model)), "</td></tr>",
            "<tr><td><strong>Best G (clusters)</strong></td><td>", gmm_cfg$best_G, "</td></tr>",
            "</tbody></table>"
          )

          # BIC plot
          bic_b64 <- plot_to_base64(quote(plot(rv$mclust_result, what = "BIC")))
          if (!is.null(bic_b64)) {
            gmm_html <- paste0(gmm_html, "<h3>BIC Plot</h3><img src='", bic_b64, "' alt='BIC Plot' />")
          }

          # Cluster report table
          report_df <- tryCatch(gmm_report_data(), error = function(e) NULL)
          if (!is.null(report_df) && nrow(report_df) > 0) {
            gmm_html <- paste0(gmm_html, "<h3>Cluster Summary (Mean &plusmn; SD)</h3>", df_to_html(report_df))
          }

          composition_df <- tryCatch(gmm_cluster_composition_data(), error = function(e) NULL)
          gmm_html <- paste0(
            gmm_html,
            processing_report_gmm_composition_html(composition_df, df_to_html)
          )

          sections <- c(sections, list(gmm_html))
        } else {
          sections <- c(sections, list("<h2>3. GMM Clustering</h2><p class='muted'>GMM was not run in this session.</p>"))
        }

        # --- UMAP ---
        umap_cfg <- if (!is.null(rv$pipeline_config)) rv$pipeline_config$umap else NULL
        sections <- c(sections, list(processing_report_umap_html(umap_cfg)))

        # --- Populations ---
        if ("population" %in% names(rv$data)) {
          pop_tbl <- as.data.frame(table(rv$data$population), stringsAsFactors = FALSE)
          names(pop_tbl) <- c("Population", "Count")
          pop_tbl <- pop_tbl[order(-pop_tbl$Count), ]
          sections <- c(sections, list(paste0(
            "<h2>5. Populations</h2>",
            df_to_html(pop_tbl)
          )))
        }

        # --- Assemble HTML ---
        css <- "
          body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; line-height: 1.5; }
          h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
          h2 { color: #2c3e50; border-bottom: 1px solid #ddd; padding-bottom: 6px; margin-top: 30px; }
          h3 { color: #555; margin-top: 18px; }
          table { border-collapse: collapse; width: 100%; margin: 12px 0 20px 0; }
          th, td { border: 1px solid #ddd; padding: 8px 12px; text-align: left; }
          th { background: #f8f9fa; font-weight: 600; }
          tr:nth-child(even) { background: #fafafa; }
          img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; margin: 10px 0; }
          .meta { color: #777; font-size: 0.9em; }
          .muted { color: #999; font-style: italic; }
          .caution { background: #fff6db; border-left: 4px solid #d89b00; padding: 10px 12px; margin: 12px 0; }
          caption { font-size: 0.95em; color: #555; }
        "

        html <- paste0(
          "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'><meta name='viewport' content='width=device-width, initial-scale=1.0'>",
          "<title>VARG-Tools Processing Report</title>",
          "<style>", css, "</style></head><body>",
          paste(sections, collapse = "\n"),
          "<hr><p class='muted' style='text-align:center;'>VARG-Tools v", APP_VERSION, "</p>",
          "</body></html>"
        )

        writeLines(html, file)
      }
    )

    # Return data for other modules
    return(reactive(rv$data))
  })
}
