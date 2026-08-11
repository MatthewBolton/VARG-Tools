# mod_visualization.R
# VARG-Tools Visualization Module

vp_score_pairs <- function(variable_pairs, predictors_with_pairs, label_complete,
                           vp_seed, vp_k, vp_n_cores) {
  library(caret)
  # caret sets this in callr subprocesses, which makes base R reject larger
  # local PSOCK clusters. Worker count is controlled explicitly below.
  Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
  library(dplyr)
  library(foreach)

  env_cap <- suppressWarnings(as.integer(Sys.getenv("VARG_MAX_CORES", unset = "")))
  worker_cap <- if (!is.na(env_cap) && env_cap > 0L) env_cap else 12L
  available <- parallel::detectCores(logical = TRUE)
  if (is.na(available) || available < 2L) {
    vp_n_cores <- 1L
  } else {
    vp_n_cores <- max(1L, min(as.integer(vp_n_cores), available - 1L, worker_cap))
  }

  set.seed(vp_seed)
  cv_index <- caret::createFolds(label_complete, k = 10, returnTrain = TRUE)
  seeds <- vector(mode = "list", length = 11)
  for (i in seq_len(11)) seeds[[i]] <- sample.int(.Machine$integer.max, 1)
  train_control <- trainControl(
    method = "cv",
    number = 10,
    index = cv_index,
    seeds = seeds
  )

  score_pair <- function(pair) {
    pair_predictors <- predictors_with_pairs[, pair, drop = FALSE]
    tryCatch({
      model <- train(
        x = pair_predictors,
        y = label_complete,
        method = "knn",
        trControl = train_control,
        tuneGrid = expand.grid(k = vp_k),
        preProcess = c("center", "scale")
      )
      data.frame(
        Variable_1 = pair[1],
        Variable_2 = pair[2],
        Accuracy = model$results$Accuracy,
        Kappa = model$results$Kappa
      )
    }, error = function(e) {
      data.frame(
        Variable_1 = pair[1],
        Variable_2 = pair[2],
        Accuracy = NA_real_,
        Kappa = NA_real_
      )
    })
  }

  is_shinyapps <- nzchar(Sys.getenv("SHINY_PORT"))
  if (is_shinyapps || vp_n_cores <= 1L) {
    foreach(pair = variable_pairs, .combine = "rbind", .packages = c("caret", "dplyr")) %do% {
      score_pair(pair)
    }
  } else {
    library(doParallel)
    cl <- parallel::makeCluster(vp_n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    foreach(
      pair = variable_pairs,
      .combine = "rbind",
      .packages = c("caret", "dplyr"),
      .export = "score_pair"
    ) %dopar% {
      score_pair(pair)
    }
  }
}

mod_visualization_ui <- function(id) {
  ns <- NS(id)

  tagList(
    bslib::navset_card_underline(
      id = ns("viz_tabs"),

      # ========================================================================
      # TAB 1: Informative Projections
      # ========================================================================
      bslib::nav_panel(
        "Informative Projections",
        icon = icon("chart-line"),
        value = "projections",
        module_banner(
          goal = "Find the best variable pairs for distinguishing between groups in your data.",
          inputs = "A dataset with numeric analyte columns and a categorical group label.",
          outputs = "Ranked list of variable pairs by classification performance (Kappa).",
          why = "Not all oxides are equally diagnostic. This tool tests every possible combination of variables and ranks them by how well they separate your groups, helping you identify the most informative bivariate plots for publication or fingerprinting."
        ),
        layout_columns(
          col_widths = c(4, 8),
          varg_card(
            title = tagList(icon("cogs"), " Configuration"),
            div(
              class = "scrollable-config",
            # Data Upload
            help_box(
              "Variable Pair Finder",
              "This tool finds the best variable pairs for discriminating between groups. It tests all possible combinations of your predictor variables and ranks them by classification performance (<b>Kappa</b> score). Use this to identify diagnostic geochemical ratios or combinations for your samples."
            ),
            div(
              class = "info-box-custom",
              h4(class = "info-box-header", icon("upload"), "Data Upload"),
              div(
                # Data source selector: Upload or use processed data from Processing module
                radioGroupButtons(
                  inputId = ns("vp_data_source"),
                  label = "Data Source:",
                  choices = c("Upload File" = "upload", "Use Processed Data" = "processed"),
                  selected = "upload",
                  justified = TRUE,
                  status = "primary",
                  checkIcon = list(yes = icon("check"))
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'upload'", ns("vp_data_source")),
                  fileInput(ns("vp_file"),
                    label = div("Upload Data File (.xlsx, .xls, or .csv)", help_icon("<strong>Upload data with numeric predictors and a categorical group label.</strong><details><summary>Learn more</summary>Should contain at least 10 complete cases (rows with no missing values in the selected predictors).<br><br><b>Tip:</b> You can use output from the Processing module, or upload independent data.</details>")),
                    accept = c(".xlsx", ".xls", ".csv")
                  ),
                  uiOutput(ns("vp_sheet_select"))
                ),
                uiOutput(ns("vp_processed_hint"))
              ),
              uiOutput(ns("vp_data_summary"))
            ),

            # Variable Selection
            conditionalPanel(
              condition = sprintf("output['%s']", ns("vp_data_loaded")),
              div(
                class = "info-box-custom",
                h4(class = "info-box-header", icon("bullseye"), "Variable Selection"),
                tags$label("Group Label Column:", help_icon("<strong>Select the target group column to classify.</strong><details><summary>Learn more</summary>Select the column containing the groups you want to predict/classify (e.g., volcanic source, tephra unit, population).<br><br>This is the 'target' variable. The tool will find which variable pairs best separate these groups.<br><br><b>Tip:</b> Ensure your labels are 'clean' (not contaminated by detrital or bad analyses), otherwise results will be misleading.</details>")),
                uiOutput(ns("vp_label_select")),
                tags$label("Predictor Variables:", help_icon("<strong>Select at least two numeric predictor variables.</strong><details><summary>Learn more</summary>The tool evaluates <b>all possible pairs</b> and ranks them by how well they discriminate between your groups (via k-NN kappa).<br><br><b>Important:</b> The VP Finder uses the selected columns as supplied and retains complete, finite, non-zero rows. For major-oxide data, preprocess and impute zeros or unmeasured analytes in Processing before running this analysis.</details>")),
                uiOutput(ns("vp_predictors_select")),
                actionLink(ns("select_all_preds"), "Select All Numeric", icon = icon("check-double"))
              )
            ),

            # Analysis Options
            conditionalPanel(
              condition = sprintf("output['%s']", ns("vp_vars_selected")),
              div(
                class = "info-box-custom",
                h4(class = "info-box-header", icon("cogs"), "Analysis Options"),
                p(class = "info-box-text", "Uses complete, finite, non-zero values from the selected columns. Preprocess or impute data before running."),
                div(
                  checkboxInput(ns("vp_use_derived"),
                    div("Create derived features (ratios & sums)", help_icon("<strong>Create ratio and sum features to improve separation power.</strong><details><summary>Learn more</summary>Generates ratio (A/B) and sum (A+B) features from your selected variables.<br><br><b>Default: OFF.</b><br><br><b>Why enable it?</b> Derived features can reveal discriminating relationships that individual oxides miss (e.g., Na\u2082O+K\u2082O separates alkaline from sub-alkaline).<br><br><b>Trade-off:</b> The exhaustive search becomes much larger. Review the workload estimate before running.</details>")),
                    value = FALSE
                  ),
                  uiOutput(ns("vp_runtime_warning"))
                ),
                p(class = "info-box-text", "Generates A/B and A+B combinations from complete cases only."),
                numericInput(ns("vp_seed"), label = div("Random Seed:", help_icon("<strong>Set a seed to reproduce identical results across runs.</strong><details><summary>Learn more</summary>Sets the random number generator seed for reproducibility.<br><br><b>Default: 123.</b> Use the same seed to get identical results across runs.<br><br><b>When to change:</b> If you want to verify that results are stable across different random splits, try a few different values (e.g., 42, 7, 999). Consistent top pairs across seeds indicate robust discrimination.</details>")), value = 123, min = 1),
                selectInput(ns("vp_k"), label = div("Neighbors (k):", help_icon("<strong>How many nearby points the classifier checks before deciding which group a sample belongs to.</strong><details><summary>Learn more</summary>The Variable Pair Finder uses a <b>k-nearest neighbors (KNN)</b> classifier to test how well each pair of variables separates your groups.<br><br><b>How it works:</b> For each sample, KNN looks at the <em>k</em> closest samples in the variable-pair space and assigns the majority group label. Kappa measures how well those predictions match the true labels.<br><br><b>Low k (1\u20133):</b> Sensitive to fine-scale boundaries. Good for detecting subtle distinctions, but more affected by noise or outliers.<br><br><b>High k (5\u20137):</b> Smoother, more stable boundaries. Better for broad group separation, may miss subtle sub-groups.<br><br><b>Default: 3</b>. This is a good balance for most tephra datasets. Try k=1 if groups are very tight; try k=5\u20137 if you have noisy data or overlapping groups.</details>")),
                  choices = c(1, 3, 5, 7), selected = 3),
                uiOutput(ns("vp_cores_select"))
              )
            ),

            # Run Button
            conditionalPanel(
              condition = sprintf("output['%s']", ns("vp_vars_selected")),
              div(
                style = "text-align: center; margin-top: 20px;",
                div(
                  class = "d-flex justify-content-center gap-2",
                  actionButton(ns("vp_run"), "Find Best Variable Pairs",
                    icon = icon("play"), class = "btn-success btn-lg"
                  ),
                  shinyjs::hidden(
                    actionButton(ns("vp_cancel"), "Cancel", class = "btn-danger btn-lg", icon = icon("stop"))
                  )
                ),
                p(class = "info-box-text", "This may take several minutes")
              )
            )
            ) # end scrollable-config
          ),
          varg_card(
            title = tagList(icon("poll"), " Results"),
            conditionalPanel(
              condition = sprintf("!output['%s']", ns("vp_results_available")),
              div(
                style = "text-align: center; padding: 60px 20px; color: #6c757d;",
                icon("chart-line", style = "font-size: 48px; color: #dee2e6;"),
                h3("Results will appear here", style = "color: #adb5bd;"),
                p("Upload data and run analysis to see variable pair rankings")
              )
            ),
            conditionalPanel(
              condition = sprintf("output['%s']", ns("vp_results_available")),
              fluidRow(
                column(
                  6,
                  h4(icon("list-ol"), " Variable Pairs Ranking"),
                  help_box(
                    "Interpreting Results",
                    "Variable pairs are ranked by <b>Kappa</b> (agreement between predicted and actual groups). Higher Kappa = better separation. <b>Kappa > 0.6</b> indicates good discrimination. Click any row to view the scatter plot."
                  ),
                  p(class = "info-box-text", "Ranked by Kappa (higher = better separation). Click row to view plot."),
                  withSpinner(DTOutput(ns("vp_results_table"))),
                  div(
                    style = "margin-top: 15px;",
                    downloadButton(ns("vp_download"), "Download Results (.xlsx)", class = "btn-primary")
                  )
                ),
                column(
                  6,
                  h4(icon("braille"), " Interactive Visualization"),
                  p(class = "info-box-text", "Scatter plot of selected variable pair"),
                  uiOutput(ns("vp_hover_label_ui")),
                  withSpinner(plotlyOutput(ns("vp_scatter"), height = "500px"))
                )
              )
            )
          )
        )
      ),

      # ========================================================================
      # TAB 2: Custom Scatter Plots
      # ========================================================================
      bslib::nav_panel(
        "Custom Scatter Plots",
        icon = icon("braille"), # FA5 compatible - represents scattered points
        value = "scatter",
        module_banner(
          goal = "Create publication-quality bivariate plots with VARG26 or custom reference fields.",
          inputs = "Processed data or uploaded CSV/Excel with numeric variables and grouping columns.",
          outputs = "Customizable scatter plots exportable as PNG, PDF, or SVG at publication resolution.",
          why = "Compare your unknowns or named populations against the VARG26 regional tephra database or your own reference sets. Generate figures ready for journal submission with full control over aesthetics, symbology, and overlays."
        ),
        layout_columns(
          col_widths = c(4, 8),

          # --- Left Column: Configuration ---
          varg_card(
            title = tagList(icon("cogs"), " Configuration"),
            div(
              style = "height: calc(100vh - 220px); overflow-y: auto; padding-right: 10px;",

              # Quick Guide (collapsed)
              tags$details(
                class = "mb-3",
                tags$summary(class = "small text-muted fw-bold", style = "cursor: pointer;", icon("info-circle"), " Quick Guide"),
                div(
                  class = "mt-2 p-2 border rounded bg-light small",
                  tags$ol(
                    tags$li(tags$b("Data Sources"), ": Upload a CSV/Excel file, or pull in your processed data from Step 1."),
                    tags$li(tags$b("Axes & Template"), ": Choose your X and Y variables, set axis labels, adjust bounds, and optionally overlay a TAS classification diagram."),
                    tags$li(tags$b("Reference Overlays"), ": Add published reference datasets (e.g., VARG26) as background density contours for visual comparison."),
                    tags$li(tags$b("Plot Aesthetics"), " (below the plot): Customize point colors, shapes, sizes, and legend appearance."),
                    tags$li(tags$b("Export"), " (below the plot): Download your figure as PNG, PDF, SVG, or export the underlying data as CSV.")
                  )
                )
              ),

              bslib::accordion(
                id = ns("scatter_accordion"),
                open = c("panel_data", "panel_axes"),

                # ── Panel 1: Data Sources ──
                bslib::accordion_panel(
                  title = tagList(icon("crosshairs"), " Data Sources"),
                  value = "panel_data",
                  tags$p(class = "text-muted small mb-2", "Unknowns or samples to plot as individual points."),
              div(
                radioGroupButtons(
                  inputId = ns("cp_data_source"),
                  label = "Data Source:",
                  choices = c("Upload File" = "upload", "Use Processed Data" = "processed"),
                  selected = "upload",
                  justified = TRUE,
                  status = "primary",
                  checkIcon = list(yes = icon("check"))
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'upload'", ns("cp_data_source")),
                  fileInput(ns("cp_file"),
                    label = div("Upload User Data (.xlsx, .xls, or .csv)", help_icon("<strong>Upload CSV/Excel data with numeric plotting variables.</strong><details><summary>Learn more</summary>Should contain numeric variables for plotting and optional grouping columns (for color/shape symbology).<br><br><b>Tip:</b> You can use the 'Import from Processing' button to pull in your current processed dataset directly.</details>")),
                    accept = c(".xlsx", ".xls", ".csv")
                  ),
                  uiOutput(ns("cp_sheet_select"))
                ),
                uiOutput(ns("cp_processed_hint"))
              ),
              uiOutput(ns("cp_data_summary")),

              # Import Variable Pairs (collapsed by default)
              tags$details(
                class = "mt-3 mb-2",
                tags$summary(class = "small text-muted fw-bold", style = "cursor: pointer;", icon("exchange-alt"), " Import Variable Pairs", help_icon("<strong>Import top variable pairs and auto-compute ratio/sum columns.</strong><details><summary>Learn more</summary>Imports the best discriminating variable pairs from Tab 1 and computes their ratios and sums for plotting.<br><br><b>When to use:</b> After running Variable Pair Analysis, click here to explore the top pairs visually. The computed ratio/sum columns will appear in the X/Y axis selectors.</details>")),
                div(
                  class = "mt-2 p-2 border rounded bg-light",
                  p(class = "text-muted small", "Import top-ranked variable pairs from the Informative Projections tab."),
                  fluidRow(
                    column(6, numericInput(ns("cp_n_pairs"), "Top N pairs:", value = 1, min = 1, max = 100, step = 1)),
                    column(6, div(style = "margin-top: 25px;", actionButton(ns("cp_import_pairs"), "Import Pairs", icon = icon("download"), class = "btn-primary btn-sm")))
                  ),
                  selectInput(ns("cp_pair_suffix"),
                      label = tags$span("Use columns with suffix:", help_icon("<strong>Pick which transformed columns to use for imported pairs.</strong><details><summary>Learn more</summary>Select which data transformation to use for computing the imported pairs.<br><br><b>Imputed (_imp)</b> is recommended because it ensures zeros or missing values do not interfere with ratio calculations.<br><br><b>Raw:</b> Uses original untransformed values.<br><br><b>Normalized (_norm):</b> Uses normalized values.</details>")),
                      choices = c(
                        "None (raw)" = "raw",
                        "Normalized (_norm)" = "_norm",
                        "Imputed (_imp)" = "_imp"
                      ),
                      selected = "_imp"
                    ),
                  uiOutput(ns("cp_imported_pairs_summary"))
                )
              ),

              # Data Filtering (collapsed by default)
              conditionalPanel(
                condition = sprintf("output['%s']", ns("cp_data_loaded")),
                tags$details(
                  class = "mt-2 mb-2",
                  tags$summary(class = "small text-muted fw-bold", style = "cursor: pointer;", icon("filter"), " Filter Data"),
                  div(
                    class = "mt-2 p-2 border rounded bg-light",
                    selectInput(ns("cp_filter_col"), label = tags$span("Filter by Column:", help_icon("<strong>Select a categorical column to filter plotted points.</strong><details><summary>Learn more</summary>Select a categorical column to filter which data points appear in the plot.<br><br><b>When to use:</b> Focus on specific subsets of your data (e.g., one site, one volcanic source) without creating separate datasets.<br><br>Select 'None' to show all data.</details>")), choices = NULL),
                    pickerInput(ns("cp_filter_values"), label = tags$span("Select Values:", help_icon("<strong>Select which values to include from the chosen filter column.</strong><details><summary>Learn more</summary>Select which groups/values to include from the filter column. Deselecting items hides them from the plot.<br><br><b>Tip:</b> Use 'Select All' / 'Deselect All' buttons for quick toggling. The search box helps find specific values in long lists.</details>")), choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE))
                  )
                )
              ),

                ), # End Panel 1: Data Sources

                # ── Panel 2: Axes & Template ──
                bslib::accordion_panel(
                  title = tagList(icon("arrows-alt"), " Axes & Template"),
                  value = "panel_axes",
                  conditionalPanel(
                    condition = sprintf("output['%s']", ns("cp_data_loaded")),
                    div(
                  selectInput(ns("cp_x_var"), label = tags$span("X-axis:", help_icon("<strong>Choose the variable for the X-axis.</strong><details><summary>Learn more</summary>UMAP coordinates are auto-selected if available.<br><br><b>Tip:</b> For geochemical comparison plots, use oxide pairs (e.g., SiO\u2082 vs K\u2082O). For population overview, use UMAP coordinates.</details>")), choices = NULL)
                ),
                div(
                  selectInput(ns("cp_y_var"), label = tags$span("Y-axis:", help_icon("<strong>Choose the variable for the Y-axis.</strong><details><summary>Learn more</summary>UMAP coordinates are auto-selected if available.<br><br>Use the same variable space as X (e.g., both oxides, or both UMAP) for interpretable plots.</details>")), choices = NULL)
                ),
                div(
                  textInput(ns("cp_x_label"), label = tags$span("X-axis Label (optional):", help_icon("<strong>Set a custom X-axis label or leave blank for auto.</strong><details><summary>Learn more</summary>Optional custom label for X-axis. Leave blank to use variable name.<br><br><b>Unicode tips (click to keep open, then copy):</b><br>Subscripts: ₀₁₂₃₄₅₆₇₈₉<br>Superscripts: ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻<br><br>Examples: SiO₂, Na₂O+K₂O, Fe²⁺, Al₂O₃</details>")), placeholder = "Auto")
                ),
                div(
                  textInput(ns("cp_y_label"), label = tags$span("Y-axis Label (optional):", help_icon("<strong>Set a custom Y-axis label or leave blank for auto.</strong><details><summary>Learn more</summary>Optional custom label for Y-axis. Leave blank to use variable name.<br><br><b>Unicode tips (click to keep open, then copy):</b><br>Subscripts: ₀₁₂₃₄₅₆₇₈₉<br>Superscripts: ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻<br><br>Examples: SiO₂, Na₂O+K₂O, Fe²⁺, Al₂O₃</details>")), placeholder = "Auto")
                ),
                
                
                    # Axis sub-options (Bounds, Units, Scale, Tick Marks)
                    
                    # Sub-menu: Bounds (min/max)
                    tags$details(
                      class = "mb-2",
                      tags$summary(class = "small text-muted", style = "cursor: pointer;", "Bounds"),
                      div(
                        class = "mt-1 ps-2",
                        fluidRow(
                          column(6, numericInput(ns("cp_x_min"), "X min:", value = NULL)),
                          column(6, numericInput(ns("cp_x_max"), "X max:", value = NULL))
                        ),
                        fluidRow(
                          column(6, numericInput(ns("cp_y_min"), "Y min:", value = NULL)),
                          column(6, numericInput(ns("cp_y_max"), "Y max:", value = NULL))
                        )
                      )
                    ),
                    
                    # Sub-menu: Units (major/minor tick intervals)
                    tags$details(
                      class = "mb-2",
                      tags$summary(class = "small text-muted", style = "cursor: pointer;", "Units"),
                      div(
                        class = "mt-1 ps-2",
                        fluidRow(
                          column(6, numericInput(ns("cp_x_tick_interval"), label = tags$span("X major interval:", help_icon("<strong>Set spacing between major X-axis ticks; leave blank for auto.</strong>")), value = NA, min = 0, step = 1)),
                          column(6, numericInput(ns("cp_y_tick_interval"), label = tags$span("Y major interval:", help_icon("<strong>Set spacing between major Y-axis ticks; leave blank for auto.</strong>")), value = NA, min = 0, step = 1))
                        ),
                        numericInput(ns("cp_minor_tick_n"), label = tags$span("Minor ticks per interval:", help_icon("<strong>Set minor ticks between major ticks (default 4).</strong><details><summary>Learn more</summary>Number of minor tick marks between each major tick.<br><br>Default: 4 (divides each major interval into 5 segments, like Excel).<br><br>Only applies when minor ticks are enabled.</details>")), value = 4, min = 1, max = 20, step = 1, width = "50%")
                      )
                    ),
                    
                    # Sub-menu: Scale (log, reverse)
                    tags$details(
                      class = "mb-2",
                      tags$summary(class = "small text-muted", style = "cursor: pointer;", "Scale"),
                      div(
                        class = "mt-1 ps-2",
                        fluidRow(
                          column(6, checkboxInput(ns("cp_log_x"), label = tags$span("Log₁₀ X", help_icon("<strong>Apply log₁₀ scaling to the X-axis.</strong><details><summary>Learn more</summary><b>Default: OFF.</b><br><br><b>When to use:</b> For highly skewed data or trace element concentrations that span several orders of magnitude.<br><br>Makes patterns in low-concentration ranges visible.</details>")), value = FALSE)),
                          column(6, checkboxInput(ns("cp_log_y"), label = tags$span("Log₁₀ Y", help_icon("<strong>Apply log₁₀ scaling to the Y-axis.</strong><details><summary>Learn more</summary><b>Default: OFF.</b><br><br>Same as Log X; useful for trace elements or any variable spanning orders of magnitude.</details>")), value = FALSE))
                        ),
                        fluidRow(
                          column(6, checkboxInput(ns("cp_reverse_x"), label = tags$span("Reverse X", help_icon("<strong>Reverse the X-axis direction.</strong><details><summary>Learn more</summary>Flip the X-axis direction (high values on the left).<br><br><b>Default: OFF.</b><br><br>Some geochemical conventions plot SiO\u2082 decreasing left-to-right to align with differentiation trends.</details>")), value = FALSE)),
                          column(6, checkboxInput(ns("cp_reverse_y"), label = tags$span("Reverse Y", help_icon("<strong>Reverse the Y-axis direction.</strong><details><summary>Learn more</summary>Flip the Y-axis direction (high values at the bottom).<br><br><b>Default: OFF.</b><br><br>Common for depth profiles where depth increases downward.</details>")), value = FALSE))
                        )
                      )
                    ),
                    
                    # Sub-menu: Tick Marks
                    tags$details(
                      class = "mb-2",
                      tags$summary(class = "small text-muted", style = "cursor: pointer;", "Tick Marks"),
                      div(
                        class = "mt-1 ps-2",
                        fluidRow(
                          column(6, selectInput(ns("cp_major_tick_type"), "Major type:", choices = c("Outside" = "outside", "None" = "none"), selected = "outside")),
                          column(6, selectInput(ns("cp_minor_tick_type"), "Minor type:", choices = c("None" = "none", "Outside" = "outside"), selected = "none"))
                        )
                      )
                    )
              ),

              # TAS Diagram Template (collapsed by default)
              conditionalPanel(
                condition = sprintf("output['%s']", ns("cp_data_loaded")),
                tags$details(
                  class = "mt-3 mb-2",
                  tags$summary(class = "small text-muted fw-bold", style = "cursor: pointer;", icon("border-all"), " TAS (Total Alkali-Silica) Template"),
                  div(
                    class = "mt-2 p-2 border rounded bg-light",
                    checkboxInput(ns("cp_show_tas"), 
                      label = tags$span("Show TAS diagram overlay", help_icon("<strong>Overlay the TAS classification diagram on the plot.</strong>")),
                      value = FALSE
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s']", ns("cp_show_tas")),
                      # Column selection section
                      div(
                        tags$h6("Select TAS Columns (with optional suffixes like _norm, _imp):", class = "text-muted mt-3 mb-2"),
                        fluidRow(
                          column(4, uiOutput(ns("cp_tas_sio2_selector"))),
                          column(4, uiOutput(ns("cp_tas_na2o_selector"))),
                          column(4, uiOutput(ns("cp_tas_k2o_selector")))
                        ),
                        actionButton(ns("cp_assign_tas_axes"), "Assign as X & Y Axes", class = "btn-sm btn-primary mt-2",
                          title = "Set the selected columns as X-axis (SiO₂) and Y-axis (Na₂O+K₂O), creating the total alkali column if needed."
                        ),
                        hr()
                      ),
                      # Display options
                      checkboxInput(ns("cp_tas_fields"), 
                        label = tags$span("Show field boundaries", help_icon("<strong>Show TAS field boundary lines on the plot.</strong><details><summary>Learn more</summary>Display the TAS classification field boundaries as lines on the plot.<br><br><b>When to turn OFF:</b> If you only want labels without the boundary lines, or if the lines clutter a dense dataset.</details>")),
                        value = TRUE
                      ),
                      checkboxInput(ns("cp_tas_labels"), 
                        label = tags$span("Show field labels", help_icon("<strong>Show TAS rock-type abbreviations within each field.</strong><details><summary>Learn more</summary>Display rock type abbreviations in each TAS field (e.g., B=Basalt, BA=Basaltic Andesite, R=Rhyolite).<br><br><b>When to turn OFF:</b> If labels overlap with your data points or reference fields.</details>")),
                        value = TRUE
                      )
                    )
                  )
                )
              ), # End Panel 2: Axes & Template (TAS CP close)
              ), # End Panel 2 accordion_panel

                # ── Panel 3: Reference Overlays ──
                bslib::accordion_panel(
                  title = tagList(icon("layer-group"), " Reference Overlays"),
                  value = "panel_overlays",

              # VARG26 Reference Data
              conditionalPanel(
                condition = sprintf("output['%s']", ns("cp_data_loaded")),
                tags$h5(icon("database"), " VARG26 Reference (Density Fields)", class = "text-primary border-bottom pb-2 mb-3 mt-4"),
                tags$p(class = "text-muted small mb-2", "Reference compositions displayed as HDR contours."),
                checkboxInput(ns("cp_use_VARG26"), 
                  label = tags$span("Show VARG26 reference fields", help_icon("<strong>Overlay VARG26 HDR contours for reference comparison.</strong><details><summary>Learn more</summary>Overlay Highest Density Region (HDR) contours from the VARG26 reference tephra database.<br><br><b>Default: OFF.</b><br><br><b>When to use:</b> To compare your unknown samples against known volcanic sources. The contours show where reference populations cluster in your chosen variable space.<br><br><b>Works best with:</b> imputed oxide pairs or UMAP coordinates.</details>")),
                  value = FALSE
                ),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("cp_use_VARG26")),
                  uiOutput(ns("cp_VARG26_selector")),
                  selectInput(ns("cp_VARG26_filter_col"), label = tags$span("Group by:", help_icon("<strong>Choose grouping level for VARG26 HDR contours.</strong><details><summary>Learn more</summary>Choose how to group the VARG26 reference data for HDR display:<br><br><b>Region Label</b>: Broad geographic regions (e.g., Alaska, Kamchatka, Japan).<br><br><b>Volcanic source</b>: Individual volcanoes (e.g., Katmai, Aniakchak).<br><br><b>Tephra</b>: Specific eruption units.<br><br><b>Tip:</b> Start broad (Region) then narrow to specific sources or tephras.</details>")),
                      choices = c("Region Label", "Volcanic source", "Tephra"),
                      selected = "Volcanic source"
                    ),
                  selectInput(
                    ns("cp_hdr_method"),
                    label = tags$span("Density estimator:", help_icon("<strong>Choose how VARG26 density is estimated.</strong><details><summary>Learn more</summary><b>Kernel density estimate</b> preserves the existing smooth HDR fields. <b>Frequency polygon</b> and <b>Histogram</b> use binned density estimates while retaining HDR probability envelopes.</details>")),
                    choices = c(
                      "Kernel density estimate" = "kde",
                      "Frequency polygon" = "freqpoly",
                      "Histogram" = "histogram"
                    ),
                    selected = "kde"
                  ),
                  sliderInput(ns("cp_hdr_prob"), label = tags$span("Envelope Size:", help_icon("<strong>Set HDR probability to control envelope size.</strong><details><summary>Learn more</summary>Controls the <b>HDR probability</b> (Highest Density Region level).<br><br><b>Default: 0.6.</b><br><br>Higher values (0.8–0.95) show larger, more inclusive envelopes around reference data.<br><br>Lower values (0.3–0.5) show tighter core regions.<br><br><b>Typical range: 0.5–0.8.</b><br><br><b>When to adjust:</b> Increase if your unknowns fall just outside a reference field; decrease if fields overlap too much.</details>")), min = 0.1, max = 0.99, value = 0.6, step = 0.05),
                  conditionalPanel(
                    condition = sprintf("input['%s'] === 'kde'", ns("cp_hdr_method")),
                    sliderInput(ns("cp_hdr_adjust"), label = tags$span("Smoothness:", help_icon("<strong>Set KDE smoothness for VARG26 HDR contours.</strong><details><summary>Learn more</summary>Controls the <b>KDE bandwidth</b> (kernel density smoothing).<br><br><b>Default: 1.0.</b><br><br><b>Lower values (0.5–0.8):</b> Show finer detail and tighter contours, but can be noisy with sparse data.<br><br><b>Higher values (1.2–2.0):</b> Show broader, smoother patterns.<br><br><b>When to adjust:</b> Increase if contours look jagged or fragmented; decrease if they're too smooth and merge distinct sub-populations.</details>")), min = 0.1, max = 2.0, value = 1.0, step = 0.1),
                    numericInput(ns("cp_hdr_n"), label = tags$span("Grid resolution:", help_icon("<strong>Set the KDE evaluation-grid resolution.</strong><details><summary>Learn more</summary>Controls how finely the KDE surface is evaluated before the HDR polygon is drawn.<br><br><b>Default: 100.</b><br><br>Increase this if a KDE boundary looks visibly angular. This changes drawing resolution, not bandwidth or the underlying observations.</details>")), value = 100, min = 50, max = 300, step = 25)
                  ),
                  conditionalPanel(
                    condition = sprintf("input['%s'] !== 'kde'", ns("cp_hdr_method")),
                    checkboxInput(ns("cp_hdr_auto_bins"), label = tags$span("Choose bin count automatically", help_icon("<strong>Use ggdensity's data-driven default bin count.</strong><details><summary>Learn more</summary>Turn this off to set the number of bins along each axis yourself. More bins can show finer structure but may fragment fields for small groups.</details>")), value = TRUE),
                    conditionalPanel(
                      condition = sprintf("!input['%s']", ns("cp_hdr_auto_bins")),
                      numericInput(ns("cp_hdr_bins"), label = tags$span("Bins per axis:", help_icon("<strong>Set the frequency-polygon or histogram resolution.</strong><details><summary>Learn more</summary>The same number of bins is used on both axes.<br><br><b>Starting value: 15.</b><br><br>Increase gradually when polygons look too coarse. Very high values can create fragmented fields when a group contains few observations.</details>")), value = 15, min = 4, max = 60, step = 1)
                    )
                  )
                )
              ),

              # Custom Reference Data
              conditionalPanel(
                condition = sprintf("output['%s']", ns("cp_data_loaded")),
                tags$h5(icon("folder-open"), " Custom Reference (Density Fields)", class = "text-primary border-bottom pb-2 mb-3 mt-4"),
                tags$p(class = "text-muted small mb-2", "Your own reference data displayed as HDR contours."),
                fileInput(ns("cp_ref_file"), 
                  label = div("Upload Reference Data (.xlsx, .xls, or .csv)", help_icon("<strong>Upload custom reference data for HDR overlay.</strong><details><summary>Learn more</summary>Upload your own reference dataset to overlay as HDR contours alongside your main data.<br><br>Must have the same column names for the variables you're plotting.<br><br><b>When to use:</b> If you have a custom reference database, previously published data, or lab standards that aren't in VARG26.</details>")),
                  accept = c(".xlsx", ".xls", ".csv")
                ),
                actionButton(ns("cp_ref_use_processed"), "Import from Processing", 
                  icon = icon("file-import"), class = "btn-outline-primary btn-sm mb-2 w-100"),
                uiOutput(ns("cp_ref_sheet_select")),
                conditionalPanel(
                  condition = sprintf("output['%s']", ns("cp_ref_loaded")),
                  selectInput(ns("cp_ref_filter_col"), 
                    label = tags$span("Group by:", help_icon("<strong>Choose the grouping column for custom HDR contours.</strong><details><summary>Learn more</summary>Choose which column in your reference data provides the group labels for HDR contours.<br><br>Each unique value becomes a separate density field on the plot.<br><br><b>Example:</b> If your reference data has a 'source' column with volcano names, select it here to display each volcano as a separate contour.</details>")),
                    choices = NULL),
                  uiOutput(ns("cp_ref_selector")),
                  div(
                    style = "margin-top: 10px;",
                    selectInput(
                      ns("cp_ref_hdr_method"),
                      label = tags$span("Density estimator:", help_icon("<strong>Choose how custom-reference density is estimated.</strong>")),
                      choices = c(
                        "Kernel density estimate" = "kde",
                        "Frequency polygon" = "freqpoly",
                        "Histogram" = "histogram"
                      ),
                      selected = "kde"
                    ),
                    sliderInput(ns("cp_ref_hdr_prob"), label = tags$span("Envelope Size:", help_icon("<strong>Set HDR probability for custom reference contours.</strong>")),
                      value = 0.8, min = 0.1, max = 0.95, step = 0.05
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] === 'kde'", ns("cp_ref_hdr_method")),
                      sliderInput(ns("cp_ref_hdr_adjust"), label = tags$span("Smoothness:", help_icon("<strong>Set KDE smoothness for custom reference contours.</strong>")),
                        value = 1.0, min = 0.1, max = 3, step = 0.1
                      ),
                      numericInput(ns("cp_ref_hdr_n"), label = tags$span("Grid resolution:", help_icon("<strong>Set the KDE evaluation-grid resolution.</strong>")), value = 100, min = 50, max = 300, step = 25)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] !== 'kde'", ns("cp_ref_hdr_method")),
                      checkboxInput(ns("cp_ref_hdr_auto_bins"), label = tags$span("Choose bin count automatically", help_icon("<strong>Use ggdensity's data-driven default bin count.</strong>")), value = TRUE),
                      conditionalPanel(
                        condition = sprintf("!input['%s']", ns("cp_ref_hdr_auto_bins")),
                        numericInput(ns("cp_ref_hdr_bins"), label = tags$span("Bins per axis:", help_icon("<strong>Set the frequency-polygon or histogram resolution.</strong>")), value = 15, min = 4, max = 60, step = 1)
                      )
                    )
                  )
                )
              ),

              # HDR Labels
              conditionalPanel(
                condition = sprintf(
                  "input['%s'] && output['%s'] || output['%s']",
                  ns("cp_use_VARG26"), ns("cp_VARG26_selected_any"), ns("cp_ref_selected_any")
                ),
                tags$h5(icon("tag"), " HDR Labels", class = "text-primary border-bottom pb-2 mb-3 mt-4"),
                checkboxInput(ns("cp_label_hdr"), label = tags$span("Label HDR centers", help_icon("<strong>Label each HDR field at its density peak.</strong><details><summary>Learn more</summary>Add text labels at the density peak of each reference field.<br><br><b>Default: ON.</b><br><br>Helps identify which contour belongs to which group.<br><br><b>When to turn OFF:</b> If labels overlap or clutter a busy plot.</details>")), value = TRUE),
                checkboxInput(ns("cp_hide_hdr_legend"), label = tags$span("Hide HDR from legend", help_icon("<strong>Hide HDR entries from the legend.</strong><details><summary>Learn more</summary>Removes HDR contour entries from the legend to reduce clutter (the contours still appear on the plot).<br><br><b>Default: ON.</b><br><br><b>When to turn OFF:</b> If you need the legend to identify which contour color belongs to which reference group.</details>")), value = TRUE),
                sliderInput(ns("cp_label_size"), label = tags$span("Label size:", help_icon("<strong>Adjust the font size of HDR labels.</strong><details><summary>Learn more</summary><b>When to adjust:</b> Increase for presentations or if labels are hard to read; decrease if labels overlap with data points.</details>")),
                    min = 2, max = 6, value = 3, step = 0.5
                  )
              ),

              # HDR Color Customization
              conditionalPanel(
                condition = sprintf(
                  "input['%s'] && output['%s'] || output['%s']",
                  ns("cp_use_VARG26"), ns("cp_VARG26_selected_any"), ns("cp_ref_selected_any")
                ),
                tags$h5(icon("palette"), " HDR Color Customization", class = "text-primary border-bottom pb-2 mb-3 mt-4"),
                checkboxInput(ns("cp_use_custom_hdr_colors"),
                    label = tagList(icon("palette"), " Use Custom HDR Colors", help_icon("<strong>Use custom colors for HDR reference groups.</strong><details><summary>Learn more</summary>Override the default HDR fill colors for VARG26 and custom reference groups.<br><br><b>Default: OFF</b> (automatic palette).<br><br><b>When to turn ON:</b> If you need specific colors for publication consistency, colorblind accessibility, or to match a co-author's figure style.</details>")),
                    value = FALSE
                  ),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("cp_use_custom_hdr_colors")),
                  uiOutput(ns("cp_hdr_color_panel"))
                )
              ),
                ) # End Panel 3: Reference Overlays
              ) # End accordion
            ) # End scrollable div
          ),

          # --- Right Column: Plot & Settings ---
          div(
            varg_card(
              title = tagList(icon("braille"), " Interactive Plot"),
              withSpinner(plotOutput(ns("cp_main_plot"), height = "600px"))
            ),
            layout_columns(
              col_widths = c(12, 12),
              varg_card(
                title = tagList(icon("download"), " Export Settings"),
                fluidRow(
                  column(3, numericInput(ns("cp_export_width"), label = div("Width (cm):", help_icon("<strong>Set exported plot width in centimeters.</strong><details><summary>Learn more</summary>Width of the exported plot in centimeters.<br><br><b>Default: 18 cm</b> (~7 inches, standard single-column journal width).<br><br><b>Common sizes:</b> 8–9 cm for half-column, 18 cm for full-column, 25–30 cm for full-page or poster figures.</details>")), value = 18, min = 5, max = 50, step = 1)),
                  column(3, numericInput(ns("cp_export_height"), label = div("Height (cm):", help_icon("<strong>Set exported plot height in centimeters.</strong><details><summary>Learn more</summary>Height of the exported plot in centimeters.<br><br><b>Default: 14 cm.</b><br><br>Adjust to match your desired aspect ratio.<br><br>Square plots (equal W×H) work well for bivariate scatter plots.</details>")), value = 14, min = 5, max = 50, step = 1)),
                  column(3, selectInput(ns("cp_export_dpi"), label = div("DPI:", help_icon("<strong>Set export resolution in dots per inch (DPI).</strong><details><summary>Learn more</summary>Dots Per Inch for the exported image.<br><br><b>Default: 300</b> (standard for print/journal submission).<br><br><b>150:</b> Faster, smaller files for drafts.<br><br><b>600:</b> High resolution for large-format printing or detailed figures.<br><br><b>Note:</b> Higher DPI = larger file size.</details>")), choices = c("150 (Draft)" = 150, "300 (Print)" = 300, "600 (High Res)" = 600), selected = 300)),
                  column(
                    3, br(),
                    fluidRow(
                      column(3, downloadButton(ns("cp_download_png"), "PNG", class = "btn-primary btn-sm")),
                      column(3, downloadButton(ns("cp_download_pdf"), "PDF", class = "btn-primary btn-sm")),
                      column(3, downloadButton(ns("cp_download_svg"), "SVG", class = "btn-primary btn-sm")),
                      column(3, downloadButton(ns("cp_download_data"), "CSV", class = "btn-info btn-sm"))
                    )
                  )
                )
              ),
              bslib::accordion(
                id = ns("scatter_aesthetics"),
                open = FALSE,
                bslib::accordion_panel(
                  title = tagList(icon("palette"), " Plot Aesthetics"),
                  value = "aesthetics",
                # --- Plot Theme ---
                selectInput(ns("cp_theme"), label = tags$span("Plot theme:", help_icon("<strong>Choose the ggplot theme for this plot.</strong><details><summary>Learn more</summary><b>bw</b>: White background with gridlines (default).<br><br><b>Classic</b>: Clean axes, no gridlines; publication ready.<br><br><b>Minimal</b>: No axis lines, light gridlines.<br><br><b>Base</b>: R base-graphics look with serif fonts.<br><br><b>Tufte</b>: Minimal ink, in the style of Edward Tufte.</details>")),
                  choices = c("Black & White" = "bw", "Classic" = "classic", "Minimal" = "minimal", "Light" = "light", "Dark" = "dark", "Linedraw" = "linedraw", "Gray" = "gray", "Base" = "base", "Tufte" = "tufte", "Few" = "few", "Void (no axes)" = "void"),
                  selected = "bw"
                ),
                fluidRow(
                  column(6, pickerInput(ns("cp_color_var"), label = tags$span("Color by:", help_icon("<strong>Choose which column controls point color.</strong><details><summary>Learn more</summary>Choose a categorical column to color-code your data points.<br><br><b>'Fixed Color':</b> All points same color (useful for overlaying on reference fields).<br><br><b>population, gmm_cluster, sample_id:</b> Distinguish groups visually.<br><br><b>Tip:</b> Color by population first to verify assignments, then switch to other variables to check for patterns.</details>")), choices = c("Fixed Color"), selected = "Fixed Color", options = list(container = "body"))),
                  column(6, pickerInput(ns("cp_shape_var"), label = tags$span("Shape by:", help_icon("<strong>Choose which column controls point shape.</strong><details><summary>Learn more</summary>Choose a categorical column for marker shapes.<br><br><b>When to use:</b> When you want to encode a second grouping variable alongside color (e.g., color by population, shape by site).<br><br><b>Limit:</b> ~6 distinct shapes before recycling, so this works best with few categories.</details>")), choices = c("Fixed Shape"), selected = "Fixed Shape", options = list(container = "body")))
                ),
                uiOutput(ns("cp_legend_title_ui")),
                checkboxInput(ns("cp_hide_legend"),
                  label = tags$span("Hide legend entirely", help_icon("<strong>Hide the plot legend entirely.</strong><details><summary>Learn more</summary>Removes the entire legend from the plot.<br><br><b>Default: OFF.</b><br><br><b>When to use:</b> For clean publication figures where the legend is redundant (e.g., you describe symbology in the figure caption), or when the legend overlaps important data.</details>")),
                  value = FALSE
                ),
                uiOutput(ns("cp_continuous_palette_ui")),

                # --- Group Customization Panel ---
                div(
                  style = "margin-top: 10px;",
                  checkboxInput(ns("cp_use_custom_symbology"),
                    label = tags$span(tagList(icon("palette"), " Use Custom Symbology"), help_icon("<strong>Use custom colors and shapes for each group.</strong><details><summary>Learn more</summary>Override default colors and shapes for each group.<br><br><b>Default: OFF</b> (automatic palette).<br><br><b>When to turn ON:</b> For publication figures where you need exact control over symbology, such as matching colors to a co-author's scheme or ensuring colorblind accessibility.</details>")),
                    value = FALSE
                  ),
                  uiOutput(ns("cp_custom_symbology_panel"))
                ),
                splitLayout(
                  cellWidths = c("50%", "50%"),
                  sliderInput(ns("cp_point_size"), label = tags$span("Point Size:", help_icon("<strong>Set the diameter of plotted points.</strong><details><summary>Learn more</summary>Controls the diameter of data points.<br><br><b>Default: 3.</b><br><br><b>Larger (4–6):</b> Better for presentations and sparse datasets.<br><br><b>Smaller (1–2):</b> Better for dense datasets with many overlapping points.<br><br><b>Tip:</b> Pair with transparency to reveal density patterns.</details>")), min = 1, max = 8, value = 3, step = 0.5),
                  sliderInput(ns("cp_point_alpha"), label = tags$span("Transparency:", help_icon("<strong>Set point transparency (opacity).</strong><details><summary>Learn more</summary>Controls point opacity.<br><br><b>Default: 0.8.</b><br><br><b>Lower values (0.3–0.5):</b> Reveal overlapping clusters, so denser areas appear darker.<br><br><b>Higher values (0.8–1.0):</b> Points are more visible individually.<br><br><b>When to adjust:</b> If you see a solid blob of color, reduce transparency to see internal structure.</details>")), min = 0.1, max = 1, value = 0.8, step = 0.1)
                ),
                sliderInput(ns("cp_stroke_width"), label = tags$span("Stroke Width:", help_icon("<strong>Set the outline thickness of points.</strong><details><summary>Learn more</summary>Controls the outline thickness of data points.<br><br><b>Default: 0.5.</b><br><br><b>Thicker (1–2):</b> Makes hollow symbols (circles, triangles, crosses) easier to see.<br><br><b>Thinner (0–0.3):</b> Reduces visual weight for filled symbols.<br><br><b>Tip:</b> Increase stroke for shapes 0–6 and 21–25 (which have visible outlines) or for crosses/X's that rely entirely on stroke.</details>")), min = 0, max = 3, value = 0.5, step = 0.25),
                # Software Citation Option
                div(
                  style = "margin-top: 10px; border-top: 1px solid #dee2e6; padding-top: 10px;",
                  checkboxInput(ns("cp_include_citation"),
                    label = tags$span(tagList(icon("quote-right"), " Include software citation on plot"), help_icon("<strong>Add software citation text to exported plots.</strong><details><summary>Learn more</summary>When checked, adds a citation annotation (app version and VARG UMAP version) to the bottom-right corner of exported plots.<br><br><b>Recommended</b> for any figure that may be published or shared, to ensure reproducibility and proper attribution.</details>")),
                    value = TRUE
                  )
                )
                )
              ),
              # Data Citations Panel
              varg_card(
                title = tagList(icon("book"), " Data Citations", help_icon("<strong>View citations for displayed VARG26 reference data.</strong><details><summary>Learn more</summary>Shows data source citations for the currently displayed VARG26 reference data.<br><br><b>Copy these</b> to reference the original publications in your work.<br><br>Citations update dynamically based on which reference groups are currently shown on the plot.</details>")),
                div(
                  style = "max-height: 200px; overflow-y: auto;",
                  uiOutput(ns("cp_data_citations"))
                )
              )
            )
          )
        )
      ),

      # ========================================================================
      # TAB 3: Stratigraphic Correlation
      # ========================================================================
      bslib::nav_panel(
        "Stratigraphic Correlation",
        icon = icon("layer-group"),
        value = "correlation",
        module_banner(
          goal = "Visually correlate depth-series between two cores and warp target depths to a reference scale.",
          inputs = "CSV/Excel with CoreID, Sample name, Depth, and a correlation variable (e.g., UMAP_1D).",
          outputs = "Tie-point table, warped depth data (CSV), and tie-point JSON for the Chronology module.",
          why = "When you have tephra data from two cores and want to establish stratigraphic links between them, for example by matching tephras across sites to build a regional tephrochronological framework. UMAP_1D makes an ideal correlation variable because it captures the full multivariate geochemical signature in a single dimension."
        ),
        layout_columns(
          col_widths = c(4, 8),

          # --- LEFT COLUMN: CONTROLS ---
          div(
            class = "sc-controls-column",

            # 1. Configuration Card
            varg_card(
              title = "1. Configuration",

              # Instructional Help Box
              help_box(
                title = "Stratigraphic Alignment",
                content = "Align two stratigraphic sequences (cores) using tie-points.
                <ul>
                  <li><b>Reference Core:</b> The master sequence that remains fixed.</li>
                  <li><b>Target Core:</b> The sequence to be warped/stretched to match the Reference.</li>
                  <li><b>Tie Points:</b> Identify equivalent layers in both cores to define the alignment.</li>
                </ul>
                <hr style='margin: 10px 0;'>
                <b>Recommended Workflow:</b><br>
                <ol style='padding-left: 20px; margin-bottom: 0;'>
                  <li>Set up the tie point interface by selecting your X variable (Corr Var)</li>
                  <li>Assign tie points manually or using 'Click to Add' functionality</li>
                  <li>Click 'Apply Warping Model' to align the target core</li>
                  <li>Inspect the Warped Alignment and Warp Fit for implausible compression, expansion, or extrapolation</li>
                  <li>Repeat with secondary correlation variables if needed</li>
                  <li>Proceed to Chronology module's 'Model Linker' to import tie points</li>
                </ol>"
              ),

              radioGroupButtons(
                inputId = ns("sc_data_source"),
                label = "Data Source:",
                choices = c("Upload File" = "upload", "Use Processed Data" = "processed"),
                selected = "upload",
                justified = TRUE,
                status = "primary",
                size = "sm",
                checkIcon = list(yes = icon("check"))
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'upload'", ns("sc_data_source")),
                fileInput(ns("sc_fileUpload"),
                  label = div("Upload CSV/XLS/XLSX", help_icon("<strong>Upload stratigraphic data from CSV or Excel.</strong><details><summary>Learn more</summary>Upload your stratigraphic data as CSV or Excel.<br><br>Required columns: <b>CoreID</b> (identifies each core), <b>Sample</b> (sample name), <b>Depth or Age</b> (vertical position), and the <b>correlation variable</b> (e.g., UMAP_1D, a geochemical ratio).<br><br>Each row represents one measurement at a specific depth.</details>")),
                  accept = c(".csv", ".xlsx", ".xls")
                ),
                uiOutput(ns("sc_sheet_select"))
              ),
              uiOutput(ns("sc_processed_hint")),

              # Column Mapping
              conditionalPanel(
                condition = sprintf("output['%s']", ns("sc_data_loaded")),
                tags$hr(class = "my-2"),
                p(class = "text-muted small mb-2", "Map your data columns:"),
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(ns("sc_col_site"), "CoreID:", choices = NULL),
                  selectInput(ns("sc_col_sample"), "Sample:", choices = NULL)
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(ns("sc_col_z"), "Depth/Age:", choices = NULL),
                  selectInput(ns("sc_col_x"), "Corr Var:", choices = NULL)
                ),
                uiOutput(ns("sc_umap_hint")),
                actionButton(ns("sc_validate_map"), "Validate Columns",
                  icon = icon("check-circle"),
                  class = "btn-primary w-100 mt-2"
                )
              ),

              # Core Selection
              conditionalPanel(
                condition = sprintf("output['%s']", ns("sc_mapping_validated")),
                tags$hr(class = "my-2"),
                p(
                  class = "text-muted small mb-2",
                  tags$strong("Reference:"), " align to. ",
                  tags$strong("Target:"), " warp/align."
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(ns("sc_reference_core"), "Reference:", choices = NULL),
                  selectInput(ns("sc_target_core"), "Target:", choices = NULL)
                ),
                p(class = "text-muted small mb-1", tags$strong("Vertical direction:")),
                layout_columns(
                  col_widths = c(6, 6),
                  selectInput(
                    ns("sc_reference_direction"), "Reference values:",
                    choices = c(
                      "Increase downward (depth/BP)" = "down",
                      "Increase upward (CE/elevation)" = "up"
                    ),
                    selected = "down"
                  ),
                  selectInput(
                    ns("sc_target_direction"), "Target values:",
                    choices = c(
                      "Increase downward (depth/BP)" = "down",
                      "Increase upward (CE/elevation)" = "up"
                    ),
                    selected = "down"
                  )
                ),
                # Display options
                tags$hr(class = "my-2"),
                p(class = "text-muted small mb-1", tags$strong("Display Options:")),
                checkboxInput(ns("sc_show_points"), label = tags$span("Show point cloud", help_icon("<strong>Show raw points over interpolated core lines.</strong><details><summary>Learn more</summary>Overlay raw data points on top of the interpolated line plot.<br><br><b>Default: OFF.</b><br><br><b>When to turn ON:</b> To inspect where actual measurements were taken vs. interpolated values, and to identify outliers or gaps in the data.</details>")), value = FALSE),
                numericInput(ns("sc_target_offset"), label = tags$span("Target X Offset:", help_icon("<strong>Shift target-core X values by a constant offset.</strong><details><summary>Learn more</summary>Shift the Target core's X-values horizontally by a constant offset.<br><br><b>Default: 0.</b><br><br><b>When to use:</b> If two cores have different background levels for the correlation variable, applying an offset can help visually align their signals before using the warping tool.</details>")), value = 0, step = 0.1),
                layout_columns(
                  col_widths = c(6, 6),
                  sliderInput(ns("sc_line_width"), "Core line width:", min = 0.25, max = 4, value = 2, step = 0.25),
                  sliderInput(ns("sc_point_size"), "Core point size:", min = 2, max = 12, value = 8, step = 1)
                ),
                shinyWidgets::radioGroupButtons(
                  ns("sc_plot_height"),
                  label = "Plot height:",
                  choices = c("1×" = "standard", "1.5×" = "tall", "2×" = "double"),
                  selected = "standard",
                  justified = TRUE,
                  size = "sm"
                ),
                # Jitter options
                layout_columns(
                  col_widths = c(6, 6),
                  numericInput(ns("sc_jitter_x"), label = tags$span("X Jitter:", help_icon("<strong>Add horizontal jitter to overlapping points.</strong><details><summary>Learn more</summary>Add random horizontal jitter to data points.<br><br><b>Default: 0.</b><br><br><b>When to use:</b> When multiple samples plot at the exact same X value, making them stack on top of each other.<br><br>A small jitter (0.01–0.05) spreads them apart for visibility.</details>")), value = 0, min = 0, step = 0.01),
                  numericInput(ns("sc_jitter_z"), label = tags$span("Z Jitter:", help_icon("<strong>Add vertical jitter to overlapping points.</strong><details><summary>Learn more</summary>Add random vertical jitter to data points.<br><br><b>Default: 0.</b><br><br><b>When to use:</b> When multiple samples share the exact same depth/age value, making them overlap vertically.<br><br>A small jitter (0.1–0.5 depth units) spreads them for visibility.</details>")), value = 0, min = 0, step = 0.1)
                )
              )
            ),

            # 2. Tie Point Manager Card
            conditionalPanel(
              condition = sprintf("output['%s']", ns("sc_cores_selected")),
              varg_card(
                title = div(
                  class = "d-flex justify-content-between align-items-center",
                  span("2. Tie Points"),
                  div(
                    actionButton(ns("sc_clickToAdd"), "", icon = icon("crosshairs"), class = "btn-sm btn-outline-success", title = "Click to Add"),
                    actionButton(ns("sc_addTiePoint"), "", icon = icon("plus"), class = "btn-sm btn-outline-primary", title = "Add"),
                    actionButton(ns("sc_removeTiePoint"), "", icon = icon("trash"), class = "btn-sm btn-outline-danger", title = "Remove"),
                    actionButton(ns("sc_selectAll"), "", icon = icon("check-double"), class = "btn-sm btn-outline-secondary", title = "Select All"),
                    actionButton(ns("sc_clearSelection"), "", icon = icon("times"), class = "btn-sm btn-outline-secondary", title = "Clear")
                  )
                ),

                # Interaction Instructions (Task 5)
                help_box(
                    title = "Linking Samples",
                    content = "Create correlation links in two ways:
                    <ol>
                        <li><b>Manual Entry:</b> Type sample names directly into the input fields below.</li>
                        <li><b>Interactive Selection:</b> Click the target reticle button above, then click a sample on the Reference plot and a sample on the Target plot to link them.</li>
                    </ol>"
                ),

                # Click mode status indicator
                uiOutput(ns("sc_click_mode_status")),

                # Header Row
                div(
                  class = "row small fw-bold border-bottom pb-1 mb-2",
                  div(class = "col-1", "Sel"),
                  div(class = "col-3", "Ref Sample"),
                  div(class = "col-3", "Target Sample"),
                  div(class = "col-3", "Custom Name"),
                  div(class = "col-2 text-center", "Use")
                ),

                # Tie Point List
                uiOutput(ns("sc_tiepoint_list_ui"))
              )
            ),

            # 3. Actions Card
            conditionalPanel(
              condition = sprintf("output['%s']", ns("sc_cores_selected")),
              varg_card(
                title = "3. Apply & Export",
                actionButton(ns("sc_applyWarp"), "Apply Warping Model",
                  icon = icon("wand-magic-sparkles"),
                  class = "btn-success w-100"
                ),
                conditionalPanel(
                  condition = sprintf("output['%s']", ns("sc_warp_applied")),
                  tags$hr(class = "my-2"),
                  layout_columns(
                    col_widths = c(6, 6),
                    downloadButton(ns("sc_downloadWarped"), "Warped Data", class = "btn-sm btn-outline-primary w-100"),
                    downloadButton(ns("sc_downloadTiePoints"), "Tie Points", class = "btn-sm btn-outline-primary w-100")
                  ),
                  actionButton(ns("sc_copyTiePoints"), "Copy JSON for Chronology",
                    icon = icon("copy"),
                    class = "btn-sm btn-outline-info w-100 mt-2"
                  )
                )
              )
            )
          ),

          # --- RIGHT COLUMN: PLOTS ---
          varg_card(
            title = NULL,
            full_screen = TRUE,
            
            # Quick-access tie point controls above plots
            conditionalPanel(
              condition = sprintf("output['%s']", ns("sc_cores_selected")),
              div(
                class = "d-flex align-items-center gap-2 px-3 py-2 bg-light border-bottom",
                span(class = "small fw-bold text-muted", "Quick Tie Points:"),
                actionButton(ns("sc_clickToAdd_top"), "", icon = icon("crosshairs"), class = "btn-sm btn-outline-success", title = "Click to Add Tie Point"),
                span(class = "small text-muted ms-2", icon("info-circle"), " Click crosshairs, then click Reference plot + Target plot to link.")
              ),
              uiOutput(ns("sc_click_mode_status_top"))
            ),
            
            bslib::navset_underline(
              id = ns("sc_plot_tabs"),
              bslib::nav_panel(
                "Initial Alignment",
                uiOutput(ns("sc_initial_plot_ui"))
              ),
              bslib::nav_panel(
                "Affine Preview",
                div(
                  class = "px-3 pt-2 border-bottom bg-light",
                  layout_columns(
                    col_widths = c(5, 3, 3, 1),
                    shinyWidgets::radioGroupButtons(
                      ns("sc_affine_mode"),
                      label = "Target coordinate preview:",
                      choices = c("Match ranges" = "auto", "Manual" = "manual"),
                      selected = "auto",
                      justified = TRUE,
                      size = "sm"
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'manual'", ns("sc_affine_mode")),
                      numericInput(ns("sc_affine_scale"), "Target Z scale magnitude:", value = 1, min = 0.001, step = 0.01)
                    ),
                    conditionalPanel(
                      condition = sprintf("input['%s'] == 'manual'", ns("sc_affine_mode")),
                      numericInput(ns("sc_affine_shift"), "Target Z shift:", value = 0, step = 1)
                    ),
                    actionButton(
                      ns("sc_affine_reset"), "", icon = icon("rotate-left"),
                      class = "btn-sm btn-outline-secondary mt-4",
                      title = "Reset affine preview"
                    )
                  ),
                  uiOutput(ns("sc_affine_status"))
                ),
                uiOutput(ns("sc_affine_plot_ui"))
              ),
              bslib::nav_panel(
                "Warped Alignment",
                uiOutput(ns("sc_warped_plot_ui")),
                verbatimTextOutput(ns("sc_rmseOutput"))
              ),
              bslib::nav_panel(
                "Warp Fit",
                plotlyOutput(ns("sc_warpFitPlot"), height = "75vh")
              ),
              bslib::nav_panel(
                "Anchor Check",
                div(
                  class = "p-3",
                  h5(icon("check-double"), " Exact-anchor Verification"),
                  p(class = "text-muted small", "The mapping passes exactly through every active tie, so these residuals should be zero apart from numerical precision. This verifies the fitted anchors; it does not test whether a tie is geologically correct."),
                  plotlyOutput(ns("sc_residualPlot"), height = "40vh"),
                  tags$hr(),
                  h6("Anchor Details:"),
                  tableOutput(ns("sc_residualTable"))
                )
              )
            )
          )
        )
      )
    )
  )
}

mod_visualization_server <- function(id, processed_data = NULL, global_rv = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    viz_pkgs <- c("ggrepel", "ggdensity", "ggthemes")
    missing_viz_pkgs <- viz_pkgs[!vapply(viz_pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_viz_pkgs) > 0) {
      stop(paste0("Missing visualization package(s): ", paste(missing_viz_pkgs, collapse = ", ")))
    }

    hdr_method_spec <- function(method = "kde", adjust = 1, bins = NULL, auto_bins = TRUE) {
      method <- if (is.null(method) || !nzchar(method)) "kde" else method
      adjust <- suppressWarnings(as.numeric(adjust))
      if (length(adjust) == 0L || is.na(adjust) || adjust <= 0) adjust <- 1
      auto_bins <- if (is.null(auto_bins)) TRUE else isTRUE(auto_bins)
      bins <- suppressWarnings(as.integer(bins))
      if (auto_bins || length(bins) == 0L || is.na(bins)) {
        bins <- NULL
      } else {
        bins <- max(4L, min(bins, 60L))
      }
      switch(
        method,
        kde = ggdensity::method_kde(adjust = adjust),
        freqpoly = ggdensity::method_freqpoly(bins = bins),
        histogram = ggdensity::method_histogram(bins = bins),
        stop("Unsupported HDR density estimator: ", method)
      )
    }

    hdr_grid_resolution <- function(value = 100) {
      value <- suppressWarnings(as.integer(value))
      if (length(value) == 0L || is.na(value)) value <- 100L
      max(50L, min(value, 300L))
    }

    prepare_hdr_groups <- function(data, x_col, y_col, group_col, source_label) {
      usable <- stats::complete.cases(data[, c(x_col, y_col, group_col), drop = FALSE]) &
        is.finite(data[[x_col]]) & is.finite(data[[y_col]])
      excluded_rows <- sum(!usable)
      data <- data[usable, , drop = FALSE]

      if (excluded_rows > 0L) {
        showNotification(
          paste0(source_label, ": excluded ", excluded_rows,
                 " row(s) with missing or non-finite axis/group values."),
          type = "warning", duration = 8
        )
      }

      if (nrow(data) == 0L) {
        showNotification(paste(source_label, "has no complete reference rows for these axes."), type = "warning")
        return(NULL)
      }

      group_checks <- data %>%
        group_by(.data[[group_col]]) %>%
        summarise(
          .hdr_n = n(),
          .hdr_x_span = diff(range(.data[[x_col]], na.rm = TRUE)),
          .hdr_y_span = diff(range(.data[[y_col]], na.rm = TRUE)),
          .groups = "drop"
        )
      invalid_groups <- group_checks[[group_col]][
        group_checks$.hdr_n < 3L |
          !is.finite(group_checks$.hdr_x_span) | group_checks$.hdr_x_span <= 0 |
          !is.finite(group_checks$.hdr_y_span) | group_checks$.hdr_y_span <= 0
      ]
      if (length(invalid_groups) > 0L) {
        preview <- paste(utils::head(invalid_groups, 4L), collapse = ", ")
        suffix <- if (length(invalid_groups) > 4L) ", ..." else ""
        showNotification(
          paste0(source_label, ": skipped ", length(invalid_groups),
                 " group(s) with fewer than 3 usable observations or no spread on an axis (",
                 preview, suffix, ")."),
          type = "warning", duration = 8
        )
        data <- data[!(data[[group_col]] %in% invalid_groups), , drop = FALSE]
      }

      if (nrow(data) == 0L) {
        showNotification(paste(source_label, "has no groups large enough for an HDR overlay."), type = "warning")
        return(NULL)
      }

      data
    }

    # Load Helper Scripts
    # IMPORTANT: sc_helpers.R must be present in the app directory for Tab 3
    source("R/functions/sc_helpers.R", local = TRUE)
    source("R/functions/plot_palettes.R", local = TRUE)

    safe_vp_cores <- function(requested = NULL, cap = NULL) {
      env_cap <- suppressWarnings(as.integer(Sys.getenv("VARG_MAX_CORES", unset = "")))
      if (is.null(cap)) {
        cap <- if (!is.na(env_cap) && env_cap > 0L) env_cap else 12L
      }

      available <- parallel::detectCores(logical = TRUE)
      if (is.na(available) || available < 2L) {
        return(1L)
      }

      max_usable <- max(1L, min(available - 1L, cap))
      requested <- suppressWarnings(as.integer(requested))
      if (length(requested) == 0L || is.na(requested)) {
        physical <- parallel::detectCores(logical = FALSE)
        requested <- if (is.na(physical)) max_usable else physical
      }

      max(1L, min(requested, max_usable))
    }

    # -------------------------------------------------------------------------
    # INTERNAL HELPERS (For Tabs 1 & 2)
    # -------------------------------------------------------------------------
    transmute2 <- function(df, .x, .y, fun, ...) {
      transmute(df, "{{ .x }}_{{ .y }}" := fun({{ .x }}, {{ .y }}, ...))
    }

    transmute_pairwise <- function(df, fun, ..., associative = TRUE) {
      var_pairs <- t(combn(names(df), 2)) %>%
        as_tibble(.name_repair = "minimal") %>%
        setNames(c(".x", ".y"))
      if (!associative) {
        var_pairs <- bind_rows(var_pairs, rename(var_pairs, .y = .x, .x = .y))
      }
      dplyr::mutate(var_pairs, across(everything(), rlang::syms)) %>%
        purrr::pmap_dfc(transmute2, df = df, fun = fun, ...)
    }

    # KDE peak finding for HDR label placement
    find_kde_peak <- function(x, y) {
      if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
        return(tibble(peak_x = mean(x, na.rm = TRUE), peak_y = mean(y, na.rm = TRUE)))
      }
      tryCatch(
        {
          kde_result <- MASS::kde2d(x, y, n = 50)
          max_idx <- which(kde_result$z == max(kde_result$z), arr.ind = TRUE)
          if (nrow(max_idx) > 1) max_idx <- max_idx[1, , drop = FALSE]
          tibble(peak_x = kde_result$x[max_idx[1, 1]], peak_y = kde_result$y[max_idx[1, 2]])
        },
        error = function(e) {
          return(tibble(peak_x = mean(x, na.rm = TRUE), peak_y = mean(y, na.rm = TRUE)))
        }
      )
    }

    # -------------------------------------------------------------------------
    # REACTIVE VALUES
    # -------------------------------------------------------------------------
    rv <- reactiveValues(
      # Tab 1: Informative Projections
      vp_data = NULL,
      vp_ranked_pairs = NULL,
      vp_plot_data = NULL,

      # Tab 2: Custom Scatter Plots
      cp_data = NULL,
      cp_VARG26_data = NULL,
      cp_ref_data = NULL,
      cp_plot_obj = NULL,
      cp_imported_pairs = NULL,
      # Group customization overrides
      cp_group_color_overrides = NULL, # Named list: group -> hex color
      cp_group_shape_overrides = NULL, # Named list: group -> pch value
      cp_current_palette = "Set1", # Currently selected RColorBrewer palette

      # Tab 3: Stratigraphic Correlation (Optimized)
      sc_data = NULL,
      sc_column_mappings = NULL,
      sc_tiepoints = data.frame(
        id = integer(), ref_sample = character(), ref_z = numeric(),
        target_sample = character(), target_z = numeric(),
        custom_name = character(), use_in_warp = logical(), stringsAsFactors = FALSE
      ),
      sc_selected_tiepoints = integer(), # Vector of selected tie point IDs (supports multiple selection)
      sc_removed_ids = integer(), # Track IDs that have been removed (to avoid stale input caching)
      sc_warped_data = NULL,
      sc_warp_result = NULL,
      sc_output_object = NULL,
      # TRIGGER: Used to force re-render of list only when rows added/removed
      sc_structure_trigger = 0,
      # Click-to-add tie point mode: "none", "ref", "target"
      sc_click_mode = "none",
      sc_pending_ref_sample = NULL,
      sc_pending_ref_z = NULL
    )

    # Background process tracking for Variable Pair Finder
    vp_bg_process <- reactiveVal(NULL)
    vp_bg_context <- reactiveValues(notif_id = NULL, token = NULL)
    vp_generation <- reactiveVal(0L)

    heavy_job_lock <- NULL

    begin_heavy_job <- function(label, min_workers = 1L, max_workers = min_workers) {
      if (!is.null(global_rv)) {
        active <- isolate(global_rv$heavy_job)
        if (!is.null(active) && nzchar(active)) {
          showNotification(
            paste0("Another analysis is already running: ", active, ". Cancel it or wait for it to finish before starting a new heavy job."),
            type = "warning",
            duration = 8
          )
          return(FALSE)
        }
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
      if (!is.null(global_rv)) global_rv$heavy_job <- label
      TRUE
    }

    finish_heavy_job <- function() {
      if (!is.null(global_rv)) {
        global_rv$heavy_job <- NULL
      }
      if (!is.null(heavy_job_lock)) {
        heavy_job_limiter_release(heavy_job_lock)
        heavy_job_lock <<- NULL
      }
    }

    reset_vp_results <- function(cancel_running = TRUE) {
      proc <- tryCatch(isolate(vp_bg_process()), error = function(e) NULL)
      had_vp_job <- !is.null(proc) || !is.null(heavy_job_lock)
      if (isTRUE(cancel_running) && !is.null(proc)) {
        try({
          if (proc$is_alive()) {
            proc$kill()
            proc$wait(timeout = 2000)
          }
        }, silent = TRUE)
        vp_bg_process(NULL)
      }
      if (had_vp_job) finish_heavy_job()
      shinyjs::hide("vp_cancel")
      shinyjs::enable("vp_run")
      if (!is.null(isolate(vp_bg_context$notif_id))) {
        try(removeNotification(isolate(vp_bg_context$notif_id)), silent = TRUE)
      }
      rv$vp_ranked_pairs <- NULL
      rv$vp_plot_data <- NULL
      vp_generation(isolate(vp_generation()) + 1L)
      invisible(NULL)
    }

    replace_vp_data <- function(data) {
      reset_vp_results(cancel_running = TRUE)
      rv$vp_data <- data
      invisible(NULL)
    }

    current_vp_context <- function() {
      list(
        label = input$vp_label,
        predictors = input$vp_predictors,
        seed = input$vp_seed,
        k = input$vp_k,
        use_derived = isTRUE(input$vp_use_derived)
      )
    }

    capture_vp_token <- function() {
      varg_make_data_token(
        isolate(rv$vp_data),
        isolate(vp_generation()),
        current_vp_context()
      )
    }

    vp_token_is_current <- function(token) {
      varg_data_token_matches(
        token,
        isolate(rv$vp_data),
        isolate(vp_generation()),
        current_vp_context()
      )
    }

    reset_cp_derived_state <- function(clear_styles = FALSE) {
      rv$cp_plot_obj <- NULL
      rv$cp_imported_pairs <- NULL
      if (isTRUE(clear_styles)) {
        rv$cp_group_color_overrides <- NULL
        rv$cp_group_shape_overrides <- NULL
      }
      invisible(NULL)
    }

    replace_cp_data <- function(data) {
      reset_cp_derived_state(clear_styles = TRUE)
      rv$cp_data <- data
      invisible(NULL)
    }

    replace_cp_reference_data <- function(data) {
      reset_cp_derived_state(clear_styles = FALSE)
      rv$cp_ref_data <- data
      invisible(NULL)
    }

    empty_sc_tiepoints <- function() {
      data.frame(
        id = integer(), ref_sample = character(), ref_z = numeric(),
        target_sample = character(), target_z = numeric(),
        custom_name = character(), use_in_warp = logical(), stringsAsFactors = FALSE
      )
    }

    reset_sc_derived_state <- function() {
      rv$sc_column_mappings <- NULL
      rv$sc_tiepoints <- empty_sc_tiepoints()
      rv$sc_selected_tiepoints <- integer()
      rv$sc_removed_ids <- integer()
      rv$sc_warped_data <- NULL
      rv$sc_warp_result <- NULL
      rv$sc_output_object <- NULL
      rv$sc_click_mode <- "none"
      rv$sc_pending_ref_sample <- NULL
      rv$sc_pending_ref_z <- NULL
      rv$sc_structure_trigger <- isolate(rv$sc_structure_trigger) + 1L
      try(sc_prev_x_var(NULL), silent = TRUE)
      try(sc_prev_z_var(NULL), silent = TRUE)
      try(sc_prev_ref_core(NULL), silent = TRUE)
      try(sc_prev_target_core(NULL), silent = TRUE)
      invisible(NULL)
    }

    replace_sc_data <- function(data) {
      reset_sc_derived_state()
      rv$sc_data <- data
      invisible(NULL)
    }

    reset_visualization_state <- function() {
      reset_vp_results(cancel_running = TRUE)
      rv$vp_data <- NULL
      reset_cp_derived_state(clear_styles = TRUE)
      rv$cp_data <- NULL
      rv$cp_ref_data <- NULL
      reset_sc_derived_state()
      rv$sc_data <- NULL
      invisible(NULL)
    }
    if (!is.null(global_rv)) global_rv$reset_visualization_state <- reset_visualization_state

    session$onSessionEnded(function() {
      proc <- tryCatch(isolate(vp_bg_process()), error = function(e) NULL)
      if (!is.null(proc)) {
        try({
          if (proc$is_alive()) {
            proc$kill()
            proc$wait(timeout = 2000)
          }
        }, silent = TRUE)
      }
      try(finish_heavy_job(), silent = TRUE)
    })

    # =========================================================================
    # TAB 1: INFORMATIVE PROJECTIONS
    # =========================================================================

    # --- UI Conditions ---
    output$vp_data_loaded <- reactive({
      !is.null(rv$vp_data)
    })
    output$vp_vars_selected <- reactive({
      !is.null(input$vp_label) && !is.null(input$vp_predictors) && length(input$vp_predictors) >= 2
    })
    output$vp_results_available <- reactive({
      !is.null(rv$vp_ranked_pairs)
    })
    output$vp_runtime_warning <- renderUI({
      predictor_count <- length(input$vp_predictors)
      if (predictor_count < 2L) return(NULL)

      derived <- isTRUE(input$vp_use_derived)
      feature_count <- if (derived) {
        predictor_count + predictor_count * (predictor_count - 1L) +
          choose(predictor_count, 2L)
      } else {
        predictor_count
      }
      inverse_ratio_pairs <- if (derived) choose(predictor_count, 2L) else 0
      pair_count <- choose(feature_count, 2L) - inverse_ratio_pairs
      resample_count <- pair_count * 10
      base_pair_count <- choose(predictor_count, 2L)
      large_run <- pair_count >= 1000

      div(
        class = if (large_run) "alert alert-warning" else "alert alert-info",
        style = "margin-top: 8px; padding: 10px; font-size: 13px;",
        icon(if (large_run) "exclamation-triangle" else "calculator"),
        span(
          style = "margin-left: 8px;",
          "About ", tags$b(format(pair_count, big.mark = ",")),
          " candidate variable-pair models; approximately ",
          tags$b(format(resample_count, big.mark = ",")),
          " cross-validation fits.",
          if (large_run) {
            paste0(
              " This is a large exhaustive run. Turning off derived features would reduce it to ",
              format(base_pair_count, big.mark = ","), " models."
            )
          }
        )
      )
    })
    outputOptions(output, "vp_data_loaded", suspendWhenHidden = FALSE)
    outputOptions(output, "vp_vars_selected", suspendWhenHidden = FALSE)
    outputOptions(output, "vp_results_available", suspendWhenHidden = FALSE)

    # --- Sheet Selection UI for Tab 1 (VP) ---
    output$vp_sheet_select <- renderUI({
      req(input$vp_file)
      file_ext <- tolower(tools::file_ext(input$vp_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        sheets <- readxl::excel_sheets(input$vp_file$datapath)
        selectInput(ns("vp_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
      }
    })

    # --- Data Loading ---
    observeEvent(input$vp_file, {
      req(input$vp_file)
      file_ext <- tolower(tools::file_ext(input$vp_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        # Wait for sheet selection
        return()
      }
      tryCatch(
        {
          replace_vp_data(read.csv(input$vp_file$datapath, stringsAsFactors = FALSE))
          showNotification("Data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading file:", e$message), type = "error")
      )
    })

    # Allow using processed data from Processing module for VP tab
    output$vp_processed_hint <- renderUI({
      if (!is.null(processed_data) && is.function(processed_data)) {
        pd <- tryCatch(processed_data(), error = function(e) NULL)
        if (!is.null(pd) && NROW(pd) > 0) {
          tagList(span(style = "color:#2e7d32; font-size:12px;", "Processed data available. Select 'Use Processed Data' to import."))
        } else {
          tagList(span(style = "color:#888; font-size:12px;", "Processed data not available in this session."))
        }
      } else {
        tagList(span(style = "color:#888; font-size:12px;", "Processed data not available in this session."))
      }
    })

    observeEvent(input$vp_data_source, {
      if (input$vp_data_source == "processed") {
        if (!is.null(processed_data) && is.function(processed_data)) {
          pd <- tryCatch(processed_data(), error = function(e) NULL)
          if (!is.null(pd) && NROW(pd) > 0) {
            replace_vp_data(pd)
            showNotification("VP: using processed data from Processing module", type = "message")
          } else {
            showNotification("Processed data is empty or unavailable", type = "warning")
          }
        } else {
          showNotification("No processed data available in this session", type = "warning")
        }
      }
    })

    # Update VP data if processed_data reactive changes and VP tab has chosen processed
    if (!is.null(processed_data) && is.function(processed_data)) {
      observeEvent(processed_data(),
        {
          if (!is.null(input$vp_data_source) && input$vp_data_source == "processed") {
            replace_vp_data(processed_data())
          }
        },
        ignoreNULL = TRUE
      )
    }

    observeEvent(input$vp_sheet, {
      req(input$vp_file, input$vp_sheet)
      tryCatch(
        {
          replace_vp_data(readxl::read_excel(input$vp_file$datapath, sheet = input$vp_sheet))
          showNotification("Data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading file:", e$message), type = "error")
      )
    })

    output$vp_data_summary <- renderUI({
      req(rv$vp_data)
      div(class = "info-box-success", h5(icon("check-circle"), "Loaded"), p(paste(nrow(rv$vp_data), "rows ×", ncol(rv$vp_data), "columns")))
    })

    # --- Variable Selection ---
    output$vp_label_select <- renderUI({
      req(rv$vp_data)
      char_cols <- names(rv$vp_data)[sapply(rv$vp_data, function(x) is.character(x) || is.factor(x))]
      selectInput(ns("vp_label"), NULL, choices = char_cols)
    })

    output$vp_predictors_select <- renderUI({
      req(rv$vp_data)
      numeric_cols <- names(rv$vp_data)[sapply(rv$vp_data, is.numeric)]
      selectizeInput(
        ns("vp_predictors"),
        NULL,
        choices = numeric_cols,
        selected = NULL,
        multiple = TRUE,
        options = list(
          plugins = list("remove_button"),
          placeholder = "Type to search; click × to remove"
        )
      )
    })

    output$vp_cores_select <- renderUI({
      if (!USE_BG_PROCESSES) {
        # On server: always sequential, hide cores input
        return(NULL)
      }
      physical_cores <- detectCores(logical = FALSE)
      default_cores <- safe_vp_cores(physical_cores)
      max_cores <- safe_vp_cores(detectCores(logical = TRUE))
      if (identical(tolower(Sys.getenv("VARG_DEPLOYMENT_MODE", unset = "")), "hosted")) {
        min_cores <- min(4L, max_cores)
        return(div(
          class = "text-muted small",
          icon("microchip"),
          paste0(
            "Hosted VP analysis receives ", min_cores, "-", max_cores,
            " CPU workers from the shared burst pool when capacity is available."
          )
        ))
      }
      numericInput(ns("vp_n_cores"), label = div(paste("CPU Cores (", physical_cores, "detected; ", max_cores, "available)"), help_icon("<strong>Set CPU cores for parallel pair evaluation.</strong><details><summary>Learn more</summary>Number of parallel CPU cores to use for pair evaluation.<br><br>More cores = faster analysis, up to the point where startup overhead dominates.<br><br><b>Default: auto-detected safe maximum.</b><br><br><b>Tip:</b> The app leaves one logical core free and caps this background analysis for local and hosted stability.<br><br><b>Running on a shared server?</b> Be conservative; using too many cores may slow down other users.</details>")), value = default_cores, min = 1, max = max_cores)
    })

    # Select All Predictors Shortcut
    observeEvent(input$select_all_preds, {
      req(rv$vp_data)
      numeric_cols <- names(rv$vp_data)[sapply(rv$vp_data, is.numeric)]
      updateSelectizeInput(session, "vp_predictors", selected = numeric_cols)
    })

    # --- Analysis Logic ---
    observeEvent(input$vp_run, {
      req(rv$vp_data, input$vp_label, input$vp_predictors, input$vp_seed)
      if (length(input$vp_predictors) < 2) {
        showNotification("Please select at least 2 predictor variables", type = "error")
        return()
      }

      vp_requested_cores <- safe_vp_cores(input$vp_n_cores)
      hosted_mode <- identical(tolower(Sys.getenv("VARG_DEPLOYMENT_MODE", unset = "")), "hosted")
      vp_min_cores <- if (hosted_mode) min(4L, vp_requested_cores) else 1L
      if (!begin_heavy_job(
        "Variable Pair Finder",
        min_workers = vp_min_cores,
        max_workers = vp_requested_cores
      )) return()

      rv$vp_ranked_pairs <- NULL
      rv$vp_plot_data <- NULL
      showNotification("Variable pair analysis started in background...", type = "message", duration = 5)

      tryCatch(
        {
          df_initial <- rv$vp_data %>% dplyr::select(dplyr::all_of(c(input$vp_label, input$vp_predictors)))
          
          # Filter out NA and "Unassigned" labels before analysis
          raw_label <- df_initial[[input$vp_label]]
          valid_label <- !is.na(raw_label) & trimws(as.character(raw_label)) != "" & 
                         tolower(trimws(as.character(raw_label))) != "unassigned"
          df_initial <- df_initial[valid_label, , drop = FALSE]
          
          if (nrow(df_initial) < 10) stop("Not enough labelled rows after excluding NA and Unassigned (need at least 10).")
          
          label <- droplevels(as.factor(make.names(df_initial[[input$vp_label]])))
          predictors_initial <- as.data.frame(df_initial %>% dplyr::select(-dplyr::all_of(input$vp_label)))
          set.seed(input$vp_seed)

          is_valid <- function(x) {
            !is.na(x) & is.finite(x) & x != 0
          }
          valid_rows <- apply(predictors_initial, 1, function(r) all(is_valid(r)))
          if (sum(valid_rows) < 10) stop("Not enough valid cases (need at least 10).")

          predictors_complete <- predictors_initial[valid_rows, , drop = FALSE]
          label_complete <- droplevels(label[valid_rows])
          
          # Ensure at least 2 classes remain
          if (nlevels(label_complete) < 2) stop("Need at least 2 distinct groups for classification. Check your label column.")

          if (input$vp_use_derived) {
            ratios <- transmute_pairwise(df = predictors_complete, fun = `/`, associative = FALSE) %>% rename_with(~ paste0(.x, "_ratio"))
            sums <- transmute_pairwise(df = predictors_complete, fun = `+`, associative = TRUE) %>% rename_with(~ paste0(.x, "_sum"))
            predictors_with_pairs <- cbind(predictors_complete, ratios, sums)
          } else {
            predictors_with_pairs <- predictors_complete
          }

          are_inverse_ratios <- function(p1, p2) {
            if (!grepl("_ratio$", p1) || !grepl("_ratio$", p2)) {
              return(FALSE)
            }
            p1_parts <- strsplit(sub("_ratio$", "", p1), "_")[[1]]
            if (length(p1_parts) == 2) {
              return(paste(p1_parts[2], p1_parts[1], sep = "_") == sub("_ratio$", "", p2))
            }
            return(FALSE)
          }

          all_pairs <- combn(names(predictors_with_pairs), 2, simplify = FALSE)
          variable_pairs <- Filter(function(p) !are_inverse_ratios(p[1], p[2]), all_pairs)

          id <- showNotification(
            if (USE_BG_PROCESSES) "Variable pair analysis running in background. You can cancel anytime."
            else "Variable pair analysis running. This may take a while...",
            type = "message", duration = NULL, closeButton = FALSE
          )
          vp_seed <- input$vp_seed
          vp_k <- as.integer(input$vp_k)
          vp_n_cores <- heavy_job_limiter_workers(heavy_job_lock, default = 1L)

          # Include identifier columns for hover labels
          # Carry all non-numeric columns from the original data for use in hover
          id_cols <- names(rv$vp_data)[sapply(rv$vp_data, function(x) is.character(x) || is.factor(x))]
          id_cols <- setdiff(id_cols, input$vp_label)  # Exclude the label col (already included)
          id_data <- rv$vp_data[valid_label, , drop = FALSE][valid_rows, id_cols, drop = FALSE]
          rv$vp_plot_data <- cbind(predictors_with_pairs, label = label_complete, id_data)

          if (USE_BG_PROCESSES) {
            # --- Background path: cancellable via callr ---
            vp_bg_context$notif_id <- id
            vp_bg_context$token <- capture_vp_token()
            shinyjs::show("vp_cancel")
            shinyjs::disable("vp_run")

            vp_bg_process(callr::r_bg(
              func = vp_score_pairs,
              args = list(
                variable_pairs = variable_pairs,
                predictors_with_pairs = predictors_with_pairs,
                label_complete = label_complete,
                vp_seed = vp_seed,
                vp_k = vp_k,
                vp_n_cores = vp_n_cores
              ),
              supervise = TRUE
            ))
          } else {
            # --- Synchronous path: no extra process, lower memory ---
            shinyjs::disable("vp_run")
            results_list <- withProgress(message = "Running variable pair analysis...", value = 0.5, {
              tryCatch(
                vp_score_pairs(variable_pairs, predictors_with_pairs, label_complete, vp_seed, vp_k, vp_n_cores),
                error = function(e) {
                  showNotification(paste("Analysis failed:", e$message), type = "error", duration = 10)
                  NULL
                }
              )
            })

            shinyjs::enable("vp_run")
            finish_heavy_job()
            removeNotification(id)

            if (!is.null(results_list)) {
              rv$vp_ranked_pairs <- results_list %>%
                dplyr::arrange(dplyr::desc(Kappa)) %>%
                dplyr::mutate(
                  Kappa_Quality = dplyr::case_when(
                    is.na(Kappa) ~ "Unknown",
                    Kappa > 0.8 ~ "Excellent",
                    Kappa >= 0.6 ~ "Good",
                    Kappa >= 0.4 ~ "Moderate",
                    TRUE ~ "Poor"
                  )
                ) %>%
                dplyr::mutate(dplyr::across(c(Accuracy, Kappa), ~ round(.x, 4)))
              showNotification("Analysis completed successfully!", type = "message")
            }
          }
        },
        error = function(e) {
          shinyjs::hide("vp_cancel")
          shinyjs::enable("vp_run")
          finish_heavy_job()
          try(removeNotification(vp_bg_context$notif_id), silent = TRUE)
          showNotification(paste("Preparation failed:", e$message), type = "error")
        }
      )
    })

    # Cancel Variable Pair Finder
    observeEvent(input$vp_cancel, {
      proc <- vp_bg_process()
      if (!is.null(proc) && proc$is_alive()) {
        proc$kill()
        try(proc$wait(timeout = 2000), silent = TRUE)
        vp_bg_process(NULL)
        shinyjs::hide("vp_cancel")
        shinyjs::enable("vp_run")
        finish_heavy_job()
        removeNotification(vp_bg_context$notif_id)
        showNotification("Variable pair analysis cancelled.", type = "warning", duration = 5)
      }
    })

    # Poll for VP Finder completion
    observe({
      proc <- vp_bg_process()
      req(proc)
      if (proc$is_alive()) {
        invalidateLater(1000)
        return()
      }

      # Process finished
      vp_bg_process(NULL)
      shinyjs::hide("vp_cancel")
      shinyjs::enable("vp_run")
      finish_heavy_job()
      removeNotification(vp_bg_context$notif_id)

      results_list <- tryCatch(proc$get_result(), error = function(e) {
        showNotification(paste("Analysis failed:", e$message), type = "error", duration = 10)
        return(NULL)
      })

      if (!is.null(results_list)) {
        if (!vp_token_is_current(vp_bg_context$token)) {
          showNotification(
            "Variable-pair analysis finished, but its input data or settings changed. The stale result was discarded.",
            type = "warning",
            duration = 8
          )
          return()
        }
        rv$vp_ranked_pairs <- results_list %>%
          arrange(desc(Kappa)) %>%
          mutate(
            Kappa_Quality = dplyr::case_when(
              is.na(Kappa) ~ "Unknown",
              Kappa > 0.8 ~ "Excellent",
              Kappa >= 0.6 ~ "Good",
              Kappa >= 0.4 ~ "Moderate",
              TRUE ~ "Poor"
            )
          ) %>%
          mutate(across(c(Accuracy, Kappa), ~ round(.x, 4)))
        showNotification("Analysis completed successfully!", type = "message")
      }
    })

    # Hover label selector for VP scatter plot
    output$vp_hover_label_ui <- renderUI({
      req(rv$vp_plot_data)
      # Get available identifier columns (non-numeric + non-label, non-predictor)
      all_cols <- names(rv$vp_plot_data)
      # Exclude label and numeric predictor columns - keep only string/identifier columns
      id_cols <- all_cols[sapply(rv$vp_plot_data, function(x) is.character(x) || is.factor(x))]
      id_cols <- setdiff(id_cols, "label")
      if (length(id_cols) == 0) return(NULL)
      # Smart default: prefer sample_point, sample_id, Samplepop, UID
      default <- intersect(c("sample_point", "sample_id", "Samplepop", "UID"), id_cols)
      selected <- if (length(default) > 0) default[1] else id_cols[1]
      selectInput(ns("vp_hover_label"),
        label = "Point identifier on hover:",
        choices = c("(none)" = "", id_cols),
        selected = selected
      )
    })

    output$vp_results_table <- renderDT({
      req(rv$vp_ranked_pairs)
      datatable(rv$vp_ranked_pairs, selection = "single", options = list(pageLength = 15, scrollX = TRUE, scrollY = "300px"), rownames = FALSE) %>%
        formatStyle("Kappa", background = styleColorBar(range(rv$vp_ranked_pairs$Kappa), "lightblue"))
    })

    output$vp_scatter <- renderPlotly({
      req(rv$vp_plot_data, rv$vp_ranked_pairs)
      selected_row <- if (length(input$vp_results_table_rows_selected)) input$vp_results_table_rows_selected else 1
      pair_info <- rv$vp_ranked_pairs[selected_row, ]
      var1 <- pair_info$Variable_1
      var2 <- pair_info$Variable_2

      # Build hover text with optional point identifier
      hover_col <- input$vp_hover_label
      pd <- rv$vp_plot_data
      if (!is.null(hover_col) && nzchar(hover_col) && hover_col %in% names(pd)) {
        hover_text <- paste0(hover_col, ": ", pd[[hover_col]],
                            "<br>Group: ", pd$label,
                            "<br>", var1, ": ", round(pd[[var1]], 3),
                            "<br>", var2, ": ", round(pd[[var2]], 3))
      } else {
        hover_text <- paste0("Group: ", pd$label,
                            "<br>", var1, ": ", round(pd[[var1]], 3),
                            "<br>", var2, ": ", round(pd[[var2]], 3))
      }

      plot_ly(
        data = pd, x = ~ get(var1), y = ~ get(var2), color = ~label, type = "scatter", mode = "markers",
        text = hover_text, hoverinfo = "text"
      ) %>%
        layout(title = paste0("<b>", var1, " vs ", var2, "</b>"), xaxis = list(title = var1), yaxis = list(title = var2))
    })

    output$vp_download <- downloadHandler(
      filename = function() paste0("knn_pairs_", format(Sys.time(), "%Y%m%d"), ".xlsx"),
      content = function(file) writexl::write_xlsx(rv$vp_ranked_pairs, file)
    )

    # =========================================================================
    # TAB 2: CUSTOM SCATTER PLOTS
    # =========================================================================

    output$cp_data_loaded <- reactive({
      !is.null(rv$cp_data)
    })
    outputOptions(output, "cp_data_loaded", suspendWhenHidden = FALSE)
    output$cp_VARG26_selected_any <- reactive({
      !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0
    })
    outputOptions(output, "cp_VARG26_selected_any", suspendWhenHidden = FALSE)
    output$cp_ref_loaded <- reactive({
      !is.null(rv$cp_ref_data)
    })
    outputOptions(output, "cp_ref_loaded", suspendWhenHidden = FALSE)
    output$cp_ref_selected_any <- reactive({
      !is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0
    })
    outputOptions(output, "cp_ref_selected_any", suspendWhenHidden = FALSE)

    normalize_varg26_viz <- function(VARG26) {
      if (is.null(VARG26) || !NROW(VARG26) || !"Data citation" %in% names(VARG26)) {
        return(VARG26)
      }

      VARG26$viz_visible <- TRUE
      sample_id <- if ("Sample_ID" %in% names(VARG26)) {
        trimws(as.character(VARG26[["Sample_ID"]]))
      } else {
        rep("", NROW(VARG26))
      }
      spm26_payne_qr <- VARG26[["Data citation"]] == "As yet undisclosed data from Alaska" &
        grepl("^SPM26$", sample_id, ignore.case = TRUE)

      # Keep the published family rows that currently sit behind informal file labels.
      # Anything not yet pinned to a published family stays hidden from the app view.
      if ("Tephra" %in% names(VARG26)) {
        tephra <- trimws(as.character(VARG26[["Tephra"]]))
        wr_rows <- VARG26[["Data citation"]] == "WR tephra master list_Jensen.xls"
        wr_unpublished <- grepl("^WRUN[5-6]$", tephra, ignore.case = FALSE)
        wr_unmapped <- is.na(VARG26[["Tephra"]]) | tephra == ""
        VARG26$viz_visible[wr_rows & (wr_unpublished | wr_unmapped)] <- FALSE
      }

      VARG26$viz_visible[VARG26[["Data citation"]] %in% c(
        "As yet undisclosed data from Alaska",
        "Unpublished internal data",
        "Excluded from counts: secondary standard",
        "Jensen et al. 2011?",
        "Reyce tephra data_dec20"
      )] <- FALSE
      VARG26$viz_visible[spm26_payne_qr] <- TRUE

      if ("Sample_ID" %in% names(VARG26)) {
        is_secondary_standard <- grepl("^(BHVO|ID3506)", sample_id, ignore.case = TRUE)
        VARG26$viz_visible[is_secondary_standard] <- FALSE
      }

      citation_short_map <- c(
        "WR tephra master list_Jensen.xls" = "Turner et al. 2013",
        "Westgate database" = "Kaufman et al. 2001",
        "Associated tephra geochem_VT_List.xls" = "Jensen 2013 PhD thesis",
        "Payne and Blackford (source papers in Payne data.xls)" = "Payne & Blackford 2008",
        "U1436 and U1435 data" = "Schindlbeck et al. 2018",
        "U1438 data" = "IODP Expedition 351 Site U1438 report",
        "J.B. Hunt and Y.M.R. Najman, Ms 186SR-107" = "Hunt & Najman (2003)",
        "Jensen et al. 2011?" = "Jensen et al. 2011",
        "VT majors_final_Jensen2011.xls" = "Jensen et al. 2011"
      )

      citation_full_map <- c(
        "WR tephra master list_Jensen.xls" = "Turner, D.G., Ward, B.C., Bond, J.D., Jensen, B.J.L., Froese, D.G., Telka, A.M., Zazula, G.D., Bigelow, N.H. (2013). Middle to Late Pleistocene ice extents, tephrochronology and paleoenvironments of the White River area, southwest Yukon. Quaternary Science Reviews 75, 59-77. DOI: 10.1016/j.quascirev.2013.05.011",
        "Westgate database" = "Kaufman, D.S., Manley, W.F., Wolfe, A.P., Hu, F.S., Preece, S.J., Westgate, J.A., Forman, S.L. (2001). The last interglacial to glacial transition, Togiak Bay, southwestern Alaska. Quaternary Research 55, 190-202. DOI: 10.1006/qres.2001.2214",
        "Associated tephra geochem_VT_List.xls" = "Jensen, B.J.L. (2013). Tephrostratigraphy and paleoenvironments of the late Quaternary in eastern Beringia. PhD thesis, University of Alberta. https://www.collectionscanada.gc.ca/obj/thesescanada/vol2/002/NR91384.pdf",
        "Payne and Blackford (source papers in Payne data.xls)" = "Payne, R.J. & Blackford, J.J. (2008). Distal volcanic impacts on peatlands: palaeoecological evidence from Alaska. Quaternary Science Reviews 27, 2012-2030. DOI: 10.1016/j.quascirev.2008.08.002",
        "U1436 and U1435 data" = "Schindlbeck, J.C., et al. (2018). One million years tephra record at IODP Sites U1436 and U1437: insights into explosive volcanism from the Japan and Izu arcs. DOI: 10.1111/iar.12244",
        "U1438 data" = "IODP Expedition 351. Site U1438 report. DOI: 10.14379/iodp.proc.351.103.2015",
        "J.B. Hunt and Y.M.R. Najman, Ms 186SR-107" = "Hunt, J.B. & Najman, Y.M.R. (2003). Ms 186SR-107. In: Proceedings of the Ocean Drilling Program, Scientific Results, Volume 186.",
        "Jensen et al. 2011?" = "Jensen, B.J.L., Preece, S.J., Lamothe, M., Pearce, N.J.G., Froese, D.G., Westgate, J.A., Schaefer, J., Begét, J. (2011). The Variegated (VT) tephra: a new regional marker for middle to late marine isotope stage 5 across Yukon and Alaska. Quaternary International 246, 312-323. DOI: 10.1016/j.quaint.2011.06.028",
        "VT majors_final_Jensen2011.xls" = "Jensen, B.J.L., Preece, S.J., Lamothe, M., Pearce, N.J.G., Froese, D.G., Westgate, J.A., Schaefer, J., Begét, J. (2011). The Variegated (VT) tephra: a new regional marker for middle to late marine isotope stage 5 across Yukon and Alaska. Quaternary International 246, 312-323. DOI: 10.1016/j.quaint.2011.06.028"
      )

      for (raw in names(citation_short_map)) {
        idx <- VARG26[["Data citation"]] == raw
        if (any(idx, na.rm = TRUE)) {
          VARG26[idx, "Data citation"] <- citation_short_map[[raw]]
          if ("Full citation" %in% names(VARG26) && !is.na(citation_full_map[[raw]]) && nchar(citation_full_map[[raw]]) > 0) {
            VARG26[idx, "Full citation"] <- citation_full_map[[raw]]
          }
        }
      }

      if (any(spm26_payne_qr, na.rm = TRUE)) {
        VARG26[spm26_payne_qr, "Data citation"] <- "Payne et al. 2008"
        if ("Full citation" %in% names(VARG26)) {
          VARG26[spm26_payne_qr, "Full citation"] <- "Payne, R.J., Blackford, J.J., van der Plicht, J. (2008). Using cryptotephras to extend regional tephrochronologies: an example from southeast Alaska and implications for hazard assessment. Quaternary Research 69(1), 42-55. DOI: 10.1016/j.yqres.2007.10.007"
        }
      }

      VARG26 <- VARG26[VARG26$viz_visible %in% TRUE, , drop = FALSE]
      VARG26$viz_visible <- NULL
      VARG26
    }

    load_VARG26 <- reactive({
      if (!is.null(rv$cp_VARG26_data)) {
        return(rv$cp_VARG26_data)
      }
      VARG26_path <- Sys.getenv("VARG26_REFERENCE_FILE", unset = "VARG26_Public_Reference_2026-07-13.csv")
      if (!file.exists(VARG26_path)) {
        VARG26_path <- "VARG-Tools_Processed_Data_2025-11-21 _VARG26 FINAL.csv"
      }
      if (!file.exists(VARG26_path)) {
        showNotification("VARG26 reference file not found", type = "warning")
        return(NULL)
      }
      tryCatch(
        {
          VARG26 <- read.csv(VARG26_path, stringsAsFactors = FALSE, check.names = FALSE)
          VARG26 <- normalize_varg26_viz(VARG26)
          rv$cp_VARG26_data <- VARG26
          VARG26
        },
        error = function(e) NULL
      )
    })

    # --- Sheet Selection UI for Tab 2 (CP) ---
    output$cp_sheet_select <- renderUI({
      req(input$cp_file)
      file_ext <- tolower(tools::file_ext(input$cp_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        sheets <- readxl::excel_sheets(input$cp_file$datapath)
        selectInput(ns("cp_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
      }
    })

    observeEvent(input$cp_file, {
      req(input$cp_file)
      file_ext <- tolower(tools::file_ext(input$cp_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        # Wait for sheet selection
        return()
      }
      tryCatch(
        {
          replace_cp_data(read.csv(input$cp_file$datapath, stringsAsFactors = FALSE, check.names = FALSE))
          showNotification("User data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading file:", e$message), type = "error")
      )
    })

    # CP processed data hint and wiring
    output$cp_processed_hint <- renderUI({
      if (!is.null(processed_data) && is.function(processed_data)) {
        pd <- tryCatch(processed_data(), error = function(e) NULL)
        if (!is.null(pd) && NROW(pd) > 0) {
          tagList(span(style = "color:#2e7d32; font-size:12px;", "Processed data available. Select 'Use Processed Data' to import."))
        } else {
          tags$p(class = "text-muted small fst-italic", 
            "No processed data available. Complete ", 
            tags$b("Step 1 (Processing)"), 
            " first, or upload a CSV/Excel file above.")
        }
      } else {
        tags$p(class = "text-muted small fst-italic", 
          "No processed data available. Complete ", 
          tags$b("Step 1 (Processing)"), 
          " first, or upload a CSV/Excel file above.")
      }
    })

    observeEvent(input$cp_data_source, {
      if (input$cp_data_source == "processed") {
        if (!is.null(processed_data) && is.function(processed_data)) {
          pd <- tryCatch(processed_data(), error = function(e) NULL)
          if (!is.null(pd) && NROW(pd) > 0) {
            replace_cp_data(pd)
            showNotification("CP: using processed data from Processing module", type = "message")
          } else {
            showNotification("Processed data is empty or unavailable", type = "warning")
          }
        } else {
          showNotification("No processed data available in this session", type = "warning")
        }
      }
    })

    if (!is.null(processed_data) && is.function(processed_data)) {
      observeEvent(processed_data(),
        {
          if (!is.null(input$cp_data_source) && input$cp_data_source == "processed") {
            replace_cp_data(processed_data())
          }
        },
        ignoreNULL = TRUE
      )
    }

    observeEvent(input$cp_sheet, {
      req(input$cp_file, input$cp_sheet)
      tryCatch(
        {
          replace_cp_data(readxl::read_excel(input$cp_file$datapath, sheet = input$cp_sheet))
          showNotification("User data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading file:", e$message), type = "error")
      )
    })

    output$cp_data_summary <- renderUI({
      req(rv$cp_data)
      div(style = "background:#d4edda;padding:8px;", h5("Loaded"), p(paste(nrow(rv$cp_data), "rows")), p(class = "small text-muted mb-0", "Select X and Y axes to start plotting."))
    })

    # Import Pairs Logic
    compute_variable <- function(var_name, data, suffix) {
      if (grepl("_ratio$", var_name)) {
        parts <- strsplit(sub("_ratio$", "", var_name), "_")[[1]]
        if (length(parts) >= 2) {
          v1 <- paste0(parts[1], suffix)
          v2 <- paste0(paste(parts[-1], collapse = "_"), suffix)
          if (v1 %in% names(data) && v2 %in% names(data)) {
            return(list(values = data[[v1]] / data[[v2]], col_name = var_name))
          }
        }
      }
      if (grepl("_sum$", var_name)) {
        parts <- strsplit(sub("_sum$", "", var_name), "_")[[1]]
        if (length(parts) >= 2) {
          v1 <- paste0(parts[1], suffix)
          v2 <- paste0(paste(parts[-1], collapse = "_"), suffix)
          if (v1 %in% names(data) && v2 %in% names(data)) {
            return(list(values = data[[v1]] + data[[v2]], col_name = var_name))
          }
        }
      }
      return(NULL) # Base var logic simplified for brevity
    }

    observeEvent(input$cp_import_pairs, {
      if (is.null(rv$vp_ranked_pairs) || is.null(rv$cp_data)) {
        showNotification("Run Tab 1 analysis and upload Tab 2 data first.", type = "warning")
        return()
      }
      tryCatch(
        {
          top_pairs <- rv$vp_ranked_pairs[1:min(input$cp_n_pairs, nrow(rv$vp_ranked_pairs)), ]
          suffix <- input$cp_pair_suffix
          # Convert "raw" to empty string for no suffix
          if (!is.null(suffix) && suffix == "raw") suffix <- ""
          computed_cols <- list()
          for (i in 1:nrow(top_pairs)) {
            v1 <- compute_variable(top_pairs$Variable_1[i], rv$cp_data, suffix)
            v2 <- compute_variable(top_pairs$Variable_2[i], rv$cp_data, suffix)
            if (!is.null(v1)) computed_cols[[v1$col_name]] <- v1$values
            if (!is.null(v2)) computed_cols[[v2$col_name]] <- v2$values
          }
          if (length(computed_cols) > 0) {
            rv$cp_data <- cbind(rv$cp_data, as.data.frame(computed_cols))

            # ALSO APPLY TO REFERENCE DATA IF LOADED
            if (!is.null(rv$cp_ref_data)) {
              ref_computed_cols <- list()
              for (i in 1:nrow(top_pairs)) {
                v1 <- compute_variable(top_pairs$Variable_1[i], rv$cp_ref_data, suffix)
                v2 <- compute_variable(top_pairs$Variable_2[i], rv$cp_ref_data, suffix)
                if (!is.null(v1)) ref_computed_cols[[v1$col_name]] <- v1$values
                if (!is.null(v2)) ref_computed_cols[[v2$col_name]] <- v2$values
              }
              if (length(ref_computed_cols) > 0) {
                rv$cp_ref_data <- cbind(rv$cp_ref_data, as.data.frame(ref_computed_cols))
              }
            }

            rv$cp_imported_pairs <- top_pairs
            showNotification(paste("Imported columns for", nrow(top_pairs), "pairs."), type = "message")
          } else {
            # No columns could be computed - warn user about missing suffix columns
            suffix_label <- if (suffix == "") "raw (no suffix)" else suffix
            showNotification(
              paste0("No columns found with '", suffix_label, "' suffix. ",
                     "Your data may not have this transformation applied. ",
                     "Try 'None (raw)' or a different suffix option."),
              type = "warning",
              duration = 8
            )
          }
        },
        error = function(e) showNotification(e$message, type = "error")
      )
    })

    output$cp_imported_pairs_summary <- renderUI({
      if (is.null(rv$cp_imported_pairs)) {
        return(p("No pairs imported yet", style = "color:#999"))
      }
      div(style = "background:#e8f5e9;padding:5px;", p("Pairs imported"))
    })

    # Axes & Aesthetics UI
    # Axes & Aesthetics Update Logic
    observeEvent(rv$cp_data, {
      req(rv$cp_data)
      numeric_cols <- names(rv$cp_data)[sapply(rv$cp_data, is.numeric)]
      all_cols <- names(rv$cp_data)

      # Auto-guess X and Y
      # Priority: UMAP > SiO2/K2O > First numeric

      # Guess X
      default_x <- numeric_cols[1]
      if (any(grepl("UMAP.*1$", numeric_cols))) {
        default_x <- grep("UMAP.*1$", numeric_cols, value = TRUE)[1]
      } else if (any(grepl("SiO2", numeric_cols, ignore.case = TRUE))) {
        default_x <- grep("SiO2", numeric_cols, ignore.case = TRUE, value = TRUE)[1]
      }

      # Guess Y
      default_y <- if (length(numeric_cols) > 1) numeric_cols[2] else numeric_cols[1]
      if (any(grepl("UMAP.*2$", numeric_cols))) {
        default_y <- grep("UMAP.*2$", numeric_cols, value = TRUE)[1]
      } else if (any(grepl("K2O", numeric_cols, ignore.case = TRUE)) && !grepl("K2O", default_x, ignore.case = TRUE)) {
        default_y <- grep("K2O", numeric_cols, ignore.case = TRUE, value = TRUE)[1]
      } else if (any(grepl("Na2O", numeric_cols, ignore.case = TRUE)) && !grepl("Na2O", default_x, ignore.case = TRUE)) {
        # TAS diagram fallback (SiO2 vs Na2O+K2O usually, but here just Na2O as proxy if K2O missing)
        default_y <- grep("Na2O", numeric_cols, ignore.case = TRUE, value = TRUE)[1]
      }

      updateSelectInput(session, "cp_x_var", choices = numeric_cols, selected = default_x)
      updateSelectInput(session, "cp_y_var", choices = numeric_cols, selected = default_y)

      # Reset color/shape to safe defaults when new data loaded
      # Check if current selections exist in new data
      current_color <- input$cp_color_var
      current_shape <- input$cp_shape_var
      
      color_valid <- is.null(current_color) || current_color == "Fixed Color" || current_color %in% all_cols
      shape_valid <- is.null(current_shape) || current_shape == "Fixed Shape" || current_shape %in% all_cols
      
      # For color: allow all columns
      if (!color_valid) {
        updatePickerInput(session, "cp_color_var", choices = c("Fixed Color", all_cols), selected = "Fixed Color")
      } else {
        updatePickerInput(session, "cp_color_var", choices = c("Fixed Color", all_cols), selected = current_color)
      }

      # For shape: exclude numeric columns (shapes don't work with continuous data)
      non_numeric_cols <- names(rv$cp_data)[!sapply(rv$cp_data, is.numeric)]
      if (!shape_valid || !(current_shape %in% c("Fixed Shape", non_numeric_cols))) {
        updatePickerInput(session, "cp_shape_var", choices = c("Fixed Shape", non_numeric_cols), selected = "Fixed Shape")
      } else {
        updatePickerInput(session, "cp_shape_var", choices = c("Fixed Shape", non_numeric_cols), selected = current_shape)
      }

      # Update Filter Column Choices
      # Prefer character/factor columns; prepend "(None)" so filter defaults to off
      char_cols <- names(rv$cp_data)[sapply(rv$cp_data, function(x) is.character(x) || is.factor(x))]
      if (length(char_cols) == 0) char_cols <- names(rv$cp_data) # Fallback to all if no char cols
      char_cols <- c("(None)", char_cols)
      updateSelectInput(session, "cp_filter_col", choices = char_cols, selected = "(None)")
    })

    # Update Filter Values when Filter Column Changes
    observeEvent(c(input$cp_filter_col, rv$cp_data), {
      req(rv$cp_data, input$cp_filter_col)
      if (input$cp_filter_col %in% names(rv$cp_data)) {
        vals <- unique(rv$cp_data[[input$cp_filter_col]])
        vals <- vals[!is.na(vals)] # Exclude NA
        vals <- sort(vals)
        updatePickerInput(session, "cp_filter_values", choices = vals, selected = vals)
      }
    })

    # Reference Data UI
    output$cp_VARG26_selector <- renderUI({
      req(input$cp_use_VARG26)
      data <- load_VARG26()
      req(data)
      filter_col <- input$cp_VARG26_filter_col
      req(filter_col, filter_col %in% names(data))
      
      col_vals <- as.character(data[[filter_col]])
      counts_table <- table(col_vals, useNA = "no")
      counts_df <- data.frame(
        group = names(counts_table),
        n = as.integer(counts_table),
        stringsAsFactors = FALSE
      )
      counts_df <- counts_df[order(-counts_df$n), ]
      
      choices <- setNames(counts_df$group, paste0(counts_df$group, " (n=", counts_df$n, ")"))
      pickerInput(ns("cp_VARG26_groups"), "Select Groups:", choices = choices, multiple = TRUE,
        options = list(`actions-box` = TRUE, `live-search` = TRUE, `selected-text-format` = "count > 3",
          `count-selected-text` = "{0} groups selected"))
    })

    # TAS Column Selectors UI
    output$cp_tas_sio2_selector <- renderUI({
      req(rv$cp_data, input$cp_show_tas)
      col_names <- names(rv$cp_data)
      # Find all columns that could be SiO2 (starts with SiO2)
      sio2_cols <- col_names[grepl("^SiO2", col_names, ignore.case = TRUE)]
      if (length(sio2_cols) == 0) {
        return(selectInput(ns("cp_tas_sio2_col"), "SiO₂ Column:", choices = col_names, selected = NULL))
      }
      selectInput(ns("cp_tas_sio2_col"), "SiO₂ Column:", choices = sio2_cols, selected = sio2_cols[1])
    })

    output$cp_tas_na2o_selector <- renderUI({
      req(rv$cp_data, input$cp_show_tas)
      col_names <- names(rv$cp_data)
      # Find all columns that could be Na2O (starts with Na2O)
      na2o_cols <- col_names[grepl("^Na2O", col_names, ignore.case = TRUE)]
      if (length(na2o_cols) == 0) {
        return(selectInput(ns("cp_tas_na2o_col"), "Na₂O Column:", choices = col_names, selected = NULL))
      }
      selectInput(ns("cp_tas_na2o_col"), "Na₂O Column:", choices = na2o_cols, selected = na2o_cols[1])
    })

    output$cp_tas_k2o_selector <- renderUI({
      req(rv$cp_data, input$cp_show_tas)
      col_names <- names(rv$cp_data)
      # Find all columns that could be K2O (starts with K2O)
      k2o_cols <- col_names[grepl("^K2O", col_names, ignore.case = TRUE)]
      if (length(k2o_cols) == 0) {
        return(selectInput(ns("cp_tas_k2o_col"), "K₂O Column:", choices = col_names, selected = NULL))
      }
      selectInput(ns("cp_tas_k2o_col"), "K₂O Column:", choices = k2o_cols, selected = k2o_cols[1])
    })

    # Continuous Color Palette UI (shows only when numeric column selected for color)
    output$cp_continuous_palette_ui <- renderUI({
      req(rv$cp_data, input$cp_color_var)
      if (input$cp_color_var != "Fixed Color" && is.numeric(rv$cp_data[[input$cp_color_var]])) {
        div(
          style = "margin-top: 10px;",
          selectInput(ns("cp_continuous_palette"), "Color Palette:",
            choices = c(
              "Viridis" = "viridis",
              "Plasma" = "plasma",
              "Inferno" = "inferno",
              "Magma" = "magma",
              "Cividis" = "cividis",
              "Blue-Red" = "RdBu",
              "Yellow-Orange-Red" = "YlOrRd",
              "Blue-Green" = "BuGn",
              "Purple-Blue" = "PuBu"
            ),
            selected = "viridis"
          )
        )
      }
    })

    # Legend Title UI - dynamic based on whether color and shape use same variable
    output$cp_legend_title_ui <- renderUI({
      color_var <- input$cp_color_var
      shape_var <- input$cp_shape_var

      color_active <- !is.null(color_var) && color_var != "Fixed Color"
      shape_active <- !is.null(shape_var) && shape_var != "Fixed Shape"

      # Don't show if both are fixed
      if (!color_active && !shape_active) {
        return(NULL)
      }

      # Check if using same variable for both (unified legend)
      same_var <- color_active && shape_active && color_var == shape_var

      if (same_var) {
        # Single unified legend title input
        div(
          style = "max-width: 50%;",
          textInput(ns("cp_legend_title"), "Legend Title:", placeholder = paste0("e.g., ", color_var))
        )
      } else if (color_active && shape_active) {
        # Both active but different - show both inputs side by side
        splitLayout(
          cellWidths = c("50%", "50%"),
          textInput(ns("cp_legend_color_title"), "Color Legend Title:", placeholder = paste0("e.g., ", color_var)),
          textInput(ns("cp_legend_shape_title"), "Shape Legend Title:", placeholder = paste0("e.g., ", shape_var))
        )
      } else if (color_active) {
        # Only color active
        div(
          style = "max-width: 50%;",
          textInput(ns("cp_legend_color_title"), "Color Legend Title:", placeholder = paste0("e.g., ", color_var))
        )
      } else {
        # Only shape active
        div(
          style = "max-width: 50%;",
          textInput(ns("cp_legend_shape_title"), "Shape Legend Title:", placeholder = paste0("e.g.", shape_var))
        )
      }
    })


    # Custom symbology panel - renders when at least one discrete aesthetic is active
    output$cp_custom_symbology_panel <- renderUI({
      if (!isTRUE(input$cp_use_custom_symbology)) {
        return(NULL)
      }

      req(rv$cp_data)
      color_var <- input$cp_color_var
      shape_var <- input$cp_shape_var

      color_in_data <- !is.null(color_var) && color_var %in% names(rv$cp_data)
      shape_in_data <- !is.null(shape_var) && shape_var %in% names(rv$cp_data)

      color_enabled <- color_in_data && color_var != "Fixed Color" && !is.numeric(rv$cp_data[[color_var]])
      shape_enabled <- shape_in_data && shape_var != "Fixed Shape"

      color_singleton <- color_enabled && length(unique(na.omit(rv$cp_data[[color_var]]))) == 1
      shape_singleton <- shape_enabled && length(unique(na.omit(rv$cp_data[[shape_var]]))) == 1

      show_color_tools <- color_enabled || color_singleton || (!color_in_data && color_var == "Fixed Color")
      show_shape_tools <- shape_enabled || shape_singleton || (!shape_in_data && shape_var == "Fixed Shape")

      if (!show_color_tools && !show_shape_tools) {
        return(div(class = "text-muted small", "Select a categorical Color or Shape to customize."))
      }

      palettes <- get_brewer_qualitative_palettes()
      current_palette <- input$cp_palette_select
      if (is.null(current_palette)) current_palette <- "Set1"

      tagList(
        div(
          style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; background-color: #f8f9fa;",
          if (color_enabled) {
            selectInput(ns("cp_palette_select"), "Color Palette:",
              choices = palettes, selected = current_palette
            )
          },
          uiOutput(ns("cp_group_controls_inner"))
        )
      )
    })

    # Inner UI for group-by-group controls (regenerates when palette, vars, or reset changes)
    output$cp_group_controls_inner <- renderUI({
      req(input$cp_use_custom_symbology, rv$cp_data)
      # Force re-render on reset button clicks
      input$cp_reset_overrides

      color_var <- input$cp_color_var
      shape_var <- input$cp_shape_var

      color_in_data <- !is.null(color_var) && color_var %in% names(rv$cp_data)
      shape_in_data <- !is.null(shape_var) && shape_var %in% names(rv$cp_data)

      color_enabled <- color_in_data && color_var != "Fixed Color" && !is.numeric(rv$cp_data[[color_var]])
      shape_enabled <- shape_in_data && shape_var != "Fixed Shape"

      color_singleton <- color_enabled && length(unique(na.omit(rv$cp_data[[color_var]]))) == 1
      shape_singleton <- shape_enabled && length(unique(na.omit(rv$cp_data[[shape_var]]))) == 1

      show_color_tools <- color_enabled || color_singleton || (!color_in_data && color_var == "Fixed Color")
      show_shape_tools <- shape_enabled || shape_singleton || (!shape_in_data && shape_var == "Fixed Shape")

      if (!show_color_tools && !show_shape_tools) {
        return(NULL)
      }

      same_var <- color_enabled && shape_enabled && identical(color_var, shape_var)

      current_palette <- input$cp_palette_select
      if (is.null(current_palette)) current_palette <- "Set1"

      shape_choices <- get_shape_choices()

      # Helper to render a table of groups with colors and/or shapes
      render_group_table <- function(groups, show_color = TRUE, show_shape = TRUE) {
        if (length(groups) == 0) {
          return(NULL)
        }

        colors <- if (show_color) generate_group_colors(groups, current_palette) else NULL
        shapes <- if (show_shape) generate_group_shapes(groups) else NULL

        header_cols <- tagList(
          div(class = "col-5", "Group"),
          if (show_color) div(class = if (show_shape) "col-4" else "col-7", "Color"),
          if (show_shape) div(class = if (show_color) "col-3" else "col-7", "Shape")
        )

        row_ui <- lapply(seq_along(groups), function(i) {
          grp <- groups[i]
          grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
          grp_color <- if (show_color && !is.null(colors[[grp]])) colors[[grp]] else "#ccc"
          input_id <- ns(paste0("cp_grp_color_", grp_id))

          div(
            class = "row py-1 border-bottom",
            style = "margin-left: 0; margin-right: 0; align-items: center;",
            div(class = "col-5 small text-truncate", title = grp,
              if (show_color) tags$span(
                class = "cp-swatch",
                `data-for` = input_id,
                style = sprintf("display:inline-block;width:12px;height:12px;border-radius:2px;background:%s;border:1px solid #999;margin-right:4px;vertical-align:middle;cursor:pointer;", grp_color)
              ),
              grp
            ),
            if (show_color) {
              div(
                class = if (show_shape) "col-4" else "col-7",
                colourpicker::colourInput(
                  inputId = input_id,
                  label = NULL,
                  value = colors[[grp]],
                  showColour = "both",
                  palette = "square"
                )
              )
            },
            if (show_shape) {
              div(
                class = if (show_color) "col-3" else "col-7",
                selectInput(
                  inputId = ns(paste0("cp_grp_shape_", grp_id)),
                  label = NULL,
                  choices = shape_choices,
                  selected = as.character(shapes[[grp]])
                )
              )
            }
          )
        })

        tagList(
          div(
            class = "row small fw-bold border-bottom pb-1 mb-2",
            style = "margin-left: 0; margin-right: 0;",
            header_cols
          ),
          tags$style(HTML("
            .cp-custom-scroll-container { 
              max-height: 200px; 
              overflow-y: auto; 
              overflow-x: hidden;
            }
            .cp-custom-scroll-container .colourpicker-input-container { 
              overflow: visible !important; 
            }
            .cp-custom-scroll-container .shiny-input-container { 
              overflow: visible !important; 
            }
          ")),
          tags$script(HTML("
            $(function() {
              // Click swatch to open its associated colourpicker
              $(document).off('click.cpSwatch').on('click.cpSwatch', '.cp-swatch', function() {
                var forId = $(this).data('for');
                if (forId) {
                  var $el = $(document.getElementById(forId));
                  if ($el.length) {
                    // colourpicker binds to the input; trigger focus then click
                    $el.trigger('focus').trigger('click');
                  }
                }
              });
              // Sync swatch color when colourpicker value changes
              $(document).off('change.cpSwatch').on('change.cpSwatch', '.colourpicker-input', function() {
                var id = $(this).attr('id');
                var val = $(this).val();
                if (id && val) {
                  var hex = val.startsWith('#') ? val : '#' + val;
                  $('.cp-swatch[data-for=\"' + id + '\"]').css('background', hex);
                }
              });
            });
          ")),
          div(
            class = "cp-custom-scroll-container",
            row_ui
          )
        )
      }

      color_groups <- if (color_enabled) sort(unique(as.character(rv$cp_data[[color_var]]))) else character()
      shape_groups <- if (shape_enabled) sort(unique(as.character(rv$cp_data[[shape_var]]))) else character()

      color_section <- NULL
      shape_section <- NULL
      fixed_color_section <- NULL
      fixed_shape_section <- NULL

      if (same_var) {
        color_section <- div(
          class = "mb-3",
          tags$p(class = "text-muted small mb-1", "Click any color swatch to change it."),
          render_group_table(color_groups, show_color = TRUE, show_shape = TRUE)
        )
      } else {
        if (color_enabled && length(color_groups) > 0) {
          color_section <- div(
            class = "mb-3",
            tags$h6(sprintf("Color overrides (%s)", color_var), class = "fw-bold"),
            tags$p(class = "text-muted small mb-1", "Click any color swatch to change it."),
            render_group_table(color_groups, show_color = TRUE, show_shape = FALSE)
          )
        }

        if (shape_enabled && length(shape_groups) > 0) {
          shape_section <- div(
            class = "mb-3",
            tags$h6(sprintf("Shape overrides (%s)", shape_var), class = "fw-bold"),
            render_group_table(shape_groups, show_color = FALSE, show_shape = TRUE)
          )
        }
      }

      if (!color_enabled && show_color_tools) {
        cur_fixed <- if (!is.null(input$cp_fixed_color_override)) input$cp_fixed_color_override else "#e74c3c"
        fixed_input_id <- ns("cp_fixed_color_override")
        fixed_color_section <- div(
          class = "mb-3",
          tags$h6("Fixed color override", class = "fw-bold"),
          tags$p(class = "text-muted small mb-1", "Click the swatch or hex code to change the color."),
          div(
            style = "display:flex;align-items:center;gap:8px;",
            tags$span(
              class = "cp-swatch",
              `data-for` = fixed_input_id,
              style = sprintf("display:inline-block;width:24px;height:24px;border-radius:4px;background:%s;border:1px solid #999;flex-shrink:0;cursor:pointer;", cur_fixed)
            ),
            colourpicker::colourInput(
              inputId = fixed_input_id,
              label = NULL,
              value = cur_fixed,
              showColour = "both",
              palette = "square"
            )
          )
        )
      }

      if (!shape_enabled && show_shape_tools) {
        fixed_shape_section <- div(
          class = "mb-3",
          tags$h6("Fixed shape override", class = "fw-bold"),
          selectInput(
            inputId = ns("cp_fixed_shape_override"),
            label = NULL,
            choices = shape_choices,
            selected = if (!is.null(input$cp_fixed_shape_override)) input$cp_fixed_shape_override else "16"
          )
        )
      }

      tagList(
        color_section,
        shape_section,
        fixed_color_section,
        fixed_shape_section,
        div(
          style = "margin-top: 10px; text-align: right;",
          actionButton(ns("cp_reset_overrides"), "Reset to Defaults", class = "btn-sm btn-outline-secondary")
        )
      )
    })

    # --- HDR Color Customization Panel ---
    output$cp_hdr_color_panel <- renderUI({
      req(input$cp_use_custom_hdr_colors)
      
      VARG26_groups <- character()
      ref_groups <- character()
      
      # Get VARG26 groups if enabled
      if (input$cp_use_VARG26 && !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0) {
        VARG26_groups <- input$cp_VARG26_groups
      }
      
      # Get custom reference groups if enabled
      if (!is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0) {
        ref_groups <- input$cp_ref_groups
      }
      
      if (length(VARG26_groups) == 0 && length(ref_groups) == 0) {
        return(div(class = "text-muted small", "Select HDR groups to customize colors."))
      }
      
      # Generate default colors
      all_hdr_groups <- c(VARG26_groups, ref_groups)
      n_groups <- length(all_hdr_groups)
      
      # Use a palette similar to the discrete color palette logic
      if (n_groups <= 8) {
        hdr_palette <- RColorBrewer::brewer.pal(min(n_groups, 8), "Dark2")
      } else if (n_groups <= 20) {
        tableau20 <- c(
          "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
          "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
          "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
          "#17BECF", "#9EDAE5"
        )
        hdr_palette <- tableau20[1:n_groups]
      } else {
        seed_colors <- c(
          "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
          "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
          "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
          "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
        )
        hdr_palette <- colorRampPalette(seed_colors)(n_groups)
      }
      
      hdr_colors <- setNames(hdr_palette, all_hdr_groups)
      
      div(
        style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; background-color: #f8f9fa;",
        
        # VARG26 section
        if (length(VARG26_groups) > 0) {
          tagList(
            tags$h6("VARG26 HDR Colors", class = "fw-bold mb-2"),
            lapply(VARG26_groups, function(grp) {
              grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
              div(
                class = "row py-1 align-items-center",
                style = "margin-left: 0; margin-right: 0;",
                div(class = "col-7 small", grp),
                div(
                  class = "col-5",
                  colourpicker::colourInput(
                    inputId = ns(paste0("cp_hdr_VARG26_color_", grp_id)),
                    label = NULL,
                    value = hdr_colors[[grp]],
                    showColour = "background",
                    palette = "square"
                  )
                )
              )
            }),
            if (length(ref_groups) > 0) tags$hr()
          )
        },
        
        # Custom reference section
        if (length(ref_groups) > 0) {
          tagList(
            tags$h6("Custom Reference HDR Colors", class = "fw-bold mb-2"),
            lapply(ref_groups, function(grp) {
              grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
              div(
                class = "row py-1 align-items-center",
                style = "margin-left: 0; margin-right: 0;",
                div(class = "col-7 small", grp),
                div(
                  class = "col-5",
                  colourpicker::colourInput(
                    inputId = ns(paste0("cp_hdr_ref_color_", grp_id)),
                    label = NULL,
                    value = hdr_colors[[grp]],
                    showColour = "background",
                    palette = "square"
                  )
                )
              )
            })
          )
        },
        
        # Reset button
        div(
          style = "margin-top: 10px; text-align: right;",
          actionButton(ns("cp_reset_hdr_colors"), "Reset HDR Colors", class = "btn-sm btn-outline-secondary")
        )
      )
    })

    # --- Sheet Selection UI for Custom Reference Data ---
    output$cp_ref_sheet_select <- renderUI({
      req(input$cp_ref_file)
      file_ext <- tolower(tools::file_ext(input$cp_ref_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        sheets <- readxl::excel_sheets(input$cp_ref_file$datapath)
        selectInput(ns("cp_ref_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
      }
    })

    observeEvent(input$cp_ref_file, {
      req(input$cp_ref_file)
      file_ext <- tolower(tools::file_ext(input$cp_ref_file$name))
      if (file_ext %in% c("xlsx", "xls")) {
        # Wait for sheet selection
        return()
      }
      tryCatch(
        {
          replace_cp_reference_data(read.csv(input$cp_ref_file$datapath, stringsAsFactors = FALSE, check.names = FALSE))
          updateSelectInput(session, "cp_ref_filter_col", choices = names(rv$cp_ref_data))

          # Check if pairs have been imported and apply to new reference data
          if (!is.null(rv$cp_imported_pairs)) {
            top_pairs <- rv$cp_imported_pairs
            suffix <- input$cp_pair_suffix
            if (!is.null(suffix) && suffix == "raw") suffix <- ""

            ref_computed_cols <- list()
            for (i in 1:nrow(top_pairs)) {
              v1 <- compute_variable(top_pairs$Variable_1[i], rv$cp_ref_data, suffix)
              v2 <- compute_variable(top_pairs$Variable_2[i], rv$cp_ref_data, suffix)
              if (!is.null(v1)) ref_computed_cols[[v1$col_name]] <- v1$values
              if (!is.null(v2)) ref_computed_cols[[v2$col_name]] <- v2$values
            }
            if (length(ref_computed_cols) > 0) {
              rv$cp_ref_data <- cbind(rv$cp_ref_data, as.data.frame(ref_computed_cols))
              showNotification("Applied derived variables to reference data.", type = "message")
            }
          }

          showNotification("Reference data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading reference file:", e$message), type = "error")
      )
    })

    observeEvent(input$cp_ref_sheet, {
      req(input$cp_ref_file, input$cp_ref_sheet)
      tryCatch(
        {
          replace_cp_reference_data(readxl::read_excel(input$cp_ref_file$datapath, sheet = input$cp_ref_sheet))
          updateSelectInput(session, "cp_ref_filter_col", choices = names(rv$cp_ref_data))

          # Check if pairs have been imported and apply to new reference data
          if (!is.null(rv$cp_imported_pairs)) {
            top_pairs <- rv$cp_imported_pairs
            suffix <- input$cp_pair_suffix
            if (!is.null(suffix) && suffix == "raw") suffix <- ""

            ref_computed_cols <- list()
            for (i in 1:nrow(top_pairs)) {
              v1 <- compute_variable(top_pairs$Variable_1[i], rv$cp_ref_data, suffix)
              v2 <- compute_variable(top_pairs$Variable_2[i], rv$cp_ref_data, suffix)
              if (!is.null(v1)) ref_computed_cols[[v1$col_name]] <- v1$values
              if (!is.null(v2)) ref_computed_cols[[v2$col_name]] <- v2$values
            }
            if (length(ref_computed_cols) > 0) {
              rv$cp_ref_data <- cbind(rv$cp_ref_data, as.data.frame(ref_computed_cols))
              showNotification("Applied derived variables to reference data.", type = "message")
            }
          }

          showNotification("Reference data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading reference file:", e$message), type = "error")
      )
    })

    # Import processed data as custom reference
    observeEvent(input$cp_ref_use_processed, {
      if (is.null(processed_data) || !is.function(processed_data)) {
        showNotification("No processed data available in this session.", type = "warning")
        return()
      }
      pd <- tryCatch(processed_data(), error = function(e) NULL)
      if (is.null(pd) || NROW(pd) == 0) {
        showNotification("Processed data is empty or unavailable.", type = "warning")
        return()
      }
      replace_cp_reference_data(pd)
      updateSelectInput(session, "cp_ref_filter_col", choices = names(rv$cp_ref_data))

      # Apply imported pairs if available
      if (!is.null(rv$cp_imported_pairs)) {
        top_pairs <- rv$cp_imported_pairs
        suffix <- input$cp_pair_suffix
        if (!is.null(suffix) && suffix == "raw") suffix <- ""

        ref_computed_cols <- list()
        for (i in 1:nrow(top_pairs)) {
          v1 <- compute_variable(top_pairs$Variable_1[i], rv$cp_ref_data, suffix)
          v2 <- compute_variable(top_pairs$Variable_2[i], rv$cp_ref_data, suffix)
          if (!is.null(v1)) ref_computed_cols[[v1$col_name]] <- v1$values
          if (!is.null(v2)) ref_computed_cols[[v2$col_name]] <- v2$values
        }
        if (length(ref_computed_cols) > 0) {
          rv$cp_ref_data <- cbind(rv$cp_ref_data, as.data.frame(ref_computed_cols))
          showNotification("Applied derived variables to reference data.", type = "message")
        }
      }

      showNotification(paste0("Processed data imported as reference (", NROW(pd), " rows)."), type = "message")
    })

    output$cp_ref_selector <- renderUI({
      req(rv$cp_ref_data, input$cp_ref_filter_col)
      col_vals <- as.character(rv$cp_ref_data[[input$cp_ref_filter_col]])
      counts_table <- table(col_vals, useNA = "no")
      counts_df <- data.frame(
        group = names(counts_table),
        n = as.integer(counts_table),
        stringsAsFactors = FALSE
      )
      counts_df <- counts_df[order(-counts_df$n), ]
      choices <- setNames(counts_df$group, paste0(counts_df$group, " (n=", counts_df$n, ")"))
      pickerInput(ns("cp_ref_groups"), "Select Groups:", choices = choices, multiple = TRUE,
        options = list(`actions-box` = TRUE, `live-search` = TRUE, `selected-text-format` = "count > 3",
          `count-selected-text` = "{0} groups selected"))
    })

    # --- Group Customization Observers ---

    # Reset button for custom symbology - triggers UI regeneration (palette change will reload defaults)
    observeEvent(input$cp_reset_overrides, {
      # Clear any stored overrides (not used anymore, but kept for safety)
      rv$cp_group_color_overrides <- NULL
      rv$cp_group_shape_overrides <- NULL
      # Reset palette selection if present; renderUI will refresh inputs
      if (!is.null(input$cp_palette_select)) {
        updateSelectInput(session, "cp_palette_select", selected = "Set1")
      }
      if (!is.null(input$cp_fixed_color_override)) {
        colourpicker::updateColourInput(session, "cp_fixed_color_override", value = "#e74c3c")
      }
      if (!is.null(input$cp_fixed_shape_override)) {
        updateSelectInput(session, "cp_fixed_shape_override", selected = "16")
      }
      showNotification("Colors and shapes reset to palette defaults.", type = "message")
    })

    # Reset button for HDR colors - forces panel to regenerate with default colors
    observeEvent(input$cp_reset_hdr_colors, {
      # Force re-render by toggling the checkbox off and back on
      updateCheckboxInput(session, "cp_use_custom_hdr_colors", value = FALSE)
      Sys.sleep(0.1)
      updateCheckboxInput(session, "cp_use_custom_hdr_colors", value = TRUE)
      showNotification("HDR colors reset to defaults.", type = "message")
    })

    # TAS Auto-Detection - when TAS overlay is activated, detect and set appropriate columns
    observeEvent(input$cp_show_tas, {
      if (input$cp_show_tas && !is.null(rv$cp_data)) {
        data <- rv$cp_data
        col_names <- names(data)
        
        # Find SiO2 column (starts with SiO2)
        sio2_col <- col_names[grepl("^SiO2", col_names, ignore.case = TRUE)]
        sio2_col <- if (length(sio2_col) > 0) sio2_col[1] else NULL
        
        # Find Na2O column
        na2o_col <- col_names[grepl("^Na2O", col_names, ignore.case = TRUE)]
        na2o_col <- if (length(na2o_col) > 0) na2o_col[1] else NULL
        
        # Find K2O column
        k2o_col <- col_names[grepl("^K2O", col_names, ignore.case = TRUE)]
        k2o_col <- if (length(k2o_col) > 0) k2o_col[1] else NULL
        
        # Check if we need to create total alkali column
        alkali_col <- col_names[grepl("^(Na2O.*K2O|K2O.*Na2O|Total.*Alk|Alk)", col_names, ignore.case = TRUE)]
        alkali_col <- if (length(alkali_col) > 0) alkali_col[1] else NULL
        
        # Create total alkali if we have both Na2O and K2O but no alkali column
        if (!is.null(na2o_col) && !is.null(k2o_col) && is.null(alkali_col)) {
          tryCatch({
            rv$cp_data <- rv$cp_data %>%
              mutate(Na2O_K2O_Total = .data[[na2o_col]] + .data[[k2o_col]])
            alkali_col <- "Na2O_K2O_Total"
            
            # Force UI update with new column
            numeric_cols <- names(rv$cp_data)[sapply(rv$cp_data, is.numeric)]
            updateSelectInput(session, "cp_y_var", choices = numeric_cols, selected = alkali_col)
            
            showNotification("Created 'Na2O_K2O_Total' column for TAS diagram.", type = "message", duration = 4)
          }, error = function(e) {
            showNotification(paste("Could not create total alkali column:", e$message), type = "warning")
          })
        }
        
        # Update X and Y axes if appropriate columns found
        if (!is.null(sio2_col)) {
          updateSelectInput(session, "cp_x_var", selected = sio2_col)
        }
        if (!is.null(alkali_col)) {
          updateSelectInput(session, "cp_y_var", selected = alkali_col)
        }
        
        # Notify user of what was detected
        if (!is.null(sio2_col) && !is.null(alkali_col)) {
          showNotification(
            sprintf("TAS axes auto-detected: X=%s, Y=%s", sio2_col, alkali_col),
            type = "message", duration = 5
          )
        } else if (is.null(sio2_col) || is.null(alkali_col)) {
          showNotification(
            "Could not auto-detect all TAS columns. Please set X-axis to SiO₂ and Y-axis to Na₂O+K₂O manually.",
            type = "warning", duration = 7
          )
        }
      }
    })

    # TAS Assign Axes Button - user explicitly selects which columns to use
    observeEvent(input$cp_assign_tas_axes, {
      req(rv$cp_data, input$cp_tas_sio2_col, input$cp_tas_na2o_col, input$cp_tas_k2o_col)
      
      sio2_col <- input$cp_tas_sio2_col
      na2o_col <- input$cp_tas_na2o_col
      k2o_col <- input$cp_tas_k2o_col
      
      # Check for existing alkali total column (various patterns)
      alkali_patterns <- c("^(Na2O.*K2O|K2O.*Na2O|Total.*Alk|Alk)", 
                          paste0("^", sub("^Na2O", "", na2o_col), ".*", sub("^K2O", "", k2o_col), "$"))
      col_names <- names(rv$cp_data)
      alkali_col <- col_names[grepl(alkali_patterns[1], col_names, ignore.case = TRUE)]
      alkali_col <- if (length(alkali_col) > 0) alkali_col[1] else NULL
      
      # Create total alkali column if both Na2O and K2O provided
      if (!is.null(na2o_col) && !is.null(k2o_col)) {
        if (is.null(alkali_col)) {
          tryCatch({
            rv$cp_data <- rv$cp_data %>%
              mutate(Na2O_K2O_Total = .data[[na2o_col]] + .data[[k2o_col]])
            alkali_col <- "Na2O_K2O_Total"
            showNotification("Created 'Na2O_K2O_Total' column.", type = "message", duration = 3)
          }, error = function(e) {
            showNotification(paste("Error creating alkali column:", e$message), type = "error")
            return()
          })
        }
      }
      
      # Update axes
      numeric_cols <- names(rv$cp_data)[sapply(rv$cp_data, is.numeric)]
      updateSelectInput(session, "cp_x_var", choices = numeric_cols, selected = sio2_col)
      if (!is.null(alkali_col)) {
        updateSelectInput(session, "cp_y_var", choices = numeric_cols, selected = alkali_col)
      }
      
      showNotification(
        sprintf("TAS axes assigned: X=%s, Y=%s", sio2_col, if (!is.null(alkali_col)) alkali_col else ""),
        type = "message", duration = 4
      )
    })


    # Plot Logic - Fully Reactive
    observe({
      req(rv$cp_data, input$cp_x_var, input$cp_y_var)

      tryCatch(
        {
          # Base plot
          p <- ggplot()

          # Use rv$cp_data as the base filtered data
          filtered_data <- rv$cp_data

          # Apply Filter
          if (!is.null(input$cp_filter_col) && !is.null(input$cp_filter_values) && input$cp_filter_col %in% names(filtered_data)) {
            filtered_data <- filtered_data %>%
              filter(.data[[input$cp_filter_col]] %in% input$cp_filter_values)
          }

          # Prepare VARG26 data and labels (but don't add labels yet)
          VARG26_label_data <- NULL
          if (input$cp_use_VARG26 && !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0) {
            VARG26_data <- load_VARG26()
            req(VARG26_data)

            # Check max groups
            if (length(input$cp_VARG26_groups) > 15) {
              showNotification("Warning: More than 15 groups selected. Performance may be slow.",
                type = "warning", duration = 5
              )
            }

            # Filter VARG26 for selected groups
            VARG26_filtered <- VARG26_data %>%
              filter(.data[[input$cp_VARG26_filter_col]] %in% input$cp_VARG26_groups)

            # Map user column names to VARG26 column names
            VARG26_x_col <- input$cp_x_var
            VARG26_y_col <- input$cp_y_var

            # Try to map if columns don't exist directly
            # VARG26 reference file uses VARG25 column names (e.g., UMAP_VARG25_2D_1)
            # while processed data uses VARG26 names (e.g., UMAP_VARG26_2D_1)
            if (!(VARG26_x_col %in% names(VARG26_filtered))) {
              # First try mapping VARG26 → VARG25 for UMAP columns
              VARG26_x_col <- gsub("VARG26", "VARG25", VARG26_x_col)
            }
            if (!(VARG26_x_col %in% names(VARG26_filtered))) {
              # Fallback: strip VARG prefix entirely
              VARG26_x_col <- gsub("_?VARG25_?|_?VARG26_?", "_", input$cp_x_var)
              VARG26_x_col <- gsub("^_|_$", "", VARG26_x_col)
            }

            if (!(VARG26_y_col %in% names(VARG26_filtered))) {
              # First try mapping VARG26 → VARG25 for UMAP columns
              VARG26_y_col <- gsub("VARG26", "VARG25", VARG26_y_col)
            }
            if (!(VARG26_y_col %in% names(VARG26_filtered))) {
              # Fallback: strip VARG prefix entirely
              VARG26_y_col <- gsub("_?VARG25_?|_?VARG26_?", "_", input$cp_y_var)
              VARG26_y_col <- gsub("^_|_$", "", VARG26_y_col)
            }

            # Check if mapped columns exist in VARG26
            if (VARG26_x_col %in% names(VARG26_filtered) && VARG26_y_col %in% names(VARG26_filtered)) {
              VARG26_filtered <- prepare_hdr_groups(
                VARG26_filtered, VARG26_x_col, VARG26_y_col,
                input$cp_VARG26_filter_col, "VARG26"
              )

              if (!is.null(VARG26_filtered)) {

              # Calculate expanded limits for KDE calculation
              VARG26_x_range <- range(VARG26_filtered[[VARG26_x_col]], na.rm = TRUE)
              VARG26_y_range <- range(VARG26_filtered[[VARG26_y_col]], na.rm = TRUE)
              VARG26_x_span <- diff(VARG26_x_range)
              VARG26_y_span <- diff(VARG26_y_range)

              # Expand by 30% on each side for KDE calculation
              VARG26_xlim <- c(
                VARG26_x_range[1] - VARG26_x_span * 0.3,
                VARG26_x_range[2] + VARG26_x_span * 0.3
              )
              VARG26_ylim <- c(
                VARG26_y_range[1] - VARG26_y_span * 0.3,
                VARG26_y_range[2] + VARG26_y_span * 0.3
              )

              # Add HDR contours (background)
              p <- p +
                ggdensity::geom_hdr(
                  data = VARG26_filtered,
                  aes(
                    x = .data[[VARG26_x_col]], y = .data[[VARG26_y_col]],
                    fill = .data[[input$cp_VARG26_filter_col]]
                  ),
                  xlim = VARG26_xlim,
                  ylim = VARG26_ylim,
                  probs = input$cp_hdr_prob,
                  method = hdr_method_spec(
                    input$cp_hdr_method, input$cp_hdr_adjust,
                    input$cp_hdr_bins, input$cp_hdr_auto_bins
                  ),
                  n = hdr_grid_resolution(input$cp_hdr_n),
                  alpha = 0.3,
                  show.legend = !input$cp_hide_hdr_legend
                )

              # Prepare HDR center labels (but don't add yet)
              if (input$cp_label_hdr) {
                VARG26_label_data <- VARG26_filtered %>%
                  group_by(!!sym(input$cp_VARG26_filter_col)) %>%
                  reframe(find_kde_peak(!!sym(VARG26_x_col), !!sym(VARG26_y_col))) %>%
                  ungroup()

                colnames(VARG26_label_data) <- c(input$cp_VARG26_filter_col, "x_label", "y_label")
              }
              }
            } else {
              showNotification("Selected axes not found in VARG26 data. HDR contours skipped.", type = "warning")
            }
          }

          # Prepare Custom Reference data and labels (but don't add labels yet)
          ref_label_data <- NULL
          if (!is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0) {
            req(rv$cp_ref_data)

            # Filter Reference for selected groups
            ref_filtered <- rv$cp_ref_data %>%
              filter(.data[[input$cp_ref_filter_col]] %in% input$cp_ref_groups)

            ref_filtered <- prepare_hdr_groups(
              ref_filtered, input$cp_x_var, input$cp_y_var,
              input$cp_ref_filter_col, "Custom reference"
            )

            if (!is.null(ref_filtered)) {

            # Calculate expanded limits for KDE calculation
            ref_x_range <- range(ref_filtered[[input$cp_x_var]], na.rm = TRUE)
            ref_y_range <- range(ref_filtered[[input$cp_y_var]], na.rm = TRUE)
            ref_x_span <- diff(ref_x_range)
            ref_y_span <- diff(ref_y_range)

            # Expand by 30% on each side for KDE calculation
            ref_xlim <- c(
              ref_x_range[1] - ref_x_span * 0.3,
              ref_x_range[2] + ref_x_span * 0.3
            )
            ref_ylim <- c(
              ref_y_range[1] - ref_y_span * 0.3,
              ref_y_range[2] + ref_y_span * 0.3
            )

            # Add HDR contours (background)
            p <- p +
              ggdensity::geom_hdr(
                data = ref_filtered,
                aes(
                  x = .data[[input$cp_x_var]], y = .data[[input$cp_y_var]],
                  fill = .data[[input$cp_ref_filter_col]]
                ),
                xlim = ref_xlim,
                ylim = ref_ylim,
                probs = input$cp_ref_hdr_prob,
                method = hdr_method_spec(
                  input$cp_ref_hdr_method, input$cp_ref_hdr_adjust,
                  input$cp_ref_hdr_bins, input$cp_ref_hdr_auto_bins
                ),
                n = hdr_grid_resolution(input$cp_ref_hdr_n),
                alpha = 0.3,
                show.legend = !input$cp_hide_hdr_legend
              )

            # Prepare HDR center labels (but don't add yet)
            if (input$cp_label_hdr) {
              ref_label_data <- ref_filtered %>%
                group_by(!!sym(input$cp_ref_filter_col)) %>%
                reframe(find_kde_peak(!!sym(input$cp_x_var), !!sym(input$cp_y_var))) %>%
                ungroup()

              colnames(ref_label_data) <- c(input$cp_ref_filter_col, "x_label", "y_label")
            }
            }
          }

          # Add user data points with dynamic aesthetics
          req(nrow(filtered_data) > 0)

          geom_args <- list(
            data = filtered_data,
            size = input$cp_point_size,
            alpha = input$cp_point_alpha,
            stroke = input$cp_stroke_width
          )

          # Build aesthetics mapping
          mapping_args <- list(x = rlang::sym(input$cp_x_var), y = rlang::sym(input$cp_y_var))

          if (input$cp_color_var != "Fixed Color") {
            # mapping_args$fill <- rlang::sym(input$cp_color_var) # Removed to avoid conflict with HDR fill
            mapping_args$colour <- rlang::sym(input$cp_color_var)
          }

          if (input$cp_shape_var != "Fixed Shape") {
            mapping_args$shape <- rlang::sym(input$cp_shape_var)
          }

          geom_args$mapping <- do.call(aes, mapping_args)

          # Add fixed aesthetics if not mapped
          if (input$cp_color_var == "Fixed Color") {
            fixed_col <- if (isTRUE(input$cp_use_custom_symbology) && !is.null(input$cp_fixed_color_override)) input$cp_fixed_color_override else "#e74c3c"
            geom_args$colour <- fixed_col
          }

          if (input$cp_shape_var == "Fixed Shape") {
            fixed_shape <- if (isTRUE(input$cp_use_custom_symbology) && !is.null(input$cp_fixed_shape_override)) as.integer(input$cp_fixed_shape_override) else 16
            geom_args$shape <- fixed_shape
          }


          p <- p + do.call(geom_point, geom_args)

          # NOW add labels (after points, so they avoid the data)
          # Prepare phantom points from user data to help ggrepel avoid them
          user_phantom <- data.frame(
            x = filtered_data[[input$cp_x_var]],
            y = filtered_data[[input$cp_y_var]],
            label = "",
            fill_group = NA_character_,
            stringsAsFactors = FALSE
          )

          # Add VARG26 labels if prepared
          if (!is.null(VARG26_label_data)) {
            # Combine label data with phantom user points
            filter_col <- input$cp_VARG26_filter_col
            VARG26_combined <- data.frame(
              x = c(VARG26_label_data$x_label, user_phantom$x),
              y = c(VARG26_label_data$y_label, user_phantom$y),
              label = c(VARG26_label_data[[filter_col]], user_phantom$label),
              fill_group = c(VARG26_label_data[[filter_col]], user_phantom$fill_group),
              stringsAsFactors = FALSE
            )

            p <- p +
              ggrepel::geom_label_repel(
                data = VARG26_combined,
                aes(x = x, y = y, label = label, fill = fill_group),
                alpha = 0.9,
                size = input$cp_label_size,
                fontface = "bold",
                box.padding = 1.0,
                point.padding = 3.0,
                min.segment.length = 0,
                max.overlaps = Inf,
                force = 5,
                force_pull = 0.2,
                show.legend = FALSE,
                na.rm = TRUE # Ignore phantom points with NA fill
              )
          }

          # Add custom reference labels if prepared
          if (!is.null(ref_label_data)) {
            # Combine label data with phantom user points
            ref_filter_col <- input$cp_ref_filter_col
            ref_combined <- data.frame(
              x = c(ref_label_data$x_label, user_phantom$x),
              y = c(ref_label_data$y_label, user_phantom$y),
              label = c(ref_label_data[[ref_filter_col]], user_phantom$label),
              fill_group = c(ref_label_data[[ref_filter_col]], user_phantom$fill_group),
              stringsAsFactors = FALSE
            )

            p <- p +
              ggrepel::geom_label_repel(
                data = ref_combined,
                aes(x = x, y = y, label = label, fill = fill_group),
                alpha = 0.9,
                size = input$cp_label_size,
                fontface = "bold",
                box.padding = 1.0,
                point.padding = 3.0,
                min.segment.length = 0,
                max.overlaps = Inf,
                force = 5,
                force_pull = 0.2,
                show.legend = FALSE,
                na.rm = TRUE # Ignore phantom points with NA fill
              )
          }

          # Add robust scales for many categories OR continuous scales for numeric data
          if (input$cp_color_var != "Fixed Color") {
            # Check if the color variable is numeric
            is_numeric_color <- is.numeric(filtered_data[[input$cp_color_var]])
            
            # Check if HDR regions are being displayed (they use discrete fill aesthetic)
            has_hdr_regions <- (isTRUE(input$cp_use_VARG26) && !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0) ||
                               (!is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0)

            if (is_numeric_color) {
              # CONTINUOUS COLOR SCALE
              palette_choice <- if (!is.null(input$cp_continuous_palette)) input$cp_continuous_palette else "viridis"

              # Apply viridis-style palettes
              if (palette_choice %in% c("viridis", "plasma", "inferno", "magma", "cividis")) {
                p <- p + scale_color_viridis_c(option = palette_choice)
                
                # Only apply continuous fill scale if NO HDR regions are displayed
                # HDR regions use discrete fill aesthetic which conflicts with continuous scales
                if (!has_hdr_regions) {
                  p <- p + scale_fill_viridis_c(option = palette_choice)
                }
              } else {
                # Apply Brewer sequential palettes
                p <- p + scale_color_distiller(palette = palette_choice, direction = 1)
                
                # Only apply continuous fill scale if NO HDR regions are displayed
                if (!has_hdr_regions) {
                  p <- p + scale_fill_distiller(palette = palette_choice, direction = 1)
                }
              }
            } else {
              # DISCRETE COLOR SCALE
              # Check if custom symbology is enabled - if so, read directly from dynamic inputs
              use_custom <- isTRUE(input$cp_use_custom_symbology)

              if (use_custom) {
                # Read colors directly from the dynamic color inputs
                groups <- sort(unique(as.character(filtered_data[[input$cp_color_var]])))
                custom_colors <- list()

                for (grp in groups) {
                  grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
                  color_val <- input[[paste0("cp_grp_color_", grp_id)]]
                  if (!is.null(color_val)) {
                    custom_colors[[grp]] <- color_val
                  }
                }

                if (length(custom_colors) > 0) {
                  p <- p + scale_color_manual(values = unlist(custom_colors))
                  
                  # Only apply fill scale if NOT using custom HDR colors
                  if (!isTRUE(input$cp_use_custom_hdr_colors)) {
                    p <- p + scale_fill_manual(values = unlist(custom_colors))
                  }
                }
              } else {
                # Use default automatic palettes
                all_color_values <- unique(filtered_data[[input$cp_color_var]])

                # Add VARG26 values if enabled
                if (input$cp_use_VARG26 && !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0) {
                  all_color_values <- c(all_color_values, input$cp_VARG26_groups)
                }

                # Add custom reference values if enabled
                if (!is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0) {
                  all_color_values <- c(all_color_values, input$cp_ref_groups)
                }

                n_colors <- length(unique(all_color_values))

                if (n_colors <= 8) {
                  p <- p + scale_color_brewer(palette = "Dark2")
                  
                  if (!isTRUE(input$cp_use_custom_hdr_colors)) {
                    p <- p + scale_fill_brewer(palette = "Dark2")
                  }
                } else if (n_colors <= 20) {
                  # Tableau 20 equivalent manual palette
                  tableau20 <- c(
                    "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A",
                    "#D62728", "#FF9896", "#9467BD", "#C5B0D5", "#8C564B", "#C49C94",
                    "#E377C2", "#F7B6D2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D",
                    "#17BECF", "#9EDAE5"
                  )
                  p <- p + scale_color_manual(values = tableau20)
                  
                  if (!isTRUE(input$cp_use_custom_hdr_colors)) {
                    p <- p + scale_fill_manual(values = tableau20)
                  }
                } else {
                  # Use extended palette for many colors
                  seed_colors <- c(
                    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
                    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
                    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
                  )
                  extended_palette <- colorRampPalette(seed_colors)(n_colors)
                  p <- p + scale_color_manual(values = extended_palette)
                  
                  if (!isTRUE(input$cp_use_custom_hdr_colors)) {
                    p <- p + scale_fill_manual(values = extended_palette)
                  }
                }
              }
            }
          }
          
          # Apply custom HDR colors if enabled
          if (isTRUE(input$cp_use_custom_hdr_colors)) {
            hdr_fill_colors <- list()
            
            # Collect VARG26 HDR colors
            if (input$cp_use_VARG26 && !is.null(input$cp_VARG26_groups) && length(input$cp_VARG26_groups) > 0) {
              for (grp in input$cp_VARG26_groups) {
                grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
                color_val <- input[[paste0("cp_hdr_VARG26_color_", grp_id)]]
                if (!is.null(color_val)) {
                  hdr_fill_colors[[grp]] <- color_val
                }
              }
            }
            
            # Collect custom reference HDR colors
            if (!is.null(input$cp_ref_groups) && length(input$cp_ref_groups) > 0) {
              for (grp in input$cp_ref_groups) {
                grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
                color_val <- input[[paste0("cp_hdr_ref_color_", grp_id)]]
                if (!is.null(color_val)) {
                  hdr_fill_colors[[grp]] <- color_val
                }
              }
            }
            
            # Apply custom HDR fill colors
            if (length(hdr_fill_colors) > 0) {
              p <- p + scale_fill_manual(values = unlist(hdr_fill_colors))
            }
          }

          if (input$cp_shape_var != "Fixed Shape") {
            # Check if custom symbology is enabled - if so, read directly from dynamic inputs
            use_custom <- isTRUE(input$cp_use_custom_symbology)

            if (use_custom) {
              # Read shapes directly from the dynamic shape inputs
              groups <- sort(unique(as.character(filtered_data[[input$cp_shape_var]])))
              custom_shapes <- list()

              for (grp in groups) {
                grp_id <- gsub("[^a-zA-Z0-9]", "_", grp)
                shape_val <- input[[paste0("cp_grp_shape_", grp_id)]]
                if (!is.null(shape_val)) {
                  custom_shapes[[grp]] <- as.integer(shape_val)
                }
              }

              if (length(custom_shapes) > 0) {
                p <- p + scale_shape_manual(values = unlist(custom_shapes))
              }
            } else {
              n_shapes <- length(unique(filtered_data[[input$cp_shape_var]]))
              # Curated distinct shapes
              distinct_shapes <- c(16, 15, 18, 17, 3, 4, 8, 10, 12, 13, 14)

              if (n_shapes > length(distinct_shapes)) {
                needed_shapes <- rep(distinct_shapes, length.out = n_shapes)
                p <- p + scale_shape_manual(values = needed_shapes)
              } else {
                p <- p + scale_shape_manual(values = distinct_shapes[1:n_shapes])
              }
            }
          }

          # Determine axis labels (use custom if provided, otherwise use variable names)
          # Convert Unicode subscripts/superscripts to plotmath for PDF compatibility
          x_label <- if (!is.null(input$cp_x_label) && nchar(trimws(input$cp_x_label)) > 0) {
            unicode_to_plotmath(input$cp_x_label)
          } else {
            input$cp_x_var
          }

          y_label <- if (!is.null(input$cp_y_label) && nchar(trimws(input$cp_y_label)) > 0) {
            unicode_to_plotmath(input$cp_y_label)
          } else {
            input$cp_y_var
          }

          # Determine if color and shape use the same variable (unified legend)
          same_legend_var <- input$cp_color_var != "Fixed Color" &&
            input$cp_shape_var != "Fixed Shape" &&
            input$cp_color_var == input$cp_shape_var

          # Determine legend titles (use custom if provided)
          if (same_legend_var) {
            # Unified legend - use single title input
            unified_title <- if (!is.null(input$cp_legend_title) && nchar(trimws(input$cp_legend_title)) > 0) {
              input$cp_legend_title
            } else {
              NULL # Use default (variable name)
            }
            color_legend <- unified_title
            shape_legend <- unified_title
          } else {
            # Separate legends
            color_legend <- if (!is.null(input$cp_legend_color_title) && nchar(trimws(input$cp_legend_color_title)) > 0) {
              input$cp_legend_color_title
            } else {
              NULL # Use default (variable name)
            }

            shape_legend <- if (!is.null(input$cp_legend_shape_title) && nchar(trimws(input$cp_legend_shape_title)) > 0) {
              input$cp_legend_shape_title
            } else {
              NULL # Use default (variable name)
            }
          }

          p <- p + labs(
            x = x_label,
            y = y_label,
            color = color_legend,
            fill = color_legend,
            shape = shape_legend
          )
          
          # Apply selected theme
          theme_fn <- switch(input$cp_theme %||% "bw",
            "bw" = theme_bw,
            "classic" = theme_classic,
            "minimal" = theme_minimal,
            "light" = theme_light,
            "dark" = theme_dark,
            "linedraw" = theme_linedraw,
            "gray" = theme_gray,
            "void" = theme_void,
            "base" = ggthemes::theme_base,
            "tufte" = ggthemes::theme_tufte,
            "few" = ggthemes::theme_few,
            theme_bw
          )
          
          p <- p +
            theme_fn(base_size = 14) +
            theme(
              legend.position = if (isTRUE(input$cp_hide_legend)) "none" else "right",
              plot.margin = margin(20, 20, 20, 20, "pt") # Add margins for HDR contours
            )

          # Set axis scales with optional log transform and reverse
          x_trans <- "identity"
          y_trans <- "identity"
          if (isTRUE(input$cp_log_x) && isTRUE(input$cp_reverse_x)) {
            x_trans <- scales::compose_trans("log10", "reverse")
          } else if (isTRUE(input$cp_log_x)) {
            x_trans <- "log10"
          } else if (isTRUE(input$cp_reverse_x)) {
            x_trans <- "reverse"
          }
          if (isTRUE(input$cp_log_y) && isTRUE(input$cp_reverse_y)) {
            y_trans <- scales::compose_trans("log10", "reverse")
          } else if (isTRUE(input$cp_log_y)) {
            y_trans <- "log10"
          } else if (isTRUE(input$cp_reverse_y)) {
            y_trans <- "reverse"
          }
          
          # Build axis scale arguments
          show_minor <- (input$cp_minor_tick_type %||% "none") == "outside"
          show_major <- (input$cp_major_tick_type %||% "outside") == "outside"
          x_breaks <- waiver()
          y_breaks <- waiver()
          x_minor_breaks <- if (show_minor) waiver() else NULL
          y_minor_breaks <- if (show_minor) waiver() else NULL
          minor_n <- if (show_minor) max(1, input$cp_minor_tick_n %||% 4) else 0
          
          # Per-side manual-bound flags
          x_has_min <- !is.null(input$cp_x_min) && !is.na(input$cp_x_min)
          x_has_max <- !is.null(input$cp_x_max) && !is.na(input$cp_x_max)
          y_has_min <- !is.null(input$cp_y_min) && !is.na(input$cp_y_min)
          y_has_max <- !is.null(input$cp_y_max) && !is.na(input$cp_y_max)
          
          # Determine effective axis range (manual limits override data range)
          x_data_range <- range(filtered_data[[input$cp_x_var]], na.rm = TRUE)
          y_data_range <- range(filtered_data[[input$cp_y_var]], na.rm = TRUE)
          x_eff_lo <- if (x_has_min) input$cp_x_min else x_data_range[1]
          x_eff_hi <- if (x_has_max) input$cp_x_max else x_data_range[2]
          y_eff_lo <- if (y_has_min) input$cp_y_min else y_data_range[1]
          y_eff_hi <- if (y_has_max) input$cp_y_max else y_data_range[2]
          
          # When TAS is active, extend auto-bounds to cover the standard TAS range
          # so that breaks and ticks are generated for the full visible area
          if (isTRUE(input$cp_show_tas)) {
            if (!x_has_min) x_eff_lo <- min(x_eff_lo, 40)
            if (!x_has_max) x_eff_hi <- max(x_eff_hi, 80)
            if (!y_has_min) y_eff_lo <- min(y_eff_lo, 0)
            if (!y_has_max) y_eff_hi <- max(y_eff_hi, 15)
          }
          
          x_interval <- input$cp_x_tick_interval
          y_interval <- input$cp_y_tick_interval
          
          if (!is.null(x_interval) && !is.na(x_interval) && x_interval > 0) {
            x_breaks <- seq(floor(x_eff_lo / x_interval) * x_interval,
                           ceiling(x_eff_hi / x_interval) * x_interval,
                           by = x_interval)
            # Compute minor breaks between each pair of major breaks
            if (minor_n > 0 && length(x_breaks) >= 2) {
              x_minor_step <- x_interval / (minor_n + 1)
              x_minor_breaks <- seq(min(x_breaks), max(x_breaks), by = x_minor_step)
              x_minor_breaks <- setdiff(x_minor_breaks, x_breaks)
            }
          } else if (x_eff_lo != x_data_range[1] || x_eff_hi != x_data_range[2]) {
            # Manual limits set but no custom interval: extend breaks to cover visible range
            x_breaks <- pretty(c(x_eff_lo, x_eff_hi))
          }
          if (!is.null(y_interval) && !is.na(y_interval) && y_interval > 0) {
            y_breaks <- seq(floor(y_eff_lo / y_interval) * y_interval,
                           ceiling(y_eff_hi / y_interval) * y_interval,
                           by = y_interval)
            if (minor_n > 0 && length(y_breaks) >= 2) {
              y_minor_step <- y_interval / (minor_n + 1)
              y_minor_breaks <- seq(min(y_breaks), max(y_breaks), by = y_minor_step)
              y_minor_breaks <- setdiff(y_minor_breaks, y_breaks)
            }
          } else if (y_eff_lo != y_data_range[1] || y_eff_hi != y_data_range[2]) {
            y_breaks <- pretty(c(y_eff_lo, y_eff_hi))
          }
          
          # Hide major ticks if set to "none"
          if (!show_major) {
            x_breaks <- NULL
            y_breaks <- NULL
          }
          
          # Per-axis expand: zero padding on sides with manual limits, small padding on auto sides
          x_expand <- expansion(mult = c(if (x_has_min) 0 else 0.02, if (x_has_max) 0 else 0.02))
          y_expand <- expansion(mult = c(if (y_has_min) 0 else 0.02, if (y_has_max) 0 else 0.02))
          
          p <- p +
            scale_x_continuous(
              expand = x_expand, trans = x_trans,
              breaks = x_breaks,
              minor_breaks = x_minor_breaks,
              guide = if (show_minor) guide_axis(minor.ticks = TRUE) else waiver()
            ) +
            scale_y_continuous(
              expand = y_expand, trans = y_trans,
              breaks = y_breaks,
              minor_breaks = y_minor_breaks,
              guide = if (show_minor) guide_axis(minor.ticks = TRUE) else waiver()
            )
          
          # Style minor ticks if enabled
          if (show_minor) {
            p <- p + theme(
              axis.minor.ticks.length = rel(0.5)
            )
          }
          
          # Hide major tick marks if "none"
          if (!show_major) {
            p <- p + theme(
              axis.ticks = element_blank(),
              axis.text = element_blank()
            )
          }

          # Apply manual axis limits (from Bounds sub-menu)
          xlim_vals <- c(input$cp_x_min, input$cp_x_max)
          ylim_vals <- c(input$cp_y_min, input$cp_y_max)
          has_limits <- !all(is.na(xlim_vals)) || !all(is.na(ylim_vals))
          
          if (has_limits) {
            p <- p + coord_cartesian(
              xlim = if (!all(is.na(xlim_vals))) xlim_vals else NULL,
              ylim = if (!all(is.na(ylim_vals))) ylim_vals else NULL,
              clip = "on"
            )
          } else {
            p <- p + coord_cartesian(clip = "off")
          }
          
          # Add TAS diagram overlay AFTER axis limits to prevent auto-scaling
          if (isTRUE(input$cp_show_tas)) {
            tryCatch({
              tas_pkg <- paste0("geo", "chem")
              if (requireNamespace(tas_pkg, quietly = TRUE)) {
                tas_fields <- isTRUE(input$cp_tas_fields)
                tas_labels <- isTRUE(input$cp_tas_labels)
                tas_diagram <- getExportedValue(tas_pkg, "tas_diagram")
                p <- p + tas_diagram(fields = tas_fields, labels = tas_labels)
              } else {
                showNotification("Install 'geochem' package to use TAS diagram overlay.", type = "warning", duration = 5)
              }
            }, error = function(e) {
              showNotification(paste("TAS overlay error:", e$message), type = "warning", duration = 5)
            })
            
            # Reapply custom axis labels to override TAS defaults (if user provided custom labels)
            # Convert Unicode subscripts/superscripts to plotmath for PDF compatibility
            if (!is.null(input$cp_x_label) && nchar(trimws(input$cp_x_label)) > 0) {
              p <- p + labs(x = unicode_to_plotmath(input$cp_x_label))
            }
            if (!is.null(input$cp_y_label) && nchar(trimws(input$cp_y_label)) > 0) {
              p <- p + labs(y = unicode_to_plotmath(input$cp_y_label))
            }
          }
          
          # Add software citation annotation if enabled
          if (isTRUE(input$cp_include_citation)) {
            # Get app version from global environment
            app_ver <- get0("APP_VERSION", envir = .GlobalEnv, inherits = TRUE, ifnotfound = "")
            
            # Check if specific VARG UMAP coordinates are being displayed to include VARG version
            # Only show VARG26 UMAP annotation for official VARG UMAP variables, not generic UMAP coordinates
            x_var <- input$cp_x_var
            y_var <- input$cp_y_var
            varg_umap_pattern <- "^UMAP_VARG25_(2D_[12]|1D)$"
            using_varg_umap <- grepl(varg_umap_pattern, x_var, ignore.case = FALSE) || 
                               grepl(varg_umap_pattern, y_var, ignore.case = FALSE)
            
            if (using_varg_umap) {
              citation_text <- paste0("Plot generated with VARG-Tools v", app_ver, " (VARG26 UMAP) | Bolton & Jensen, in prep.")
            } else {
              citation_text <- paste0("Plot generated with VARG-Tools v", app_ver, " | Bolton & Jensen, in prep.")
            }
            
            p <- p + annotate(
              "text",
              x = Inf,
              y = -Inf,
              label = citation_text,
              hjust = 1.02,
              vjust = -0.5,
              size = 3.5,
              color = "grey40",
              fontface = "italic"
            )
          }

          rv$cp_plot_obj <- p
        },
        error = function(e) {
          showNotification(paste("Plot error:", e$message), type = "error", duration = 10)
        }
      )
    })

    output$cp_main_plot <- renderPlot({
      if (is.null(rv$cp_data)) {
        plot.new()
        text(0.5, 0.5, "Upload data or use Processed Data\nto generate a scatter plot", 
             cex = 1.3, col = "gray60", font = 2)
        return()
      }
      req(rv$cp_plot_obj)
      print(rv$cp_plot_obj)
    })

    # Data Citations output - shows citations from VARG26 reference data
    output$cp_data_citations <- renderUI({
      # Only show citations when VARG26 is enabled with groups selected
      if (!isTRUE(input$cp_use_VARG26) || is.null(input$cp_VARG26_groups) || length(input$cp_VARG26_groups) == 0) {
        return(tags$p(class = "text-muted small", "Enable VARG26 reference data and select groups to see data citations."))
      }
      
      VARG26_data <- load_VARG26()
      if (is.null(VARG26_data)) {
        return(tags$p(class = "text-muted small", "VARG26 data not available."))
      }
      
      # Check if Data citation column exists
      if (!"Data citation" %in% names(VARG26_data)) {
        return(tags$p(class = "text-muted small", "Data citation column not found in VARG26 data."))
      }
      
      # Filter to selected groups
      filter_col <- input$cp_VARG26_filter_col
      if (is.null(filter_col) || !filter_col %in% names(VARG26_data)) {
        return(tags$p(class = "text-muted small", "Filter column not found."))
      }
      
      VARG26_filtered <- VARG26_data[VARG26_data[[filter_col]] %in% input$cp_VARG26_groups, ]
      
      # Get unique short citations, excluding NA and empty strings
      short_citations <- unique(VARG26_filtered[["Data citation"]])
      short_citations <- short_citations[!is.na(short_citations) & nchar(trimws(short_citations)) > 0]
      short_citations <- sort(short_citations)
      
      # Get unique full citations if available
      has_full_citations <- "Full citation" %in% names(VARG26_data)
      full_citations_map <- list()
      if (has_full_citations) {
        # Create mapping from short to full citation
        for (short in short_citations) {
          matching_full <- unique(VARG26_filtered[VARG26_filtered[["Data citation"]] == short, "Full citation"])
          matching_full <- matching_full[!is.na(matching_full) & nchar(trimws(matching_full)) > 0]
          if (length(matching_full) > 0) {
            full_citations_map[[short]] <- matching_full[1]
          }
        }
      }
      
      if (length(short_citations) == 0) {
        return(tags$p(class = "text-muted small", "No data citations found for selected groups."))
      }
      
      # Create formatted list with short citations, and full citations as expandable
      tagList(
        tags$p(class = "small text-muted mb-2", 
               paste0(length(short_citations), " unique data source", if(length(short_citations) > 1) "s" else "", ":")),
        tags$ul(
          class = "small list-unstyled",
          style = "margin-bottom: 0;",
          lapply(short_citations, function(cit) {
            full_cit <- full_citations_map[[cit]]
            if (!is.null(full_cit) && full_cit != cit && nchar(full_cit) > nchar(cit) + 10) {
              # Has meaningful full citation - show expandable details
              tags$li(
                style = "padding: 4px 0; border-bottom: 1px solid #eee;",
                tags$details(
                  tags$summary(style = "cursor: pointer; font-weight: 500;", cit),
                  tags$p(style = "margin: 4px 0 0 12px; font-size: 0.9em; color: #555;", full_cit)
                )
              )
            } else {
              # No full citation or same as short - just show short
              tags$li(style = "padding: 2px 0; border-bottom: 1px solid #eee;", cit)
            }
          })
        )
      )
    })

    # Download handlers for Custom Scatter Plot
    output$cp_download_png <- downloadHandler(
      filename = function() {
        paste0("scatter_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        req(rv$cp_plot_obj)
        width_cm <- if (!is.null(input$cp_export_width)) input$cp_export_width else 18
        height_cm <- if (!is.null(input$cp_export_height)) input$cp_export_height else 14
        dpi <- as.numeric(if (!is.null(input$cp_export_dpi)) input$cp_export_dpi else 300)

        ggsave(
          filename = file,
          plot = rv$cp_plot_obj,
          device = "png",
          width = width_cm,
          height = height_cm,
          units = "cm",
          dpi = dpi,
          bg = "white"
        )
      }
    )

    output$cp_download_pdf <- downloadHandler(
      filename = function() {
        paste0("scatter_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        req(rv$cp_plot_obj)
        width_cm <- if (!is.null(input$cp_export_width)) input$cp_export_width else 18
        height_cm <- if (!is.null(input$cp_export_height)) input$cp_export_height else 14
        
        # Check for unicode subscripts/superscripts in axis labels
        unicode_pattern <- "[\u2080-\u2089\u2070\u00B9\u00B2\u00B3\u2074-\u2079\u207A\u207B]"
        x_label <- input$cp_x_label
        y_label <- input$cp_y_label
        has_unicode <- ((!is.null(x_label) && grepl(unicode_pattern, x_label)) ||
                        (!is.null(y_label) && grepl(unicode_pattern, y_label)))
        
        if (has_unicode) {
          showNotification(
            HTML("<b>Note:</b> Your axis labels contain Unicode subscripts/superscripts (e.g., ₂, ²⁺). 
                  These may not render correctly in all PDF viewers. Please check the output carefully. 
                  Consider exporting as SVG or PNG for more reliable rendering, or edit the PDF manually if needed."),
            type = "warning",
            duration = 10
          )
        }

        ggsave(
          filename = file,
          plot = rv$cp_plot_obj,
          device = "pdf",
          width = width_cm,
          height = height_cm,
          units = "cm"
        )
        # Note: Unicode subscripts/superscripts are converted to plotmath expressions
        # via unicode_to_plotmath() function, but rendering quality varies by PDF viewer
      }
    )

    output$cp_download_svg <- downloadHandler(
      filename = function() {
        paste0("scatter_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg")
      },
      content = function(file) {
        req(rv$cp_plot_obj)
        width_cm <- if (!is.null(input$cp_export_width)) input$cp_export_width else 18
        height_cm <- if (!is.null(input$cp_export_height)) input$cp_export_height else 14

        ggsave(
          filename = file,
          plot = rv$cp_plot_obj,
          device = "svg",
          width = width_cm,
          height = height_cm,
          units = "cm"
        )
      }
    )

    output$cp_download_data <- downloadHandler(
      filename = function() {
        paste0("plot_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        write.csv(rv$cp_data, file, row.names = FALSE)
      }
    )


    # =========================================================================
    # TAB 3: STRATIGRAPHIC CORRELATION (OPTIMIZED)
    # =========================================================================

    # --- UI Conditions ---
    output$sc_data_loaded <- reactive({
      !is.null(rv$sc_data)
    })
    output$sc_mapping_validated <- reactive({
      !is.null(rv$sc_column_mappings)
    })
    output$sc_cores_selected <- reactive({
      !is.null(rv$sc_column_mappings) && !is.null(input$sc_reference_core)
    })
    output$sc_has_tiepoints <- reactive({
      nrow(rv$sc_tiepoints) >= 2
    })
    output$sc_warp_applied <- reactive({
      !is.null(rv$sc_warped_data)
    })

    outputOptions(output, "sc_data_loaded", suspendWhenHidden = FALSE)
    outputOptions(output, "sc_mapping_validated", suspendWhenHidden = FALSE)
    outputOptions(output, "sc_cores_selected", suspendWhenHidden = FALSE)
    outputOptions(output, "sc_has_tiepoints", suspendWhenHidden = FALSE)
    outputOptions(output, "sc_warp_applied", suspendWhenHidden = FALSE)

    # --- Sheet Selection UI for Tab 3 (SC) ---
    output$sc_sheet_select <- renderUI({
      req(input$sc_fileUpload)
      file_ext <- tolower(tools::file_ext(input$sc_fileUpload$name))
      if (file_ext %in% c("xlsx", "xls")) {
        sheets <- readxl::excel_sheets(input$sc_fileUpload$datapath)
        selectInput(ns("sc_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
      }
    })

    # --- 1. Upload & Mapping ---
    observeEvent(input$sc_fileUpload, {
      req(input$sc_fileUpload)
      file_ext <- tolower(tools::file_ext(input$sc_fileUpload$name))
      if (file_ext %in% c("xlsx", "xls")) {
        # Wait for sheet selection
        return()
      }
      data <- sc_load_file(
        input$sc_fileUpload$datapath,
        file_name = input$sc_fileUpload$name
      )
      if (is.null(data)) {
        return()
      }
      replace_sc_data(data)

      # Smart column detection
      cols <- names(data)
      nums <- names(data)[sapply(data, is.numeric)]
      detected <- sc_detect_columns(data)

      updateSelectInput(session, "sc_col_site", choices = cols, selected = detected$site_id)
      updateSelectInput(session, "sc_col_sample", choices = cols, selected = detected$sample_name)
      updateSelectInput(session, "sc_col_z", choices = nums, selected = detected$z_var)
      updateSelectInput(session, "sc_col_x", choices = nums, selected = detected$x_var)
    })

    # SC processed-data hint and wiring
    output$sc_processed_hint <- renderUI({
      if (is.null(processed_data)) {
        tagList(span(style = "color:#888; font-size:12px;", "Processed data not available in this session."))
      } else {
        tagList(span(style = "color:#2e7d32; font-size:12px;", "Processed data available. Select 'Use Processed Data' to import."))
      }
    })

    observeEvent(input$sc_data_source, {
      if (input$sc_data_source == "processed") {
        if (!is.null(processed_data) && is.function(processed_data)) {
          pd <- tryCatch(processed_data(), error = function(e) NULL)
          if (!is.null(pd)) {
            replace_sc_data(pd)
            showNotification("SC: using processed data from Processing module", type = "message")
            # Smart column detection
            cols <- names(pd)
            nums <- names(pd)[sapply(pd, is.numeric)]
            detected <- sc_detect_columns(pd)

            updateSelectInput(session, "sc_col_site", choices = cols, selected = detected$site_id)
            updateSelectInput(session, "sc_col_sample", choices = cols, selected = detected$sample_name)
            updateSelectInput(session, "sc_col_z", choices = nums, selected = detected$z_var)
            updateSelectInput(session, "sc_col_x", choices = nums, selected = detected$x_var)
          } else {
            showNotification("Processed data is empty or unavailable", type = "warning")
          }
        } else {
          showNotification("No processed data available in this session", type = "warning")
        }
      }
    })

    if (!is.null(processed_data) && is.function(processed_data)) {
      observeEvent(processed_data(),
        {
          if (!is.null(input$sc_data_source) && input$sc_data_source == "processed") {
            pd <- tryCatch(processed_data(), error = function(e) NULL)
            if (!is.null(pd)) {
              replace_sc_data(pd)
            }
          }
        },
        ignoreNULL = TRUE
      )
    }

    observeEvent(input$sc_sheet, {
      req(input$sc_fileUpload, input$sc_sheet)
      tryCatch(
        {
          data <- readxl::read_excel(input$sc_fileUpload$datapath, sheet = input$sc_sheet)
          if (is.null(data) || nrow(data) == 0) {
            showNotification("Sheet is empty or could not be read", type = "error")
            return()
          }
          replace_sc_data(data)

          # Smart column detection
          cols <- names(data)
          nums <- names(data)[sapply(data, is.numeric)]
          detected <- sc_detect_columns(data)

          updateSelectInput(session, "sc_col_site", choices = cols, selected = detected$site_id)
          updateSelectInput(session, "sc_col_sample", choices = cols, selected = detected$sample_name)
          updateSelectInput(session, "sc_col_z", choices = nums, selected = detected$z_var)
          updateSelectInput(session, "sc_col_x", choices = nums, selected = detected$x_var)
          showNotification("Data loaded successfully", type = "message")
        },
        error = function(e) showNotification(paste("Error loading sheet:", e$message), type = "error")
      )
    })

    # Track previous column mappings for X-variable persistence
    sc_prev_x_var <- reactiveVal(NULL)
    sc_prev_z_var <- reactiveVal(NULL)
    
    observeEvent(input$sc_validate_map, {
      req(rv$sc_data)
      mappings <- list(site_id = input$sc_col_site, sample_name = input$sc_col_sample, z_var = input$sc_col_z, x_var = input$sc_col_x)
      res <- sc_validate_columns(rv$sc_data, mappings)
      if (res$valid) {
        # Check if only X-variable changed (tiepoints should persist)
        prev_x <- sc_prev_x_var()
        prev_z <- sc_prev_z_var()
        x_var_changed <- !is.null(prev_x) && prev_x != input$sc_col_x
        z_var_changed <- !is.null(prev_z) && prev_z != input$sc_col_z
        
        rv$sc_data <- sc_normalize_identifier_columns(rv$sc_data, mappings)
        rv$sc_column_mappings <- mappings
        cores <- unique(rv$sc_data[[mappings$site_id]])
        updateSelectInput(session, "sc_reference_core", choices = cores, selected = cores[1])
        updateSelectInput(session, "sc_target_core", choices = setdiff(cores, cores[1]))
        
        # Handle Z-variable change - must clear tie points as z values are now invalid
        if (z_var_changed && nrow(rv$sc_tiepoints) > 0) {
          showNotification(
            paste("Z-axis variable changed to", input$sc_col_z, "- tie points cleared (z values no longer valid)."), 
            type = "warning", 
            duration = 6
          )
          # Clear tie points since z values are based on old column
          rv$sc_tiepoints <- data.frame(
            id = integer(), ref_sample = character(), ref_z = numeric(),
            target_sample = character(), target_z = numeric(),
            custom_name = character(), use_in_warp = logical(), stringsAsFactors = FALSE
          )
          rv$sc_structure_trigger <- rv$sc_structure_trigger + 1
          rv$sc_warped_data <- NULL
          rv$sc_warp_result <- NULL
          rv$sc_selected_tiepoints <- integer()
        } else if (x_var_changed && nrow(rv$sc_tiepoints) > 0) {
          # X-variable changed but tiepoints preserved - notify user
          showNotification(
            paste("Correlation variable changed to", input$sc_col_x, "- tie points preserved."), 
            type = "message", 
            duration = 4
          )
          # Clear warped data since it needs recalculation with new X variable
          rv$sc_warped_data <- NULL
          rv$sc_warp_result <- NULL
        } else {
          showNotification("Mapped.", type = "message")
        }
        
        # Update tracked variables
        sc_prev_x_var(input$sc_col_x)
        sc_prev_z_var(input$sc_col_z)
      } else {
        showNotification(res$message, type = "error")
      }
    })

    output$sc_umap_hint <- renderUI({
      req(rv$sc_data)
      if (!("UMAP_1D" %in% names(rv$sc_data))) {
        tags$div(
          class = "alert alert-info",
          style = "padding: 8px 10px; margin-top: 8px; margin-bottom: 8px;",
          icon("info-circle"),
          " UMAP_1D is the recommended correlation variable. To create it, run UMAP 1D projection in the Processing module (Step 4)."
        )
      }
    })

    # Track previous core selections to detect actual changes vs re-validation
    sc_prev_ref_core <- reactiveVal(NULL)
    sc_prev_target_core <- reactiveVal(NULL)
    
    # When reference core changes, clear tie points (only if it's a real change)
    observeEvent(input$sc_reference_core, {
      req(rv$sc_data, rv$sc_column_mappings)
      
      current_ref <- input$sc_reference_core
      prev_ref <- sc_prev_ref_core()
      
      cores <- unique(rv$sc_data[[rv$sc_column_mappings$site_id]])
      updateSelectInput(session, "sc_target_core", choices = setdiff(cores, current_ref))

      # Only clear tie points if reference core actually changed (not just re-validated)
      if (!is.null(prev_ref) && prev_ref != current_ref) {
        if (nrow(rv$sc_tiepoints) > 0) {
          showNotification("Tie points cleared due to core selection change.", type = "warning")
          rv$sc_tiepoints <- data.frame(
            id = integer(), ref_sample = character(), ref_z = numeric(),
            target_sample = character(), target_z = numeric(),
            custom_name = character(), use_in_warp = logical(), stringsAsFactors = FALSE
          )
          rv$sc_structure_trigger <- rv$sc_structure_trigger + 1
        }
        rv$sc_warped_data <- NULL
        rv$sc_warp_result <- NULL
        rv$sc_selected_tiepoints <- integer()
        # Reset click mode
        rv$sc_click_mode <- "none"
        rv$sc_pending_ref_sample <- NULL
        rv$sc_pending_ref_z <- NULL
      }
      
      # Update tracked value
      sc_prev_ref_core(current_ref)
    })

    # When target core changes, also clear tie points (only if actual change)
    observeEvent(input$sc_target_core,
      {
        req(rv$sc_data, rv$sc_column_mappings, input$sc_reference_core)
        
        current_target <- input$sc_target_core
        prev_target <- sc_prev_target_core()

        # Only clear tie points if target core actually changed (not just re-validated)
        if (!is.null(prev_target) && prev_target != current_target) {
          if (nrow(rv$sc_tiepoints) > 0) {
            showNotification("Tie points cleared due to core selection change.", type = "warning")
            rv$sc_tiepoints <- data.frame(
              id = integer(), ref_sample = character(), ref_z = numeric(),
              target_sample = character(), target_z = numeric(),
              custom_name = character(), use_in_warp = logical(), stringsAsFactors = FALSE
            )
          }
          rv$sc_structure_trigger <- rv$sc_structure_trigger + 1
          rv$sc_warped_data <- NULL
          rv$sc_warp_result <- NULL
          rv$sc_selected_tiepoints <- integer()
          # Reset click mode
          rv$sc_click_mode <- "none"
          rv$sc_pending_ref_sample <- NULL
          rv$sc_pending_ref_z <- NULL
        }
        
        # Update tracked value
        sc_prev_target_core(current_target)
      },
      ignoreInit = TRUE
    )

    # --- 2. Tie Point Manager (Optimized) ---

    # Add Row
    observeEvent(input$sc_addTiePoint, {
      # Use max of all IDs ever used (current + removed) to avoid reusing IDs
      all_known_ids <- c(rv$sc_tiepoints$id, rv$sc_removed_ids)
      new_id <- if (length(all_known_ids) == 0) 1 else max(all_known_ids) + 1
      rv$sc_tiepoints <- rbind(rv$sc_tiepoints, sc_create_tiepoint(new_id))
      # NEW: Do NOT auto-select - leave it unselected by default
      # Increment trigger to force UI rebuild only on Add/Remove
      rv$sc_structure_trigger <- rv$sc_structure_trigger + 1
    })

    # Remove Row
    observeEvent(input$sc_removeTiePoint, {
      req(length(rv$sc_selected_tiepoints) > 0)
      # Track removed IDs to prevent stale input caching issues
      rv$sc_removed_ids <- c(rv$sc_removed_ids, rv$sc_selected_tiepoints)
      # Remove all selected tie points
      rv$sc_tiepoints <- rv$sc_tiepoints[!rv$sc_tiepoints$id %in% rv$sc_selected_tiepoints, ]
      rv$sc_selected_tiepoints <- integer() # Clear selection
      rv$sc_structure_trigger <- rv$sc_structure_trigger + 1
    })

    # Select All
    observeEvent(input$sc_selectAll, {
      req(nrow(rv$sc_tiepoints) > 0)
      # Select all tie point IDs
      rv$sc_selected_tiepoints <- rv$sc_tiepoints$id
      # Update all checkboxes to checked
      for (id in rv$sc_tiepoints$id) {
        shinyWidgets::updatePrettyCheckbox(session, paste0("tp_sel_", id), value = TRUE)
      }
    })

    # Clear Selection
    observeEvent(input$sc_clearSelection, {
      req(nrow(rv$sc_tiepoints) > 0)
      # Unselect all tie points
      rv$sc_selected_tiepoints <- integer()
      # Update all checkboxes to unchecked
      for (id in rv$sc_tiepoints$id) {
        shinyWidgets::updatePrettyCheckbox(session, paste0("tp_sel_", id), value = FALSE)
      }
    })

    # --- Click-to-Add Tie Point Mode ---

    # Toggle click mode (both original and duplicate button)
    observeEvent(input$sc_clickToAdd, {
      if (rv$sc_click_mode == "none") {
        rv$sc_click_mode <- "ref"
        rv$sc_pending_ref_sample <- NULL
        rv$sc_pending_ref_z <- NULL
      } else {
        rv$sc_click_mode <- "none"
        rv$sc_pending_ref_sample <- NULL
        rv$sc_pending_ref_z <- NULL
      }
    })

    # Duplicate button near plots - same behavior
    observeEvent(input$sc_clickToAdd_top, {
      if (rv$sc_click_mode == "none") {
        rv$sc_click_mode <- "ref"
        rv$sc_pending_ref_sample <- NULL
        rv$sc_pending_ref_z <- NULL
      } else {
        rv$sc_click_mode <- "none"
        rv$sc_pending_ref_sample <- NULL
        rv$sc_pending_ref_z <- NULL
      }
    })



    # Click mode status indicator (sidebar)
    output$sc_click_mode_status <- renderUI({
      if (rv$sc_click_mode == "none") {
        return(NULL)
      } else if (rv$sc_click_mode == "ref") {
        div(
          class = "alert alert-info py-1 px-2 mb-2 small",
          icon("crosshairs"), " Click a ", tags$strong("Reference"), " point on the plot..."
        )
      } else if (rv$sc_click_mode == "target") {
        div(
          class = "alert alert-warning py-1 px-2 mb-2 small",
          icon("crosshairs"), " Ref: ", tags$strong(rv$sc_pending_ref_sample),
          ". Now click a ", tags$strong("Target"), " point..."
        )
      }
    })

    # Click mode status indicator (above plots)
    output$sc_click_mode_status_top <- renderUI({
      if (rv$sc_click_mode == "none") {
        return(NULL)
      } else if (rv$sc_click_mode == "ref") {
        div(
          class = "alert alert-info py-1 px-2 mb-1 small",
          icon("crosshairs"), " Click a ", tags$strong("Reference"), " point on the plot..."
        )
      } else if (rv$sc_click_mode == "target") {
        div(
          class = "alert alert-warning py-1 px-2 mb-1 small",
          icon("crosshairs"), " Ref: ", tags$strong(rv$sc_pending_ref_sample),
          ". Now click a ", tags$strong("Target"), " point..."
        )
      }
    })

    # --- Helper function for handling plot clicks ---
    handle_plot_click <- function(click) {
      req(rv$sc_click_mode != "none")
      req(click)

      resolved_click <- tryCatch(
        sc_resolve_plot_click(
          click = click,
          data = rv$sc_data,
          mappings = rv$sc_column_mappings,
          ref_core = input$sc_reference_core,
          target_core = input$sc_target_core,
          show_points = isTRUE(input$sc_show_points)
        ),
        error = function(e) {
          showNotification(e$message, type = "warning")
          NULL
        }
      )
      if (is.null(resolved_click)) {
        return(invisible(NULL))
      }

      sample_name <- resolved_click$sample_name
      y_value <- resolved_click$z
      is_ref <- resolved_click$is_ref

      if (rv$sc_click_mode == "ref") {
        if (is_ref) {
          # Store reference selection
          rv$sc_pending_ref_sample <- sample_name
          rv$sc_pending_ref_z <- y_value
          rv$sc_click_mode <- "target"
          showNotification(
            paste("Reference selected:", sample_name, "- now click a Target point."),
            type = "message", duration = 3
          )
        } else {
          showNotification("Please click a Reference (blue) point first.", type = "warning")
        }
      } else if (rv$sc_click_mode == "target") {
        if (!is_ref) {
          # Create the tie point - use max of all IDs ever used (current + removed) to avoid ID reuse
          all_known_ids <- c(rv$sc_tiepoints$id, rv$sc_removed_ids)
          new_id <- if (length(all_known_ids) == 0) 1L else max(all_known_ids) + 1L
          new_row <- data.frame(
            id = new_id,
            ref_sample = rv$sc_pending_ref_sample,
            ref_z = rv$sc_pending_ref_z,
            target_sample = sample_name,
            target_z = y_value,
            custom_name = NA_character_,
            use_in_warp = TRUE,
            stringsAsFactors = FALSE
          )
          rv$sc_tiepoints <- rbind(rv$sc_tiepoints, new_row)
          rv$sc_structure_trigger <- rv$sc_structure_trigger + 1

          # Reset click mode
          rv$sc_click_mode <- "none"
          rv$sc_pending_ref_sample <- NULL
          rv$sc_pending_ref_z <- NULL

          showNotification(paste("Tie point added:", new_row$ref_sample, "↔", new_row$target_sample), type = "message")
        } else {
          showNotification("Please click a Target (orange) point.", type = "warning")
        }
      }
    }

    # Handle plot clicks from Initial Alignment plot
    observeEvent(event_data("plotly_click", source = "sc_initial", priority = "event"),
      {
        click <- event_data("plotly_click", source = "sc_initial", priority = "event")
        handle_plot_click(click)
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # Affine-preview clicks resolve through sample keys to original uploaded depths.
    observeEvent(event_data("plotly_click", source = "sc_affine", priority = "event"),
      {
        click <- event_data("plotly_click", source = "sc_affine", priority = "event")
        handle_plot_click(click)
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # Handle plot clicks from Warped Alignment plot
    observeEvent(event_data("plotly_click", source = "sc_warped", priority = "event"),
      {
        click <- event_data("plotly_click", source = "sc_warped", priority = "event")
        handle_plot_click(click)
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE
    )

    # Render List
    output$sc_tiepoint_list_ui <- renderUI({
      req(rv$sc_data, rv$sc_column_mappings)
      # Depend on BOTH structure trigger (for Add/Remove) AND selection changes (for highlighting)
      trigger <- rv$sc_structure_trigger
      selection_state <- rv$sc_selected_tiepoints # This makes the UI reactive to selection changes

      # Use isolate to grab current data for initial population
      current_tp <- isolate(rv$sc_tiepoints)

      if (nrow(current_tp) == 0) {
        return(p("No tie points yet.", style = "padding:10px; color:#999;"))
      }

      ref_samps <- rv$sc_data %>%
        filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_reference_core) %>%
        pull(!!sym(rv$sc_column_mappings$sample_name))
      tgt_samps <- rv$sc_data %>%
        filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_target_core) %>%
        pull(!!sym(rv$sc_column_mappings$sample_name))

      lapply(seq_len(nrow(current_tp)), function(i) {
        tp <- current_tp[i, ]
        id <- tp$id

        # Determine initial selection state (check if id is in the vector of selected ids)
        # Use selection_state (reactive value we created above) to check current state
        is_sel <- id %in% selection_state
        row_class <- if (is_sel) "tie-point-row selected-row" else "tie-point-row"

        div(
          class = row_class,
          # Col 1: Select (Pretty Checkbox)
          div(
            class = "tie-point-col", style = "width:10%; padding-left:10px;",
            shinyWidgets::prettyCheckbox(ns(paste0("tp_sel_", id)), paste0("TP ", id),
              value = is_sel,
              status = "primary", icon = icon("check"), shape = "curve"
            )
          ),
          # Col 2: Ref Dropdown (Selectize with body parent)
          div(
            class = "tie-point-col", style = "width:25%;",
            selectizeInput(ns(paste0("tp_ref_", id)), NULL,
              choices = c("", unique(ref_samps)), selected = tp$ref_sample, width = "100%",
              options = list(dropdownParent = "body", placeholder = "Select Ref...")
            )
          ),
          # Col 3: Target Dropdown (Selectize with body parent)
          div(
            class = "tie-point-col", style = "width:25%;",
            selectizeInput(ns(paste0("tp_tgt_", id)), NULL,
              choices = c("", unique(tgt_samps)), selected = tp$target_sample, width = "100%",
              options = list(dropdownParent = "body", placeholder = "Select Target...")
            )
          ),
          # Col 4: Custom Name (Text Input)
          div(
            class = "tie-point-col", style = "width:25%;",
            textInput(ns(paste0("tp_name_", id)), NULL,
              value = if (is.null(tp$custom_name) || length(tp$custom_name) == 0 || is.na(tp$custom_name)) "" else tp$custom_name,
              placeholder = paste0("TP_", id)
            )
          ),
          # Col 5: Use (Pretty Checkbox)
          div(
            class = "tie-point-col", style = "width:15%; text-align:center;",
            shinyWidgets::prettyCheckbox(ns(paste0("tp_use_", id)), "Use",
              value = tp$use_in_warp,
              status = "success", icon = icon("check"), shape = "curve"
            )
          )
        )
      })
    })

    # Data Sync Logic
    observe({
      req(nrow(rv$sc_tiepoints) > 0, rv$sc_column_mappings, input$sc_reference_core, input$sc_target_core)

      # Get valid sample names for current cores
      valid_ref_samples <- rv$sc_data %>%
        filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_reference_core) %>%
        pull(!!sym(rv$sc_column_mappings$sample_name)) %>%
        unique()
      valid_target_samples <- rv$sc_data %>%
        filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_target_core) %>%
        pull(!!sym(rv$sc_column_mappings$sample_name)) %>%
        unique()

      # Iterate over IDs currently in the dataframe
      for (id in rv$sc_tiepoints$id) {
        # 1. Handle Selection (Highlight) - Multiple Selection Support
        sel_id <- paste0("tp_sel_", id)
        if (!is.null(input[[sel_id]])) {
          if (input[[sel_id]]) {
            # Add to selection vector if not already there
            if (!id %in% rv$sc_selected_tiepoints) {
              rv$sc_selected_tiepoints <- c(rv$sc_selected_tiepoints, id)
            }
          } else {
            # Remove from selection vector
            rv$sc_selected_tiepoints <- rv$sc_selected_tiepoints[rv$sc_selected_tiepoints != id]
          }
        }

        # 2. Update Data Values (Without triggering UI redraw)
        idx <- which(rv$sc_tiepoints$id == id)
        if (length(idx) == 0) next # Skip if tie point was removed

        ref_id <- paste0("tp_ref_", id)
        tgt_id <- paste0("tp_tgt_", id)
        use_id <- paste0("tp_use_", id)

        # Only sync ref dropdown if input exists AND value is valid for current reference core
        if (!is.null(input[[ref_id]])) {
          curr <- input[[ref_id]]
          # Only update if curr is a valid sample for the current reference core (or empty)
          if (curr == "" || curr %in% valid_ref_samples) {
            if (!is.na(rv$sc_tiepoints$ref_sample[idx]) && rv$sc_tiepoints$ref_sample[idx] != curr) {
              rv$sc_tiepoints$ref_sample[idx] <- curr
              if (curr != "") {
                z <- rv$sc_data %>%
                  filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_reference_core, .data[[rv$sc_column_mappings$sample_name]] == curr) %>%
                  pull(!!sym(rv$sc_column_mappings$z_var)) %>%
                  head(1)
                if (length(z) > 0) rv$sc_tiepoints$ref_z[idx] <- z
              } else {
                rv$sc_tiepoints$ref_z[idx] <- NA_real_
              }
            } else if (is.na(rv$sc_tiepoints$ref_sample[idx]) && curr != "") {
              # Handle initial empty case
              rv$sc_tiepoints$ref_sample[idx] <- curr
              z <- rv$sc_data %>%
                filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_reference_core, .data[[rv$sc_column_mappings$sample_name]] == curr) %>%
                pull(!!sym(rv$sc_column_mappings$z_var)) %>%
                head(1)
              if (length(z) > 0) rv$sc_tiepoints$ref_z[idx] <- z
            }
          }
          # If curr is NOT valid for current core, don't sync - keep the stored value
        }

        # Only sync target dropdown if input exists AND value is valid for current target core
        if (!is.null(input[[tgt_id]])) {
          curr <- input[[tgt_id]]
          # Only update if curr is a valid sample for the current target core (or empty)
          if (curr == "" || curr %in% valid_target_samples) {
            if (!is.na(rv$sc_tiepoints$target_sample[idx]) && rv$sc_tiepoints$target_sample[idx] != curr) {
              rv$sc_tiepoints$target_sample[idx] <- curr
              if (curr != "") {
                z <- rv$sc_data %>%
                  filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_target_core, .data[[rv$sc_column_mappings$sample_name]] == curr) %>%
                  pull(!!sym(rv$sc_column_mappings$z_var)) %>%
                  head(1)
                if (length(z) > 0) rv$sc_tiepoints$target_z[idx] <- z
              } else {
                rv$sc_tiepoints$target_z[idx] <- NA_real_
              }
            } else if (is.na(rv$sc_tiepoints$target_sample[idx]) && curr != "") {
              rv$sc_tiepoints$target_sample[idx] <- curr
              z <- rv$sc_data %>%
                filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_target_core, .data[[rv$sc_column_mappings$sample_name]] == curr) %>%
                pull(!!sym(rv$sc_column_mappings$z_var)) %>%
                head(1)
              if (length(z) > 0) rv$sc_tiepoints$target_z[idx] <- z
            }
          }
          # If curr is NOT valid for current core, don't sync - keep the stored value
        }

        if (!is.null(input[[use_id]]) && rv$sc_tiepoints$use_in_warp[idx] != input[[use_id]]) {
          rv$sc_tiepoints$use_in_warp[idx] <- input[[use_id]]
        }
        
        # Sync custom name
        name_id <- paste0("tp_name_", id)
        if (!is.null(input[[name_id]])) {
          curr_name <- input[[name_id]]
          stored_name <- rv$sc_tiepoints$custom_name[idx]
          # Compare with NA handling
          if ((is.na(stored_name) && curr_name != "") || (!is.na(stored_name) && stored_name != curr_name)) {
            rv$sc_tiepoints$custom_name[idx] <- if (curr_name == "") NA_character_ else curr_name
          }
        }
      }
    })

    # --- 3. Warping ---
    observeEvent(input$sc_applyWarp, {
      req(rv$sc_tiepoints)
      val <- sc_validate_tiepoints(
        rv$sc_tiepoints,
        ref_direction = input$sc_reference_direction,
        target_direction = input$sc_target_direction
      )
      if (!val$valid) {
        # Use duration = NULL so crossing/duplicate errors stay until user dismisses them
        showNotification(val$message, type = "error", duration = NULL)
        return()
      }

      tryCatch(
        {
          res <- sc_calculate_warp(
            rv$sc_tiepoints,
            method = "auto",
            ref_direction = input$sc_reference_direction,
            target_direction = input$sc_target_direction,
            extrapolation = "linear"
          )
          rv$sc_warp_result <- res
          target_df <- rv$sc_data %>% filter(.data[[rv$sc_column_mappings$site_id]] == input$sc_target_core)
          rv$sc_warped_data <- sc_apply_warp(target_df, rv$sc_column_mappings$z_var, res$warp_func)

          rv$sc_output_object <- sc_create_output(rv$sc_tiepoints, res, list(
            reference_core = input$sc_reference_core, target_core = input$sc_target_core,
            x_var = rv$sc_column_mappings$x_var, z_var = rv$sc_column_mappings$z_var,
            reference_direction = input$sc_reference_direction,
            target_direction = input$sc_target_direction
          ))
          showNotification("Warping applied.", type = "message")
        },
        error = function(e) showNotification(e$message, type = "error")
      )
    })

    # --- 4. Plots & Export ---
    output$sc_initial_plot_ui <- renderUI({
      plotlyOutput(
        ns("sc_initialPlot"),
        height = sc_plot_height_css(input$sc_plot_height)
      )
    })

    affine_parameters <- reactive({
      req(rv$sc_data, rv$sc_column_mappings, input$sc_reference_core, input$sc_target_core)
      sc_affine_depth_parameters(
        rv$sc_data,
        rv$sc_column_mappings,
        input$sc_reference_core,
        input$sc_target_core,
        mode = if (!is.null(input$sc_affine_mode)) input$sc_affine_mode else "auto",
        manual_scale = if (!is.null(input$sc_affine_scale)) input$sc_affine_scale else 1,
        manual_shift = if (!is.null(input$sc_affine_shift)) input$sc_affine_shift else 0,
        ref_direction = if (!is.null(input$sc_reference_direction)) input$sc_reference_direction else "down",
        target_direction = if (!is.null(input$sc_target_direction)) input$sc_target_direction else "down"
      )
    })

    observeEvent(
      list(input$sc_reference_direction, input$sc_target_direction),
      {
        rv$sc_warped_data <- NULL
        rv$sc_warp_result <- NULL
        rv$sc_output_object <- NULL
      },
      ignoreInit = TRUE
    )

    observeEvent(input$sc_affine_reset, {
      shinyWidgets::updateRadioGroupButtons(session, "sc_affine_mode", selected = "auto")
      updateNumericInput(session, "sc_affine_scale", value = 1)
      updateNumericInput(session, "sc_affine_shift", value = 0)
    })

    output$sc_affine_status <- renderUI({
      pars <- affine_parameters()
      div(
        class = "small text-muted pb-2",
        tags$strong("Display only."),
        sprintf(
          " Target preview coordinate = %.4g x uploaded coordinate %s %.4g. Tie points retain uploaded coordinates.",
          pars$scale,
          if (pars$shift < 0) "-" else "+",
          abs(pars$shift)
        )
      )
    })

    output$sc_affine_plot_ui <- renderUI({
      plotlyOutput(
        ns("sc_affinePlot"),
        height = sc_plot_height_css(input$sc_plot_height)
      )
    })

    output$sc_warped_plot_ui <- renderUI({
      plotlyOutput(
        ns("sc_warpedPlot"),
        height = sc_plot_height_css(input$sc_plot_height)
      )
    })

    output$sc_initialPlot <- renderPlotly({
      req(rv$sc_data, rv$sc_column_mappings, input$sc_reference_core, input$sc_target_core)
      # Pass highlight_sample when in "target" mode (ref already selected)
      highlight <- if (rv$sc_click_mode == "target") rv$sc_pending_ref_sample else NULL
      tryCatch(
        {
          sc_plot_initial(
            rv$sc_data, rv$sc_column_mappings,
            input$sc_reference_core, input$sc_target_core,
            rv$sc_tiepoints, rv$sc_selected_tiepoints,
            show_points = isTRUE(input$sc_show_points),
            plot_source = "sc_initial",
            highlight_sample = highlight,
            target_offset = if (!is.null(input$sc_target_offset)) input$sc_target_offset else 0,
            line_width = if (!is.null(input$sc_line_width)) input$sc_line_width else 2,
            point_size = if (!is.null(input$sc_point_size)) input$sc_point_size else 8,
            jitter_x = if (!is.null(input$sc_jitter_x)) input$sc_jitter_x else 0,
            jitter_z = if (!is.null(input$sc_jitter_z)) input$sc_jitter_z else 0,
            reference_direction = if (!is.null(input$sc_reference_direction)) input$sc_reference_direction else "down"
          )
        },
        error = function(e) {
          # Return empty plot with error message
          plot_ly(source = "sc_initial") %>%
            layout(
              title = "Error rendering plot",
              annotations = list(
                list(
                  x = 0.5, y = 0.5, text = paste("Error:", e$message),
                  showarrow = FALSE, xref = "paper", yref = "paper"
                )
              )
            ) %>%
            event_register("plotly_click")
        }
      )
    })
    output$sc_affinePlot <- renderPlotly({
      req(rv$sc_data, rv$sc_column_mappings, input$sc_reference_core, input$sc_target_core)
      highlight <- if (rv$sc_click_mode == "target") rv$sc_pending_ref_sample else NULL
      tryCatch(
        {
          sc_plot_affine_preview(
            rv$sc_data, rv$sc_column_mappings,
            input$sc_reference_core, input$sc_target_core,
            affine_parameters = affine_parameters(),
            tiepoints = rv$sc_tiepoints,
            selected_ids = rv$sc_selected_tiepoints,
            show_points = isTRUE(input$sc_show_points),
            plot_source = "sc_affine",
            highlight_sample = highlight,
            target_offset = if (!is.null(input$sc_target_offset)) input$sc_target_offset else 0,
            line_width = if (!is.null(input$sc_line_width)) input$sc_line_width else 2,
            point_size = if (!is.null(input$sc_point_size)) input$sc_point_size else 8,
            jitter_x = if (!is.null(input$sc_jitter_x)) input$sc_jitter_x else 0,
            jitter_z = if (!is.null(input$sc_jitter_z)) input$sc_jitter_z else 0,
            reference_direction = if (!is.null(input$sc_reference_direction)) input$sc_reference_direction else "down"
          )
        },
        error = function(e) {
          plot_ly(source = "sc_affine") %>%
            layout(
              title = "Error rendering affine preview",
              annotations = list(
                list(
                  x = 0.5, y = 0.5, text = paste("Error:", e$message),
                  showarrow = FALSE, xref = "paper", yref = "paper"
                )
              )
            ) %>%
            event_register("plotly_click")
        }
      )
    })
    output$sc_warpedPlot <- renderPlotly({
      req(rv$sc_warped_data)
      # Pass highlight_sample when in "target" mode (ref already selected)
      highlight <- if (rv$sc_click_mode == "target") rv$sc_pending_ref_sample else NULL
      tryCatch(
        {
          sc_plot_warped(
            rv$sc_data, rv$sc_warped_data, rv$sc_column_mappings,
            input$sc_reference_core, input$sc_target_core,
            rv$sc_tiepoints, rv$sc_selected_tiepoints,
            show_points = isTRUE(input$sc_show_points),
            plot_source = "sc_warped",
            highlight_sample = highlight,
            target_offset = if (!is.null(input$sc_target_offset)) input$sc_target_offset else 0,
            line_width = if (!is.null(input$sc_line_width)) input$sc_line_width else 2,
            point_size = if (!is.null(input$sc_point_size)) input$sc_point_size else 8,
            jitter_x = if (!is.null(input$sc_jitter_x)) input$sc_jitter_x else 0,
            jitter_z = if (!is.null(input$sc_jitter_z)) input$sc_jitter_z else 0,
            reference_direction = if (!is.null(input$sc_reference_direction)) input$sc_reference_direction else "down"
          )
        },
        error = function(e) {
          plot_ly(source = "sc_warped") %>%
            layout(
              title = "Error rendering plot",
              annotations = list(
                list(
                  x = 0.5, y = 0.5, text = paste("Error:", e$message),
                  showarrow = FALSE, xref = "paper", yref = "paper"
                )
              )
            ) %>%
            event_register("plotly_click")
        }
      )
    })
    output$sc_warpFitPlot <- renderPlotly({
      req(rv$sc_warp_result)
      sc_plot_fit(rv$sc_warp_result)
    })
    output$sc_rmseOutput <- renderText({
      req(rv$sc_warp_result)
      paste("Anchor RMSE (verification only):", round(rv$sc_warp_result$diagnostics$rmse, 4))
    })
    
    # Exact-anchor verification plot
    output$sc_residualPlot <- renderPlotly({
      req(rv$sc_warp_result, rv$sc_warp_result$diagnostics$residuals)
      resid_df <- rv$sc_warp_result$diagnostics$residuals
      
      # Exact-knot residuals should be zero apart from numerical precision.
      plot_ly(
        data = resid_df,
        x = ~id,
        y = ~residual,
        type = "bar",
        marker = list(color = "#3498db"),
        text = ~paste0("TP", id, ": ", ref_sample, " ↔ ", target_sample, 
                       "<br>Residual: ", round(residual, 3)),
        hoverinfo = "text"
      ) %>%
        layout(
          title = "Exact-anchor Verification",
          xaxis = list(title = "Tie Point ID", dtick = 1),
          yaxis = list(title = "Reference Z - mapped anchor Z"),
          annotations = list(
            list(
              x = 0.5, y = 1.08, xref = "paper", yref = "paper",
              text = "Non-zero values indicate numerical or implementation error, not a weak geological tie.",
              showarrow = FALSE, font = list(size = 11, color = "#5f6368")
            )
          )
        )
    })
    
    # Exact-anchor verification table
    output$sc_residualTable <- renderTable({
      req(rv$sc_warp_result, rv$sc_warp_result$diagnostics$residuals)
      resid_df <- rv$sc_warp_result$diagnostics$residuals
      
      # Create display-friendly table
      display_df <- data.frame(
        `Tie Point` = paste0("TP", resid_df$id),
        `Ref Sample` = resid_df$ref_sample,
        `Target Sample` = resid_df$target_sample,
        `Ref Z` = round(resid_df$ref_z, 3),
        `Target Z` = round(resid_df$target_z, 3),
        `Predicted Z` = round(resid_df$pred_z, 3),
        `Anchor Residual` = signif(resid_df$residual, 6),
        check.names = FALSE
      )
      display_df
    }, striped = TRUE, hover = TRUE, bordered = TRUE)
    
    output$sc_downloadWarped <- downloadHandler(filename = function() {
      "warped_data.csv"
    }, content = function(file) {
      write.csv(rv$sc_warped_data, file, row.names = FALSE)
    })
    output$sc_downloadTiePoints <- downloadHandler(filename = function() {
      "tiepoints.csv"
    }, content = function(file) {
      write.csv(rv$sc_tiepoints, file, row.names = FALSE)
    })
    observeEvent(input$sc_copyTiePoints, {
      req(rv$sc_output_object)
      if (requireNamespace("clipr", quietly = TRUE)) {
        clipr::write_clip(sc_to_json(rv$sc_output_object))
        showNotification("Copied!", type = "message")
      } else {
        showNotification("Install clipr", type = "warning")
      }
    })

    # Return processed data and tie points for other modules
    return(list(
      vp_data = reactive(rv$vp_data),
      tiepoints = reactive(rv$sc_tiepoints),
      # Population/symbology styles for project save/load
      pop_styles = reactive(list(
        cp_group_color_overrides = rv$cp_group_color_overrides,
        cp_group_shape_overrides = rv$cp_group_shape_overrides,
        cp_current_palette = rv$cp_current_palette
      )),
      set_pop_styles = function(styles) {
        if (!is.null(styles)) {
          rv$cp_group_color_overrides <- styles$cp_group_color_overrides
          rv$cp_group_shape_overrides <- styles$cp_group_shape_overrides
          if (!is.null(styles$cp_current_palette)) {
            rv$cp_current_palette <- styles$cp_current_palette
          }
        }
      }
    ))
  })
}
