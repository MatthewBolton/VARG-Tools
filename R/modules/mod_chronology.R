chronology_as_finite_numeric <- function(x) {
    text <- trimws(as.character(x))
    blank <- is.na(x) | !nzchar(text)
    value <- suppressWarnings(as.numeric(text))
    valid <- !blank & is.finite(value)
    list(value = value, blank = blank, valid = valid)
}

chronology_validate_table <- function(df, model = c("age_depth", "phase"), allow_duplicates = FALSE) {
    model <- match.arg(model)
    errors <- character(0)
    if (!is.data.frame(df) || nrow(df) == 0) {
        return(list(ok = FALSE, data = df, errors = "No chronology rows were loaded."))
    }

    fields <- if (identical(model, "age_depth")) {
        list(name = "Name", age = "Age", uncertainty = "Uncertainty", type = "Type", depth = "Depth")
    } else {
        list(
            name = "Name", age = "Age.Mean", uncertainty = "Age.SD",
            type = "Date.Type", tephra = "Tephra.Name", phase = "Phase"
        )
    }
    required_columns <- unname(unlist(fields, use.names = FALSE))
    missing_columns <- setdiff(required_columns, names(df))
    if (length(missing_columns) > 0) {
        return(list(
            ok = FALSE,
            data = df,
            errors = paste("Missing required columns:", paste(missing_columns, collapse = ", "))
        ))
    }

    row_error <- function(mask, message) {
        rows <- which(mask)
        if (length(rows) > 0) {
            errors <<- c(errors, sprintf("Row(s) %s: %s", paste(rows, collapse = ", "), message))
        }
    }

    names_clean <- trimws(as.character(df[[fields$name]]))
    row_error(is.na(df[[fields$name]]) | !nzchar(names_clean), "Name is required.")
    if (!isTRUE(allow_duplicates)) {
        duplicate_names <- !is.na(names_clean) & nzchar(names_clean) &
            (duplicated(names_clean) | duplicated(names_clean, fromLast = TRUE))
        row_error(duplicate_names, "Name must be unique, or enable automatic duplicate renaming.")
    }
    df[[fields$name]] <- names_clean

    allowed_types <- c(
        "R_Date", "R_F14C", "Date", "Tephra", "Date_calBP", "Tephra_calBP",
        "Date_CE", "Tephra_CE", "CE"
    )
    if (identical(model, "age_depth")) allowed_types <- c(allowed_types, "Boundary")
    types <- trimws(as.character(df[[fields$type]]))
    row_error(is.na(df[[fields$type]]) | !nzchar(types), "Date type is required.")
    row_error(!is.na(types) & nzchar(types) & !types %in% allowed_types, "Date type is not supported.")
    df[[fields$type]] <- types

    if (identical(model, "age_depth")) {
        depth <- chronology_as_finite_numeric(df[[fields$depth]])
        row_error(!depth$valid, "Depth must be a finite number.")
        df[[fields$depth]] <- depth$value
    } else {
        tephra_names <- trimws(as.character(df[[fields$tephra]]))
        phases <- trimws(as.character(df[[fields$phase]]))
        row_error(is.na(df[[fields$tephra]]) | !nzchar(tephra_names), "Tephra name is required.")
        row_error(is.na(df[[fields$phase]]) | !nzchar(phases), "Phase is required.")
        row_error(!is.na(phases) & nzchar(phases) & !phases %in% c("Before", "During", "After"),
                  "Phase must be Before, During, or After.")
        df[[fields$tephra]] <- tephra_names
        df[[fields$phase]] <- phases
    }

    age <- chronology_as_finite_numeric(df[[fields$age]])
    uncertainty <- chronology_as_finite_numeric(df[[fields$uncertainty]])
    required_age_types <- c(
        "R_Date", "R_F14C", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE", "CE"
    )
    requires_age <- types %in% required_age_types
    generic_event <- types %in% c("Date", "Tephra")
    boundary <- types == "Boundary"

    row_error(requires_age & !age$valid, "Age must be a finite number for this date type.")
    row_error(requires_age & !uncertainty$valid, "Uncertainty must be a finite number for this date type.")
    row_error(generic_event & xor(age$blank, uncertainty$blank),
              "Generic Date/Tephra age and uncertainty must either both be blank or both be finite.")
    row_error(generic_event & !age$blank & !age$valid, "Age must be finite when supplied.")
    row_error(generic_event & !uncertainty$blank & !uncertainty$valid,
              "Uncertainty must be finite when supplied.")
    row_error(boundary & (!age$blank | !uncertainty$blank),
              "Boundary rows must leave age and uncertainty blank.")
    row_error(!uncertainty$blank & uncertainty$valid & uncertainty$value < 0,
              "Uncertainty cannot be negative.")
    row_error(types == "R_F14C" & age$valid & age$value <= 0, "R_F14C age must be positive.")
    df[[fields$age]] <- age$value
    df[[fields$uncertainty]] <- uncertainty$value

    unc_field <- intersect(c("Unc_Type", "Unc.Type"), names(df))
    if (length(unc_field) > 0) {
        unc_values <- toupper(trimws(as.character(df[[unc_field[[1]]]])))
        unc_blank <- is.na(df[[unc_field[[1]]]]) | !nzchar(unc_values)
        row_error(!unc_blank & !unc_values %in% c("N", "NORMAL", "U", "UNIFORM"),
                  "Uncertainty type must be N/Normal or U/Uniform.")
    }

    outlier_field <- intersect(c("Outlier.Type", "Outlier.type", "Outlier_Type", "Outlier_type"), names(df))
    if (length(outlier_field) > 0) {
        outlier_values <- trimws(as.character(df[[outlier_field[[1]]]]))
        outlier_blank <- is.na(df[[outlier_field[[1]]]]) | !nzchar(outlier_values)
        valid_outliers <- c("General", "Charcoal", "T_Outlier", "SSimple", "RSimple", "TSimple", "RScaled")
        row_error(!outlier_blank & !tolower(outlier_values) %in% tolower(valid_outliers),
                  "Outlier type is not supported.")
    }

    probability_field <- intersect(
        c("Outlier.Prob", "Outlier.prob.", "Outlier_Prob", "Outlier_prob", "Manual.outlier.prob."),
        names(df)
    )
    if (length(probability_field) > 0) {
        probability <- chronology_as_finite_numeric(df[[probability_field[[1]]]])
        row_error(!probability$blank & !probability$valid, "Outlier probability must be finite.")
        row_error(probability$valid & (probability$value < 0 | probability$value > 1),
                  "Outlier probability must be between 0 and 1.")
        df[[probability_field[[1]]]] <- probability$value
    }

    list(ok = length(errors) == 0, data = df, errors = unique(errors))
}

chronology_validate_interp_rate <- function(value) {
    parsed <- chronology_as_finite_numeric(value)
    if (length(parsed$valid) != 1 || !parsed$valid || parsed$value <= 0) {
        return(list(
            ok = FALSE,
            value = NULL,
            error = "Interpolation rate must be a finite number greater than zero."
        ))
    }
    list(ok = TRUE, value = parsed$value, error = NULL)
}

chronology_normalize_outlier_type <- function(value) {
    if (length(value) != 1 || is.na(value)) return(NA_character_)
    value <- trimws(as.character(value))
    if (!nzchar(value)) return(NA_character_)

    valid <- c("General", "Charcoal", "T_Outlier", "SSimple", "RSimple", "TSimple", "RScaled")
    match_index <- match(tolower(value), tolower(valid))
    normalized <- if (is.na(match_index)) "General" else valid[[match_index]]
    if (identical(normalized, "T_Outlier")) "General" else normalized
}

chronology_outlier_model_types <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(character(0))

    outlier_column <- intersect(
        c("Outlier.Type", "Outlier.type", "Outlier_Type", "Outlier_type"),
        names(df)
    )
    outlier_values <- if (length(outlier_column) > 0) {
        df[[outlier_column[[1]]]]
    } else {
        rep(NA_character_, nrow(df))
    }
    normalized <- vapply(
        outlier_values,
        chronology_normalize_outlier_type,
        character(1),
        USE.NAMES = FALSE
    )
    model_types <- normalized[!is.na(normalized)]

    type_column <- intersect(c("Type", "Date.Type"), names(df))
    if (length(type_column) > 0) {
        date_types <- trimws(as.character(df[[type_column[[1]]]]))
        uses_default_general <- date_types %in% c("R_Date", "R_F14C", "Date", "CE") &
            is.na(normalized)
        if (any(uses_default_general)) model_types <- c(model_types, "General")
    }

    unique(model_types)
}

chronology_select_curve <- function(type, age, primary_curve, use_bomb, bomb_curve_name = "Bomb21NH1") {
    standard_curves <- c("IntCal20", "Marine20", "SHCal20")
    bomb_curves <- c("Bomb21NH1", "Bomb21NH2", "Bomb21NH3", "Bomb21SH12", "Bomb21SH3")
    if (length(primary_curve) != 1 || is.na(primary_curve) || !primary_curve %in% standard_curves) {
        stop("Primary calibration curve is not supported.")
    }
    if (identical(primary_curve, "Marine20") || !isTRUE(use_bomb)) return(primary_curve)
    if (length(bomb_curve_name) != 1 || is.na(bomb_curve_name) || !bomb_curve_name %in% bomb_curves) {
        stop("Bomb calibration curve is not supported.")
    }

    type <- as.character(type[[1]])
    age <- suppressWarnings(as.numeric(age[[1]]))
    if (!is.finite(age)) return(primary_curve)
    if (identical(type, "R_Date") && age < 1000) return(bomb_curve_name)
    if (identical(type, "R_F14C") && age > 1) return(bomb_curve_name)
    primary_curve
}

chronology_validate_delta_r <- function(curve, value, uncertainty) {
    if (!identical(curve, "Marine20")) {
        return(list(ok = TRUE, value = NULL, uncertainty = NULL, error = NULL))
    }
    parsed_value <- chronology_as_finite_numeric(value)
    parsed_uncertainty <- chronology_as_finite_numeric(uncertainty)
    if (length(parsed_value$valid) != 1 || !parsed_value$valid ||
        length(parsed_uncertainty$valid) != 1 || !parsed_uncertainty$valid) {
        return(list(ok = FALSE, value = NULL, uncertainty = NULL,
                    error = "Marine20 requires finite Delta_R value and uncertainty inputs."))
    }
    if (parsed_uncertainty$value < 0) {
        return(list(ok = FALSE, value = NULL, uncertainty = NULL,
                    error = "Marine20 Delta_R uncertainty cannot be negative."))
    }
    list(ok = TRUE, value = parsed_value$value, uncertainty = parsed_uncertainty$value, error = NULL)
}

chronology_format_delta_r <- function(curve, value, uncertainty, indent = "   ") {
    checked <- chronology_validate_delta_r(curve, value, uncertainty)
    if (!checked$ok) stop(checked$error)
    if (!identical(curve, "Marine20")) return(character(0))
    sprintf(
        '%sDelta_R("LocalMarine",%s,%s);',
        indent,
        format(checked$value, scientific = FALSE, trim = TRUE, digits = 15),
        format(checked$uncertainty, scientific = FALSE, trim = TRUE, digits = 15)
    )
}

chronology_validate_oxcal_code <- function(code) {
    if (length(code) == 0 || anyNA(code)) {
        return(list(ok = FALSE, error = "Generated OxCal code is empty or contains missing lines."))
    }
    collapsed <- paste(code, collapse = "\n")
    invalid_token <- grepl(
        "(^|[^A-Za-z0-9_.])(?:NA|NaN|-?Inf)(?=$|[^A-Za-z0-9_.])",
        collapsed,
        perl = TRUE
    )
    if (invalid_token) {
        return(list(ok = FALSE, error = "Generated OxCal code contains an invalid numeric token."))
    }
    list(ok = TRUE, error = NULL)
}

mod_chronology_ui <- function(id) {
    ns <- NS(id)
    
    tagList(
        bslib::navset_card_underline(
            id = ns("chron_tabs"),
            title = "Chronology Tools",
            
            # ====================================================================
            # TAB 1: Age-Depth Model
            # ====================================================================
            bslib::nav_panel(
                title = "Age-Depth Model",
                icon = icon("layer-group"),
                
                module_banner(
                  goal = "Generate OxCal P_Sequence code for continuous age-depth modeling.",
                  inputs = "CSV/Excel with Name, Age, Uncertainty, Depth, and Date Type columns.",
                  outputs = "OxCal code ready to paste into the online OxCal program (c14.arch.ox.ac.uk).",
                  why = "Building an age-depth model lets you estimate the age of any depth in your core. This is critical for dating tephras found between radiocarbon-dated horizons. The P_Sequence model uses Bayesian statistics to produce calibrated age ranges with proper uncertainty."
                ),
                layout_columns(
                    col_widths = c(4, 8),
                    
                    # --- LEFT: Configuration ---
                    div(
                        varg_card(
                            title = tagList(icon("upload"), " Data Upload"),
                            # Instructional Help Box
                            help_box(
                                title = "Quick Reference",
                                content = "Upload a CSV or Excel file with columns: <b>Name</b>, <b>Age</b>, <b>Uncertainty</b>, <b>Depth</b>, and <b>Type</b>.
                                <ul class='small mb-1'>
                                    <li><b>Download the template</b> below for the correct format with dropdown validation.</li>
                                    <li>Hover over any <b>(?)</b> icon for field-specific guidance.</li>
                                </ul>
                                <details class='small'><summary style='cursor:pointer; color:#0d6efd;'>Supported date types &amp; outlier models</summary>
                                <ul class='mt-1 mb-0'>
                                    <li><b>R_Date:</b> Radiocarbon (¹⁴C BP)</li>
                                    <li><b>R_F14C:</b> Fraction modern carbon</li>
                                    <li><b>Date / Tephra:</b> Calendar age (calBP, CE, or blank)</li>
                                    <li><b>Boundary:</b> Stratigraphic marker (no age)</li>
                                </ul>
                                <p class='mt-1 mb-0'>Outlier models: General, Charcoal, SSimple, RSimple, TSimple, RScaled.</p>
                                <p class='mt-1 mb-0'>See the <b>User Guide → Chronology Module</b> for full details.</p>
                                </details>"
                            ),
                            fileInput(ns("file_age_depth"), 
                                label = div("Upload Age-Depth File (.xlsx or .csv)", 
                                    help_icon("<strong>Upload a correctly formatted age-depth CSV/XLSX file.</strong><details><summary>Learn more</summary>Upload a CSV or Excel file with columns:<br><b>Name</b> (unique date label), <b>Age</b> (radiocarbon or calendar years), <b>Uncertainty</b> (±1\u03c3), <b>Depth</b> (position in core), <b>Type</b> (R_Date, R_F14C, Date, Tephra, Date_CE, Tephra_CE, or Boundary).<br><br>Optional: <b>Unc_Type</b> (N for normal or U for uniform), <b>Outlier_Type</b> (General, Charcoal, etc.).<br><br><b>Tip:</b> Download the template for the correct format.</details>")),
                                accept = c(".xlsx", ".csv")
                            ),
                            uiOutput(ns("age_depth_sheet_select")),
                            uiOutput(ns("age_depth_column_check")),
                            downloadLink(ns("download_template_ad"), 
                                label = "Download Template (.xlsx)", 
                                icon = icon("file-excel")
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("cogs"), " Model Settings"),
                            textInput(ns("core_name"), 
                                label = div("Core Name", help_icon("<strong>Name that will appear as the P_Sequence label in OxCal. Use a descriptive, unique identifier (e.g., 'Skilak_Core2', 'WL1_2024'). Spaces are fine.</strong>")), 
                                value = "My Core"
                            ),
                            selectInput(ns("curve_name"),
                                label = div("Primary Curve", help_icon("<strong>Choose the calibration curve matching sample setting.</strong><details><summary>Learn more</summary>The radiocarbon calibration curve to use.<br><br><b>IntCal20</b>: Northern Hemisphere terrestrial (Reimer et al. 2020). Use this for most NH land-based samples.<br><b>Marine20</b>: Marine samples (Heaton et al. 2020). This requires a \u0394R correction for local reservoir effects (you will need to specify this manually, e.g., using Reimer's database).<br><b>SHCal20</b>: Southern Hemisphere terrestrial (Hogg et al. 2020).<br><br><b>Why it matters:</b> Using the wrong curve introduces systematic age offsets of decades to centuries.</details>")),
                choices = c("IntCal20" = "IntCal20", "Marine20" = "Marine20", "SHCal20" = "SHCal20"),
                                selected = "IntCal20"
                            ),
                            checkboxInput(ns("use_bomb"), 
                                label = div("Use Bomb Curve", help_icon("<strong>Enable atmospheric bomb curves for modern or near-modern ¹⁴C dates.</strong><details><summary>Learn more</summary><b>Default: ON.</b><br><br><b>When to use:</b> If your core includes modern or near-modern ¹⁴C dates that may postdate 1950 CE.<br><br><b>How it works:</b> R_Dates with ¹⁴C ages &lt; 1000 BP are automatically calibrated against the selected bomb curve instead of the primary curve. Dates ≥ 1000 BP always use the primary curve.<br><br><b>Atmospheric only:</b> Bomb curves record atmospheric ¹⁴C and are only appropriate with IntCal20 or SHCal20. When using Marine20, bomb curves are not applicable and this option is ignored.<br><br>Without a bomb curve, OxCal may warn: <em>'Date may extend out of range'</em>.<br><br>Reference: Hua et al. 2021.</details>")), 
                                value = TRUE
                            ),
                            conditionalPanel(
                                condition = "input.use_bomb && input.curve_name != 'Marine20'",
                                ns = ns,
                                selectInput(ns("bomb_curve"),
                                    label = div("Post-Bomb Curve Zone", help_icon("<strong>Select the post-bomb ¹⁴C zone for your site location.</strong><details><summary>Learn more</summary>Select the geographic zone for your sampling site.<br><br>Each curve covers a different atmospheric \u00B9\u2074C region (Hua et al. 2021):<br><br><b>Bomb21NH1</b>: North of 40\u00B0N (prepended with IntCal20)<br><b>Bomb21NH2</b>: 40\u00B0N to mean summer ITCZ (prepended with IntCal20)<br><b>Bomb21NH3</b>: NH mean summer ITCZ (prepended with IntCal20)<br><b>Bomb21SH12</b>: SH outside mean summer ITCZ (prepended with SHCal20)<br><b>Bomb21SH3</b>: SH mean summer ITCZ (prepended with SHCal20)<br><br><b>Tip:</b> For most North American and European sites, use NH Zone 1.</details>")),
                                    choices = c(
                                        "NH Zone 1 (>40°N)" = "Bomb21NH1",
                                        "NH Zone 2 (40°N-ITCZ)" = "Bomb21NH2",
                                        "NH Zone 3 (ITCZ)" = "Bomb21NH3",
                                        "SH Zones 1-2" = "Bomb21SH12",
                                        "SH Zone 3" = "Bomb21SH3"
                                    ),
                                    selected = "Bomb21NH1"
                                )
                            ),
                            conditionalPanel(
                                condition = "input.curve_name == 'Marine20'",
                                ns = ns,
                                tags$div(
                                    class = "p-2 border rounded bg-light mb-2",
                                    tags$div(class = "fw-bold small mb-1", icon("water"), " Marine Reservoir Correction (\u0394R)"),
                                    tags$p(class = "small text-muted mb-2", 
                                        "Marine20 requires a local \u0394R correction. Look up your value at ",
                                        tags$a(href = "http://calib.org/marine/", target = "_blank", "calib.org/marine"),
                                        "."
                                    ),
                                    layout_columns(
                                        col_widths = c(6, 6),
                                        numericInput(ns("delta_r_value"),
                                            label = div("\u0394R Value", help_icon("<strong>Local marine reservoir offset in \u00b9\u2074C years.</strong><details><summary>Learn more</summary>\u0394R represents the difference between the global marine reservoir effect (already built into Marine20) and your local marine reservoir.<br><br><b>Positive values:</b> Local reservoir is older than global average.<br><b>Negative values:</b> Local reservoir is younger.<br><b>Zero:</b> Use global Marine20 without local correction.<br><br>Look up values at <a href='http://calib.org/marine/' target='_blank'>calib.org/marine</a>.<br><br>In OxCal: <code>Delta_R(\"Local\", value, error);</code></details>")),
                                            value = 0, step = 1
                                        ),
                                        numericInput(ns("delta_r_error"),
                                            label = div("\u0394R Error (\u00b11\u03c3)", help_icon("<strong>Uncertainty on the \u0394R value, in \u00b9\u2074C years.</strong>")),
                                            value = 0, step = 1
                                        )
                                    )
                                )
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("ruler-vertical"), " Depth Settings"),
                            layout_columns(
                                col_widths = c(6, 6),
                                numericInput(ns("top_depth"), 
                                    label = div("Top Depth", help_icon("<strong>Set an optional top boundary depth for P_Sequence.</strong><details><summary>Learn more</summary>Optional. Defines the depth for OxCal's top Boundary in your P_Sequence (OxCal requires boundaries at both ends).<br><br><b>When to use:</b> Set to 0 for surface cores where the top is known to be modern, or to the depth of your uppermost sample.<br><br>If left blank, the boundary defaults to the depth of your shallowest date.</details>")), 
                                    value = NA, step = 0.1
                                ),
                                numericInput(ns("bottom_depth"), 
                                    label = div("Bottom Depth", help_icon("<strong>Set an optional bottom boundary depth for P_Sequence.</strong><details><summary>Learn more</summary>Optional. Defines the depth for OxCal's bottom Boundary in your P_Sequence (OxCal requires boundaries at both ends).<br><br><b>When to use:</b> Set to the maximum depth of your dated section.<br><br>If left blank, the boundary defaults to the depth of your deepest date.</details>")), 
                                    value = NA, step = 0.1
                                )
                            ),
                            selectInput(ns("depth_unit"), 
                                label = div("Depth Unit", help_icon("<strong>Match depth units to the correct k0 value.</strong><details><summary>Learn more</summary>Sets the k0 parameter in OxCal's P_Sequence.<br><br><b>cm (k0=1):</b> Use when depths are in centimeters.<br><b>m (k0=100):</b> Use when depths are in meters.<br><br><b>Why it matters:</b> k0 controls the prior on accumulation rate variability. Mismatching units and k0 can make your age model over-constrained (too rigid) or under-constrained (too flexible).<br><br><b>Tip:</b> Always match this to the units in your depth column.</details>")), 
                                choices = c("cm", "m"), selected = "cm"
                            ),
                            layout_columns(
                                col_widths = c(10, 2),
                                numericInput(ns("interp_rate"), 
                                    label = div("Interpolation Rate", help_icon("<strong>Control interpolation density versus runtime.</strong><details><summary>Learn more</summary>How many interpolation points per depth unit.<br><br><b>Default: 0.05</b> per cm (one point every 20 cm).<br><br><b>When to adjust:</b> Smaller values (e.g., 0.02) = fewer segments = faster but coarser. Larger values (e.g., 0.1) = finer resolution but slower.<br><br>This mostly affects the appearance of the resulting uncertainty ribbon on the age-depth plot.<br><br>The auto-calculate button estimates a sensible rate from your depth range and date count.</details>")), 
                                    value = 0.05, min = 0.000001, step = 0.01
                                ),
                                div(
                                    style = "padding-top: 32px;",
                                    actionButton(ns("update_interp_rate"), "", 
                                        icon = icon("sync"), 
                                        title = "Auto-calculate rate",
                                        class = "btn-outline-secondary btn-sm"
                                    )
                                )
                            ),
                            uiOutput(ns("interp_rate_equivalent")),
                            uiOutput(ns("interp_warning")),
                            checkboxInput(ns("auto_rename_ad"), 
                                label = div("Auto-rename duplicate names", help_icon("<strong>Auto-suffix duplicate date names for OxCal compatibility.</strong><details><summary>Learn more</summary>Automatically suffix duplicate date names with _01, _02, etc.<br><br><b>Why?</b> OxCal requires unique names for every dated element.<br><br>We leave this optional so that it's a conscious choice. If your data has duplicates, you should be made aware.<br><br><b>Turn ON</b> if you have multiple dates with the same name and want automatic resolution.</details>")),
                                value = FALSE
                            ),
                            actionButton(ns("gen_age_depth"), "Generate OxCal Code", 
                                class = "btn-primary w-100 mt-2", 
                                icon = icon("code")
                            )
                        )
                    ),
                    
                    # --- RIGHT: Output ---
                    varg_card(
                        title = tagList(icon("code"), " Generated OxCal Code"),
                        full_screen = TRUE,
                        p(class = "text-muted small", 
                            "Copy this code and paste it into ", 
                            tags$a(href = "https://c14.arch.ox.ac.uk/oxcal.html", target = "_blank", "OxCal"),
                            " to run your age-depth model."
                        ),
                        actionButton(ns("copy_age_depth"), "Copy to Clipboard", 
                            icon = icon("clipboard"), 
                            class = "btn-outline-primary btn-sm mb-2 varg-copy-to-clipboard",
                            `data-copy-target` = ns("out_age_depth"),
                            `data-copy-status-input` = ns("copy_age_depth_status")
                        ),
                        verbatimTextOutput(ns("out_age_depth"))
                    )
                )
            ),
            
            # ====================================================================
            # TAB 2: Phase Model
            # ====================================================================
            bslib::nav_panel(
                title = "Phase Model",
                icon = icon("sitemap"),
                
                module_banner(
                  goal = "Generate OxCal Phase/Sequence code to constrain tephra ages relative to bracketing dates.",
                  inputs = "CSV/Excel with Tephra Name, Phase (Before/During/After), Date Type, Age, and SD.",
                  outputs = "OxCal Sequence code with Phase boundaries, outlier models, and calibration curves.",
                  why = "When you have radiocarbon dates bracketing a tephra but no continuous depth-series, a Phase model constrains the tephra's age using Bayesian sequencing. This is common for archaeological sites or exposures where only a few dates frame the eruption event."
                ),
                layout_columns(
                    col_widths = c(4, 8),
                    
                    # --- LEFT: Configuration ---
                    div(
                        varg_card(
                            title = tagList(icon("upload"), " Data Upload"),
                            # Instructional Help Box
                            help_box(
                                title = "Quick Reference",
                                content = "Upload a CSV or Excel file with columns: <b>Tephra Name</b>, <b>Phase</b> (Before/During/After), <b>Name</b>, <b>Date Type</b>, <b>Age Mean</b>, <b>Age SD</b>.
                                <ul class='small mb-1'>
                                    <li><b>Download the template</b> below for the correct format.</li>
                                    <li>Same date types and outlier models as the Age-Depth tab.</li>
                                    <li>Hover over any <b>(?)</b> icon for field-specific guidance.</li>
                                </ul>
                                <p class='small mb-0'>See the <b>User Guide → Phase Model</b> for full details.</p>"
                            ),
                            fileInput(ns("file_phase"), 
                                label = div("Upload Phase Model File (.xlsx or .csv)", 
                                    help_icon("<strong>Upload a correctly formatted phase-model CSV/XLSX file.</strong><details><summary>Learn more</summary>Upload a CSV or Excel with columns:<br><b>Tephra Name</b>, <b>Phase</b> (Before/During/After), <b>Name</b> (date label), <b>Date Type</b> (R_Date, R_F14C, Date, Tephra, Date_CE, or Tephra_CE), <b>Age Mean</b>, <b>Age SD</b>.<br><br>Optional: <b>Outlier_Type/Outlier.Prob</b>, <b>Unc_Type</b> (N for normal or U for uniform).<br><br>Supported outlier models: General, Charcoal, SSimple, RSimple, TSimple, RScaled.<br><br><b>Tip:</b> Download the tephra chronology template for the correct format.</details>")),
                                accept = c(".xlsx", ".csv")
                            ),
                            uiOutput(ns("phase_sheet_select")),
                            downloadLink(ns("download_template_phase"), 
                                label = "Download Template (.xlsx)", 
                                icon = icon("file-excel")
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("sliders-h"), " Calibration Settings"),
                            selectInput(ns("phase_curve_name"),
                                label = div("Primary Curve", help_icon("<strong>Use the same curve logic as P_Sequence calibration.</strong><details><summary>Learn more</summary>Same as P_Sequence.<br><br><b>IntCal20</b>: NH terrestrial.<br><b>Marine20</b>: Marine (requires \u0394R).<br><b>SHCal20</b>: SH terrestrial.<br><br>Must match the hemisphere and environment of your dated samples.</details>")),
                                choices = c("IntCal20" = "IntCal20", "Marine20" = "Marine20", "SHCal20" = "SHCal20"),
                                selected = "IntCal20"
                            ),
                            checkboxInput(ns("phase_use_bomb"), 
                                label = div("Use Bomb Curve", help_icon("<strong>Enable atmospheric bomb curves for modern or near-modern \u00b9\u2074C dates.</strong><details><summary>Learn more</summary>Same as P_Sequence.<br><br><b>Default: ON.</b><br><br><b>Atmospheric only:</b> Bomb curves are not applicable with Marine20 and will be ignored.<br><br>Reference: Hua et al. 2021.</details>")), 
                                value = TRUE
                            ),
                            conditionalPanel(
                                condition = "input.phase_use_bomb && input.phase_curve_name != 'Marine20'",
                                ns = ns,
                                selectInput(ns("phase_bomb_curve"),
                                    label = div("Post-Bomb Curve Zone", help_icon("<strong>Select the geographic post-bomb ¹⁴C calibration zone.</strong><details><summary>Learn more</summary>Same as P_Sequence.<br><br>Select the geographic zone for your site's post-bomb \u00b9\u2074C curve (Hua et al. 2021).<br><br><b>Default: NH Zone 1</b> (most of North America and Europe).</details>")),
                                    choices = c(
                                        "NH Zone 1 (>40°N)" = "Bomb21NH1",
                                        "NH Zone 2 (40°N-ITCZ)" = "Bomb21NH2",
                                        "NH Zone 3 (ITCZ)" = "Bomb21NH3",
                                        "SH Zones 1-2" = "Bomb21SH12",
                                        "SH Zone 3" = "Bomb21SH3"
                                    ),
                                    selected = "Bomb21NH1"
                                )
                            ),
                            conditionalPanel(
                                condition = "input.phase_curve_name == 'Marine20'",
                                ns = ns,
                                tags$div(
                                    class = "p-2 border rounded bg-light mb-2",
                                    tags$div(class = "fw-bold small mb-1", icon("water"), " Marine Reservoir Correction (\u0394R)"),
                                    tags$p(class = "small text-muted mb-2", 
                                        "Marine20 requires a local \u0394R correction. Look up your value at ",
                                        tags$a(href = "http://calib.org/marine/", target = "_blank", "calib.org/marine"),
                                        "."
                                    ),
                                    layout_columns(
                                        col_widths = c(6, 6),
                                        numericInput(ns("phase_delta_r_value"),
                                            label = div("\u0394R Value", help_icon("<strong>Local marine reservoir offset in \u00b9\u2074C years. Same as P_Sequence. See calib.org/marine for values.</strong>")),
                                            value = 0, step = 1
                                        ),
                                        numericInput(ns("phase_delta_r_error"),
                                            label = div("\u0394R Error (\u00b11\u03c3)", help_icon("<strong>Uncertainty on the \u0394R value, in \u00b9\u2074C years.</strong>")),
                                            value = 0, step = 1
                                        )
                                    )
                                )
                            ),
                            checkboxInput(ns("auto_rename_phase"), 
                                label = div("Auto-rename duplicate names", help_icon("<strong>Auto-suffix duplicate names to keep OxCal identifiers unique.</strong><details><summary>Learn more</summary>Automatically suffix duplicate date names with _01, _02, etc.<br><br>Same as P_Sequence: ensures OxCal gets unique names.<br><br>We leave this optional so that it's a conscious choice. You should be aware if your data has duplicates.</details>")),
                                value = FALSE
                            ),
                            tags$hr(class = "my-2"),
                            tags$h6("Phase Boundary Settings"),
                            checkboxInput(ns("use_tau_before"), 
                                label = tags$span("Use Tau_Boundary for Before phase", help_icon("<strong>Choose Tau_Boundary when pre-tephra dates cluster near event.</strong><details><summary>Learn more</summary><b>Tau_Boundary</b>: Assumes dates cluster exponentially near the tephra (more likely close, less likely far). Use when bracketing dates were specifically collected close to the tephra.<br><br><b>Boundary (unchecked)</b>: Assumes uniform distribution of dates throughout the phase. Use when dates are randomly or evenly distributed.<br><br><b>Reference:</b> Bronk Ramsey 2009, 'Bayesian Analysis of Radiocarbon Dates', <i>Radiocarbon</i> 51(1).</details>")), 
                                value = TRUE
                            ),
                            checkboxInput(ns("use_tau_after"), 
                                label = tags$span("Use Tau_Boundary for After phase", help_icon("<strong>Choose Tau_Boundary when post-tephra dates cluster near event.</strong><details><summary>Learn more</summary>Same as above, but for the After (post-eruption) phase.<br><br><b>Tau_Boundary</b> if dates are expected to cluster close to the tephra; <b>Boundary</b> if dates are distributed throughout the phase.</details>")), 
                                value = TRUE
                            ),
                            actionButton(ns("gen_phase"), "Generate OxCal Code", 
                                class = "btn-primary w-100 mt-2", 
                                icon = icon("code")
                            )
                        )
                    ),
                    
                    # --- RIGHT: Output ---
                    varg_card(
                        title = tagList(icon("code"), " Generated OxCal Code"),
                        full_screen = TRUE,
                        p(class = "text-muted small", 
                            "Copy this code and paste it into ", 
                            tags$a(href = "https://c14.arch.ox.ac.uk/oxcal.html", target = "_blank", "OxCal"),
                            " to run your phase model."
                        ),
                        actionButton(ns("copy_phase"), "Copy to Clipboard", 
                            icon = icon("clipboard"), 
                            class = "btn-outline-primary btn-sm mb-2 varg-copy-to-clipboard",
                            `data-copy-target` = ns("out_phase"),
                            `data-copy-status-input` = ns("copy_phase_status")
                        ),
                        verbatimTextOutput(ns("out_phase"))
                    )
                )
            ),
            
            # ====================================================================
            # TAB 3: Model Linker
            # ====================================================================
            bslib::nav_panel(
                title = "Model Linker",
                icon = icon("link"),
                
                module_banner(
                  goal = "Link two independent OxCal models via shared tephra tie-points into one regional chronology.",
                  inputs = "Two OxCal model code blocks (from Tab 1 or external), plus tie-point JSON from the Visualization module.",
                  outputs = "Combined OxCal script with cross-referencing constraints ready for execution.",
                  why = "If you've correlated tephras between two sites (using the Stratigraphic Correlation tool), you can combine their independent age models into one linked regional chronology. Shared isochrons provide cross-constraints that tighten age estimates at both sites simultaneously."
                ),
                layout_columns(
                    col_widths = c(4, 8),
                    
                    # --- LEFT: Configuration ---
                    div(
                        varg_card(
                            title = tagList(icon("file-alt"), " Model A"),
                            # Instructional Help Box
                            help_box(
                                title = "Model Linker",
                                content = "Combine two separate OxCal models using shared tie-points.
                                <ul>
                                    <li><b>Purpose:</b> Create a regional chronology by linking independent site models via tephra isochrons.</li>
                                    <li><b>Workflow:</b> Paste code for Model A and Model B, then define which events are equivalent (Tie Points).</li>
                                    <li><b>Result:</b> A single OxCal script that runs both models simultaneously with cross-referencing constraints.</li>
                                </ul>"
                            ),
                            textAreaInput(ns("model_a_code"), 
                                label = div("OxCal Code", help_icon("<strong>Paste full OxCal code for Model A.</strong><details><summary>Learn more</summary>Paste the full OxCal code for Model A.<br><br>The Plot() wrapper is removed automatically if present.<br><br><b>Tip:</b> You can paste code generated from Tab 1 (P_Sequence) or Tab 2 (Tephra Chronology), or any manually written OxCal model.</details>")),
                                rows = 8, 
                                placeholder = "Paste OxCal code for Model A..."
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("file-alt"), " Model B"),
                            textAreaInput(ns("model_b_code"), 
                                label = div("OxCal Code", help_icon("<strong>Paste full OxCal code for Model B.</strong><details><summary>Learn more</summary>Paste the full OxCal code for Model B.<br><br>Same rules as Model A.<br><br><b>Tip:</b> Compare two versions of the same model (e.g., with vs. without outlier detection) or two different cores.</details>")),
                                rows = 8, 
                                placeholder = "Paste OxCal code for Model B..."
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("link"), " Tie Points"),
                            actionButton(ns("import_tiepoints"), "Import from Viz", 
                                icon = icon("download"), 
                                class = "btn-info btn-sm w-100 mb-2"
                            ),
                            textAreaInput(ns("tie_points"), 
                                label = div("Define Tie Points", help_icon("<strong>Define event tie-points linking Model A and Model B.</strong><details><summary>Learn more</summary>Define stratigraphic tie points between the two models.<br><br>Format: <b>TieName, NameInModelA, NameInModelB</b> (one per line).<br><br>The tool will correlate these dated events across models to compare age estimates.<br><br><b>Example:</b> 'Aniakchak_CFE, Aniakchak_A, Aniakchak_B'.</details>")), 
                                rows = 6, 
                                placeholder = "Tie1, NameInA, NameInB\nTie2, NameInA, NameInB"
                            ),
                            checkboxInput(ns("link_append_suffix"), 
                                label = div("Append suffixes (_A/_B)", help_icon("<strong>Add suffixes to prevent cross-model name collisions.</strong><details><summary>Learn more</summary>Add _A and _B suffixes to all date names in Model A and Model B respectively.<br><br><b>When to use:</b> If both models contain dates with the same names (common when comparing two cores with the same tephra layers).<br><br>We leave this optional so it's a conscious choice. You should verify that names are consistent across models.</details>")),
                                value = FALSE
                            ),
                            uiOutput(ns("link_warning")),
                            actionButton(ns("link_models"), "Link Models", 
                                class = "btn-primary w-100 mt-2", 
                                icon = icon("link")
                            )
                        ),
                        
                        varg_card(
                            title = tagList(icon("info-circle"), " Tie Point Format"),
                            tags$p(class = "small mb-1", "One tie point per line:"),
                            tags$code(class = "small", "TieName, NameInModelA, NameInModelB"),
                            tags$details(class = "mt-2",
                                tags$summary(tags$small("Example")),
                                tags$pre(class = "small p-2 bg-light", 
                                    "Tephra X, Core A Sample 5, Core B Sample 12\nTephra Y, Core A Sample 12, Core B Sample 25"
                                )
                            )
                        )
                    ),
                    
                    # --- RIGHT: Output ---
                    varg_card(
                        title = tagList(icon("code"), " Linked OxCal Code"),
                        full_screen = TRUE,
                        p(class = "text-muted small", 
                            "Tie points create cross-references between chronologies using OxCal's '=' syntax."
                        ),
                        actionButton(ns("copy_linked"), "Copy to Clipboard", 
                            icon = icon("clipboard"), 
                            class = "btn-outline-primary btn-sm mb-2 varg-copy-to-clipboard",
                            `data-copy-target` = ns("out_linked"),
                            `data-copy-status-input` = ns("copy_linked_status")
                        ),
                        verbatimTextOutput(ns("out_linked"))
                    )
                )
            )
        )
    )
}

mod_chronology_server <- function(id, tiepoints_reactive = NULL, global_rv = NULL) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # --- Interpolation Rate Logic ---

        calculate_interp_rate <- function(df, top = NA, bottom = NA) {
            if (is.null(df) || nrow(df) == 0) {
                return(0.05)
            }

            # Count valid dates (R_Date, R_F14C, Date, CE, etc.)
            valid_types <- c("R_Date", "R_F14C", "Date", "Tephra", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE", "CE")
            n_dates <- 0
            if ("Type" %in% names(df)) {
                n_dates <- sum(df$Type %in% valid_types, na.rm = TRUE)
            } else if ("Date.Type" %in% names(df)) {
                n_dates <- sum(df$Date.Type %in% valid_types, na.rm = TRUE)
            }

            if (n_dates == 0) {
                return(0.05)
            }

            # Determine depth range
            min_z <- NA
            max_z <- NA

            if ("Depth" %in% names(df)) {
                depths <- as.numeric(df$Depth)
                depths <- depths[!is.na(depths)]
                if (length(depths) > 0) {
                    min_z <- min(depths)
                    max_z <- max(depths)
                }
            }

            # Override with manual boundaries if provided
            if (!is.null(top) && !is.na(top)) min_z <- top
            if (!is.null(bottom) && !is.na(bottom)) max_z <- bottom

            if (is.na(min_z) || is.na(max_z) || min_z >= max_z) {
                return(0.05)
            }

            range_z <- max_z - min_z

            # Formula: 1.3 * (N_dates / Depth_Range)
            rate <- 1.3 * (n_dates / range_z)

            # Round to 3 decimal places
            rate <- round(rate, 3)

            # Safety check
            if (rate <= 0) {
                return(0.05)
            }

            return(rate)
        }

        observeEvent(input$update_interp_rate, {
            req(rv_age_depth_data())
            df <- rv_age_depth_data()
            rate <- calculate_interp_rate(df, input$top_depth, input$bottom_depth)
            updateNumericInput(session, "interp_rate", value = rate)
            showNotification(paste("Interpolation rate updated to", rate), type = "message")
        })

        output$interp_rate_equivalent <- renderUI({
            req(input$interp_rate)
            if (is.na(input$interp_rate) || input$interp_rate <= 0) {
                return(NULL)
            }

            spacing_in_selected_unit <- 1 / input$interp_rate
            spacing_cm <- if (!is.null(input$depth_unit) && identical(input$depth_unit, "m")) {
                spacing_in_selected_unit * 100
            } else {
                spacing_in_selected_unit
            }

            tags$div(
                class = "text-muted small",
                style = "margin-top: -4px; margin-bottom: 6px;",
                paste0("Equivalent resolution: 1 point every ", round(spacing_cm, 2), " cm")
            )
        })

        output$interp_warning <- renderUI({
            req(input$interp_rate, rv_age_depth_data())
            rate <- input$interp_rate
            df <- rv_age_depth_data()

            # Calculate range
            min_z <- NA
            max_z <- NA
            if ("Depth" %in% names(df)) {
                depths <- as.numeric(df$Depth)
                depths <- depths[!is.na(depths)]
                if (length(depths) > 0) {
                    min_z <- min(depths)
                    max_z <- max(depths)
                }
            }
            if (!is.na(input$top_depth)) min_z <- input$top_depth
            if (!is.na(input$bottom_depth)) max_z <- input$bottom_depth

            if (is.na(min_z) || is.na(max_z) || min_z >= max_z) {
                return(NULL)
            }

            range_z <- max_z - min_z
            n_interp <- range_z * rate

            if (n_interp > 100) {
                tags$div(
                    style = "color: red; font-weight: bold; margin-top: 5px;",
                    icon("exclamation-triangle"),
                    paste("Warning: High interpolation count (~", round(n_interp), "points). This may be very slow!")
                )
            } else {
                NULL
            }
        })

        observeEvent(input$age_depth_sheet, {
            req(input$file_age_depth, input$age_depth_sheet)
            generated_age_depth_code("")
            rv_age_depth_data(NULL)
            rv_age_depth_required_cols(NULL)
            tryCatch(
                {
                    df <- readxl::read_excel(input$file_age_depth$datapath, sheet = input$age_depth_sheet)
                    df <- normalize_age_depth_columns(df)
                    rv_age_depth_data(df)
                    update_age_depth_column_status(df)

                    # Calculate default interpolation rate
                    rate <- calculate_interp_rate(df)
                    updateNumericInput(session, "interp_rate", value = rate)

                    showNotification("Data loaded successfully", type = "message")
                },
                error = function(e) {
                    showNotification(paste("Error reading file:", e$message), type = "error")
                }
            )
        })

        # --- Linker Logic ---

        # Helper to extract names from OxCal code (R_Date("Name",...) or Date("Name",...))
        extract_oxcal_names <- function(code_text) {
            # Regex to find names inside R_Date("...", or Date("...",
            # Matches: Word chars, spaces, hyphens, etc.
            # R_Date("Name"  or  Date("Name"
            pattern <- '(?:R_Date|Date)\\s*\\(\\s*"([^"]+)"'
            matches <- gregexpr(pattern, code_text)
            names_list <- regmatches(code_text, matches)

            extracted <- c()
            for (m in names_list[[1]]) {
                # Clean up to just the name
                name <- sub('(?:R_Date|Date)\\s*\\(\\s*"([^"]+)"', "\\1", m)
                extracted <- c(extracted, name)
            }
            return(extracted)
        }

        # Reactive value to hold Age-Depth data (from file or import)
        rv_age_depth_data <- reactiveVal(NULL)
        rv_age_depth_required_cols <- reactiveVal(NULL)

        update_age_depth_column_status <- function(df) {
            required_cols <- c("Name", "Age", "Uncertainty", "Depth", "Type")
            present <- setNames(required_cols %in% names(df), required_cols)
            rv_age_depth_required_cols(present)

            missing_cols <- names(present)[!present]
            if (length(missing_cols) > 0) {
                showNotification(
                    paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
                    type = "warning",
                    duration = 8
                )
            }
        }

        output$age_depth_column_check <- renderUI({
            status <- rv_age_depth_required_cols()
            req(status)

            tags$div(
                class = "mt-2 mb-2 p-2 border rounded bg-light",
                tags$div(class = "small fw-bold mb-1", "Required Columns Check"),
                lapply(names(status), function(col_name) {
                    ok <- isTRUE(status[[col_name]])
                    tags$div(
                        class = if (ok) "text-success small" else "text-danger small",
                        icon(if (ok) "check-circle" else "exclamation-triangle"),
                        paste0(" ", col_name)
                    )
                })
            )
        })

        # Helper to normalize column names for age-depth data
        normalize_age_depth_columns <- function(df) {
            rename_col <- function(df, patterns, new_name) {
                match_idx <- which(tolower(names(df)) %in% tolower(patterns))
                if (length(match_idx) > 0) {
                    names(df)[match_idx[1]] <- new_name
                }
                return(df)
            }

            df <- rename_col(df, c("Name", "Sample Name", "Sample ID", "Lab Code"), "Name")
            df <- rename_col(df, c("Age", "Age Mean", "Age_Mean", "Age.Mean"), "Age")
            df <- rename_col(df, c("Uncertainty", "Error", "Age SD", "Age_SD", "Age.SD", "sd"), "Uncertainty")
            df <- rename_col(df, c("Depth"), "Depth")
            df <- rename_col(df, c("Type", "Date Type", "Date_Type", "Date.Type"), "Type")
            df <- rename_col(df, c("Outlier type", "Outlier Type", "Outlier.Type", "Outlier.type", "Outlier_Type", "Outlier_type"), "Outlier.Type")
            df <- rename_col(df, c("Outlier prob", "Outlier Prob", "Outlier.Prob", "Outlier.prob.", "Outlier_Prob", "Outlier_prob"), "Outlier.Prob")

            return(df)
        }

        # --- Sheet Selection UI for Age-Depth ---
        output$age_depth_sheet_select <- renderUI({
            req(input$file_age_depth)
            file_ext <- tolower(tools::file_ext(input$file_age_depth$name))
            if (file_ext == "xlsx") {
                sheets <- readxl::excel_sheets(input$file_age_depth$datapath)
                selectInput(ns("age_depth_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
            }
        })

        # --- Download Handlers ---
        output$download_template_ad <- downloadHandler(
            filename = function() {
                "VARG_AgeDepth_Template.xlsx"
            },
            content = function(file) {
                # Generate template using advanced template generation with dropdowns
                if (!exists("generate_age_depth_template")) {
                    source("generate_chronology_templates.R", local = TRUE)
                }
                generate_age_depth_template(output_file = file, verbose = FALSE)
            }
        )

        output$download_template_phase <- downloadHandler(
            filename = function() {
                "VARG_Phase_Template.xlsx"
            },
            content = function(file) {
                # Generate template using advanced template generation with dropdowns
                if (!exists("generate_phase_template")) {
                    source("generate_chronology_templates.R", local = TRUE)
                }
                generate_phase_template(output_file = file, verbose = FALSE)
            }
        )


        # Observer for file upload
        observeEvent(input$file_age_depth, {
            req(input$file_age_depth)
            generated_age_depth_code("")
            rv_age_depth_data(NULL)
            rv_age_depth_required_cols(NULL)
            file_ext <- tolower(tools::file_ext(input$file_age_depth$name))
            if (file_ext == "xlsx") {
                # Load first sheet immediately for early column validation feedback
                tryCatch(
                    {
                        sheets <- readxl::excel_sheets(input$file_age_depth$datapath)
                        if (length(sheets) > 0) {
                            df <- readxl::read_excel(input$file_age_depth$datapath, sheet = sheets[1])
                            df <- normalize_age_depth_columns(df)
                            rv_age_depth_data(df)
                            update_age_depth_column_status(df)

                            # Calculate default interpolation rate
                            rate <- calculate_interp_rate(df)
                            updateNumericInput(session, "interp_rate", value = rate)
                        }
                    },
                    error = function(e) {
                        showNotification(paste("Error reading Excel file:", e$message), type = "error")
                    }
                )
                return()
            }
            tryCatch(
                {
                    df <- read.csv(input$file_age_depth$datapath, stringsAsFactors = FALSE)
                    df <- normalize_age_depth_columns(df)
                    rv_age_depth_data(df)
                    update_age_depth_column_status(df)

                    # Calculate default interpolation rate
                    rate <- calculate_interp_rate(df)
                    updateNumericInput(session, "interp_rate", value = rate)

                    showNotification("Data loaded successfully", type = "message")
                },
                error = function(e) {
                    showNotification(paste("Error reading file:", e$message), type = "error")
                }
            )
        })

        # Observer for Tie Point Import
        observeEvent(input$import_tiepoints, {
            if (is.null(tiepoints_reactive)) {
                showNotification("Tie points data not available.", type = "warning")
                return()
            }
            tp_data <- tiepoints_reactive()
            if (is.null(tp_data) || nrow(tp_data) == 0) {
                showNotification("No tie points to import.", type = "warning")
                return()
            }

            # Format: TieName, NameInModelA (Ref), NameInModelB (Target)
            # Use custom_name if provided, otherwise default to TP_id
            lines <- apply(tp_data, 1, function(row) {
                tie_name <- if (!is.null(row["custom_name"]) && !is.na(row["custom_name"]) && row["custom_name"] != "") {
                    row["custom_name"]
                } else {
                    paste0("TP_", row["id"])
                }
                sprintf("%s, %s, %s", tie_name, row["ref_sample"], row["target_sample"])
            })
            text_block <- paste(lines, collapse = "\n")

            updateTextAreaInput(session, "tie_points", value = text_block)
            showNotification(paste("Imported", nrow(tp_data), "tie points into Model Linker."), type = "message")
        })

        # --- Helper: Curve Switching Logic ---
        # Returns a list: list(code = c(...), current_curve = "...", defined_curves = c(...))
        apply_curve_logic <- function(data, primary_curve, use_bomb, bomb_curve_name = "Bomb21NH1", current_curve = NA, defined_curves = c()) {
            code <- c()

            # Define primary curve initially if not set
            if (is.null(primary_curve) || is.na(primary_curve)) primary_curve <- "IntCal20"

            # If current_curve is NA (start of sequence), set it to primary
            if (is.na(current_curve)) {
                current_curve <- primary_curve

                # Emit initial curve definition
                if (primary_curve == "IntCal20") {
                    code <- c(code, '   Curve("IntCal20", "intcal20.14c");')
                    defined_curves <- c(defined_curves, "IntCal20")
                } else if (primary_curve == "Marine20") {
                    code <- c(code, '   Curve("Marine20", "marine20.14c");')
                    defined_curves <- c(defined_curves, "Marine20")
                } else if (primary_curve == "SHCal20") {
                    code <- c(code, '   Curve("SHCal20", "shcal20.14c");')
                    defined_curves <- c(defined_curves, "SHCal20")
                }
            }

            if (nrow(data) == 0) {
                return(list(code = code, current_curve = current_curve, defined_curves = defined_curves))
            }

            # Map bomb curve name to file
            bomb_curve_files <- list(
                "Bomb21NH1" = "bomb21nh1.14c",
                "Bomb21NH2" = "bomb21nh2.14c",
                "Bomb21NH3" = "bomb21nh3.14c",
                "Bomb21SH12" = "bomb21sh12.14c",
                "Bomb21SH3" = "bomb21sh3.14c"
            )
            bomb_curve_file <- bomb_curve_files[[bomb_curve_name]]
            if (is.null(bomb_curve_file)) bomb_curve_file <- "bomb21nh1.14c"

            for (i in seq_len(nrow(data))) {
                row <- data[i, ]

                # Safe extraction
                age_val <- NA
                if ("Age" %in% names(row)) {
                    age_val <- row$Age
                } else if ("Age.Mean" %in% names(row)) age_val <- row$Age.Mean

                age <- as.numeric(age_val) # NA if age_val is NA or non-numeric string

                type <- NA
                if ("Type" %in% names(row)) {
                    type <- row$Type
                } else if ("Date.Type" %in% names(row)) type <- row$Date.Type

                target_curve <- chronology_select_curve(
                    type = type,
                    age = age,
                    primary_curve = primary_curve,
                    use_bomb = use_bomb,
                    bomb_curve_name = bomb_curve_name
                )

                # Switch if needed - handle NAs safely
                curves_differ <- FALSE
                if (!is.na(target_curve) && !is.na(current_curve)) {
                    curves_differ <- target_curve != current_curve
                } else if (is.na(target_curve) != is.na(current_curve)) {
                    curves_differ <- TRUE
                }

                if (curves_differ) {
                    if (target_curve %in% names(bomb_curve_files)) {
                        # Switching TO Bomb Curve
                        if (target_curve %in% defined_curves) {
                            # Already defined, just reference it
                            code <- c(code, sprintf('   Curve("=%s");', target_curve))
                        } else {
                            # Define it for the first time
                            code <- c(code, sprintf('   Curve("%s", "%s");', target_curve, bomb_curve_file))
                            defined_curves <- c(defined_curves, target_curve)
                        }
                    } else {
                        # Switching TO Primary Curve (or other standard curve)
                        # We assume primary curve was defined at start, so just reference it
                        code <- c(code, sprintf('   Curve("=%s");', target_curve))
                    }
                    current_curve <- target_curve
                }

                # Format the date line
                line <- format_date_generic(row)
                code <- c(code, paste0("   ", line))
            }
            return(list(code = code, current_curve = current_curve, defined_curves = defined_curves))
        }

        # Helper function to extract unique outlier types from data and generate Outlier_Model declarations
        get_outlier_models_from_data <- function(df) {
            # Mapping of outlier types to their Outlier_Model declarations
            outlier_model_defs <- list(
                "General" = '  Outlier_Model("General", T(5), U(0,4), "t");',
                "SSimple" = '  Outlier_Model("SSimple", N(0,2), 0, "s");',
                "RSimple" = '  Outlier_Model("RSimple", N(0,100), 0, "r");',
                "TSimple" = '  Outlier_Model("TSimple", N(0,100), 0, "t");',
                "RScaled" = '  Outlier_Model("RScaled", T(5), U(0,4), "r");',
                "Charcoal" = '  Outlier_Model("Charcoal", Exp(1,-10,0), U(0,3), "t");'
            )

            outlier_types <- chronology_outlier_model_types(df)

            # Generate code for needed outlier models
            code <- c()
            for (otype in outlier_types) {
                if (otype %in% names(outlier_model_defs)) {
                    code <- c(code, outlier_model_defs[[otype]])
                }
            }
            
            return(code)
        }

        format_date_generic <- function(row) {
            # Handle both column naming conventions (Age Depth vs Phase)
            name <- if ("Name" %in% names(row)) row$Name else "Unknown"

            age_val <- NA
            if ("Age" %in% names(row)) {
                age_val <- row$Age
            } else if ("Age.Mean" %in% names(row)) age_val <- row$Age.Mean
            age <- age_val

            error_val <- NA
            if ("Uncertainty" %in% names(row)) {
                error_val <- row$Uncertainty
            } else if ("Age.SD" %in% names(row)) error_val <- row$Age.SD
            error <- error_val

            depth <- if ("Depth" %in% names(row)) row$Depth else NA

            type <- NA
            if ("Type" %in% names(row)) {
                type <- row$Type
            } else if ("Date.Type" %in% names(row)) type <- row$Date.Type

            # Uncertainty distribution type: N (normal, default) or U (uniform)
            unc_type <- "N"
            if ("Unc_Type" %in% names(row)) {
                ut <- row$Unc_Type
                if (!is.na(ut) && !is.null(ut) && toupper(trimws(ut)) %in% c("U", "UNIFORM")) {
                    unc_type <- "U"
                }
            } else if ("Unc.Type" %in% names(row)) {
                ut <- row$Unc.Type
                if (!is.na(ut) && !is.null(ut) && toupper(trimws(ut)) %in% c("U", "UNIFORM")) {
                    unc_type <- "U"
                }
            }

            # Outlier logic
            outlier_str <- "" # Default empty, only add if needed

            # Determine outlier type and probability
            otype <- NA
            if ("Outlier.Type" %in% names(row)) otype <- row$Outlier.Type
            if ("Outlier.type" %in% names(row)) otype <- row$Outlier.type # Fallback
            if ("Outlier_Type" %in% names(row)) otype <- row$Outlier_Type # Underscore variant
            if ("Outlier_type" %in% names(row)) otype <- row$Outlier_type # Underscore variant

            # Trim whitespace from outlier type
            if (!is.na(otype) && !is.null(otype)) {
                otype <- trimws(as.character(otype))
            }

            # Check if outlier info is present
            if (!is.na(otype) && nzchar(otype)) {
                # Validate outlier type - default to General if unknown
                valid_outliers <- c("General", "Charcoal", "T_Outlier", "SSimple", "RSimple", "TSimple", "RScaled")
                # Case insensitive check
                match_idx <- which(tolower(valid_outliers) == tolower(otype))
                if (length(match_idx) == 0) {
                    otype <- "General"
                } else {
                    otype <- valid_outliers[match_idx]
                }
                otype <- chronology_normalize_outlier_type(otype)

                oprob <- NA
                if ("Outlier.Prob" %in% names(row)) oprob <- row$Outlier.Prob
                if ("Outlier.prob." %in% names(row)) oprob <- row$Outlier.prob.
                if ("Outlier_Prob" %in% names(row)) oprob <- row$Outlier_Prob
                if ("Outlier_prob" %in% names(row)) oprob <- row$Outlier_prob
                if ("Manual.outlier.prob." %in% names(row)) oprob <- row$Manual.outlier.prob.

                # Defaults if probability is missing
                if (is.null(oprob) || length(oprob) == 0 || is.na(oprob) || as.character(oprob) == "") {
                    if (tolower(otype) == "charcoal") {
                        oprob <- 1
                    } else {
                        oprob <- 0.05
                    }
                }

                outlier_str <- sprintf('Outlier("%s", %s)', otype, oprob)
            } else {
                if (!is.na(type) && (type == "R_Date" || type == "R_F14C" || type == "Date" || type == "CE")) {
                    # Default to General if not specified
                    otype <- "General"
                    oprob <- 0.05
                    outlier_str <- sprintf('Outlier("%s", %s)', otype, oprob)
                }
            }

            z_str <- ""
            if (!is.null(depth) && !is.na(depth)) {
                z_str <- sprintf("z=%s;", depth)
            }

            # Helper to build the content inside braces without leading stray semicolons
            build_brace <- function(outlier_str, z_str) {
                parts <- c()
                if (nzchar(outlier_str)) parts <- c(parts, paste0(outlier_str, ";"))
                if (nzchar(z_str)) parts <- c(parts, z_str)
                if (length(parts) > 0) {
                    return(paste0("{", paste(parts, collapse = " "), "}"))
                }
                return("")
            }

            if (is.na(type)) {
                return("")
            }

            # Normalize type
            type <- trimws(type)

            # Helper function to format distribution based on unc_type
            # For U (uniform): U(age-error, age+error)
            # For N (normal): N(wrapper(age), error)
            format_dist <- function(wrapper_fn, age, error, unc_type) {
                if (unc_type == "U") {
                    from_val <- as.numeric(age) - as.numeric(error)
                    to_val <- as.numeric(age) + as.numeric(error)
                    return(sprintf('U(%s(%s), %s(%s))', wrapper_fn, from_val, wrapper_fn, to_val))
                } else {
                    return(sprintf('N(%s(%s), %s)', wrapper_fn, age, error))
                }
            }

            if (type == "R_Date") {
                return(sprintf('R_Date("%s", %s, %s)%s;', name, age, error, build_brace(outlier_str, z_str)))
            } else if (type == "R_F14C") {
                # R_F14C uses Fraction Modern Carbon - no curve needed, just F14C value
                return(sprintf('R_F14C("%s", %s, %s)%s;', name, age, error, build_brace(outlier_str, z_str)))
            } else if (type == "Boundary") {
                return(sprintf('Boundary("%s")%s;', name, build_brace("", z_str)))
            } else if (type %in% c("Date", "Tephra")) {
                # Generic Date/Tephra - default to calBP if age present
                if (is.na(age) || as.character(age) == "") {
                    return(sprintf('Date("%s")%s;', name, build_brace("", z_str)))
                } else {
                    dist_str <- format_dist("calBP", age, error, unc_type)
                    return(sprintf('Date("%s", %s)%s;', name, dist_str, build_brace(outlier_str, z_str)))
                }
            } else if (type %in% c("Date_calBP", "Tephra_calBP")) {
                dist_str <- format_dist("calBP", age, error, unc_type)
                return(sprintf('Date("%s", %s)%s;', name, dist_str, build_brace(outlier_str, z_str)))
            } else if (type %in% c("Date_CE", "Tephra_CE", "CE")) {
                dist_str <- format_dist("CE", age, error, unc_type)
                return(sprintf('Date("%s", %s)%s;', name, dist_str, build_brace(outlier_str, z_str)))
            }
            return("")
        }

        # --- Age-Depth Model Generation ---

        generated_age_depth_code <- reactiveVal("")

        observeEvent(
            list(
                input$core_name, input$curve_name, input$use_bomb, input$bomb_curve,
                input$delta_r_value, input$delta_r_error, input$top_depth,
                input$bottom_depth, input$depth_unit, input$interp_rate,
                input$auto_rename_ad
            ),
            generated_age_depth_code(""),
            ignoreInit = TRUE
        )

        observeEvent(input$gen_age_depth, {
            generated_age_depth_code("")
            if (is.null(rv_age_depth_data())) {
                showNotification("Upload age-depth data before generating OxCal code.", type = "error")
                return()
            }

            tryCatch(
                {
                    df <- rv_age_depth_data()
                    validation <- chronology_validate_table(
                        df,
                        model = "age_depth",
                        allow_duplicates = isTRUE(input$auto_rename_ad)
                    )
                    if (!validation$ok) stop(paste(validation$errors, collapse = "\n"))
                    df <- validation$data

                    # Sort by Depth descending (Bottom to Top)
                    df <- df[order(df$Depth, decreasing = TRUE), ]

                    # Duplicate Name Handling
                    if ("Name" %in% names(df)) {
                        dups <- duplicated(df$Name) | duplicated(df$Name, fromLast = TRUE)
                        if (any(dups)) {
                            if (input$auto_rename_ad) {
                                # Rename duplicates with _01, _02 etc.
                                # Since we are ordered by Depth decreasing (Bottom to Top),
                                # we can just iterate and append count
                                counts <- table(df$Name)
                                dup_names <- names(counts)[counts > 1]

                                for (dn in dup_names) {
                                    idx <- which(df$Name == dn)
                                    # idx is ordered from bottom to top (deepest first)
                                    for (i in seq_along(idx)) {
                                        df$Name[idx[i]] <- paste0(dn, sprintf("_%02d", i))
                                    }
                                }
                                showNotification("Duplicate names auto-renamed.", type = "message")
                            } else {
                                showNotification("Warning: Duplicate names detected! OxCal requires unique names.", type = "warning")
                            }
                        }
                    }

                    core_name <- input$core_name

                    code <- c()
                    
                    # Check if F14C dates are used or if any uncertainty is < 5
                    needs_resolution <- FALSE
                    if ("Type" %in% names(df) && "Uncertainty" %in% names(df)) {
                        has_f14c <- any(!is.na(df$Type) & df$Type == "R_F14C")
                        has_low_unc <- any(!is.na(df$Uncertainty) & df$Uncertainty < 5)
                        needs_resolution <- has_f14c || has_low_unc
                    }
                    
                    # Add Options() BEFORE Plot() if needed
                    if (needs_resolution) {
                        code <- c(code, "Options(){Resolution=1;};")
                    }
                    
                    code <- c(code, "Plot() {")

                    # Generate Outlier_Model declarations from data
                    outlier_models_code <- get_outlier_models_from_data(df)
                    code <- c(code, outlier_models_code)

                    # Determine k0 based on unit
                    k0 <- if (input$depth_unit == "m") 100 else 1
                    checked_interp_rate <- chronology_validate_interp_rate(input$interp_rate)
                    if (!checked_interp_rate$ok) stop(checked_interp_rate$error)
                    interp_rate <- checked_interp_rate$value

                    code <- c(code, sprintf('  P_Sequence("%s",%s,%s,U(-2,2)) {', core_name, k0, interp_rate))

                    # Bottom Boundary
                    if (!is.na(input$bottom_depth)) {
                        code <- c(code, sprintf('   Boundary("%s bottom") { z=%s; };', core_name, input$bottom_depth))
                    } else {
                        code <- c(code, sprintf('   Boundary("%s bottom");', core_name))
                    }

                    # Validate for NA/empty fields that would produce "NA" in OxCal code
                    na_warnings <- c()
                    
                    # Types that require Age and Uncertainty
                    types_needing_age <- c("R_Date", "R_F14C", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE", "CE")
                    
                    # Name and Type are always required
                    for (col_label in c("Name", "Type")) {
                      if (col_label %in% names(df)) {
                        na_rows <- which(is.na(df[[col_label]]) | trimws(as.character(df[[col_label]])) == "")
                        if (length(na_rows) > 0) {
                          row_names <- if ("Name" %in% names(df)) {
                            paste(na.omit(df$Name[na_rows][1:min(3, length(na_rows))]), collapse = ", ")
                          } else {
                            paste("row(s)", paste(na_rows[1:min(3, length(na_rows))], collapse = ", "))
                          }
                          na_warnings <- c(na_warnings, sprintf("'%s' is empty for %s%s", 
                            col_label, row_names, if (length(na_rows) > 3) sprintf(" (+%d more)", length(na_rows) - 3) else ""))
                        }
                      }
                    }
                    
                    # Age and Uncertainty only required for types that need them
                    for (col_label in c("Age", "Uncertainty")) {
                      if (col_label %in% names(df)) {
                        # Only check rows where Type requires age/uncertainty
                        check_rows <- which(
                          (is.na(df[[col_label]]) | trimws(as.character(df[[col_label]])) == "") &
                          (!is.na(df$Type) & trimws(as.character(df$Type)) %in% types_needing_age)
                        )
                        if (length(check_rows) > 0) {
                          row_names <- if ("Name" %in% names(df)) {
                            paste(na.omit(df$Name[check_rows][1:min(3, length(check_rows))]), collapse = ", ")
                          } else {
                            paste("row(s)", paste(check_rows[1:min(3, length(check_rows))], collapse = ", "))
                          }
                          na_warnings <- c(na_warnings, sprintf("'%s' is empty for %s%s", 
                            col_label, row_names, if (length(check_rows) > 3) sprintf(" (+%d more)", length(check_rows) - 3) else ""))
                        }
                      }
                    }
                    
                    if (length(na_warnings) > 0) {
                      showNotification(
                        HTML(paste0(
                          "<b>Warning: Empty fields detected!</b> These will appear as 'NA' in the OxCal code:<br>",
                          paste("• ", na_warnings, collapse = "<br>")
                        )),
                        type = "warning", duration = 15
                      )
                    }

                    # Apply curve logic and dates, with Combine blocks for same-depth dates
                    # Group rows by depth to detect same-depth dates
                    combined_depths <- c()
                    if ("Depth" %in% names(df)) {
                      depth_counts <- table(df$Depth[!is.na(df$Depth)])
                      dup_depths <- as.numeric(names(depth_counts[depth_counts > 1]))
                      
                      if (length(dup_depths) > 0) {
                        combined_depths <- dup_depths
                        # Notify user about same-depth combining
                        depth_list <- paste(dup_depths, collapse = ", ")
                        showNotification(
                          HTML(paste0(
                            "<b>Same-depth dates detected</b> at depth(s): ", depth_list, ".<br><br>",
                            "These are automatically wrapped in <code>Combine</code> blocks, treating them as ",
                            "independent chronological information for the same context.<br><br>",
                            "<b>Note:</b> This may not be appropriate if you dated the same object multiple times. ",
                            "In that case, <code>R_Combine</code> may be more suitable. ",
                            "See the <a href='https://c14.arch.ox.ac.uk/oxcalhelp/hlp_commands.html' target='_blank'>OxCal command reference</a>."
                          )),
                          type = "warning", duration = 30
                        )
                      }
                    }
                    
                    # Process rows with depth grouping
                    current_curve <- NA
                    defined_curves <- c()
                    processed_depths <- c()
                    
                    # Emit initial curve definition
                    primary_curve <- input$curve_name
                    if (is.null(primary_curve) || is.na(primary_curve)) primary_curve <- "IntCal20"
                    curve_files <- list("IntCal20" = "intcal20.14c", "Marine20" = "marine20.14c", "SHCal20" = "shcal20.14c")
                    if (primary_curve %in% names(curve_files)) {
                      code <- c(code, sprintf('   Curve("%s", "%s");', primary_curve, curve_files[[primary_curve]]))
                      defined_curves <- c(defined_curves, primary_curve)
                    }
                    
                    # Insert Delta_R for Marine20
                    if (primary_curve == "Marine20") {
                      code <- c(code, chronology_format_delta_r(
                        primary_curve,
                        input$delta_r_value,
                        input$delta_r_error,
                        indent = "   "
                      ))
                    }
                    
                    current_curve <- primary_curve
                    
                    for (i in seq_len(nrow(df))) {
                      row <- df[i, ]
                      row_depth <- if ("Depth" %in% names(row)) row$Depth else NA
                      
                      # Skip if this depth was already processed as part of a Combine group
                      if (!is.na(row_depth) && row_depth %in% processed_depths) next
                      
                      if (!is.na(row_depth) && row_depth %in% combined_depths) {
                        # Same-depth group — wrap in Combine block
                        group_rows <- df[!is.na(df$Depth) & df$Depth == row_depth, ]
                        
                        code <- c(code, sprintf('   Combine("%s depth") {', row_depth))
                        
                        # Generate each date in the group WITHOUT z= (z goes on Combine)
                        for (j in seq_len(nrow(group_rows))) {
                          grow <- group_rows[j, ]
                          # Temporarily remove Depth so format_date_generic doesn't add z=
                          grow$Depth <- NA
                          
                          # Handle curve switching for this row
                          type_val <- NA
                          if ("Type" %in% names(grow)) type_val <- grow$Type
                          age_val <- NA
                          if ("Age" %in% names(grow)) age_val <- as.numeric(grow$Age)
                          
                          target_curve <- chronology_select_curve(
                            type = type_val,
                            age = age_val,
                            primary_curve = input$curve_name,
                            use_bomb = input$use_bomb,
                            bomb_curve_name = input$bomb_curve
                          )
                          
                          # Emit curve switch if needed
                          if (!is.na(target_curve) && target_curve != current_curve) {
                            bomb_curve_files <- list(
                              "Bomb21NH1" = "bomb21nh1.14c", "Bomb21NH2" = "bomb21nh2.14c",
                              "Bomb21NH3" = "bomb21nh3.14c", "Bomb21SH12" = "bomb21sh12.14c",
                              "Bomb21SH3" = "bomb21sh3.14c"
                            )
                            if (target_curve %in% names(bomb_curve_files) && !(target_curve %in% defined_curves)) {
                              code <- c(code, sprintf('    Curve("%s", "%s");', target_curve, bomb_curve_files[[target_curve]]))
                              defined_curves <- c(defined_curves, target_curve)
                            } else {
                              code <- c(code, sprintf('    Curve("=%s");', target_curve))
                            }
                            current_curve <- target_curve
                          }
                          
                          line <- format_date_generic(grow)
                          code <- c(code, paste0("    ", line))
                        }
                        
                        code <- c(code, sprintf('    z=%s;', row_depth))
                        code <- c(code, '   };')
                        
                        processed_depths <- c(processed_depths, row_depth)
                      } else {
                        # Single date at this depth — use apply_curve_logic for one row
                        single_df <- df[i, , drop = FALSE]
                        res <- apply_curve_logic(single_df, input$curve_name, input$use_bomb, input$bomb_curve,
                                                 current_curve, defined_curves)
                        code <- c(code, res$code)
                        current_curve <- res$current_curve
                        defined_curves <- res$defined_curves
                      }
                    }

                    # Top Boundary
                    if (!is.na(input$top_depth)) {
                        code <- c(code, sprintf('   Boundary("%s top") { z=%s; };', core_name, input$top_depth))
                    } else {
                        code <- c(code, sprintf('   Boundary("%s top");', core_name))
                    }
                    code <- c(code, "  };")
                    code <- c(code, "};")

                    code_check <- chronology_validate_oxcal_code(code)
                    if (!code_check$ok) stop(code_check$error)
                    generated_age_depth_code(paste(code, collapse = "\n"))
                },
                error = function(e) {
                    showNotification(paste("Error generating code:", e$message), type = "error")
                }
            )
        })

        output$out_age_depth <- renderText({
            generated_age_depth_code()
        })

        observeEvent(input$copy_age_depth_status, {
            status <- input$copy_age_depth_status
            if (isTRUE(status$ok)) {
                showNotification("Code copied to clipboard!", type = "message", duration = 2)
            } else if (identical(status$reason, "empty")) {
                showNotification("No valid age-depth code is available to copy.", type = "warning")
            } else {
                showNotification(
                    "Browser clipboard access was blocked. Select the generated code and copy it manually.",
                    type = "error",
                    duration = 6
                )
            }
        })

        # --- Phase Model Generation ---

        # Reactive value to hold Phase Model data
        rv_phase_data <- reactiveVal(NULL)

        # Helper to normalize column names
        normalize_phase_columns <- function(df) {
            # Standardize column names to: Tephra.Name, Phase, Name, Date.Type, Age.Mean, Age.SD, Outlier.Type, Outlier.Prob

            rename_col <- function(df, patterns, new_name) {
                # Find first matching column
                match_idx <- which(tolower(names(df)) %in% tolower(patterns))
                if (length(match_idx) > 0) {
                    names(df)[match_idx[1]] <- new_name
                }
                return(df)
            }

            df <- rename_col(df, c("Tephra Name", "Tephra_Name", "Tephra.Name"), "Tephra.Name")
            df <- rename_col(df, c("Phase"), "Phase")
            df <- rename_col(df, c("Name", "Sample Name", "Sample ID", "Lab Code"), "Name")
            df <- rename_col(df, c("Date Type", "Date_Type", "Date.Type"), "Date.Type")
            df <- rename_col(df, c("Age Mean", "Age_Mean", "Age.Mean", "Age"), "Age.Mean")
            df <- rename_col(df, c("Age SD", "Age_SD", "Age.SD", "Uncertainty", "Error", "sd"), "Age.SD")
            df <- rename_col(df, c("Outlier type", "Outlier Type", "Outlier.Type", "Outlier_type", "Outlier_Type"), "Outlier.Type")
            df <- rename_col(df, c("Outlier prob", "Outlier Prob", "Outlier.Prob", "Outlier_Prob", "Outlier_prob"), "Outlier.Prob")

            return(df)
        }

        # --- Sheet Selection UI for Phase Model ---
        output$phase_sheet_select <- renderUI({
            req(input$file_phase)
            file_ext <- tolower(tools::file_ext(input$file_phase$name))
            if (file_ext == "xlsx") {
                sheets <- readxl::excel_sheets(input$file_phase$datapath)
                selectInput(ns("phase_sheet"), "Select Sheet:", choices = sheets, selected = sheets[1])
            }
        })

        # Observer for Phase Model file upload
        observeEvent(input$file_phase, {
            req(input$file_phase)
            generated_phase_code("")
            rv_phase_data(NULL)
            file_ext <- tolower(tools::file_ext(input$file_phase$name))
            if (file_ext == "xlsx") {
                # Force read of first sheet immediately to ensure data is loaded
                # even if sheet selection doesn't trigger a change
                tryCatch(
                    {
                        sheets <- readxl::excel_sheets(input$file_phase$datapath)
                        if (length(sheets) > 0) {
                            df <- readxl::read_excel(input$file_phase$datapath, sheet = sheets[1])
                            df <- normalize_phase_columns(df)
                            rv_phase_data(df)
                            showNotification(paste("Loaded sheet:", sheets[1]), type = "message")
                        }
                    },
                    error = function(e) {
                        showNotification(paste("Error reading Excel file:", e$message), type = "error")
                    }
                )
                return()
            }
            tryCatch(
                {
                    df <- read.csv(input$file_phase$datapath, stringsAsFactors = FALSE)
                    df <- normalize_phase_columns(df)
                    rv_phase_data(df)
                    showNotification("Data loaded successfully", type = "message")
                },
                error = function(e) {
                    showNotification(paste("Error reading file:", e$message), type = "error")
                }
            )
        })

        observeEvent(input$phase_sheet, {
            req(input$file_phase, input$phase_sheet)
            generated_phase_code("")
            rv_phase_data(NULL)
            tryCatch(
                {
                    df <- readxl::read_excel(input$file_phase$datapath, sheet = input$phase_sheet)
                    df <- normalize_phase_columns(df)
                    rv_phase_data(df)
                    showNotification("Data loaded successfully", type = "message")
                },
                error = function(e) {
                    showNotification(paste("Error reading file:", e$message), type = "error")
                }
            )
        })

        generated_phase_code <- reactiveVal("")

        observeEvent(
            list(
                input$phase_curve_name, input$phase_use_bomb, input$phase_bomb_curve,
                input$phase_delta_r_value, input$phase_delta_r_error,
                input$auto_rename_phase, input$use_tau_before, input$use_tau_after
            ),
            generated_phase_code(""),
            ignoreInit = TRUE
        )

        observeEvent(input$gen_phase, {
            generated_phase_code("")
            if (is.null(rv_phase_data())) {
                showNotification("Upload phase-model data before generating OxCal code.", type = "error")
                return()
            }
            tryCatch(
                {
                    df <- rv_phase_data()
                    validation <- chronology_validate_table(
                        df,
                        model = "phase",
                        allow_duplicates = isTRUE(input$auto_rename_phase)
                    )
                    if (!validation$ok) stop(paste(validation$errors, collapse = "\n"))
                    df <- validation$data

                    # Duplicate Name Handling
                    if ("Name" %in% names(df)) {
                        dups <- duplicated(df$Name) | duplicated(df$Name, fromLast = TRUE)
                        if (any(dups)) {
                            if (input$auto_rename_phase) {
                                # Rename duplicates
                                counts <- table(df$Name)
                                dup_names <- names(counts)[counts > 1]
                                for (dn in dup_names) {
                                    idx <- which(df$Name == dn)
                                    for (i in seq_along(idx)) {
                                        df$Name[idx[i]] <- paste0(dn, sprintf("_%02d", i))
                                    }
                                }
                                showNotification("Duplicate names auto-renamed.", type = "message")
                            } else {
                                showNotification("Warning: Duplicate names detected! OxCal requires unique names.", type = "warning")
                            }
                        }
                    }

                    code <- c()
                    
                    # Validate for NA/empty fields that would produce "NA" in OxCal code
                    na_warnings <- c()
                    
                    # Types that require Age and SD
                    types_needing_age <- c("R_Date", "R_F14C", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE", "CE")
                    
                    # Always-required columns
                    always_required <- list("Name" = "Name", "Tephra Name" = "Tephra.Name", 
                                            "Phase" = "Phase", "Date Type" = "Date.Type")
                    for (col_label in names(always_required)) {
                      col_name <- always_required[[col_label]]
                      if (col_name %in% names(df)) {
                        na_rows <- which(is.na(df[[col_name]]) | trimws(as.character(df[[col_name]])) == "")
                        if (length(na_rows) > 0) {
                          row_ids <- if ("Name" %in% names(df)) {
                            paste(na.omit(df$Name[na_rows][1:min(3, length(na_rows))]), collapse = ", ")
                          } else {
                            paste("row(s)", paste(na_rows[1:min(3, length(na_rows))], collapse = ", "))
                          }
                          na_warnings <- c(na_warnings, sprintf("'%s' is empty for %s%s", 
                            col_label, row_ids, if (length(na_rows) > 3) sprintf(" (+%d more)", length(na_rows) - 3) else ""))
                        }
                      }
                    }
                    
                    # Age and SD only required for types that need them
                    age_sd_cols <- list("Age" = "Age.Mean", "SD" = "Age.SD")
                    type_col <- if ("Date.Type" %in% names(df)) "Date.Type" else NULL
                    for (col_label in names(age_sd_cols)) {
                      col_name <- age_sd_cols[[col_label]]
                      if (col_name %in% names(df) && !is.null(type_col)) {
                        check_rows <- which(
                          (is.na(df[[col_name]]) | trimws(as.character(df[[col_name]])) == "") &
                          (!is.na(df[[type_col]]) & trimws(as.character(df[[type_col]])) %in% types_needing_age)
                        )
                        if (length(check_rows) > 0) {
                          row_ids <- if ("Name" %in% names(df)) {
                            paste(na.omit(df$Name[check_rows][1:min(3, length(check_rows))]), collapse = ", ")
                          } else {
                            paste("row(s)", paste(check_rows[1:min(3, length(check_rows))], collapse = ", "))
                          }
                          na_warnings <- c(na_warnings, sprintf("'%s' is empty for %s%s", 
                            col_label, row_ids, if (length(check_rows) > 3) sprintf(" (+%d more)", length(check_rows) - 3) else ""))
                        }
                      }
                    }
                    
                    if (length(na_warnings) > 0) {
                      showNotification(
                        HTML(paste0(
                          "<b>Warning: Empty fields detected!</b> These will appear as 'NA' in the OxCal code:<br>",
                          paste("&bull; ", na_warnings, collapse = "<br>")
                        )),
                        type = "warning", duration = 15
                      )
                    }

                    # Check if F14C dates are used or if any uncertainty is < 5
                    needs_resolution <- FALSE
                    if ("Date.Type" %in% names(df) && "Age.SD" %in% names(df)) {
                        has_f14c <- any(!is.na(df$Date.Type) & df$Date.Type == "R_F14C")
                        has_low_unc <- any(!is.na(df$Age.SD) & df$Age.SD < 5)
                        needs_resolution <- has_f14c || has_low_unc
                    }
                    
                    # Add Options() BEFORE Plot() if needed
                    if (needs_resolution) {
                        code <- c(code, "Options(){Resolution=1;};")
                    }
                    
                    code <- c(code, "Plot() {")

                    # Generate Outlier_Model declarations from data
                    outlier_models_code <- get_outlier_models_from_data(df)
                    code <- c(code, outlier_models_code)

                    code <- c(code, '  Sequence("") {')
                    
                    # Use Tau_Boundary or Boundary for Before phase based on checkbox
                    start_boundary <- if (input$use_tau_before) 'Tau_Boundary("Start");' else 'Boundary("Start");'
                    code <- c(code, paste0('    ', start_boundary))

                    if ("Phase" %in% names(df)) {
                        df <- df[!is.na(df$Phase) & trimws(df$Phase) != "", ]
                    } else {
                        showNotification("Warning: 'Phase' column not found!", type = "warning")
                    }

                    # We need to know the Tephra Name for the labels
                    tephra_name <- "Tephra"
                    if ("Tephra.Name" %in% names(df)) {
                        t_names <- unique(df$Tephra.Name)
                        t_names <- t_names[!is.na(t_names) & t_names != ""]
                        if (length(t_names) > 0) tephra_name <- t_names[1]
                    }

                    # Check if "During" phase has any data
                    has_during_data <- FALSE
                    if ("Phase" %in% names(df)) {
                        has_during_data <- any(df$Phase == "During")
                    }

                    # Initialize curve state
                    current_curve <- NA
                    defined_curves <- c()
                    
                    # Emit initial curve + Delta_R for Marine20 before phases
                    phase_primary <- input$phase_curve_name
                    if (is.null(phase_primary) || is.na(phase_primary)) phase_primary <- "IntCal20"
                    phase_curve_files <- list("IntCal20" = "intcal20.14c", "Marine20" = "marine20.14c", "SHCal20" = "shcal20.14c")
                    if (phase_primary %in% names(phase_curve_files)) {
                        code <- c(code, sprintf('    Curve("%s", "%s");', phase_primary, phase_curve_files[[phase_primary]]))
                        defined_curves <- c(defined_curves, phase_primary)
                        current_curve <- phase_primary
                    }
                    if (phase_primary == "Marine20") {
                        code <- c(code, chronology_format_delta_r(
                            phase_primary,
                            input$phase_delta_r_value,
                            input$phase_delta_r_error,
                            indent = "    "
                        ))
                    }

                    # --- Phase: Before ---
                    code <- c(code, sprintf('    Phase("Before %s") {', tephra_name))
                    if ("Phase" %in% names(df)) {
                        phase_data <- df[df$Phase == "Before", ]
                        if (nrow(phase_data) > 0) {
                            res <- apply_curve_logic(phase_data, input$phase_curve_name, input$phase_use_bomb, input$phase_bomb_curve, current_curve, defined_curves)
                            code <- c(code, paste0("   ", res$code))
                            current_curve <- res$current_curve
                            defined_curves <- res$defined_curves
                        }
                    }
                    code <- c(code, "    };")

                    # --- Phase: During / Boundary ---
                    if (has_during_data) {
                        # Scenario A: During phase has data
                        code <- c(code, '    Boundary("Tephra Start");')
                        code <- c(code, sprintf('    Phase("During %s") {', tephra_name))

                        phase_data <- df[df$Phase == "During", ]
                        if (nrow(phase_data) > 0) {
                            res <- apply_curve_logic(phase_data, input$phase_curve_name, input$phase_use_bomb, input$phase_bomb_curve, current_curve, defined_curves)
                            code <- c(code, paste0("   ", res$code))
                            current_curve <- res$current_curve
                            defined_curves <- res$defined_curves
                        }

                        # Automatically add the Tephra date placeholder
                        code <- c(code, sprintf('       Date("%s");', tephra_name))

                        code <- c(code, "    };")
                        code <- c(code, '    Boundary("Tephra End");')
                    } else {
                        # Scenario B: No During phase data -> Single Boundary
                        code <- c(code, sprintf('    Boundary("%s");', tephra_name))
                    }

                    # --- Phase: After ---
                    code <- c(code, sprintf('    Phase("After %s") {', tephra_name))
                    if ("Phase" %in% names(df)) {
                        phase_data <- df[df$Phase == "After", ]
                        if (nrow(phase_data) > 0) {
                            res <- apply_curve_logic(phase_data, input$phase_curve_name, input$phase_use_bomb, input$phase_bomb_curve, current_curve, defined_curves)
                            code <- c(code, paste0("   ", res$code))
                            current_curve <- res$current_curve
                            defined_curves <- res$defined_curves
                        }
                    }
                    code <- c(code, "    };")

                    # Use Tau_Boundary or Boundary for After phase based on checkbox
                    end_boundary <- if (input$use_tau_after) 'Tau_Boundary("End");' else 'Boundary("End");'
                    code <- c(code, paste0('    ', end_boundary))
                    code <- c(code, "  };")
                    code <- c(code, "};")

                    code_check <- chronology_validate_oxcal_code(code)
                    if (!code_check$ok) stop(code_check$error)
                    generated_phase_code(paste(code, collapse = "\n"))
                },
                error = function(e) {
                    showNotification(paste("Error generating phase code:", e$message), type = "error")
                }
            )
        })

        output$out_phase <- renderText({
            generated_phase_code()
        })

        observeEvent(input$copy_phase_status, {
            status <- input$copy_phase_status
            if (isTRUE(status$ok)) {
                showNotification("Code copied to clipboard!", type = "message", duration = 2)
            } else if (identical(status$reason, "empty")) {
                showNotification("No valid phase-model code is available to copy.", type = "warning")
            } else {
                showNotification(
                    "Browser clipboard access was blocked. Select the generated code and copy it manually.",
                    type = "error",
                    duration = 6
                )
            }
        })

        # --- Model Linking ---

        linked_code <- reactiveVal("")

        observeEvent(input$link_models, {
            req(input$model_a_code, input$model_b_code)

            # Robust parser for OxCal code
            parse_oxcal <- function(code) {
                # Helper to find balanced blocks
                # returns list(remainder, extracted_blocks)
                extract_blocks <- function(txt, start_pattern, open_char, close_char) {
                    blocks <- c()
                    repeat {
                        # Find start
                        m <- regexpr(start_pattern, txt, perl = TRUE)
                        if (m == -1) break
                        
                        start_idx <- as.integer(m)
                        match_len <- attr(m, "match.length")
                        
                        # Scan for balanced close_char
                        # Start scanning after the match
                        scan_pos <- start_idx + match_len
                        chars <- strsplit(substring(txt, scan_pos), "")[[1]]
                        
                        depth <- 1 
                        
                        found_end <- FALSE
                        end_rel_idx <- 0
                        
                        for (i in seq_along(chars)) {
                            if (chars[i] == open_char) {
                                depth <- depth + 1
                            } else if (chars[i] == close_char) {
                                depth <- depth - 1
                            }
                            
                            if (depth == 0) {
                                end_rel_idx <- i
                                found_end <- TRUE
                                break
                            }
                        }
                        
                        if (found_end) {
                            end_idx <- scan_pos + end_rel_idx - 1
                            
                            # Check for trailing semicolon
                            full_end_idx <- end_idx
                            remainder_check <- substring(txt, end_idx + 1)
                            if (grepl("^\\s*;", remainder_check)) {
                                semi_match <- regexpr("^\\s*;", remainder_check)
                                full_end_idx <- end_idx + attr(semi_match, "match.length")
                            }
                            
                            # Extract block
                            block <- substring(txt, start_idx, full_end_idx)
                            blocks <- c(blocks, block)
                            
                            # Remove from txt
                            txt <- paste0(substring(txt, 1, start_idx - 1), substring(txt, full_end_idx + 1))
                        } else {
                            # Unbalanced or malformed, stop to avoid infinite loop
                            break
                        }
                    }
                    list(txt = txt, blocks = blocks)
                }
                
                # 1. Extract Options
                # Pattern includes the opening {
                opt_res <- extract_blocks(code, "Options\\s*\\(\\s*\\)\\s*\\{", "{", "}")
                code <- opt_res$txt
                options_blocks <- opt_res$blocks
                
                # 2. Extract Outlier_Model
                # Pattern includes the opening (
                out_res <- extract_blocks(code, "Outlier_Model\\s*\\(", "(", ")")
                code <- out_res$txt
                outlier_models <- out_res$blocks
                
                # 3. Handle Plot() wrappers
                # We want to strip the wrapper but keep the content.
                # We can use extract_blocks to find them, then strip the wrapper string from the block.
                plot_res <- extract_blocks(code, "Plot\\s*\\(\\s*\\)\\s*\\{", "{", "}")
                
                # The remainder in plot_res$txt is code that was NOT in a Plot() block.
                # The blocks in plot_res$blocks are the Plot() blocks.
                
                # Strip wrappers from Plot blocks
                stripped_content <- c()
                for (p_block in plot_res$blocks) {
                    # Remove "Plot(){" prefix (flexible whitespace)
                    # We know it matches the pattern, but let's be precise
                    # Find the first {
                    first_brace <- regexpr("{", p_block, fixed = TRUE)
                    
                    # Content starts after first {
                    content_start <- first_brace + 1
                    
                    # Content ends before last }
                    # Find last }
                    # Since we extracted it carefully, it should be near the end.
                    # Check if it ends with ;
                    trimmed_block <- trimws(p_block, which = "right")
                    if (endsWith(trimmed_block, ";")) {
                        trimmed_block <- sub(";$", "", trimmed_block)
                    }
                    trimmed_block <- trimws(trimmed_block, which = "right")
                    # Now it should end with }
                    if (endsWith(trimmed_block, "}")) {
                         # We need to extract the content inside the braces
                         # We can use substring on the original p_block
                         # But we need indices.
                         # Let's just use regex to strip start and end
                         inner <- sub("^Plot\\s*\\(\\s*\\)\\s*\\{", "", p_block)
                         inner <- sub("\\}\\s*;?\\s*$", "", inner)
                         stripped_content <- c(stripped_content, inner)
                    } else {
                        # Fallback
                        stripped_content <- c(stripped_content, p_block)
                    }
                }
                
                # Combine non-plot code and stripped plot content
                # Usually non-plot code is just whitespace if everything is in Plot
                final_content <- paste(c(plot_res$txt, stripped_content), collapse = "\n")
                
                # Clean up excessive whitespace/newlines
                final_content <- gsub("\n{3,}", "\n\n", final_content)
                final_content <- trimws(final_content)
                
                list(
                    content = final_content,
                    options = options_blocks,
                    outliers = outlier_models
                )
            }
            
            # Function to deduplicate Curve calls - first occurrence gets full definition,
            # subsequent occurrences use the equals notation
            deduplicate_curves <- function(txt) {
                # Find all Curve("Name", "File") patterns
                # We need to be careful not to match Curve("=Name")
                # Regex: Curve\s*\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)
                
                pattern <- 'Curve\\s*\\(\\s*"([^"]+)"\\s*,\\s*"([^"]+)"\\s*\\)'
                
                # We iterate through matches
                matches <- gregexpr(pattern, txt, perl = TRUE)
                if (matches[[1]][1] == -1) return(txt)
                
                match_data <- regmatches(txt, matches)[[1]]
                match_indices <- as.integer(matches[[1]])
                match_lengths <- attr(matches[[1]], "match.length")
                
                seen_curves <- c()
                    
                # Let's identify which matches are duplicates
                replacements <- list() # list of (start, end, replacement_string)
                
                for (i in seq_along(match_data)) {
                    m_text <- match_data[i]
                    # Extract name
                    # We can use sub to get the name
                    # Curve("Name", "File") -> Name
                    curve_name <- sub('^Curve\\s*\\(\\s*"([^"]+)".*$', "\\1", m_text)
                    
                    if (curve_name %in% seen_curves) {
                        # Duplicate!
                        replacements[[length(replacements) + 1]] <- list(
                            start = match_indices[i],
                            end = match_indices[i] + match_lengths[i] - 1,
                            new_text = paste0('Curve("=', curve_name, '")')
                        )
                    } else {
                        seen_curves <- c(seen_curves, curve_name)
                    }
                }
                
                # Apply replacements from back to front
                if (length(replacements) > 0) {
                    for (i in rev(seq_along(replacements))) {
                        rep <- replacements[[i]]
                        txt <- paste0(
                            substring(txt, 1, rep$start - 1),
                            rep$new_text,
                            substring(txt, rep$end + 1)
                        )
                    }
                }
                
                return(txt)
            }
            
            # Function to deduplicate Outlier_Model definitions
            deduplicate_outlier_models <- function(models) {
                seen <- list()
                unique_models <- c()
                
                for (model in models) {
                    # Extract model name (first quoted string)
                    name_match <- regmatches(model, regexec('"([^"]+)"', model, perl = TRUE))[[1]]
                    if (length(name_match) >= 2) {
                        model_name <- name_match[2]
                        if (!model_name %in% names(seen)) {
                            seen[[model_name]] <- TRUE
                            # Ensure model ends with semicolon
                            model <- trimws(model)
                            if (!grepl(";$", model)) {
                                model <- paste0(model, ";")
                            }
                            unique_models <- c(unique_models, model)
                        }
                    }
                }
                return(unique_models)
            }

            result_a <- parse_oxcal(input$model_a_code)
            result_b <- parse_oxcal(input$model_b_code)
            
            code_a <- result_a$content
            code_b <- result_b$content
            
            # Collect and deduplicate Options() blocks
            all_options <- c(result_a$options, result_b$options)
            
            # Collect and deduplicate Outlier_Model definitions
            all_outlier_models <- c(result_a$outliers, result_b$outliers)
            unique_outlier_models <- deduplicate_outlier_models(all_outlier_models)

            ties_txt <- input$tie_points

            # Collect all tie point names
            tie_names <- c()

            if (nzchar(ties_txt)) {
                lines <- strsplit(ties_txt, "\n")[[1]]
                for (line in lines) {
                    # Skip empty lines
                    if (nzchar(trimws(line))) {
                        # Try tab-separated first (Excel default), then comma-separated
                        parts <- strsplit(line, "\t")[[1]]
                        if (length(parts) < 3) {
                            # Fall back to comma-separated
                            parts <- strsplit(line, ",")[[1]]
                        }

                        if (length(parts) >= 3) {
                            tie_name <- trimws(parts[1])
                            syn_a <- trimws(parts[2])
                            syn_b <- trimws(parts[3])

                            # Store tie name
                            if (!tie_name %in% tie_names) {
                                tie_names <- c(tie_names, tie_name)
                            }

                            # Replace "SynA" with "=TieName"
                            code_a <- gsub(paste0('"', syn_a, '"'), paste0('"=', tie_name, '"'), code_a, fixed = TRUE)
                            code_b <- gsub(paste0('"', syn_b, '"'), paste0('"=', tie_name, '"'), code_b, fixed = TRUE)
                        }
                    }
                }
            }

            # Build base tie point definitions
            tie_definitions <- ""
            if (length(tie_names) > 0) {
                # Reverse order so they appear in descending order (as in the example)
                tie_defs <- sapply(rev(tie_names), function(tn) {
                    paste0('  Date("', tn, '");')
                })
                tie_definitions <- paste0(paste(tie_defs, collapse = "\n"), "\n")
            }
            
            # Combine the model code
            combined_code <- paste0(code_a, "\n", code_b)
            
            # Deduplicate curves across the combined code
            combined_code <- deduplicate_curves(combined_code)
            
            # Build Options block (merge and deduplicate settings)
            options_block <- ""
            if (length(all_options) > 0) {
                # Extract individual option settings and deduplicate
                option_settings <- c()
                for (opt_block in all_options) {
                    # Extract content between Options(){ and }
                    inner <- gsub("^\\s*Options\\s*\\(\\s*\\)\\s*\\{\\s*", "", opt_block)
                    inner <- gsub("\\s*\\}\\s*;?\\s*$", "", inner)
                    # Split by semicolons or newlines to get individual settings
                    settings <- strsplit(inner, "[;\n]+")[[1]]
                    settings <- trimws(settings)
                    settings <- settings[nzchar(settings)]
                    option_settings <- c(option_settings, settings)
                }
                # Deduplicate based on setting name (before =)
                unique_settings <- list()
                for (setting in option_settings) {
                    setting_name <- trimws(strsplit(setting, "=")[[1]][1])
                    unique_settings[[setting_name]] <- setting
                }
                if (length(unique_settings) > 0) {
                    options_block <- paste0(
                        "Options()\n{\n  ",
                        paste(unlist(unique_settings), collapse = ";\n  "),
                        ";\n};\n"
                    )
                }
            }
            
            # Build Outlier_Model block
            outlier_block <- ""
            if (length(unique_outlier_models) > 0) {
                outlier_block <- paste0(paste(unique_outlier_models, collapse = "\n"), "\n")
            }

            final <- paste0(
                options_block,
                "Plot()\n",
                "{\n",
                outlier_block,
                tie_definitions,
                combined_code, "\n",
                "};"
            )
            
            # Clean up excessive blank lines in final output
            final <- gsub("\n{3,}", "\n\n", final)

            linked_code(final)
        })

        output$out_linked <- renderText({
            linked_code()
        })

        observeEvent(input$copy_linked_status, {
            status <- input$copy_linked_status
            if (isTRUE(status$ok)) {
                showNotification("Code copied to clipboard!", type = "message", duration = 2)
            } else if (identical(status$reason, "empty")) {
                showNotification("No valid linked-model code is available to copy.", type = "warning")
            } else {
                showNotification(
                    "Browser clipboard access was blocked. Select the generated code and copy it manually.",
                    type = "error",
                    duration = 6
                )
            }
        })

        observeEvent(
            list(input$model_a_code, input$model_b_code, input$tie_points),
            linked_code(""),
            ignoreInit = TRUE
        )

        if (!is.null(global_rv)) {
            global_rv$reset_chronology_state <- function() {
                rv_age_depth_data(NULL)
                rv_age_depth_required_cols(NULL)
                generated_age_depth_code("")
                rv_phase_data(NULL)
                generated_phase_code("")
                linked_code("")
                invisible(NULL)
            }
        }
    })
}
