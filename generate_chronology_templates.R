# generate_chronology_templates.R
# Advanced template generation for Chronology module with data validation

library(openxlsx)

# ==============================================================================
# AGE-DEPTH MODEL TEMPLATE
# ==============================================================================

generate_age_depth_template <- function(output_file, verbose = FALSE) {
  if (verbose) cat("Generating Age-Depth Model Template...\n")
  
  wb <- createWorkbook()
  
  # Valid categories for dropdowns
  VALID_DATE_TYPES <- c("R_Date", "R_F14C", "Date", "Tephra", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE", "Boundary")
  VALID_UNC_TYPES <- c("N", "U")
  VALID_OUTLIER_TYPES <- c("General", "Charcoal", "SSimple", "RSimple", "TSimple", "RScaled", "")
  
  # --- SHEET 1: Data Entry ---
  addWorksheet(wb, "Data")
  
  headers <- c("Name", "Age", "Uncertainty", "Depth", "Type", "Unc_Type", "Outlier_Type")
  
  # Example data
  data_df <- data.frame(
    Name = c(
      "GCR-12345",
      "MSMB-14C",
      "Unknown-Tephra",
      "VARG-2025",
      "Dated Tephra-A",
      "Hist-Event",
      "Modern-F14C",
      "Surface"
    ),
    Age = c(8500, 6400, NA, 5000, 4255, 1912, 1.0235, 2026),
    Uncertainty = c(75, 60, NA, 50, 100, 1, 0.0025, 1),
    Depth = c(210.0, 175.0, 145.0, 125.5, 95.0, 20.0, 5.0, 0.0),
    Type = c("R_Date", "R_Date", "Tephra", "R_Date", "Tephra", "Date_CE", "R_F14C", "Date_CE"),
    Unc_Type = c("N", "N", "N", "N", "N", "N", "N", "U"),
    Outlier_Type = c("General", "Charcoal", "", "General", "General", "TSimple", "General", ""),
    stringsAsFactors = FALSE
  )
  
  writeData(wb, "Data", data_df)
  
  # Formatting
  header_style <- createStyle(
    fontSize = 11, textDecoration = "bold",
    fgFill = "#3c8dbc", fontColour = "white",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight", borderColour = "black"
  )
  addStyle(wb, "Data", header_style, rows = 1, cols = 1:7, gridExpand = TRUE)
  freezePane(wb, "Data", firstRow = TRUE)
  setColWidths(wb, "Data", cols = 1:7, widths = c(18, 12, 12, 12, 16, 12, 16))
  
  # Data validation - Type dropdown (allow custom input)
  dataValidation(wb, "Data",
    col = 5, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$A$2:$A$10",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # Data validation - Unc_Type dropdown
  dataValidation(wb, "Data",
    col = 6, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$B$2:$B$3",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # Data validation - Outlier_Type dropdown (allow custom input)
  dataValidation(wb, "Data",
    col = 7, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$C$2:$C$7",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # --- SHEET 2: Valid Categories ---
  addWorksheet(wb, "Valid_Categories")
  
  max_length <- max(length(VALID_DATE_TYPES), length(VALID_UNC_TYPES), length(VALID_OUTLIER_TYPES))
  valid_data <- data.frame(
    Type = c(VALID_DATE_TYPES, rep("", max_length - length(VALID_DATE_TYPES))),
    Unc_Type = c(VALID_UNC_TYPES, rep("", max_length - length(VALID_UNC_TYPES))),
    Outlier_Type = c(VALID_OUTLIER_TYPES, rep("", max_length - length(VALID_OUTLIER_TYPES)))
  )
  writeData(wb, "Valid_Categories", valid_data)
  
  header_style2 <- createStyle(fontSize = 11, textDecoration = "bold", fgFill = "#d4edda")
  addStyle(wb, "Valid_Categories", header_style2, rows = 1, cols = 1:3, gridExpand = TRUE)
  setColWidths(wb, "Valid_Categories", cols = 1:3, widths = 20)
  
  # --- SHEET 3: Instructions ---
  addWorksheet(wb, "Instructions")
  
  inst_text <- c(
    "AGE-DEPTH MODEL TEMPLATE",
    "",
    "═══════════════════════════════════════════════════════════════════════════",
    "COLUMN DEFINITIONS",
    "═══════════════════════════════════════════════════════════════════════════",
    "",
    "Name: Unique identifier for the sample/event.",
    "",
    "Age: Age value in years. Leave BLANK to let the model estimate ('solve for x').",
    "  • R_Date: Radiocarbon years BP (e.g., 3500)",
    "  • R_F14C: Fraction Modern Carbon for post-1950 samples (e.g., 1.0235)",
    "  • Date_CE/Tephra_CE: Calendar year CE/AD (use negative for BCE, e.g., -44 = 44 BCE)",
    "  • Date_calBP/Tephra_calBP: Calibrated years BP (e.g., 2500)",
    "  • Date/Tephra: Generic dated event",
    "",
    "Uncertainty: Error/Standard Deviation of the age (1-sigma).",
    "  • Leave blank if unknown or for Boundaries",
    "  • For Uniform distribution (U), this is the half-range",
    "",
    "Depth: Stratigraphic depth in consistent units (e.g., cm).",
    "",
    "Type: Date type - use dropdowns or type your own:",
    "  • R_Date: Radiocarbon date (BP)",
    "  • R_F14C: Fraction Modern Carbon (for post-1950 samples)",
    "  • Date: Generic known-age date",
    "  • Tephra: Generic tephra event",
    "  • Date_calBP / Tephra_calBP: Explicit calibrated BP",
    "  • Date_CE / Tephra_CE: Explicit CE/AD calendar date",
    "  • Boundary: Stratigraphic boundary",
    "  • Custom values accepted",
    "",
    "Unc_Type: Uncertainty distribution type:",
    "  • N: Normal distribution (default, most common)",
    "  • U: Uniform distribution (Age ± Uncertainty defines range)",
    "",
    "Outlier_Type: Outlier detection model (optional):",
    "  • General: T-distribution (recommended default for most dates)",
    "  • Charcoal: Exponential older bias (for old-wood effect)",
    "  • SSimple: Shift model",
    "  • RSimple: Ratio model",
    "  • TSimple: Time model",
    "  • RScaled: Scaled ratio model",
    "  • Leave blank for non-radiocarbon dates or to disable outlier detection",
    "  • Custom values accepted",
    "",
    "═══════════════════════════════════════════════════════════════════════════",
    "USAGE TIPS",
    "═══════════════════════════════════════════════════════════════════════════",
    "",
    "• Delete example rows before uploading your data",
    "• Dropdowns provide common options but you can type custom values",
    "• For unknown ages (solve for x): Leave Age and Uncertainty blank",
    "• Depth values must be numeric and in consistent units",
    "• For historical dates: Use Date_CE type with negative values for BCE",
    "• For modern samples: Use R_F14C type with Fraction Modern Carbon values",
    "",
    "═══════════════════════════════════════════════════════════════════════════"
  )
  
  instructions <- data.frame(Content = inst_text)
  writeData(wb, "Instructions", instructions, colNames = FALSE)
  
  title_style <- createStyle(fontSize = 14, textDecoration = "bold", fgFill = "#3c8dbc", fontColour = "white")
  addStyle(wb, "Instructions", title_style, rows = 1, cols = 1)
  setColWidths(wb, "Instructions", cols = 1, widths = 100)
  
  # Save
  if (verbose) cat("  Saving to", output_file, "...\n")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  if (verbose) cat("✓ Age-Depth template created successfully!\n")
  
  return(invisible(TRUE))
}

# ==============================================================================
# PHASE MODEL TEMPLATE
# ==============================================================================

generate_phase_template <- function(output_file, verbose = FALSE) {
  if (verbose) cat("Generating Phase Model Template...\n")
  
  wb <- createWorkbook()
  
  # Valid categories for dropdowns
  VALID_PHASES <- c("Before", "During", "After")
  VALID_DATE_TYPES <- c("R_Date", "R_F14C", "Date", "Tephra", "Date_calBP", "Tephra_calBP", "Date_CE", "Tephra_CE")
  VALID_UNC_TYPES <- c("N", "U")
  VALID_OUTLIER_TYPES <- c("General", "Charcoal", "SSimple", "RSimple", "TSimple", "RScaled", "")
  
  # --- SHEET 1: Data Entry ---
  addWorksheet(wb, "Data")
  
  headers <- c("Tephra_Name", "Phase", "Name", "Date_Type", "Age_Mean", "Age_SD", "Unc_Type", "Outlier_Type")
  
  # Example data
  data_df <- data.frame(
    Tephra_Name = c("My Tephra", "My Tephra", "My Tephra", "My Tephra", "My Tephra", "My Tephra"),
    Phase = c("Before", "During", "After", "Before", "After", "After"),
    Name = c("Sample-1", "Tephra independent date", "Sample-2", "Hist-Ref", "Old-Ref", "Modern-Sample"),
    Date_Type = c("R_Date", "Date", "R_Date", "Date_CE", "Date_calBP", "R_F14C"),
    Age_Mean = c(3500, 2670, 2200, -44, 2500, 1.0125),
    Age_SD = c(40, 50, 35, 10, 50, 0.0018),
    Unc_Type = c("N", "N", "N", "N", "U", "N"),
    Outlier_Type = c("General", "General", "Charcoal", "TSimple", "SSimple", "General"),
    stringsAsFactors = FALSE
  )
  
  writeData(wb, "Data", data_df)
  
  # Formatting
  header_style <- createStyle(
    fontSize = 11, textDecoration = "bold",
    fgFill = "#3c8dbc", fontColour = "white",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight", borderColour = "black"
  )
  addStyle(wb, "Data", header_style, rows = 1, cols = 1:8, gridExpand = TRUE)
  freezePane(wb, "Data", firstRow = TRUE)
  setColWidths(wb, "Data", cols = 1:8, widths = c(18, 12, 20, 16, 12, 12, 12, 16))
  
  # Data validation - Phase dropdown (strict - only Before/During/After)
  # Note: some environments ship older openxlsx versions that do not support
  # errorTitle/error arguments; keep arguments minimal for compatibility.
  dataValidation(wb, "Data",
    col = 2, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$A$2:$A$4",
    showErrorMsg = TRUE, allowBlank = FALSE
  )
  
  # Data validation - Date_Type dropdown (allow custom input)
  dataValidation(wb, "Data",
    col = 4, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$B$2:$B$9",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # Data validation - Unc_Type dropdown
  dataValidation(wb, "Data",
    col = 7, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$C$2:$C$3",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # Data validation - Outlier_Type dropdown (allow custom input)
  dataValidation(wb, "Data",
    col = 8, rows = 2:1000, type = "list",
    value = "'Valid_Categories'!$D$2:$D$7",
    showErrorMsg = FALSE, allowBlank = TRUE
  )
  
  # --- SHEET 2: Valid Categories ---
  addWorksheet(wb, "Valid_Categories")
  
  max_length <- max(length(VALID_PHASES), length(VALID_DATE_TYPES), length(VALID_UNC_TYPES), length(VALID_OUTLIER_TYPES))
  valid_data <- data.frame(
    Phase = c(VALID_PHASES, rep("", max_length - length(VALID_PHASES))),
    Date_Type = c(VALID_DATE_TYPES, rep("", max_length - length(VALID_DATE_TYPES))),
    Unc_Type = c(VALID_UNC_TYPES, rep("", max_length - length(VALID_UNC_TYPES))),
    Outlier_Type = c(VALID_OUTLIER_TYPES, rep("", max_length - length(VALID_OUTLIER_TYPES)))
  )
  writeData(wb, "Valid_Categories", valid_data)
  
  header_style2 <- createStyle(fontSize = 11, textDecoration = "bold", fgFill = "#d4edda")
  addStyle(wb, "Valid_Categories", header_style2, rows = 1, cols = 1:4, gridExpand = TRUE)
  setColWidths(wb, "Valid_Categories", cols = 1:4, widths = 20)
  
  # --- SHEET 3: Instructions ---
  addWorksheet(wb, "Instructions")
  
  inst_text <- c(
    "PHASE MODEL TEMPLATE",
    "",
    "═══════════════════════════════════════════════════════════════════════════",
    "COLUMN DEFINITIONS",
    "═══════════════════════════════════════════════════════════════════════════",
    "",
    "Tephra_Name: Name of the event defining the phases (can be same for all rows).",
    "",
    "Phase: Temporal relationship to the tephra event (REQUIRED - use dropdown):",
    "  • Before: Dates that predate the event",
    "  • During: Direct dates on the tephra itself (use Date type, not R_Date)",
    "  • After: Dates that postdate the event",
    "",
    "Name: Unique identifier for this specific date/sample.",
    "",
    "Date_Type: Type of date - use dropdowns or type your own:",
    "  • R_Date: Radiocarbon date (BP)",
    "  • R_F14C: Fraction Modern Carbon (for post-1950 samples)",
    "  • Date: Generic known-age date (e.g., independent tephra correlation)",
    "  • Tephra: Generic tephra event",
    "  • Date_calBP / Tephra_calBP: Explicit calibrated BP",
    "  • Date_CE / Tephra_CE: Explicit CE/AD calendar date",
    "  • Custom values accepted",
    "",
    "Age_Mean: Age value in years. Leave BLANK to let model estimate ('solve for x').",
    "  • R_Date: Radiocarbon years BP (e.g., 3500)",
    "  • R_F14C: Fraction Modern Carbon (e.g., 1.0125)",
    "  • Date_CE: Calendar year CE/AD (negative for BCE, e.g., -44 = 44 BCE)",
    "  • Date_calBP: Calibrated years BP",
    "",
    "Age_SD: Uncertainty (1-sigma) of the age.",
    "  • For Uniform distribution (U), this is the half-range",
    "",
    "Unc_Type: Uncertainty distribution type:",
    "  • N: Normal distribution (default, most common)",
    "  • U: Uniform distribution (Age_Mean ± Age_SD defines range)",
    "",
    "Outlier_Type: Outlier detection model (optional):",
    "  • General: T-distribution (recommended default)",
    "  • Charcoal: Exponential older bias (for old-wood effect)",
    "  • SSimple: Shift model",
    "  • RSimple: Ratio model",
    "  • TSimple: Time model",
    "  • RScaled: Scaled ratio model",
    "  • Leave blank to disable outlier detection",
    "  • Custom values accepted",
    "",
    "═══════════════════════════════════════════════════════════════════════════",
    "USAGE TIPS",
    "═══════════════════════════════════════════════════════════════════════════",
    "",
    "• Delete example rows before uploading your data",
    "• Phase MUST be exactly 'Before', 'During', or 'After' (case-sensitive)",
    "• For independent age constraints on the tephra (e.g., from correlation to",
    "  a dated eruption), use Phase='During' with Date_Type='Date'",
    "• Most dropdowns allow custom values except Phase (strict validation)",
    "• For unknown ages: Leave Age_Mean and Age_SD blank",
    "",
    "═══════════════════════════════════════════════════════════════════════════",
    "EXAMPLE USAGE",
    "═══════════════════════════════════════════════════════════════════════════",
    "",
    "Scenario: Dating 'Mazama Tephra' using dates before and after it",
    "",
    "Row 1: Tephra_Name='Mazama', Phase='Before', Name='Sample-1', Date_Type='R_Date',",
    "       Age_Mean=8500, Age_SD=75",
    "",
    "Row 2: Tephra_Name='Mazama', Phase='During', Name='Independent correlation',",
    "       Date_Type='Date', Age_Mean=7627, Age_SD=150",
    "       (This is an independent age estimate from correlation)",
    "",
    "Row 3: Tephra_Name='Mazama', Phase='After', Name='Sample-2', Date_Type='R_Date',",
    "       Age_Mean=6400, Age_SD=60",
    "",
    "═══════════════════════════════════════════════════════════════════════════"
  )
  
  instructions <- data.frame(Content = inst_text)
  writeData(wb, "Instructions", instructions, colNames = FALSE)
  
  title_style <- createStyle(fontSize = 14, textDecoration = "bold", fgFill = "#3c8dbc", fontColour = "white")
  addStyle(wb, "Instructions", title_style, rows = 1, cols = 1)
  setColWidths(wb, "Instructions", cols = 1, widths = 100)
  
  # Save
  if (verbose) cat("  Saving to", output_file, "...\n")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  if (verbose) cat("✓ Phase template created successfully!\n")
  
  return(invisible(TRUE))
}
