# global.R
# VARG-Tools: Global Configuration and Shared Resources

# --- Global App Version ---
APP_VERSION <- "0.9.8"

# --- Load Libraries ---
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(dplyr)
library(readxl)
library(DT)
library(plotly)
library(mclust)
library(robCompositions)
library(uwot)
library(RColorBrewer)
library(Polychrome)
library(rlang)
library(shinycssloaders)
library(sdcMicro)
library(future)
library(promises)
library(shinyWidgets)
library(shinyjs)

# Visualization Module dependencies
library(caret) # k-NN classification and kappa
library(ks) # KDE operations
library(MASS) # KDE for density peaks
library(robustbase) # Explicit: ensures ltsReg is on search path for impCoda's get()
library(foreach) # Parallel looping
library(doParallel) # Parallel backend
library(writexl) # Excel export
library(mgcv) # GAM for axis warping
library(tidyr) # Data tidying
library(purrr) # Functional programming
library(jsonlite) # JSON export for tie points
library(bslib) # Modern layout components
library(httr) # HTTP requests for Google Apps Script
library(digest) # Hashing for submission IDs
library(zip) # Zip archive creation for project bundles
library(markdown) # Markdown to HTML conversion for user guide
library(openxlsx) # Excel workbook creation for advanced templates
library(callr) # Background process execution with cancellation support

# Deployment shim: keep xfun in the manifest so sdcMicro's build-time dependency resolves on shinyapps.io.
if (FALSE) {
  xfun::attr(list(), "names")
}

# --- Global Options ---

# Source the VARG26 template generator — hidden with Contribute module
# if (file.exists("generate_varg26_template.R")) {
#   source("generate_varg26_template.R")
# }
options(shiny.maxRequestSize = 500 * 1024^2) # 500 MB limit

# --- Global Constants ---
VARG26_OXIDES <- c("SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O")

# Create a fixed color palette
P50 <- Polychrome::createPalette(50, c("#e58606", "#5d69b1", "#52bca3", "#99c945", "#cc61b0", "#24796c"))
names(P50) <- NULL

VALID_SHAPES <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "cross", "x", "star", "hexagram")

# --- Load Pre-trained Models (Lazy Loading) ---
# We check for existence here, but load them inside the module to avoid startup delays or errors if missing
MODEL_PATH_1D <- "varg26_model_1d.uwot"
MODEL_PATH_2D <- "varg26_model_2d.uwot"

check_models <- function() {
  missing <- c()
  if (!file.exists(MODEL_PATH_1D)) missing <- c(missing, MODEL_PATH_1D)
  if (!file.exists(MODEL_PATH_2D)) missing <- c(missing, MODEL_PATH_2D)

  if (length(missing) > 0) {
    warning(paste("Missing model files:", paste(missing, collapse = ", ")))
    return(FALSE)
  }
  return(TRUE)
}

# --- Async Configuration ---
# Detect if running on shinyapps.io (cloud) vs local development
# shinyapps.io sets SHINY_PORT environment variable
is_shinyapps <- nzchar(Sys.getenv("SHINY_PORT"))

# Background processes (callr::r_bg) spawn a separate R session which doubles
# memory usage. On memory-constrained servers (e.g., shinyapps.io free tier at
# 1GB), this causes OOM crashes. Disable on such servers; enable for local use
# or VMs with sufficient RAM. Override with VARG_BG_PROCESSES env var if needed.
USE_BG_PROCESSES <- if (nzchar(Sys.getenv("VARG_BG_PROCESSES"))) {
  tolower(trimws(Sys.getenv("VARG_BG_PROCESSES"))) %in% c("true", "1", "yes", "on")
} else {
  !is_shinyapps
}
message(sprintf("Background processes: %s", if (USE_BG_PROCESSES) "ENABLED" else "DISABLED (synchronous fallback)"))

if (is_shinyapps) {
  # shinyapps.io does NOT support multisession/multicore parallelization
  # Use sequential mode (slower but compatible)
  plan(sequential)
  message("Running on shinyapps.io - using sequential processing")
} else {
  # Local development: use 1 background worker for async tasks
  # This keeps the UI responsive while long tasks run in background
  plan(multisession, workers = 1)
  message("Running locally - using multisession with 1 worker")
}

# --- Data Normalization Helper ---
normalize_geochem_data <- function(df) {
  # Map common synonyms to VARG26 internal standard (FeO)
  # We assume FeOT/FeOt represents Total Iron calculated as FeO
  col_map <- c(
    "FeOT" = "FeO",
    "FeOt" = "FeO",
    "FeO_T" = "FeO",
    "feot" = "FeO",
    "Fe_Total" = "FeO"
  )

  current_names <- names(df)

  # Only rename if target "FeO" implies duplication or isn't present
  # Strategy: If "FeO" is missing but a synonym exists, rename it.
  if (!"FeO" %in% current_names) {
    for (syn in names(col_map)) {
      if (syn %in% current_names) {
        names(df)[names(df) == syn] <- col_map[syn]
        # Stop after first match to avoid confusion if multiple synonyms exist
        break
      }
    }
  }

  return(df)
}

# --- Template Generation Helper ---
generate_template <- function(type) {
  if (type == "geochem") {
    # Comprehensive VARG26 Workflow Template
    # Includes all columns needed from data entry through SC and Chronology modules

    data_df <- data.frame(
      # Sample identification
      Sample_ID = c("ABC-001", "ABC-001", "ABC-002", "ABC-002", "ABC-003"),
      Sample_point = c(1, 2, 1, 2, 1),

      # Stratigraphic context
      Site = c("Core A", "Core A", "Core A", "Core A", "Core B"),
      Depth = c(10.5, 10.5, 25.0, 25.0, 105.2),
      Depth_unit = c("cm", "cm", "cm", "cm", "cm"),

      # Chronology (optional - for age modeling)
      Age = c(NA, NA, 1200, 1200, 5400),
      Age_sigma = c(NA, NA, 50, 50, 120),
      Age_unit = c(NA, NA, "cal BP", "cal BP", "cal BP"),

      # VARG26 Major oxides (wt%) - required for UMAP projection
      SiO2 = c(75.4, 75.1, 74.8, 74.5, 68.2),
      TiO2 = c(0.25, 0.27, 0.28, 0.30, 0.65),
      Al2O3 = c(12.5, 12.6, 12.8, 12.7, 14.5),
      FeOT = c(1.8, 1.85, 1.9, 1.88, 3.5),
      MnO = c(0.05, 0.06, 0.06, 0.07, 0.12),
      MgO = c(0.15, 0.16, 0.18, 0.17, 0.85),
      CaO = c(1.2, 1.22, 1.3, 1.28, 2.5),
      Na2O = c(4.5, 4.48, 4.4, 4.42, 3.8),
      K2O = c(3.8, 3.82, 3.9, 3.88, 2.5),

      # Additional oxides (optional)
      P2O5 = c(0.02, 0.02, 0.03, 0.03, 0.08),
      Cl = c(0.18, 0.17, 0.20, 0.19, 0.12),
      stringsAsFactors = FALSE
    )

    instr_df <- data.frame(
      Column = c(
        "Sample_ID", "Sample_point",
        "Site", "Depth", "Depth_unit",
        "Age", "Age_sigma", "Age_unit",
        "SiO2-K2O", "P2O5, Cl, etc."
      ),
      Required = c(
        "Yes", "Optional",
        "Recommended", "Recommended", "Recommended",
        "Optional", "Optional", "Optional",
        "Yes (for VARG26)", "Optional"
      ),
      Description = c(
        "Unique identifier for the physical sample (e.g., tephra layer, glass population).",
        "Point number for multiple analyses of the same sample (e.g., different glass shards). Use sequential numbers: 1, 2, 3...",
        "Core ID, site name, or location identifier. Used for Stratigraphic Correlation.",
        "Depth below surface or other stratigraphic position. Used for SC and Chronology.",
        "Unit of depth measurement: cm, m, mm.",
        "Age of sample if known (from radiocarbon, tephra correlation, etc.).",
        "Uncertainty (1-sigma) of the age estimate.",
        "Age unit: ka, cal BP, Ma, CE, BCE.",
        "Major oxide concentrations in weight percent (wt%). These 9 oxides are required for VARG26 UMAP projection: SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O.",
        "Additional analytes can be included and will be preserved through the workflow."
      ),
      stringsAsFactors = FALSE
    )

    notes_df <- data.frame(
      Topic = c(
        "VARG26 Compatibility",
        "Multiple Analyses",
        "Missing Values",
        "Additional Columns",
        "Stratigraphic Correlation",
        "Chronology Module"
      ),
      Details = c(
        "To use the VARG26 pretrained UMAP, you MUST have all 9 major oxides: SiO2, TiO2, Al2O3, FeOT, MnO, MgO, CaO, Na2O, K2O. Note: Use FeOT (total iron as FeO), not Fe2O3.",
        "Each row = one analysis point. If you measured 5 points on sample ABC-001, you'll have 5 rows with Sample_ID='ABC-001' and Sample_point=1,2,3,4,5.",
        "Leave cells empty or use NA for missing values. The app will handle imputation if needed.",
        "You can add ANY additional columns (trace elements, other metadata). They will be preserved and available for visualization and export.",
        "For Stratigraphic Correlation, ensure you have Site/CoreID, Sample_ID, Depth, and a correlation variable (e.g., one of your UMAP coordinates).",
        "For age modeling in the Chronology module, include Age, Age_sigma, and Age_unit for dated samples."
      ),
      stringsAsFactors = FALSE
    )

    return(list(Data = data_df, Instructions = instr_df, Notes = notes_df))
  } else if (type == "chron_age_depth") {
    # 2. Age-Depth Template
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
      Outlier_Type = c("General", "Charcoal", NA, "General", "General", "TSimple", "General", NA),
      stringsAsFactors = FALSE
    )

    instr_df <- data.frame(
      Column = c("Name", "Age", "Uncertainty", "Depth", "Type", "Unc_Type", "Outlier_Type"),
      Description = c(
        "Unique identifier for the sample.",
        "Age value. Leave BLANK to let the model estimate the age ('solve for x'). For CE dates, use negative values for BCE. For R_F14C, use Fraction Modern Carbon (e.g., 1.0235).",
        "Error/Standard Deviation of the age. Leave blank if unknown or for Boundaries. For U type, this is the half-range.",
        "Depth in consistent units (e.g., cm).",
        "Accepted types: 'R_Date' (BP), 'R_F14C' (Fraction Modern Carbon for post-1950), 'Date' or 'Tephra' (Generic), 'Date_calBP' or 'Tephra_calBP' (Explicit calBP), 'Date_CE' or 'Tephra_CE' (Explicit CE/AD), 'Boundary'.",
        "Uncertainty distribution: 'N' (Normal, default) or 'U' (Uniform). For U, use Age±Uncertainty range.",
        "Outlier model (optional): 'General' (T-distribution), 'Charcoal' (Exp older bias), 'SSimple' (Shift), 'RSimple' (Ratio), 'TSimple' (Time), 'RScaled' (Scaled ratio). Leave blank for non-radiocarbon dates or to use default (General for R_Date/R_F14C)."
      ),
      stringsAsFactors = FALSE
    )

    notes_df <- data.frame(
      Topic = c(
        "Date Types",
        "Outlier Models",
        "Uncertainty Types"
      ),
      Details = c(
        "R_Date: Radiocarbon date (BP). R_F14C: Fraction Modern Carbon (for post-1950 samples, e.g., 1.0235). Date_CE: Calendar date (CE/AD, use negative for BCE). Date_calBP: Calibrated years BP. Tephra/Date: Generic dated events.",
        "General: T-distribution (recommended default for most dates). Charcoal: Exponential older bias (for old-wood effect). SSimple: Shift model. RSimple: Ratio model. TSimple: Time model.",
        "N: Normal distribution (default, most common). U: Uniform distribution (Age ± Uncertainty defines range, useful for historical references with known bounds)."
      ),
      stringsAsFactors = FALSE
    )

    return(list(Data = data_df, Instructions = instr_df, Notes = notes_df))
  } else if (type == "chron_phase") {
    # 3. Phase Model Template
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

    instr_df <- data.frame(
      Column = c("Tephra Name", "Phase", "Date Type", "Age Mean", "Age SD", "Unc_Type", "Outlier Type"),
      Description = c(
        "Name of the event defining the phases.",
        "Must be exactly: 'Before', 'During', or 'After'.",
        "Type of date: 'R_Date', 'R_F14C' (Fraction Modern Carbon for post-1950), 'Date', 'Tephra', 'Date_calBP', 'Date_CE', etc.",
        "Age value. Leave BLANK to let the model estimate the age. For CE, negative = BCE. For R_F14C, use F14C ratio.",
        "Uncertainty (1-sigma) of the age. For U distribution, this is the half-range.",
        "Uncertainty distribution: 'N' (Normal, default) or 'U' (Uniform). For U, Age_Mean ± Age_SD defines the range.",
        "Outlier model: 'General', 'Charcoal', 'SSimple', 'RSimple', 'TSimple', 'RScaled'. Unknown types default to 'General'."
      ),
      stringsAsFactors = FALSE
    )

    notes_df <- data.frame(
      Topic = c(
        "Phase Structure",
        "Date Types",
        "Outlier Models",
        "Independent Tephra Dates"
      ),
      Details = c(
        "Before: Dates that predate the event. During: Direct dates on the tephra (use Date type, not R_Date). After: Dates that postdate the event.",
        "R_Date: Radiocarbon date (BP). R_F14C: Fraction Modern Carbon (post-1950). Date_CE: Calendar date (CE/AD, negative for BCE). Date: Generic known-age date (e.g., independent tephra correlation). Date_calBP: Calibrated years BP.",
        "General: T-distribution (recommended default). Charcoal: Exponential older bias. SSimple: Shift model. RSimple: Ratio model. TSimple: Time model.",
        "For independent age constraints on the tephra itself (e.g., from correlation to a dated eruption), use Phase='During' with Date_Type='Date' and the known age."
      ),
      stringsAsFactors = FALSE
    )

    return(list(Data = data_df, Instructions = instr_df, Notes = notes_df))
  }
  return(NULL)
}
