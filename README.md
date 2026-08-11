# VARG-Tools v0.9.8

**VARG-Tools** is an R Shiny toolkit for tephrochronology and volcanic-glass geochemical analysis. It supports the workflow from data preparation through clustering, visualization, stratigraphic correlation, and chronology-code generation.

## Key Features

### 1. Data Processing
*   **Imputation & Transformation**: Handles missing values and applies Isometric Log-Ratio (ILR) transformations for compositional data analysis.
*   **GMM Clustering**: Automated Gaussian Mixture Model clustering with BIC optimization to identify geochemical populations.
*   **Dimensionality Reduction**: Generate 1D and 2D UMAP projections to explore data structure.
*   **Population Definition**: Interactive tools to assign samples to populations based on clustering or manual selection.
*   **Advanced Mode**: Flexible column selection for mixing compositional and non-compositional variables.

### 2. Visualization
*   **Custom Scatter Plots**: Publication-ready plots with customizable aesthetics, dimensions, and export options (PDF/PNG).
*   **Stratigraphic Correlation**: Interactive tools to align records using tie points, exact-anchor monotonic piecewise-linear warping, and display-only offset/affine previews.

### 3. Chronology
*   **Age-Depth Model Generation**: Produces auditable OxCal code for age-depth and phase models; OxCal is run separately.
*   **Tie-Point Management**: Link tephra layers across sites to build regional chronologies.

## Deployment & Installation

### Online (Recommended)
VARG-Tools is available through the project portal at <https://matthewbolton.github.io/VARG_Tools/>. The portal links to the hosted review application and current access information.

### Local Distributions
Versioned local downloads are distributed through the project portal and the [GitHub releases page](https://github.com/MatthewBolton/VARG-Tools/releases).

- **Windows standalone:** includes its own private R 4.4.2 runtime and package library. Users extract the ZIP and run `Launch VARG-Tools.cmd`; no R or RStudio installation is required.
- **R-user bundle:** requires R 4.4 or newer; RStudio is optional. Users run `Install Dependencies.cmd` once on Windows (or `Rscript app/install_dependencies.R` on any supported platform), then launch with `Run VARG-Tools.cmd` or `Rscript app/run_vargtools.R`.

Both local distributions bind only to `127.0.0.1` and process uploaded data on the local computer. The Windows standalone uses only its bundled R runtime and package library. The R-user launcher places its bundle-local library first while retaining base, recommended, and existing compatible R libraries as fallbacks.

The curated application source is available at <https://github.com/MatthewBolton/VARG-Tools>. Scientific outputs remain conditional on the input data and user decisions and should be checked before interpretation.

## Documentation

* **User guide:** Open **User Guide** from the app header for the first-run walkthrough, task-based workflows, scientific checkpoints, and troubleshooting.
* **Tooltips:** Hover over a `?` icon for contextual help on an individual control.
* **Templates and practice data:** Download the current files from Home or the relevant Chronology workflow so column names match the installed version.

## License

VARG-Tools is released under the MIT License. See `LICENSE`.
