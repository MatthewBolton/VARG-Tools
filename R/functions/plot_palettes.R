# plot_palettes.R
# Helper functions for color palettes and shape definitions
# Used by mod_visualization.R for group customization

# =============================================================================
# RColorBrewer Palette Helpers
# =============================================================================

#' Get available qualitative RColorBrewer palettes suitable for discrete groups
get_brewer_qualitative_palettes <- function() {
    c(
        "Set1" = "Set1",
        "Set2" = "Set2",
        "Set3" = "Set3",
        "Paired" = "Paired",
        "Dark2" = "Dark2",
        "Accent" = "Accent",
        "Pastel1" = "Pastel1",
        "Pastel2" = "Pastel2"
    )
}

#' Get colors from a RColorBrewer palette, with fallback interpolation for many groups
#' @param palette_name Name of RColorBrewer palette
#' @param n Number of colors needed
#' @return Character vector of hex colors
get_brewer_colors <- function(palette_name, n) {
    # Get max colors for this palette
    max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]

    if (n <= max_colors) {
        # Direct palette extraction (minimum 3 for brewer.pal)
        colors <- RColorBrewer::brewer.pal(max(3, n), palette_name)
        return(colors[1:n])
    } else {
        # Interpolate for more colors
        base_colors <- RColorBrewer::brewer.pal(max_colors, palette_name)
        return(colorRampPalette(base_colors)(n))
    }
}

#' Generate default color assignments for groups using a palette
#' @param groups Character vector of group names
#' @param palette_name RColorBrewer palette name
#' @return Named list: list(group_name = hex_color)
generate_group_colors <- function(groups, palette_name = "Set1") {
    n <- length(groups)
    colors <- get_brewer_colors(palette_name, n)
    stats::setNames(as.list(colors), groups)
}

# =============================================================================
# Shape Definitions
# =============================================================================

#' Get available shapes with friendly names
#' @return Named vector: friendly_name = pch_value
get_available_shapes <- function() {
    c(
        "Circle (solid)" = 16,
        "Square (solid)" = 15,
        "Diamond (solid)" = 18,
        "Triangle up (solid)" = 17,
        "Triangle down (solid)" = 25,
        "Circle (hollow)" = 1,
        "Square (hollow)" = 0,
        "Diamond (hollow)" = 5,
        "Triangle up (hollow)" = 2,
        "Triangle down (hollow)" = 6,
        "Plus" = 3,
        "Cross" = 4,
        "Star" = 8,
        "Circle (filled)" = 21,
        "Square (filled)" = 22,
        "Diamond (filled)" = 23,
        "Triangle (filled)" = 24
    )
}

#' Get shape choices for selectInput (name = value format)
get_shape_choices <- function() {
    shapes <- get_available_shapes()
    stats::setNames(as.character(shapes), names(shapes))
}

#' Generate default shape assignments for groups
#' @param groups Character vector of group names
#' @return Named list: list(group_name = pch_value)
generate_group_shapes <- function(groups) {
    # Use a curated set of distinct solid shapes first, then cycle
    distinct_shapes <- c(16, 15, 18, 17, 3, 4, 8, 25, 6)
    n <- length(groups)

    # Cycle through shapes if needed
    shape_values <- rep(distinct_shapes, length.out = n)
    stats::setNames(as.list(shape_values), groups)
}

# =============================================================================
# Preset Tephra Palette (ashplotR-inspired, colorblind-friendly)
# =============================================================================

#' Get a curated tephra palette with 35 accessible colors
get_tephra_palette <- function() {
    c(
        "#209A24", "#A769EE", "#B4C428", "#e30303", "#000000",
        "#EF139D", "#3A700F", "#5288FF", "#E64701", "#7F9AF5",
        "#F36D20", "#0B3075", "#C18522", "#9634A6", "#1D7D46",
        "#80ddb8", "#7CBC92", "#B10E89", "#384B27", "#F896FD",
        "#C6A766", "#20628E", "#9C0323", "#C3AFE0", "#572D24",
        "#FF79A2", "#3E6D6F", "#E18C63", "#6F2B40", "#DF9FB4",
        "#291700", "#7e005b", "#00597a", "#7b003a", "#470600"
    )
}
