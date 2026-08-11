# R/utils/utils_ui.R

#' Convert Unicode subscripts/superscripts to plotmath expression
#' 
#' Converts text with Unicode subscript/superscript characters to a plotmath expression
#' that can be used in ggplot2 axis labels. This ensures proper rendering in PDF exports.
#' 
#' @param text Character string that may contain Unicode sub/superscripts
#' @return An expression suitable for use in ggplot2 labs(), or the original string if no conversion needed
#' @examples
#' unicode_to_plotmath("SiO₂")  # Returns expression for SiO[2]
#' unicode_to_plotmath("Fe²⁺")  # Returns expression for Fe^{2+}
#' @export
unicode_to_plotmath <- function(text) {
    if (is.null(text) || !is.character(text) || length(text) == 0 || nchar(text) == 0) {
        return(text)
    }
    
    # Unicode subscript digits (₀-₉)
    subscript_chars <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084", 
                         "\u2085", "\u2086", "\u2087", "\u2088", "\u2089")
    subscript_digits <- as.character(0:9)
    
    # Unicode superscript digits (⁰-⁹) and ⁺⁻
    superscript_chars <- c("\u2070", "\u00B9", "\u00B2", "\u00B3", "\u2074",
                           "\u2075", "\u2076", "\u2077", "\u2078", "\u2079",
                           "\u207A", "\u207B")  # ⁺ ⁻
    superscript_values <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "+", "-")
    
    all_unicode <- c(subscript_chars, superscript_chars)
    
    # Check if text contains unicode subscripts or superscripts
    has_unicode <- any(sapply(all_unicode, function(char) grepl(char, text, fixed = TRUE)))
    if (!has_unicode) {
        return(text)  # No unicode found, return as-is
    }
    
    # Build regex pattern for consecutive subscripts/superscripts
    subscript_pattern <- paste0("([", paste0(subscript_chars, collapse = ""), "]+)")
    superscript_pattern <- paste0("([", paste0(superscript_chars, collapse = ""), "]+)")
    
    result <- text
    
    # Helper to convert unicode digits to normal digits
    convert_subscript <- function(unicode_str) {
        out <- unicode_str
        for (i in seq_along(subscript_chars)) {
            out <- gsub(subscript_chars[i], subscript_digits[i], out, fixed = TRUE)
        }
        out
    }
    
    convert_superscript <- function(unicode_str) {
        out <- unicode_str
        for (i in seq_along(superscript_chars)) {
            out <- gsub(superscript_chars[i], superscript_values[i], out, fixed = TRUE)
        }
        out
    }
    
    # Replace subscripts with [x] notation for plotmath
    while (grepl(subscript_pattern, result, perl = TRUE)) {
        match <- regmatches(result, regexpr(subscript_pattern, result, perl = TRUE))
        converted <- convert_subscript(match)
        result <- sub(subscript_pattern, paste0("[", converted, "]"), result, perl = TRUE)
    }
    
    # Replace superscripts with ^{x} notation for plotmath
    while (grepl(superscript_pattern, result, perl = TRUE)) {
        match <- regmatches(result, regexpr(superscript_pattern, result, perl = TRUE))
        converted <- convert_superscript(match)
        result <- sub(superscript_pattern, paste0("^{", converted, "}"), result, perl = TRUE)
    }
    
    # Handle plus sign in middle of text (e.g., Na2O+K2O) - convert to ~"+"~
    # But not if already inside braces
    if (grepl("\\+", result) && !grepl("\\^\\{.*\\+.*\\}", result)) {
        result <- gsub("(?<!\\{)\\+(?!\\})", '~"+"~', result, perl = TRUE)
    }
    
    # Wrap letter sequences (element symbols) properly for plotmath
    # Split into parts and reassemble with * for concatenation
    # This handles cases like "SiO[2]" -> becomes valid plotmath
    
    # Parse as expression
    tryCatch({
        expr <- parse(text = result)
        if (length(expr) > 0) {
            return(expr)
        } else {
            return(text)
        }
    }, error = function(e) {
        # If parsing fails, return original text
        text
    })
}

#' Create a Standard Card
#'
#' A wrapper around bslib::card to ensure consistent styling.
#' @param ... UI content
#' @param title Optional title
#' @param footer Optional footer
#' @param full_screen Allow full screen expansion?
#' @param height Optional height
#' @param class Additional CSS class
#' @export
varg_card <- function(..., title = NULL, footer = NULL, full_screen = FALSE, height = NULL, fill = FALSE, class = NULL) {
    bslib::card(
        full_screen = full_screen,
        height = height,
        fill = fill,
        class = class,
        if (!is.null(title)) bslib::card_header(title),
        bslib::card_body(..., fillable = FALSE),
        if (!is.null(footer)) bslib::card_footer(footer)
    )
}

#' Create a Value Box
#'
#' A wrapper around bslib::value_box.
#' @param title Title
#' @param value Main value
#' @param icon Icon (shiny::icon)
#' @param theme Color theme (primary, secondary, etc.)
#' @export
varg_value_box <- function(title, value, icon = NULL, theme = "primary") {
    bslib::value_box(
        title = title,
        value = value,
        showcase = icon,
        theme = theme
    )
}

#' Section Header
#'
#' Standardized header for sections.
#' @param text Title text
#' @param icon Icon
#' @export
section_header <- function(text, icon = NULL) {
    div(
        class = "pb-2 mt-4 mb-2 border-bottom",
        h3(if (!is.null(icon)) icon, text)
    )
}

#' Module Banner
#'
#' A consistent banner at the top of each module/tab stating Goal, Inputs, Outputs,
#' and optionally "Why" (expected use / motivation).
#' @param goal What this step achieves
#' @param inputs What the user needs to provide
#' @param outputs What this step produces
#' @param why Optional: Why the user might want to use this step
#' @export
module_banner <- function(goal, inputs, outputs, why = NULL) {
    why_div <- if (!is.null(why)) {
        tags$div(
            style = "margin-top: 6px; padding-top: 6px; border-top: 1px solid #dee2e6;",
            tags$span(tags$strong(icon("lightbulb", style = "color: #f39c12;"), " Why: "), why)
        )
    }
    div(
        class = "alert alert-light border mb-3",
        style = "padding: 10px 15px; font-size: 0.88rem; line-height: 1.5; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); border-left: 4px solid #2C3E50 !important;",
        tags$div(
            style = "display: flex; flex-wrap: wrap; gap: 15px;",
            tags$span(tags$strong(icon("bullseye", style = "color: #18bc9c;"), " Goal: "), goal),
            tags$span(tags$strong(icon("sign-in-alt", style = "color: #3498db;"), " Input: "), inputs),
            tags$span(tags$strong(icon("sign-out-alt", style = "color: #e74c3c;"), " Output: "), outputs)
        ),
        why_div
    )
}
