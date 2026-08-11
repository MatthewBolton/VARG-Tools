# tooltip_helpers.R
# Shared tooltip utility functions for VARG-Tools

# Simple tooltip helper using title attribute
tooltip <- function(text, tooltip_text) {
  tags$span(
    text,
    tags$i(
      class = "fa fa-question-circle",
      style = "color: #3c8dbc; margin-left: 5px; cursor: help;",
      title = tooltip_text,
      `data-toggle` = "tooltip",
      `data-html` = "true"
    )
  )
}

# Help text box for longer explanations
help_box <- function(title, content, icon_name = "info-circle") {
  div(
    class = "alert alert-info",
    style = "margin-top: 10px; padding: 0; background-color: #f8f9fa; border-left: 4px solid #3c8dbc; border-radius: 4px; overflow: hidden;",
    h5(
      style = "margin: 0; padding: 12px 15px; background: linear-gradient(135deg, #3c8dbc 0%, #2980b9 100%); color: #ffffff; font-weight: 600;",
      icon(icon_name),
      " ",
      title
    ),
    div(
      style = "margin: 0; padding: 15px; color: #000000; line-height: 1.6; font-size: 14px;",
      HTML(content)
    )
  )
}

# Inline help icon - FIXED to use bslib::popover for copyable content with hover + click
help_icon <- function(tooltip_text) {
  bslib::popover(
    trigger = tags$i(class = "fa fa-question-circle", style = "color: #3c8dbc; cursor: help; margin-left: 5px; font-size: 14px;", onclick = "event.stopPropagation(); event.preventDefault();"),
    HTML(tooltip_text),
    placement = "right",
    options = list(trigger = "hover click")
  )
}
