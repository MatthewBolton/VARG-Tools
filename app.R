# app.R
source("R/utils/package_compatibility.R")
varg_assert_package_compatibility()
source("global.R")
source("R/utils/tooltip_helpers.R")
source("R/utils/utils_ui.R")
source("R/utils/heavy_job_limiter.R")
source("R/functions/sc_helpers.R") # Stratigraphic correlation helpers
source("R/functions/fct_analysis.R") # Analysis functions (GMM, UMAP)
# source("R/functions/validate_varg26_submission.R") # VARG26 submission validation — hidden with Contribute module
source("R/modules/mod_processing.R")
source("R/modules/mod_visualization.R")
source("R/modules/mod_chronology.R")
# NOTE: The Contribute module is hidden pending development of the curation
# and citation framework. Code preserved for future reactivation.
# source("R/modules/mod_contribution.R")
source("R/modules/mod_project_manager.R")
library(shinyWidgets)
library(shiny) # ensure shiny namespace is attached for global APP_VERSION usage

guide_file <- "docs/VARG-Tools_User_Guide.md"
guide_asset_dir <- NULL

if (dir.exists("docs/guide")) {
  guide_asset_dir <- normalizePath("docs/guide", winslash = "/", mustWork = TRUE)
  addResourcePath("guide-assets", guide_asset_dir)
}

rewrite_guide_markdown_images <- function(markdown_text, transform_path) {
  pattern <- "!\\[([^\\]]*)\\]\\((guide-assets/[^)]+)\\)"
  matches <- gregexpr(pattern, markdown_text, perl = TRUE)

  if (matches[[1]][1] == -1) {
    return(markdown_text)
  }

  image_refs <- regmatches(markdown_text, matches)[[1]]
  replacements <- vapply(image_refs, function(image_ref) {
    alt_text <- sub(pattern, "\\1", image_ref, perl = TRUE)
    relative_path <- sub(pattern, "\\2", image_ref, perl = TRUE)
    sprintf("![%s](%s)", alt_text, transform_path(relative_path))
  }, character(1))

  regmatches(markdown_text, matches) <- list(replacements)
  markdown_text
}

rewrite_guide_html_images <- function(html_text, transform_path) {
  pattern <- '(<img[^>]+src=["\'])([^"\']*guide-assets/[^"\']+)(["\'][^>]*>)'
  matches <- gregexpr(pattern, html_text, perl = TRUE)

  if (matches[[1]][1] == -1) {
    return(html_text)
  }

  image_tags <- regmatches(html_text, matches)[[1]]
  replacements <- vapply(image_tags, function(image_tag) {
    prefix <- sub(pattern, "\\1", image_tag, perl = TRUE)
    src <- sub(pattern, "\\2", image_tag, perl = TRUE)
    suffix <- sub(pattern, "\\3", image_tag, perl = TRUE)
    paste0(prefix, transform_path(src), suffix)
  }, character(1))

  regmatches(html_text, matches) <- list(replacements)
  html_text
}

guide_image_resource_path <- function(relative_path) {
  relative_path
}

guide_image_data_uri <- function(relative_path) {
  if (is.null(guide_asset_dir)) {
    stop("Guide asset directory is not available.", call. = FALSE)
  }

  asset_path <- file.path(guide_asset_dir, sub("^/?guide-assets/", "", relative_path))
  if (!file.exists(asset_path)) {
    stop(sprintf("Guide asset not found: %s", relative_path), call. = FALSE)
  }

  extension <- tolower(tools::file_ext(asset_path))
  mime_type <- switch(
    extension,
    png = "image/png",
    jpg = "image/jpeg",
    jpeg = "image/jpeg",
    svg = "image/svg+xml",
    gif = "image/gif",
    webp = "image/webp",
    "application/octet-stream"
  )

  asset_size <- file.info(asset_path)$size
  asset_raw <- readBin(asset_path, what = "raw", n = asset_size)
  asset_b64 <- jsonlite::base64_enc(asset_raw)
  asset_b64 <- gsub("[\r\n]", "", asset_b64)
  paste0("data:", mime_type, ";base64,", asset_b64)
}

render_guide_markdown_fragment <- function(markdown_text) {
  if (requireNamespace("commonmark", quietly = TRUE)) {
    return(commonmark::markdown_html(
      markdown_text,
      extensions = c("table", "autolink", "strikethrough")
    ))
  }

  external_prefix <- "https://guide-assets.local/"
  safe_markdown <- gsub(
    "\\((/)?guide-assets/",
    paste0("(", external_prefix),
    markdown_text,
    perl = TRUE
  )
  safe_html <- markdown::markdownToHTML(
    text = safe_markdown,
    fragment.only = TRUE,
    options = ""
  )
  gsub(external_prefix, "guide-assets/", safe_html, fixed = TRUE)
}

render_guide_html <- function(embed_images = FALSE) {
  if (!file.exists(guide_file)) {
    stop(sprintf("User guide file not found: %s", guide_file), call. = FALSE)
  }

  guide_content <- paste(readLines(guide_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  if (!embed_images) {
    guide_content <- rewrite_guide_markdown_images(guide_content, guide_image_resource_path)
    return(render_guide_markdown_fragment(guide_content))
  }

  guide_html <- render_guide_markdown_fragment(guide_content)
  rewrite_guide_html_images(guide_html, guide_image_data_uri)
}

guide_modal_styles <- "
  .modal-dialog.modal-xl {
    max-width: min(96vw, 1500px);
    width: min(96vw, 1500px);
  }
  .guide-layout {
    display: flex;
    gap: 0;
    height: 78vh;
    min-height: 560px;
  }
  .guide-toc {
    width: 260px;
    min-width: 260px;
    max-width: 260px;
    overflow-y: auto;
    border-right: 1px solid #dee2e6;
    padding: 14px 14px 18px;
    font-size: 0.84rem;
    background: #f8f9fa;
  }
  .guide-toc a {
    display: block;
    color: #2c3e50;
    text-decoration: none;
    border-left: 2px solid transparent;
    padding: 0.25rem 0 0.25rem 0.55rem;
    line-height: 1.35;
  }
  .guide-toc a:hover,
  .guide-toc a.active {
    color: #18bc9c;
    border-left-color: #18bc9c;
  }
  .guide-toc a.toc-h2 { padding-left: 1.1rem; }
  .guide-toc a.toc-h3 {
    padding-left: 1.8rem;
    font-size: 0.78rem;
    color: #666;
  }
  .guide-content {
    flex: 1;
    min-width: 0;
    overflow-y: auto;
    overflow-x: hidden;
    overflow-wrap: anywhere;
    padding: 14px 28px 24px;
  }
  .guide-content > p,
  .guide-content > ul,
  .guide-content > ol,
  .guide-export-content > p,
  .guide-export-content > ul,
  .guide-export-content > ol {
    max-width: 76ch;
  }
  .guide-content h1,
  .guide-content h2,
  .guide-content h3 {
    line-height: 1.2;
  }
  .guide-content h1 {
    font-size: 1.65rem;
    margin-top: 0.5rem;
  }
  .guide-content h2 {
    font-size: 1.28rem;
    margin-top: 1.8rem;
  }
  .guide-content h3 {
    font-size: 1.08rem;
    margin-top: 1.1rem;
  }
  .guide-meta {
    display: flex;
    flex-wrap: wrap;
    gap: 0.25rem 1rem;
  }
  .guide-meta span {
    white-space: nowrap;
  }
  .guide-content img,
  .guide-export img,
  .modal-body img {
    display: block;
    max-width: 100%;
    width: auto;
    height: auto;
    max-height: 60vh;
    margin: 1rem auto;
    object-fit: contain;
  }
  .guide-content figure,
  .guide-export figure {
    margin: 1.25rem 0;
  }
  .guide-content pre,
  .guide-export pre {
    white-space: pre-wrap;
    word-break: break-word;
  }
  .guide-content table,
  .guide-export table,
  .modal-body table {
    width: 100%;
    margin-top: 1.2rem;
    margin-bottom: 1.5rem;
    border-collapse: collapse;
    border: none;
    font-size: 0.92rem;
    line-height: 1.45;
  }
  .guide-content table th,
  .guide-export table th,
  .modal-body table th {
    background-color: #e9eff1;
    color: #17252a;
    padding: 0.7rem 0.75rem;
    text-align: left;
    font-weight: 600;
    border: none;
    border-bottom: 1px solid #bfcacd;
    vertical-align: top;
  }
  .guide-content table td,
  .guide-export table td,
  .modal-body table td {
    padding: 0.7rem 0.75rem;
    border: none;
    border-bottom: 1px solid #d8dfe1;
    vertical-align: top;
  }
  .guide-content table tbody tr:nth-child(even),
  .guide-export table tbody tr:nth-child(even),
  .modal-body table tbody tr:nth-child(even) {
    background-color: rgba(233, 239, 241, 0.38);
  }
  .guide-content table tbody tr:hover,
  .guide-export table tbody tr:hover,
  .modal-body table tbody tr:hover {
    background-color: #f1f5f6;
  }
  .guide-content table p,
  .guide-export table p,
  .modal-body table p {
    margin-bottom: 0.5rem;
  }
  .guide-content table p:last-child,
  .guide-export table p:last-child,
  .modal-body table p:last-child {
    margin-bottom: 0;
  }
  @media (max-width: 1199px) {
    .modal-dialog.modal-xl {
      max-width: calc(100vw - 24px);
      width: calc(100vw - 24px);
      margin: 12px auto;
    }
    .modal-dialog.modal-xl > .modal-content,
    .modal-dialog.modal-xl .modal-body {
      min-width: 0;
      max-width: 100%;
    }
    .guide-layout {
      flex-direction: column;
      height: auto;
      min-height: 0;
    }
    .guide-toc {
      width: 100%;
      min-width: 0;
      max-width: none;
      max-height: 22vh;
      border-right: none;
      border-bottom: 1px solid #dee2e6;
    }
    .guide-toc a {
      overflow-wrap: anywhere;
    }
    .guide-content {
      width: 100%;
      max-height: 56vh;
      padding: 14px 18px 20px;
      box-sizing: border-box;
    }
    .guide-content h1 {
      font-size: 1.4rem;
    }
    .guide-content table,
    .guide-content pre {
      display: block;
      max-width: 100%;
      overflow-x: auto;
    }
  }
"

guide_modal_script <- function(content_id, toc_id, heading_prefix) {
  sprintf(
    "
      setTimeout(function() {
        var content = document.getElementById('%s');
        var toc = document.getElementById('%s');
        if (!content || !toc) return;
        toc.innerHTML = '';
        var headings = content.querySelectorAll('h1');
        headings.forEach(function(h, i) {
          if (!h.id) h.id = '%s' + i;
          var a = document.createElement('a');
          a.href = '#' + h.id;
          a.textContent = h.textContent;
          a.className = 'toc-' + h.tagName.toLowerCase();
          a.onclick = function(e) {
            e.preventDefault();
            h.scrollIntoView({behavior: 'smooth', block: 'start'});
          };
          toc.appendChild(a);
        });
      }, 200);
    ",
    content_id,
    toc_id,
    heading_prefix
  )
}

show_guide_modal <- function(download_id, content_id, toc_id, heading_prefix) {
  guide_html <- HTML(render_guide_html(embed_images = FALSE))

  showModal(modalDialog(
    title = "VARG-Tools User Guide",
    size = "xl",
    easyClose = TRUE,
    tags$head(tags$style(HTML(guide_modal_styles))),
    div(
      class = "guide-layout",
      div(class = "guide-toc", id = toc_id),
      div(class = "guide-content", id = content_id, guide_html)
    ),
    tags$script(HTML(guide_modal_script(content_id, toc_id, heading_prefix))),
    footer = tagList(
      downloadButton(download_id, "Download as HTML", class = "btn-primary"),
      modalButton("Close")
    )
  ))
}

build_guide_export_document <- function() {
  guide_body <- render_guide_html(embed_images = TRUE)

  paste0(
    '<!DOCTYPE html>\n',
    '<html>\n',
    '<head>\n',
    '  <meta charset="utf-8">\n',
    '  <meta name="viewport" content="width=device-width, initial-scale=1">\n',
    '  <title>VARG-Tools User Guide</title>\n',
    '  <style>\n',
    '    body {\n',
    '      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;\n',
    '      line-height: 1.6;\n',
    '      margin: 0;\n',
    '      padding: 0;\n',
    '      background-color: #ffffff;\n',
    '      color: #333;\n',
    '    }\n',
    '    .guide-export-shell {\n',
    '      display: flex;\n',
    '      min-height: 100vh;\n',
    '    }\n',
    '    .guide-export-toc {\n',
    '      width: 280px;\n',
    '      min-width: 280px;\n',
    '      max-width: 280px;\n',
    '      overflow-y: auto;\n',
    '      border-right: 1px solid #dee2e6;\n',
    '      padding: 18px 14px 24px;\n',
    '      background: #f8f9fa;\n',
    '      position: sticky;\n',
    '      top: 0;\n',
    '      height: 100vh;\n',
    '      box-sizing: border-box;\n',
    '    }\n',
    '    .guide-export-shell.is-collapsed .guide-export-toc {\n',
    '      width: 56px;\n',
    '      min-width: 56px;\n',
    '      max-width: 56px;\n',
    '      padding: 14px 8px;\n',
    '      overflow: hidden;\n',
    '    }\n',
    '    .guide-export-toggle {\n',
    '      display: inline-flex;\n',
    '      align-items: center;\n',
    '      justify-content: center;\n',
    '      width: 100%;\n',
    '      margin: 0 0 0.75rem;\n',
    '      padding: 0.45rem 0.4rem;\n',
    '      border: 1px solid #cfd8dc;\n',
    '      border-radius: 0.375rem;\n',
    '      background: #ffffff;\n',
    '      color: #2c3e50;\n',
    '      font-size: 0.85rem;\n',
    '      font-weight: 600;\n',
    '      cursor: pointer;\n',
    '    }\n',
    '    .guide-export-toggle:hover {\n',
    '      background: #eef3f6;\n',
    '    }\n',
    '    .guide-export-toc .toc-title {\n',
    '      font-size: 0.95rem;\n',
    '      font-weight: 700;\n',
    '      margin: 0 0 0.75rem;\n',
    '      color: #2c3e50;\n',
    '    }\n',
    '    .guide-export-toc a {\n',
    '      display: block;\n',
    '      color: #2c3e50;\n',
    '      text-decoration: none;\n',
    '      border-left: 2px solid transparent;\n',
    '      padding: 0.25rem 0 0.25rem 0.55rem;\n',
    '      line-height: 1.35;\n',
    '      margin: 0.08rem 0;\n',
    '    }\n',
    '    .guide-export-toc a:hover,\n',
    '    .guide-export-toc a.active {\n',
    '      color: #18bc9c;\n',
    '      border-left-color: #18bc9c;\n',
    '    }\n',
    '    .guide-export-toc a.toc-h2 { padding-left: 1.1rem; }\n',
    '    .guide-export-toc a.toc-h3 {\n',
    '      padding-left: 1.8rem;\n',
    '      font-size: 0.78rem;\n',
    '      color: #666;\n',
    '    }\n',
    '    .guide-export-shell.is-collapsed .guide-export-toc .toc-title,\n',
    '    .guide-export-shell.is-collapsed .guide-export-toc a {\n',
    '      display: none;\n',
    '    }\n',
    '    .guide-export-shell.is-collapsed .guide-export-toggle {\n',
    '      writing-mode: vertical-rl;\n',
    '      transform: rotate(180deg);\n',
    '      width: 100%;\n',
    '      min-height: 180px;\n',
    '      margin: 0;\n',
    '    }\n',
    '    .guide-export-content {\n',
    '      flex: 1;\n',
    '      min-width: 0;\n',
    '      padding: 24px 32px 32px;\n',
    '      max-width: 1280px;\n',
    '      margin: 0 auto;\n',
    '      box-sizing: border-box;\n',
    '    }\n',
    '    .guide-export-content h1,\n',
    '    .guide-export-content h2,\n',
    '    .guide-export-content h3 {\n',
    '      line-height: 1.2;\n',
    '    }\n',
    '    .guide-export-content h1 { font-size: 1.9rem; margin-top: 0.4rem; }\n',
    '    .guide-export-content h2 { font-size: 1.4rem; margin-top: 1.5rem; }\n',
    '    .guide-export-content h3 { font-size: 1.1rem; margin-top: 1.1rem; }\n',
    '    a { color: #0d6efd; text-decoration: none; }\n',
    '    a:hover { text-decoration: underline; }\n',
    '    code { background-color: #f8f9fa; padding: 2px 6px; border-radius: 3px; font-family: "Courier New", monospace; }\n',
    '    pre { background-color: #f8f9fa; border-radius: 0.375rem; padding: 1rem; overflow-x: auto; }\n',
    '    pre code { background-color: transparent; padding: 0; }\n',
    '    blockquote { border-left: 1px solid #cbd3d5; padding-left: 1rem; margin-left: 0; color: #5f6b70; font-style: italic; }\n',
    '    hr { border: none; border-top: 1px solid #dee2e6; margin: 2rem 0; }\n',
    '    .guide-export-content img,\n',
    '    .guide-export-content figure {\n',
    '      display: block;\n',
    '      max-width: 100%;\n',
    '      width: auto;\n',
    '      height: auto;\n',
    '      max-height: 60vh;\n',
    '      margin: 1rem auto;\n',
    '      object-fit: contain;\n',
    '    }\n',
    '    .guide-export-content figure { margin: 1.25rem 0; }\n',
    '    .guide-export-content table {\n',
    '      width: 100%;\n',
    '      margin-top: 1.5rem;\n',
    '      margin-bottom: 1.5rem;\n',
    '      border-collapse: collapse;\n',
    '      border: 1px solid #dee2e6;\n',
    '      font-size: 0.95rem;\n',
    '      line-height: 1.6;\n',
    '    }\n',
    '    .guide-export-content table th {\n',
    '      background-color: #2c3e50;\n',
    '      color: white;\n',
    '      padding: 0.875rem;\n',
    '      text-align: left;\n',
    '      font-weight: 600;\n',
    '      border: 1px solid #dee2e6;\n',
    '      vertical-align: top;\n',
    '    }\n',
    '    .guide-export-content table td {\n',
    '      padding: 0.75rem;\n',
    '      border: 1px solid #dee2e6;\n',
    '      vertical-align: top;\n',
    '    }\n',
    '    .guide-export-content table tbody tr:nth-child(even) { background-color: #f8f9fa; }\n',
    '    .guide-export-content table tbody tr:hover { background-color: #e9ecef; }\n',
    '    .guide-export-content table p { margin-bottom: 0.5rem; }\n',
    '    .guide-export-content table p:last-child { margin-bottom: 0; }\n',
           guide_modal_styles,
    '    .guide-layout { display: block; }\n',
    '    .guide-toc { display: none; }\n',
    '    .guide-content { padding: 0; overflow: visible; }\n',
    '  </style>\n',
    '</head>\n',
    '<body>\n',
    '  <div class="guide-export-shell">\n',
    '    <nav class="guide-export-toc" id="guideExportToc">\n',
    '      <button type="button" class="guide-export-toggle" id="guideExportToggle" aria-expanded="true" aria-controls="guideExportContent">Collapse contents</button>\n',
    '      <div class="toc-title">Contents</div>\n',
    '    </nav>\n',
    '    <main class="guide-export-content guide-export" id="guideExportContent">\n',
         guide_body,
    '    </main>\n',
    '  </div>\n',
    '  <script>\n',
    '    (function() {\n',
    '      var content = document.getElementById("guideExportContent");\n',
    '      var toc = document.getElementById("guideExportToc");\n',
    '      var toggle = document.getElementById("guideExportToggle");\n',
    '      if (!content || !toc || !toggle) return;\n',
    '      var headings = content.querySelectorAll("h1");\n',
    '      headings.forEach(function(h, i) {\n',
    '        if (!h.id) h.id = "guide-export-heading-" + i;\n',
    '        var a = document.createElement("a");\n',
    '        a.href = "#" + h.id;\n',
    '        a.textContent = h.textContent;\n',
    '        a.className = "toc-" + h.tagName.toLowerCase();\n',
    '        a.addEventListener("click", function(e) {\n',
    '          e.preventDefault();\n',
    '          h.scrollIntoView({behavior: "smooth", block: "start"});\n',
    '        });\n',
    '        toc.appendChild(a);\n',
    '      });\n',
    '      toggle.addEventListener("click", function() {\n',
    '        var collapsed = content.parentElement.classList.toggle("is-collapsed");\n',
    '        toggle.setAttribute("aria-expanded", collapsed ? "false" : "true");\n',
    '        toggle.textContent = collapsed ? "Expand contents" : "Collapse contents";\n',
    '      });\n',
    '    })();\n',
    '  </script>\n',
    '</body>\n',
    '</html>\n'
  )
}

# APP_VERSION is defined in global.R

deployment_mode <- tolower(Sys.getenv("VARG_DEPLOYMENT_MODE", unset = "local"))
is_hosted_deployment <- identical(deployment_mode, "hosted")

# Define Theme
varg_theme <- bslib::bs_theme(
  version = 5,
  preset = "flatly", # Clean, modern look
  primary = "#2C3E50",
  secondary = "#95a5a6",
  success = "#18bc9c",
  info = "#3498db",
  warning = "#f39c12",
  danger = "#e74c3c",
  base_font = "Arial",
  heading_font = "Arial"
)

ui <- bslib::page_sidebar(
  fillable = FALSE,
  title = tags$div(
    style = "display: flex; align-items: center; gap: 10px; width: 100%;",
    tags$img(src = "varg-tools-symbol.svg", height = "48px", style = "vertical-align: middle; background: white; border-radius: 50%; padding: 3px;"),
    paste0("VARG-Tools ", APP_VERSION),
    tags$div(style = "flex-grow: 1;"),
    actionButton("header_user_guide", 
      label = tagList(icon("book-open"), " User Guide"),
      class = "btn btn-outline-light",
      style = "font-size: 1.1rem; padding: 8px 20px; font-weight: 500;"
    )
  ),
  theme = varg_theme,
  shinyjs::useShinyjs(),
  sidebar = bslib::sidebar(
    width = 280,
    open = "desktop",  # Collapsible on desktop, closed on mobile

    # Project Manager and Close Application at top
    div(
      class = "mb-3",
      h5("Save / Load", class = "mb-2 text-muted", style = "font-size: 0.9rem;"),
      mod_project_manager_ui("project_mgr"),
      tags$hr(class = "my-2"),
      div(
        class = "mt-3 pt-2",
        actionButton(
          "shutdown_app",
          "Close Application",
          icon = icon("power-off"),
          class = "btn-outline-danger btn-sm w-100"
        ),
        p(class = "text-muted text-center mt-1", style = "font-size: 0.75rem;", 
          "Closes the Shiny session")
      )
    ),
    tags$hr(),

    # Navigation Pills (nav only, no content)
    h5("Workflow Navigation", class = "mb-3 text-muted", style = "font-size: 0.9rem;"),
    bslib::navset_pill_list(
      id = "workflow_nav",
      widths = c(12, 12),
      bslib::nav_panel("Home", icon = icon("home"), value = "project"),
      bslib::nav_panel("1. Processing", icon = icon("cogs"), value = "processing"),
      bslib::nav_panel("2. Visualization", icon = icon("chart-bar"), value = "viz"),
      bslib::nav_panel("3. Chronology", icon = icon("clock"), value = "chron"),
      # bslib::nav_panel("Contribute", icon = icon("upload"), value = "contribute"),  # Hidden — see note at top of file
      bslib::nav_panel("Feedback", icon = icon("comment-dots"), value = "feedback"),
      bslib::nav_panel("About", icon = icon("info-circle"), value = "about")
    )
  ),

  # Main content area - conditionally show based on selected nav
  conditionalPanel(
    condition = "input.workflow_nav == 'project'",
    layout_columns(
      col_widths = c(12, 12, 12),

      # Hero Section
      varg_card(
        div(
          class = "text-center",
          div(
            style = "display: flex; align-items: center; justify-content: center; gap: 20px; margin-bottom: 15px;",
            tags$img(src = "varg-tools-symbol.svg", height = "90px", alt = "VARG-Tools Symbol"),
            tags$img(src = "VARG-Tools_wordlogo.svg", height = "60px", alt = "VARG-Tools")
          ),
          p("A comprehensive toolkit for volcanic ash and glass geochemical analysis, visualization, and chronology",
            style = "font-size: 1.1rem; color: #666; margin-bottom: 10px;"
          )
        )
      ),

      # Quick Start Workflow
      varg_card(
        title = "Quick Start: 3-Module Workflow",
        p(class = "text-muted mb-3", "VARG-Tools is modular. Work through the full pipeline or jump to any step. Each module can accept its own data independently."),
        layout_columns(
          col_widths = c(4, 4, 4),
          div(
            class = "text-center",
            div(style = "font-size: 2rem; color: #18bc9c; margin-bottom: 10px;", icon("cogs")),
            h5("1. Processing"),
            p(class = "small text-start", style = "color: #555;",
              tags$strong("What you do:"), " Import a CSV/Excel of major-oxide analyses, apply compositional transformations (ILR/pivot coordinates), run Gaussian Mixture Modeling (GMM) to find clusters, and use UMAP to visualize your data in 2D. You can also project it onto the VARG26 reference database.", tags$br(), tags$br(),
              tags$strong("What you need:"), " A spreadsheet with rows = analyses and columns = oxide wt.% (SiO\u2082, TiO\u2082, Al\u2082O\u2083, etc.) plus optional metadata (Depth, Age, CoreID).", tags$br(), tags$br(),
              tags$strong("What you get:"), " A labeled dataset with cluster IDs, UMAP coordinates, and named population assignments ready for visualization or export."
            )
          ),
          div(
            class = "text-center",
            div(style = "font-size: 2rem; color: #3498db; margin-bottom: 10px;", icon("chart-bar")),
            h5("2. Visualization"),
            p(class = "small text-start", style = "color: #555;",
              tags$strong("What you do:"), " Find the most diagnostic oxide pairs for separating your populations (KNN pair-finder), create publication-quality scatter plots with VARG26 reference overlays, and correlate depth-series between cores using visual tie-pointing.", tags$br(), tags$br(),
              tags$strong("What you need:"), " Processed data from Step 1, or any CSV/Excel with numeric analyte columns and a grouping variable.", tags$br(), tags$br(),
              tags$strong("What you get:"), " Exportable figures (PNG/PDF/SVG), variable pair rankings, stratigraphic tie-point tables, and warped depth data."
            )
          ),
          div(
            class = "text-center",
            div(style = "font-size: 2rem; color: #e74c3c; margin-bottom: 10px;", icon("clock")),
            h5("3. Chronology"),
            p(class = "small text-start", style = "color: #555;",
              tags$strong("What you do:"), " Generate OxCal Bayesian age-modeling code from your radiocarbon dates and stratigraphic constraints, including P_Sequence (age-depth), Phase models, and linked multi-site models using tie points from Step 2.", tags$br(), tags$br(),
              tags$strong("What you need:"), " A spreadsheet with sample names, ages, uncertainties, depths, and date types (e.g., R_Date, Tephra). Optionally, tie-point JSON from the Visualization module.", tags$br(), tags$br(),
              tags$strong("What you get:"), " Ready-to-run OxCal code that you paste into ", tags$a("c14.arch.ox.ac.uk", href = "https://c14.arch.ox.ac.uk/oxcal.html", target = "_blank"), " for Bayesian calibration and age-depth modeling."
            )
          )
        ),
        div(
          class = "text-center mt-3",
          downloadButton("download_template_home", "Download Data Template", class = "btn-info", style = "margin-right: 10px;", icon = icon("file-download")),
          actionButton("goto_processing", "Start Processing", class = "btn-success", icon = icon("play"))
        )
      ),

      # Resources & Help
      layout_columns(
        col_widths = c(6, 6),
        varg_card(
          title = tagList(icon("book-open", class = "text-primary"), " Help & Resources"),
          p(class = "small mb-3", "Every input in VARG-Tools has a ", icon("question-circle", class = "text-primary"), " tooltip explaining what it does and when to use it. For a full walkthrough:"),
          div(
            class = "d-flex flex-column gap-2",
            actionButton("show_user_guide",
                         label = tagList(icon("book-open"), " View User Guide"),
                         class = "btn btn-outline-primary w-100"),
            downloadButton("download_beta_data",
                           label = tagList(icon("database"), " Download Test Data"),
                           class = "btn btn-outline-info w-100")
          )
        ),
        varg_card(
          title = tagList(icon("save", class = "text-success"), " Session Tips"),
          tags$ul(
            class = "small mb-0", style = "padding-left: 18px;",
            tags$li("Use the ", tags$strong("Save / Load"), " buttons in the sidebar to preserve your work between sessions."),
            if (is_hosted_deployment) {
              tags$li("Hosted sessions close after 1 hour without user activity. Active analyses remain open, and the 1-hour idle allowance restarts when an analysis finishes.")
            } else {
              tags$li("This local session continues until you close VARG-Tools or its console window.")
            },
            if (is_hosted_deployment) {
              tags$li("Close the app when finished (sidebar button) to free resources for others.")
            } else {
              tags$li("Use Close Application in the sidebar to stop the local session cleanly.")
            },
            tags$li("All long-running tasks (GMM, Variable Pair Finder) can be cancelled.")
          )
        )
      ),
      
      # Footer
      div(
        class = "text-center text-muted small mt-1 py-2",
        tags$span(style = "font-size: 0.75rem;", paste0("VARG-Tools ", APP_VERSION))
      )
    )
  ),
  conditionalPanel(
    condition = "input.workflow_nav == 'processing'",
    mod_processing_ui("proc")
  ),
  conditionalPanel(
    condition = "input.workflow_nav == 'viz'",
    mod_visualization_ui("viz")
  ),
  conditionalPanel(
    condition = "input.workflow_nav == 'chron'",
    mod_chronology_ui("chron")
  ),
  # conditionalPanel(
  #   condition = "input.workflow_nav == 'contribute'",
  #   mod_contribution_ui("contribute")
  # ),
  conditionalPanel(
    condition = "input.workflow_nav == 'feedback'",
    varg_card(
      title = "Help Improve VARG-Tools",
      div(
        class = "text-center py-4",
        h4(icon("comments"), "Share Your Experience"),
        p("This app is in active development. Your feedback guides improvements."),
        actionButton(
          "feedbackLink", "Submit Feedback",
          icon = icon("external-link-alt"),
          class = "btn-primary btn-lg my-3",
          onclick = "window.open('https://forms.gle/Phe3zW7HaEa2ncPc8', '_blank')"
        )
      )
    )
  ),
  conditionalPanel(
    condition = "input.workflow_nav == 'about'",
    layout_columns(
      col_widths = c(12),
      varg_card(
        title = "About VARG-Tools",
        div(
          class = "py-3",
          h4(icon("flask"), " VARG-Tools ", tags$span(class = "badge bg-primary", APP_VERSION)),
          p(class = "lead", "An Integrated Tephrochronology Workflow from UMAP Visualization to Linked Bayesian Chronologies"),
          tags$hr(),
          
          h5(icon("user"), " Author"),
          p("Matthew S. M. Bolton"),
          
          h5(icon("book"), " Citation"),
          p("If you use VARG-Tools in your research, please cite:"),
          div(class = "bg-light p-3 border rounded mb-3",
            p(class = "mb-0", style = "font-family: monospace;",
              "Bolton, M.S.M. and Jensen, B.J.L. (in prep.) An Integrated Tephrochronology Workflow from UMAP Visualization to Linked Bayesian Chronologies with VARG-Tools."
            )
          ),
          
          h5(icon("database"), " VARG Reference Data"),
          p("The embedded reference data (VARG26) includes tephra labels and UMAP coordinates for visualization. The raw geochemical data underlying these embeddings is currently unpublished."),
          
          h5(icon("quote-left"), " Recommended Citations"),
          p("Similar to OxCal, please cite the appropriate methods when using specific features:"),
          tags$ul(
            tags$li(tags$strong("UMAP embeddings:"), " Cite the VARG-Tools paper above"),
            tags$li(tags$strong("GMM clustering:"), " Cite the VARG-Tools paper above"),
            tags$li(tags$strong("OxCal chronologies:"), " Cite Bronk Ramsey (2009) and the relevant calibration curve"),
            tags$li(tags$strong("Deposition models:"), " Cite Bronk Ramsey (2008)")
          ),
          
          tags$hr(),
          h5(icon("balance-scale"), " License"),
          p("The VARG-Tools application code is open source. The embedded VARG26 coordinates are provided for visualization purposes."),
          
          p(class = "text-muted text-center mt-4", 
            paste0("VARG-Tools ", APP_VERSION, " | © ", format(Sys.Date(), "%Y"), " Matthew S. M. Bolton")
          )
        )
      )
    )
  ),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(src = "custom.js")
  ),
  # Splash Screen (keep existing)
  tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function() {
      var modalEl = document.getElementById('welcomeModal');
      if (modalEl && window.bootstrap && window.bootstrap.Modal) {
        var welcomeModal = bootstrap.Modal.getOrCreateInstance(modalEl, {
          backdrop: true,
          focus: true,
          keyboard: true
        });

        modalEl.addEventListener('hidden.bs.modal', function() {
          document.querySelectorAll('.modal-backdrop').forEach(function(backdrop) {
            backdrop.remove();
          });
          document.body.classList.remove('modal-open');
          document.body.style.removeProperty('overflow');
          document.body.style.removeProperty('padding-right');
        });

        welcomeModal.show();
      } else if (window.jQuery) {
        $('#welcomeModal').modal('show');
      }

    });
    
    // Release hosted sessions through ShinyProxy's own stop control.
    Shiny.addCustomMessageHandler('closeWindow', function(message) {
      var closedHtml = '<div style=\"box-sizing:border-box;display:flex;justify-content:center;align-items:center;min-height:100vh;padding:24px;flex-direction:column;font-family:Arial,sans-serif;background:#f5f7f8;color:#22313f;text-align:center;\"><h2 style=\"margin:0 0 12px;\">VARG-Tools session closed</h2><p style=\"margin:0 0 18px;max-width:520px;line-height:1.45;\">Your VARG-Tools session is closing. You can close this browser tab.</p></div>';

      function showClosedMessage(targetWindow) {
        targetWindow.document.body.innerHTML = closedHtml;
      }

      var pathParts = window.location.pathname.split('/');
      var isHosted = pathParts[1] === 'app_proxy' && Boolean(pathParts[2]);
      if (!isHosted) {
        showClosedMessage(window);
        Shiny.setInputValue('local_close_ready', Date.now(), {priority: 'event'});
        return;
      }

      var stopHostedSession = window.VARGTools && window.VARGTools.stopHostedSession;
      if (typeof stopHostedSession !== 'function') {
        Shiny.setInputValue(
          'hosted_close_failed',
          'Hosted session release helper was not found',
          {priority: 'event'}
        );
        return;
      }

      stopHostedSession(window.top).catch(function(e) {
        Shiny.setInputValue(
          'hosted_close_failed',
          e && e.message ? e.message : 'Unknown release error',
          {priority: 'event'}
        );
      });
    });
  ")),

  # Welcome Modal (keep existing structure but ensure it works with BS5)
  tags$div(
    class = "modal fade", id = "welcomeModal", tabindex = "-1",
    tags$div(
      class = "modal-dialog modal-lg",
      tags$div(
        class = "modal-content",
        tags$div(
          class = "modal-header bg-primary text-white",
          tags$div(
            style = "display: flex; align-items: center; gap: 12px;",
            tags$img(src = "varg-tools-symbol.svg", height = "36px", style = "filter: brightness(0) invert(1);"),
            tags$h5(class = "modal-title mb-0", "Welcome to VARG-Tools")
          ),
          tags$button(type = "button", class = "btn-close btn-close-white", `data-bs-dismiss` = "modal")
        ),
        tags$div(
          class = "modal-body",
          div(class = "text-center mb-3",
            tags$img(src = "VARG-Tools_wordlogo.svg", height = "50px", alt = "VARG-Tools")
          ),
          h5(icon("compass"), " What you can do"),
          p("Welcome! VARG-Tools is a modular toolkit for volcanic glass geochemical analysis, from raw data to age models. You can work through the complete Processing \u2192 Visualization \u2192 Chronology pipeline, or jump directly to the module you need."),
          tags$ul(
            tags$li("Processing: import, preprocess, cluster, and project your data."),
            tags$li("Visualization: build variable-pair diagnostics, custom plots, and stratigraphic correlations."),
            tags$li("Chronology: generate OxCal age-depth and phase-model code for dating workflows.")
          ),
          tags$hr(),
          h5(icon("book-open"), " Getting Help"),
          p("The ", tags$strong("User Guide"), " is available anytime from the header bar at the top of the page, or from the home page. Hover or click the ",
            tags$i(class = "fa fa-question-circle", style = "color: #3c8dbc;"),
            " icons throughout the app for contextual explanations of each feature."),
          tags$hr(),
          h5(icon("comment-dots"), " Feedback Welcome"),
          p("Researchers and other users are welcome to explore VARG-Tools. Please report problems or suggestions through the 'Feedback' tab."),
          tags$hr(),
          if (is_hosted_deployment) {
            tagList(
              h5(icon("server"), " Hosting Notes"),
              p("Hosted sessions close after 1 hour without user activity. Active analyses remain open, and the 1-hour idle allowance restarts when an analysis finishes. Use the Save/Load buttons in the sidebar to preserve important work."),
              tags$div(
                class = "alert alert-info mt-2 mb-0",
                icon("lightbulb"),
                tags$strong(" Save resources: "),
                "When you're done, please click the ",
                tags$strong("'Close Application'"),
                " button in the sidebar. This frees up active hours for other users. Thank you!"
              )
            )
          } else {
            tagList(
              h5(icon("desktop"), " Local Session"),
              p("This copy runs on your computer and has no hosted one-hour limit. Use Save / Load to preserve work between sessions."),
              tags$div(
                class = "alert alert-info mt-2 mb-0",
                icon("lightbulb"),
                tags$strong(" When finished: "),
                "Click ", tags$strong("'Close Application'"),
                " in the sidebar or close the VARG-Tools console window."
              )
            )
          }
        ),
        tags$div(
          class = "modal-footer",
          tags$button(type = "button", class = "btn btn-primary", `data-bs-dismiss` = "modal", "Get Started")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Initialize global reactive values for project state
  # We need to extract this from mod_processing to make it shared
  # For now, we'll access the module's return value, but ideally we should lift the state up.
  # However, mod_processing currently owns the data.
  # Strategy: mod_processing returns a reactive of its internal data.
  # But to LOAD data, we need to be able to push data INTO mod_processing.
  # This requires refactoring mod_processing to accept `initial_data` or expose a `update_data` function.

  # Shared reactive state contract (`global_rv`)
  #
  # Field ownership and types:
  # - `data` (data.frame | NULL) — owned by Processing module; read by Visualization/Project Manager
  # - `mclust_result` (list | NULL) — owned by Processing module (GMM fit object)
  # - `umap_mode_ran` (character | NULL) — owned by Processing module ("new", "pretrained", "loaded")
  # - `data_stale` (logical) — owned by Processing module; indicates settings changed since last transform
  # - `pipeline_config` (list | NULL) — owned by Processing module; stores preprocessing/GMM/UMAP settings for save/load
  # - `user_umap_model` (list | NULL) — owned by Processing module; trained model for projection reuse
  # - `source_filename` (character | NULL) — owned by Processing module/Project Manager for provenance
  # - `original_cols` (character vector | NULL) — owned by Processing module; original upload columns for export/report ordering
  # - `restore_trigger` (integer scalar) — owned by Project Manager; increments on load to trigger UI restoration
  # - `data_generation` (integer scalar) — increments on source replacement or successful data mutation
  # - `heavy_job` (character | NULL) — cross-module session lock for memory-heavy background work
  #
  # Cross-module rule: modules may read shared fields, but only the owning module should mutate them.
  global_rv <- reactiveValues(
    data = NULL,
    mclust_result = NULL,
    umap_mode_ran = NULL,
    data_stale = FALSE,
    pipeline_config = NULL,
    user_umap_model = NULL,
    source_filename = NULL,
    original_cols = NULL,
    restore_trigger = 0, # Counter to trigger UI restoration only on project load
    data_generation = 0L,
    heavy_job = NULL
  )

  observe({
    heavy_job <- global_rv$heavy_job
    heavy_active <- !is.null(heavy_job) && length(heavy_job) > 0 &&
      !is.na(heavy_job[[1]]) && nzchar(heavy_job[[1]])
    session$sendCustomMessage(
      "sessionHeavyJobState",
      list(
        active = heavy_active,
        label = if (heavy_active) as.character(heavy_job[[1]]) else NULL
      )
    )
  })

  # Pass this global_rv to modules
  # Note: This requires updating mod_processing to use global_rv instead of its internal rv
  # OR we sync them. Syncing is easier for a quick refactor.

  # Load Processing Module
  # We pass global_rv to it so it can read/write shared state
  processed_data <- mod_processing_server("proc", global_rv)

  # Load Visualization Module (can accept processed data)
  viz_output <- mod_visualization_server("viz", processed_data, global_rv = global_rv)

  # Create pop_styles reactive wrapper that bridges viz_output to project manager
  pop_styles <- reactiveVal(NULL)
  
  # Sync pop_styles from viz_output (for saving)
  observe({
    styles <- viz_output$pop_styles()
    if (!is.null(styles)) {
      pop_styles(styles)
    }
  })
  
  # When pop_styles changes (from load), push to viz module
  observeEvent(pop_styles(), {
    styles <- pop_styles()
    if (!is.null(styles) && !is.null(viz_output$set_pop_styles)) {
      viz_output$set_pop_styles(styles)
    }
  }, ignoreInit = TRUE)
  
  # Load Project Manager with pop_styles for save/load
  mod_project_manager_server("project_mgr", global_rv, pop_styles = pop_styles)

  # Load Chronology Module
  mod_chronology_server("chron", tiepoints_reactive = viz_output$tiepoints, global_rv = global_rv)

  # Load Contribution Module — hidden; see note at top of file
  # mod_contribution_server("contribute", app_version = APP_VERSION)

  idle_warning_id <- reactiveVal(NULL)

  observeEvent(input$session_idle_warning, {
    if (!is.null(idle_warning_id())) {
      removeNotification(idle_warning_id())
    }
    idle_warning_id(showNotification(
      "This hosted session has been inactive for 50 minutes and will close after 10 more inactive minutes. Move the pointer, type, or interact with the app to continue.",
      type = "warning",
      duration = NULL
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$session_idle_resumed, {
    if (!is.null(idle_warning_id())) {
      removeNotification(idle_warning_id())
      idle_warning_id(NULL)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$session_idle_timeout, {
    if (is_hosted_deployment) {
      session$sendCustomMessage("closeWindow", list(reason = "idle_timeout"))
    }
  }, ignoreInit = TRUE)

  # Shutdown button handler
  observeEvent(input$shutdown_app, {
    showModal(modalDialog(
      title = tagList(icon("exclamation-triangle"), " Close VARG-Tools?"),
      "Are you sure you want to close the application? Any unsaved work will be lost.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_shutdown", "Close App", class = "btn-danger", icon = icon("power-off"))
      ),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$confirm_shutdown, {
    # Remove the confirmation modal first
    removeModal()
    
    # Hosted sessions are released through ShinyProxy; local sessions report
    # back to R before stopApp() is called.
    session$sendCustomMessage("closeWindow", list())
  })

  observeEvent(input$local_close_ready, {
    stopApp()
  }, ignoreInit = TRUE)

  observeEvent(input$hosted_close_failed, {
    showNotification(
      "The hosted session could not be released automatically. Save your work, then close this browser tab and return to the VARG-Tools launch page.",
      type = "error",
      duration = NULL
    )
  }, ignoreInit = TRUE)

  # Handle Start Processing button click
  observeEvent(input$goto_processing, {
    bslib::nav_select("workflow_nav", selected = "processing")
  })

  # Header User Guide button handler (mirrors show_user_guide)
  observeEvent(input$header_user_guide, {
    if (file.exists(guide_file)) {
      show_guide_modal(
        download_id = "download_user_guide_header",
        content_id = "guideContentHeader",
        toc_id = "guideTocHeader",
        heading_prefix = "guide-hdr-h-"
      )
    } else {
      showNotification("User guide file not found.", type = "error", duration = 5)
    }
  })

  # Template download handler for Welcome Page
  output$download_template_home <- downloadHandler(
    filename = function() {
      "VARG_Tools_ProcessViz_Template.xlsx"
    },
    content = function(file) {
      sheets <- generate_template("geochem")
      writexl::write_xlsx(sheets, path = file)
    }
  )
  outputOptions(output, "download_template_home", suspendWhenHidden = FALSE)
  
  # Test data download handler
  output$download_beta_data <- downloadHandler(
    filename = function() {
      "Data_for_Testing.zip"
    },
    content = function(file) {
      zip_file <- "Data for Beta Testers.zip"
      if (file.exists(zip_file)) {
        file.copy(zip_file, file)
      } else {
        # If zip doesn't exist yet, create a placeholder message
        showNotification("Test data folder is being prepared. Please check back soon!", 
                         type = "warning", duration = 5)
        # Create an empty placeholder file
        writeLines("Test data coming soon! This folder will contain example datasets for testing.", file)
      }
    }
  )
  outputOptions(output, "download_beta_data", suspendWhenHidden = FALSE)
  
  # User Guide handlers
  observeEvent(input$show_user_guide, {
    if (file.exists(guide_file)) {
      show_guide_modal(
        download_id = "download_user_guide",
        content_id = "guideContentMain",
        toc_id = "guideTocMain",
        heading_prefix = "guide-hdr-m-"
      )
    } else {
      showNotification("User guide file not found. Please check docs/VARG-Tools_User_Guide.md", 
                       type = "error", duration = 5)
    }
  })
  
  # Download User Guide as standalone HTML
  output$download_user_guide <- downloadHandler(
    filename = function() {
      paste0("VARG-Tools_User_Guide_", Sys.Date(), ".html")
    },
    content = function(file) {
      if (file.exists(guide_file)) {
        writeLines(build_guide_export_document(), file)
      }
    }
  )
  outputOptions(output, "download_user_guide", suspendWhenHidden = FALSE)

  output$download_user_guide_header <- downloadHandler(
    filename = function() {
      paste0("VARG-Tools_User_Guide_", Sys.Date(), ".html")
    },
    content = function(file) {
      if (file.exists(guide_file)) {
        writeLines(build_guide_export_document(), file)
      }
    }
  )
  outputOptions(output, "download_user_guide_header", suspendWhenHidden = FALSE)
  
  # Download handlers removed - guide is available in-app only for now
}

# Run the application
shinyApp(ui = ui, server = server)
