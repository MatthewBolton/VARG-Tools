args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)

module_text <- paste(readLines(file.path(app_dir, "R", "modules", "mod_chronology.R"), warn = FALSE), collapse = "\n")
script_text <- paste(readLines(file.path(app_dir, "www", "custom.js"), warn = FALSE), collapse = "\n")

stopifnot(
  !grepl("clipr::write_clip", module_text, fixed = TRUE),
  lengths(regmatches(module_text, gregexpr("varg-copy-to-clipboard", module_text, fixed = TRUE))) == 3L,
  lengths(regmatches(module_text, gregexpr("data-copy-status-input", module_text, fixed = TRUE))) == 3L,
  grepl("copyTextToClipboard", script_text, fixed = TRUE),
  grepl("navigator.clipboard.writeText", script_text, fixed = TRUE),
  grepl("document.execCommand(\"copy\")", script_text, fixed = TRUE)
)

cat("chronology_clipboard_contract_ok\n")
