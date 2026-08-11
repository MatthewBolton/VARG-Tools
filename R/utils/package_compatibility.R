VARG_MINIMUM_PACKAGE_VERSIONS <- c(
  plotly = "4.12.0"
)

varg_package_compatibility_issues <- function(
    minimum_versions = VARG_MINIMUM_PACKAGE_VERSIONS,
    installed = NULL) {
  issues <- character(0)

  for (package in names(minimum_versions)) {
    minimum <- minimum_versions[[package]]
    if (is.null(installed)) {
      package_available <- requireNamespace(package, quietly = TRUE)
      installed_version <- if (package_available) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    } else {
      package_available <- package %in% rownames(installed)
      installed_version <- if (package_available) installed[package, "Version"] else NA_character_
    }

    if (!package_available) {
      issues <- c(issues, sprintf("%s is not installed", package))
      next
    }

    if (utils::compareVersion(installed_version, minimum) < 0L) {
      issues <- c(
        issues,
        sprintf("%s %s is installed; version %s or newer is required", package, installed_version, minimum)
      )
    }
  }

  issues
}

varg_assert_package_compatibility <- function(
    minimum_versions = VARG_MINIMUM_PACKAGE_VERSIONS,
    installed = NULL) {
  issues <- varg_package_compatibility_issues(minimum_versions, installed = installed)
  if (length(issues) > 0L) {
    stop(
      paste(
        "VARG-Tools package compatibility check failed:",
        paste0("- ", issues, collapse = "\n"),
        "Update the listed package(s), then restart VARG-Tools.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
