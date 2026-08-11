# R/modules/mod_project_manager.R
# Failure-atomic project save/load with versioned UMAP model persistence.

varg_normalize_project_pop_styles <- function(styles) {
    if (is.null(styles)) {
        return(list(ok = TRUE, styles = NULL, migrated = FALSE, errors = character(0)))
    }
    if (!is.list(styles) || is.null(names(styles))) {
        return(list(
            ok = FALSE,
            styles = NULL,
            migrated = FALSE,
            errors = "Project population-style metadata is invalid."
        ))
    }

    current_fields <- c("cp_group_color_overrides", "cp_group_shape_overrides")
    legacy_fields <- c("colors", "shapes")
    style_names <- names(styles)

    if (all(current_fields %in% style_names)) {
        normalized <- list(
            cp_group_color_overrides = styles[["cp_group_color_overrides"]],
            cp_group_shape_overrides = styles[["cp_group_shape_overrides"]],
            cp_current_palette = styles[["cp_current_palette"]]
        )
        return(list(ok = TRUE, styles = normalized, migrated = FALSE, errors = character(0)))
    }

    if (all(legacy_fields %in% style_names)) {
        palette <- if ("cp_current_palette" %in% style_names) {
            styles[["cp_current_palette"]]
        } else if ("palette" %in% style_names) {
            styles[["palette"]]
        } else {
            NULL
        }
        normalized <- list(
            cp_group_color_overrides = styles[["colors"]],
            cp_group_shape_overrides = styles[["shapes"]],
            cp_current_palette = palette
        )
        return(list(ok = TRUE, styles = normalized, migrated = TRUE, errors = character(0)))
    }

    list(
        ok = FALSE,
        styles = NULL,
        migrated = FALSE,
        errors = "Project population-style metadata is invalid."
    )
}

varg_validate_project_payload <- function(payload) {
    errors <- character(0)
    if (!is.list(payload)) {
        return(list(ok = FALSE, errors = "Project payload must be a list."))
    }
    if (!is.data.frame(payload$data) || nrow(payload$data) == 0) {
        errors <- c(errors, "Project payload must contain a non-empty data frame.")
    }
    if (!is.null(payload$original_cols)) {
        original_cols <- as.character(unlist(payload$original_cols, use.names = FALSE))
        if (length(original_cols) == 0 || anyNA(original_cols) || any(!nzchar(original_cols))) {
            errors <- c(errors, "Project original-column metadata is invalid.")
        }
    }
    if (!is.null(payload$pop_styles)) {
        style_result <- varg_normalize_project_pop_styles(payload$pop_styles)
        if (!style_result$ok) errors <- c(errors, style_result$errors)
    }

    has_user_umap <- isTRUE(payload$has_user_umap)
    if (has_user_umap) {
        contract <- validate_umap_projection_contract(payload$pipeline_config)
        if (!contract$valid) errors <- c(errors, contract$errors)
        if (is.null(payload$umap_mode_ran) || !identical(payload$umap_mode_ran, "new")) {
            errors <- c(errors, "Project UMAP model is inconsistent with its saved mode.")
        }
    }

    list(ok = length(errors) == 0, errors = unique(errors))
}

mod_project_manager_ui <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            class = "d-flex w-100",
            style = "gap: 6px;",
            downloadButton(ns("save_project"), "Save Project", class = "btn-default btn-sm flex-fill"),
            actionButton(ns("load_trigger"), "Load Project", icon = icon("folder-open"), class = "btn-default btn-sm flex-fill")
        )
    )
}

mod_project_manager_server <- function(id, global_rv, pop_styles = NULL) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        loaded_model_dir <- reactiveVal(NULL)

        unload_current_project_model <- function() {
            model <- isolate(global_rv$user_umap_model)
            config <- isolate(global_rv$pipeline_config)
            if (!is.null(model)) {
                unload_umap_model_collection(model, if (is.null(config)) NULL else config$umap)
            }
            old_dir <- isolate(loaded_model_dir())
            if (!is.null(old_dir) && dir.exists(old_dir)) {
                unlink(old_dir, recursive = TRUE, force = TRUE)
            }
            loaded_model_dir(NULL)
            invisible(NULL)
        }

        session$onSessionEnded(function() {
            try(unload_current_project_model(), silent = TRUE)
        })

        output$save_project <- downloadHandler(
            filename = function() {
                paste0("VARG_Project_", format(Sys.time(), "%Y%m%d_%H%M"), ".varg")
            },
            content = function(file) {
                showModal(modalDialog("Saving project...", footer = NULL))
                on.exit(removeModal())

                temp_dir <- tempfile("varg_bundle_")
                dir.create(temp_dir)
                on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

                current_model <- isolate(global_rv$user_umap_model)
                current_mode <- isolate(global_rv$umap_mode_ran)
                has_user_umap <- !is.null(current_model) && identical(current_mode, "new")
                model_files <- character(0)

                if (has_user_umap) {
                    contract <- validate_umap_projection_contract(isolate(global_rv$pipeline_config))
                    if (!contract$valid) stop(paste(contract$errors, collapse = "\n"))
                    model_files <- save_umap_model_collection(
                        current_model,
                        isolate(global_rv$pipeline_config$umap),
                        temp_dir
                    )
                }

                save_list <- list(
                    data = isolate(global_rv$data),
                    mclust_result = isolate(global_rv$mclust_result),
                    umap_mode_ran = current_mode,
                    pipeline_config = isolate(global_rv$pipeline_config),
                    source_filename = isolate(global_rv$source_filename),
                    original_cols = isolate(global_rv$original_cols),
                    data_stale = isTRUE(isolate(global_rv$data_stale)),
                    has_user_umap = has_user_umap,
                    metadata = list(
                        version = "2.0.0",
                        format = "zip_bundle",
                        date = Sys.time(),
                        app_version = APP_VERSION
                    )
                )
                if (!is.null(pop_styles)) {
                    style_result <- varg_normalize_project_pop_styles(isolate(pop_styles()))
                    if (!style_result$ok) stop(paste(style_result$errors, collapse = "\n"))
                    save_list$pop_styles <- style_result$styles
                }

                validation <- varg_validate_project_payload(save_list)
                if (!validation$ok) stop(paste(validation$errors, collapse = "\n"))

                saveRDS(save_list, file.path(temp_dir, "project.rds"))
                old_wd <- setwd(temp_dir)
                on.exit(setwd(old_wd), add = TRUE)
                zip::zip(
                    zipfile = file,
                    files = c("project.rds", model_files),
                    mode = "cherry-pick"
                )
                showNotification("Project saved successfully!", type = "message", duration = 3)
            }
        )

        observeEvent(input$load_trigger, {
            showModal(modalDialog(
                title = "Load Project",
                fileInput(ns("project_file"), "Choose .varg file", accept = ".varg"),
                footer = modalButton("Cancel"),
                easyClose = TRUE
            ))
        })

        invoke_shared_reset <- function(name) {
            callback <- isolate(global_rv[[name]])
            if (is.function(callback)) callback()
            invisible(NULL)
        }

        commit_candidate <- function(candidate, staged_models = NULL, retained_dir = NULL) {
            invoke_shared_reset("cancel_processing_jobs")
            invoke_shared_reset("reset_visualization_state")
            invoke_shared_reset("reset_chronology_state")
            unload_current_project_model()

            current_generation <- isolate(global_rv$data_generation)
            if (is.null(current_generation) || length(current_generation) != 1 || is.na(current_generation)) {
                current_generation <- 0L
            }

            global_rv$data <- candidate$data
            global_rv$mclust_result <- candidate$mclust_result
            global_rv$umap_mode_ran <- candidate$umap_mode_ran
            global_rv$pipeline_config <- candidate$pipeline_config
            global_rv$source_filename <- candidate$source_filename
            global_rv$original_cols <- candidate$original_cols
            global_rv$data_stale <- isTRUE(candidate$data_stale)
            global_rv$user_umap_model <- staged_models
            global_rv$data_generation <- as.integer(current_generation) + 1L
            if (!is.null(pop_styles) && !is.null(candidate$pop_styles)) pop_styles(candidate$pop_styles)
            loaded_model_dir(retained_dir)
            global_rv$restore_trigger <- isolate(global_rv$restore_trigger) + 1L
            invisible(NULL)
        }

        stage_zip_bundle <- function(file_path) {
            temp_dir <- tempfile("varg_load_")
            dir.create(temp_dir)
            staged_models <- NULL
            keep_dir <- FALSE
            on.exit({
                if (!keep_dir) {
                    if (!is.null(staged_models)) unload_umap_model_collection(staged_models)
                    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE, force = TRUE)
                }
            }, add = TRUE)

            zip::unzip(file_path, exdir = temp_dir)
            project_rds <- file.path(temp_dir, "project.rds")
            if (!file.exists(project_rds)) stop("Invalid project file: project.rds not found in bundle.")
            candidate <- readRDS(project_rds)
            validation <- varg_validate_project_payload(candidate)
            if (!validation$ok) stop(paste(validation$errors, collapse = "\n"))
            candidate$pop_styles <- varg_normalize_project_pop_styles(candidate$pop_styles)$styles

            if (is.null(candidate$original_cols)) candidate$original_cols <- names(candidate$data)
            if (is.null(candidate$data_stale)) candidate$data_stale <- FALSE
            if (isTRUE(candidate$has_user_umap)) {
                staged_models <- load_umap_model_collection(temp_dir, candidate$pipeline_config$umap)
                keep_dir <- TRUE
            }

            list(
                candidate = candidate,
                models = staged_models,
                retained_dir = if (keep_dir) temp_dir else NULL,
                cleanup_dir = if (keep_dir) NULL else temp_dir
            )
        }

        stage_legacy_rds <- function(file_path) {
            legacy <- readRDS(file_path)
            if (!is.data.frame(legacy$data) || nrow(legacy$data) == 0) {
                stop("Invalid legacy project file: no non-empty data frame was found.")
            }
            candidate <- list(
                data = legacy$data,
                mclust_result = legacy$mclust_result,
                umap_mode_ran = legacy$umap_mode_ran,
                pipeline_config = NULL,
                source_filename = NULL,
                original_cols = names(legacy$data),
                data_stale = FALSE,
                has_user_umap = FALSE,
                pop_styles = legacy$pop_styles
            )
            validation <- varg_validate_project_payload(candidate)
            if (!validation$ok) stop(paste(validation$errors, collapse = "\n"))
            candidate$pop_styles <- varg_normalize_project_pop_styles(candidate$pop_styles)$styles
            list(candidate = candidate, models = NULL, retained_dir = NULL, cleanup_dir = NULL)
        }

        observeEvent(input$project_file, {
            req(input$project_file)
            removeModal()
            staged <- NULL
            committed <- FALSE

            tryCatch({
                showModal(modalDialog("Loading project...", footer = NULL))
                file_path <- input$project_file$datapath
                con <- file(file_path, "rb")
                on.exit(try(close(con), silent = TRUE), add = TRUE)
                magic <- readBin(con, "raw", 2)
                close(con)
                is_zip <- identical(magic, as.raw(c(0x50, 0x4B)))

                staged <- if (is_zip) stage_zip_bundle(file_path) else stage_legacy_rds(file_path)
                commit_candidate(staged$candidate, staged$models, staged$retained_dir)
                committed <- TRUE
                if (!is.null(staged$cleanup_dir) && dir.exists(staged$cleanup_dir)) {
                    unlink(staged$cleanup_dir, recursive = TRUE, force = TRUE)
                }
                if (!is_zip) {
                    showNotification(
                        "Loaded legacy project data. Re-save it to use the current project format.",
                        type = "warning",
                        duration = 6
                    )
                }
                showNotification("Project loaded successfully!", type = "message")
            }, error = function(e) {
                if (!committed && !is.null(staged)) {
                    if (!is.null(staged$models)) unload_umap_model_collection(staged$models)
                    for (path in c(staged$retained_dir, staged$cleanup_dir)) {
                        if (!is.null(path) && dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
                    }
                }
                showNotification(paste("Error loading project:", e$message), type = "error", duration = 12)
            }, finally = {
                removeModal()
            })
        })
    })
}
