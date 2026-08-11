heavy_job_limiter_config <- function() {
  slots <- suppressWarnings(as.integer(Sys.getenv("VARG_HEAVY_SLOTS", unset = "0")))
  if (is.na(slots) || slots < 1L) {
    slots <- 2L
  }

  cpu_tokens <- suppressWarnings(as.integer(Sys.getenv("VARG_HEAVY_CPU_TOKENS", unset = "0")))
  if (is.na(cpu_tokens) || cpu_tokens < 1L) {
    available <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (is.na(available) || available < 2L) available <- 2L
    worker_cap <- suppressWarnings(as.integer(Sys.getenv("VARG_MAX_CORES", unset = "12")))
    if (is.na(worker_cap) || worker_cap < 1L) worker_cap <- 12L
    cpu_tokens <- max(1L, min(available - 1L, worker_cap))
  }

  stale_minutes <- suppressWarnings(as.numeric(Sys.getenv("VARG_HEAVY_LOCK_STALE_MINUTES", unset = "180")))
  if (is.na(stale_minutes) || stale_minutes < 1) {
    stale_minutes <- 180
  }

  list(
    lock_dir = Sys.getenv(
      "VARG_HEAVY_LOCK_DIR",
      unset = file.path(tempdir(), "vargtools-heavy-locks")
    ),
    slots = slots,
    cpu_tokens = cpu_tokens,
    stale_seconds = stale_minutes * 60
  )
}

heavy_job_limiter_token <- function(label) {
  paste(
    Sys.info()[["nodename"]],
    Sys.getpid(),
    gsub("[^A-Za-z0-9_.-]+", "_", label),
    sprintf("%.0f", as.numeric(Sys.time())),
    sep = "-"
  )
}

heavy_job_limiter_prune_stale <- function(config) {
  if (!dir.exists(config$lock_dir)) {
    return(invisible(NULL))
  }

  lock_dirs <- list.dirs(config$lock_dir, full.names = TRUE, recursive = FALSE)
  if (length(lock_dirs) == 0) {
    return(invisible(NULL))
  }

  now <- Sys.time()
  for (lock_path in lock_dirs) {
    if (basename(lock_path) == ".allocator.lock") next
    age <- difftime(now, file.info(lock_path)$mtime, units = "secs")
    if (!is.na(age) && as.numeric(age) > config$stale_seconds) {
      unlink(lock_path, recursive = TRUE, force = TRUE)
    }
  }
  invisible(NULL)
}

heavy_job_limiter_acquire_mutex <- function(config, wait_seconds = 5, stale_seconds = 30) {
  mutex_path <- file.path(config$lock_dir, ".allocator.lock")
  deadline <- Sys.time() + wait_seconds

  repeat {
    if (dir.create(mutex_path, showWarnings = FALSE)) {
      writeLines(
        c(
          paste0("pid=", Sys.getpid()),
          paste0("host=", Sys.info()[["nodename"]]),
          paste0("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
        ),
        file.path(mutex_path, "owner.txt")
      )
      return(mutex_path)
    }

    if (dir.exists(mutex_path)) {
      age <- difftime(Sys.time(), file.info(mutex_path)$mtime, units = "secs")
      if (!is.na(age) && as.numeric(age) > stale_seconds) {
        unlink(mutex_path, recursive = TRUE, force = TRUE)
        next
      }
    }

    if (Sys.time() >= deadline) return(NULL)
    Sys.sleep(0.02)
  }
}

heavy_job_limiter_release_mutex <- function(mutex_path) {
  if (!is.null(mutex_path) && dir.exists(mutex_path)) {
    unlink(mutex_path, recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

heavy_job_limiter_write_owner <- function(path, token, label, kind, index) {
  writeLines(
    c(
      paste0("token=", token),
      paste0("label=", label),
      paste0("kind=", kind),
      paste0("index=", index),
      paste0("pid=", Sys.getpid()),
      paste0("host=", Sys.info()[["nodename"]]),
      paste0("started=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    ),
    file.path(path, "owner.txt")
  )
}

heavy_job_limiter_acquire <- function(label, min_tokens = 1L, max_tokens = min_tokens) {
  config <- heavy_job_limiter_config()

  min_tokens <- suppressWarnings(as.integer(min_tokens))
  max_tokens <- suppressWarnings(as.integer(max_tokens))
  if (length(min_tokens) != 1L || is.na(min_tokens) || min_tokens < 1L) min_tokens <- 1L
  if (length(max_tokens) != 1L || is.na(max_tokens) || max_tokens < min_tokens) max_tokens <- min_tokens
  max_tokens <- min(max_tokens, config$cpu_tokens)

  if (!dir.exists(config$lock_dir)) {
    dir.create(config$lock_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!dir.exists(config$lock_dir)) {
    return(list(
      acquired = FALSE,
      message = "The server-wide heavy-analysis limiter is unavailable. Try again shortly.",
      lock = NULL,
      workers = 0L
    ))
  }

  mutex_path <- heavy_job_limiter_acquire_mutex(config)
  if (is.null(mutex_path)) {
    return(list(
      acquired = FALSE,
      message = "The server-wide resource allocator is busy. Try this analysis again shortly.",
      lock = NULL,
      workers = 0L
    ))
  }
  on.exit(heavy_job_limiter_release_mutex(mutex_path), add = TRUE)

  heavy_job_limiter_prune_stale(config)

  token <- heavy_job_limiter_token(label)
  slot_path <- NULL
  for (slot in seq_len(config$slots)) {
    candidate <- file.path(config$lock_dir, sprintf("slot_%02d.lock", slot))
    if (dir.create(candidate, showWarnings = FALSE)) {
      slot_path <- candidate
      heavy_job_limiter_write_owner(slot_path, token, label, "job", slot)
      break
    }
  }

  if (is.null(slot_path)) {
    return(list(
      acquired = FALSE,
      message = paste0(
        "The server is already running ",
        config$slots,
        " heavy analyses. You can keep using the app, but start this analysis again shortly."
      ),
      lock = NULL,
      workers = 0L
    ))
  }

  cpu_paths <- character(0)
  for (cpu in seq_len(config$cpu_tokens)) {
    if (length(cpu_paths) >= max_tokens) break
    cpu_path <- file.path(config$lock_dir, sprintf("cpu_%02d.lock", cpu))
    if (dir.create(cpu_path, showWarnings = FALSE)) {
      heavy_job_limiter_write_owner(cpu_path, token, label, "cpu", cpu)
      cpu_paths <- c(cpu_paths, cpu_path)
    }
  }

  if (length(cpu_paths) < min_tokens) {
    unlink(cpu_paths, recursive = TRUE, force = TRUE)
    unlink(slot_path, recursive = TRUE, force = TRUE)
    return(list(
      acquired = FALSE,
      message = paste0(
        "The server does not currently have the ", min_tokens,
        " CPU workers needed to start this analysis. You can keep using the app and try again shortly."
      ),
      lock = NULL,
      workers = 0L
    ))
  }

  list(
    acquired = TRUE,
    message = NULL,
    lock = list(
      path = slot_path,
      cpu_paths = cpu_paths,
      token = token,
      label = label,
      workers = length(cpu_paths)
    ),
    workers = length(cpu_paths)
  )
}

heavy_job_limiter_workers <- function(lock, default = 1L) {
  workers <- if (is.null(lock)) NULL else lock$workers
  workers <- suppressWarnings(as.integer(workers))
  if (length(workers) == 0L || is.na(workers) || workers < 1L) {
    return(as.integer(default))
  }
  workers
}

heavy_job_limiter_release <- function(lock) {
  if (is.null(lock) || is.null(lock$path)) {
    return(invisible(FALSE))
  }

  cpu_paths <- if (is.null(lock$cpu_paths)) character(0) else lock$cpu_paths
  paths <- c(cpu_paths, lock$path)
  existing_paths <- paths[dir.exists(paths)]
  if (length(existing_paths) == 0L) {
    return(invisible(FALSE))
  }

  config <- heavy_job_limiter_config()
  mutex_path <- heavy_job_limiter_acquire_mutex(config)
  if (is.null(mutex_path)) return(invisible(FALSE))
  on.exit(heavy_job_limiter_release_mutex(mutex_path), add = TRUE)

  existing_paths <- paths[dir.exists(paths)]
  if (length(existing_paths) == 0L) return(invisible(FALSE))
  expected <- paste0("token=", lock$token)
  for (path in existing_paths) {
    owner_file <- file.path(path, "owner.txt")
    if (file.exists(owner_file)) {
      owner <- readLines(owner_file, warn = FALSE)
      if (!expected %in% owner) {
        return(invisible(FALSE))
      }
    }
  }

  unlink(cpu_paths, recursive = TRUE, force = TRUE)
  unlink(lock$path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}
