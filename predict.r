#!/usr/bin/env Rscript
## predict.r
## Usage: Rscript predict.r <model.rds> <historic_data.csv> <future_data.csv> <predictions.csv>
##
## Produces a CSV with columns: time_period, location, sample_0
## using posterior mean regression coefficients + spatial effects from the trained CARBayesST model.

library(CARBayesST)
library(splines)

source("model_helpers.R")

predict_chap <- function(model_fn, historic_data_fn, future_data_fn, predictions_fn) {

  # ── Load trained model and companions ────────────────────────────────────
  cat("Loading trained model from:", model_fn, "\n")
  model <- readRDS(model_fn)

  quartiles_fn <- paste0(model_fn, ".quartiles")
  if (!file.exists(quartiles_fn)) stop("Quartiles file not found: ", quartiles_fn)
  quartiles <- readRDS(quartiles_fn)
  cat("Loaded quartiles from:", quartiles_fn, "\n")

  metadata_fn <- "output/training_metadata.rds"
  if (!file.exists(metadata_fn)) stop("Training metadata not found: ", metadata_fn)
  metadata <- readRDS(metadata_fn)
  cat("Training time span:", metadata$start_time, "-", metadata$end_time,
      " | orgunits:", metadata$n_units, "\n")

  # ── Augment future data using historic context for lag computation ───────
  # Lags must be computed over the combined historic+future window so that the
  # first future time steps have valid lag1/lag2/lag3 values.
  cat("Loading historic data from:", historic_data_fn, "\n")
  historic <- augment_harmonized_data(historic_data_fn)

  cat("Loading and augmenting future data from:", future_data_fn, "\n")
  raw_future <- read.csv(future_data_fn, stringsAsFactors = FALSE)
  # Identify the future time periods before merging
  if ("time_period" %in% names(raw_future)) {
    future_times_raw <- unique(as.integer(gsub("-", "", raw_future$time_period)))
  } else {
    future_times_raw <- unique(as.integer(gsub("-", "", raw_future$time)))
  }

  # Combine historic + future, compute lags together, then extract future rows
  combined <- augment_harmonized_data(
    csv_path = future_data_fn  # will be replaced below via rbind trick
  )
  # Rebuild by row-binding the raw CSVs and re-running augment on combined file
  hist_raw  <- read.csv(historic_data_fn, stringsAsFactors = FALSE)
  fut_raw   <- read.csv(future_data_fn,   stringsAsFactors = FALSE)
  # Align columns: fill any column present in one but not the other with NA
  all_cols <- union(names(hist_raw), names(fut_raw))
  for (col in setdiff(all_cols, names(hist_raw))) hist_raw[[col]] <- NA
  for (col in setdiff(all_cols, names(fut_raw)))  fut_raw[[col]]  <- NA
  hist_raw     <- hist_raw[, all_cols, drop = FALSE]
  fut_raw      <- fut_raw[,  all_cols, drop = FALSE]
  combined_raw <- rbind(hist_raw, fut_raw)
  tmp_combined <- tempfile(fileext = ".csv")
  write.csv(combined_raw, tmp_combined, row.names = FALSE)
  combined <- augment_harmonized_data(tmp_combined)
  file.remove(tmp_combined)

  # Keep only the future time steps
  future <- combined[combined$time %in% future_times_raw, ]

  # Keep only orgunits present in training
  trained_units <- sort(metadata$selected_orgunits)
  future <- future[future$orgunitname %in% trained_units, ]
  if (nrow(future) == 0) stop("No future rows match trained orgunits.")
  future <- future[order(future$orgunitname, future$timeid), ]
  cat("Future rows after filtering:", nrow(future), "| orgunits:", length(unique(future$orgunitname)), "\n")

  # ── Extract posterior mean parameters ────────────────────────────────────
  sr <- model$summary.results

  # Beta: all rows except known scalar hyperparameters
  named_params <- c("tau2", "rho.S", "rho1.T", "rho2.T",
                    "DIC", "WAIC", "p.d", "LMPL", "Deviance")
  beta_rows <- rownames(sr)[!rownames(sr) %in% named_params]
  beta_mean  <- sr[beta_rows, "Mean"]
  cat("Beta rows from summary.results:", paste(beta_rows, collapse=", "), "\n")

  # Posterior mean of spatial random effects (phi), one value per area
  # model$phi is area×time matrix; use column means (time-average) as area baseline
  if (!is.null(model$phi)) {
    phi_mat   <- model$phi          # K × N or N × K — check
    if (nrow(phi_mat) == metadata$n_units) {
      phi_area  <- rowMeans(phi_mat)
    } else {
      phi_area  <- colMeans(phi_mat)
    }
    names(phi_area) <- trained_units
  } else {
    # Fall back: extract from fitted values / summary not available
    cat("WARNING: model$phi not found; setting spatial effects to 0\n")
    phi_area <- setNames(rep(0, length(trained_units)), trained_units)
  }

  rho_s  <- if ("rho.S"  %in% rownames(sr)) sr["rho.S",  "Mean"] else 0
  rho1_t <- if ("rho1.T" %in% rownames(sr)) sr["rho1.T", "Mean"] else 0
  rho2_t <- if ("rho2.T" %in% rownames(sr)) sr["rho2.T", "Mean"] else 0
  tau2   <- if ("tau2"   %in% rownames(sr)) sr["tau2",   "Mean"] else 0
  cat(sprintf("  rho.S=%.3f  rho1.T=%.3f  rho2.T=%.3f  tau2=%.3f\n",
              rho_s, rho1_t, rho2_t, tau2))

  # ── Build design matrix for future data ──────────────────────────────────
  tryCatch({
    X_future <<- model.matrix(
      ~ ns(lag1_PRCP,    knots = c(quartiles$lag1_PRCP_q1,
                                   quartiles$lag1_PRCP_q2,
                                   quartiles$lag1_PRCP_q3)) +
        ns(lag2_PRCP,    knots = c(quartiles$lag2_PRCP_q1,
                                   quartiles$lag2_PRCP_q2,
                                   quartiles$lag2_PRCP_q3)) +
        temp_max +
        ns(lag1_TEMPmax, knots = c(quartiles$lag1_TEMPmax_q1,
                                   quartiles$lag1_TEMPmax_q2,
                                   quartiles$lag1_TEMPmax_q3)) +
        ns(lag2_TEMPmax, knots = c(quartiles$lag2_TEMPmax_q1,
                                   quartiles$lag2_TEMPmax_q2,
                                   quartiles$lag2_TEMPmax_q3)) +
        ns(lag3_TEMPmax, knots = c(quartiles$lag3_TEMPmax_q1,
                                   quartiles$lag3_TEMPmax_q2,
                                   quartiles$lag3_TEMPmax_q3)),
      data = future
    )
  }, error = function(e) {
    cat("ERROR building design matrix:", e$message, "\n")
    cat("Future data columns:", paste(names(future), collapse=", "), "\n")
    stop(e)
  })

  cat("Design matrix:", nrow(X_future), "rows x", ncol(X_future), "cols | beta length:", length(beta_mean), "\n")

  # Align beta to design matrix columns (intercept is in X_future col 1)
  p <- ncol(X_future)
  if (length(beta_mean) != p) {
    warning(sprintf("Beta length (%d) != design matrix cols (%d). Adjusting.", length(beta_mean), p))
    if (length(beta_mean) > p) {
      beta_mean <- beta_mean[seq_len(p)]
    } else {
      beta_mean <- c(beta_mean, rep(0, p - length(beta_mean)))
    }
  }

  # ── Compute predictions ───────────────────────────────────────────────────
  offset_future <- log(pmax(future$pop, 1))
  phi_future    <- phi_area[future$orgunitname]
  phi_future[is.na(phi_future)] <- 0

  lp             <- as.numeric(X_future %*% beta_mean) + phi_future + offset_future
  future$sample_0 <- round(exp(lp))

  cat("Prediction summary (sample_0):\n")
  print(summary(future$sample_0))

  # ── Write output ──────────────────────────────────────────────────────────
  out <- data.frame(
    time_period = future$time_period,
    location    = future$orgunitname,
    sample_0    = future$sample_0,
    stringsAsFactors = FALSE
  )

  write.csv(out, predictions_fn, row.names = FALSE)
  cat("Wrote", nrow(out), "predictions to:", predictions_fn, "\n")
  return(invisible(out))
}

# ── CLI entry point ───────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  result <- list(historic = NULL, future = NULL, output = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--historic" && i + 1 <= length(args)) {
      result$historic <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--future" && i + 1 <= length(args)) {
      result$future <- args[i + 1]; i <- i + 2
    } else if (args[i] == "--output" && i + 1 <= length(args)) {
      result$output <- args[i + 1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  result
}

parsed <- parse_args(args)
if (!is.null(parsed$historic) && !is.null(parsed$future) && !is.null(parsed$output)) {
  predict_chap(
    model_fn         = "model.rds",
    historic_data_fn = parsed$historic,
    future_data_fn   = parsed$future,
    predictions_fn   = parsed$output
  )
} else {
  cat("Usage: Rscript predict.r --historic <historic.csv> --future <future.csv> --output <predictions.csv>\n")
  quit(status = 1)
}

