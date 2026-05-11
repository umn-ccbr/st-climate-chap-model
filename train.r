#!/usr/bin/env Rscript

## train.r
## Usage: Rscript train.r data/harmonized_data.csv output_model.rds

library(CARBayesST)
library(splines)

# Source helper functions: augment_harmonized_data(), compute_lag_quartiles()
source("model_helpers.R")

# Path to adjacency matrix (orgunitname x orgunitname), adjust if needed
W_FN <- "W_orgunits_CARBayesST.rds"

train_chap <- function(csv_fn, model_fn, W_path = W_FN, config = default_model_config) {
  cat("Loading and augmenting data from:", csv_fn, "\n")
  aug <- augment_harmonized_data(csv_fn)

  # Keep rows with enough history (timeid > 4, like final_data_lagged4)
  aug4 <- subset(aug, timeid > 4)

  # Load tradeoff analysis and select optimal configuration
  cat("Loading tradeoff analysis from: output/tradeoff_orgunits_vs_months.csv\n")
  tradeoff_path <- "output/tradeoff_orgunits_vs_months.csv"
  
  if (!file.exists(tradeoff_path)) {
    stop("Tradeoff analysis file not found. Please run scripts/find_complete_timespan.py first.")
  }
  
  tradeoff <- read.csv(tradeoff_path, stringsAsFactors = FALSE)
  
  # Filter for configurations that end at the most recent month
  # Handle both TRUE/FALSE (R) and "True"/"False" (Python) formats
  tradeoff$ends_at_recent <- as.logical(ifelse(tradeoff$ends_at_recent == "True", TRUE, 
                                               ifelse(tradeoff$ends_at_recent == "False", FALSE, 
                                                      tradeoff$ends_at_recent)))
  recent_configs <- tradeoff[tradeoff$ends_at_recent == TRUE, ]
  
  # Pick the configuration with the largest total_obs
  # Prefer recent configs, but if none exist, pick from all configs
  if (nrow(recent_configs) > 0) {
    best_config <- recent_configs[which.max(recent_configs$total_obs), ]
    cat("Selected configuration ends at most recent month\n")
  } else {
    cat("WARNING: No configurations end at most recent month. Selecting best available configuration.\n")
    # Remove configs with 0 observations
    valid_configs <- tradeoff[tradeoff$total_obs > 0, ]
    if (nrow(valid_configs) == 0) {
      stop("No valid configurations found with any observations.")
    }
    best_config <- valid_configs[which.max(valid_configs$total_obs), ]
  }
  
  cat("Selected optimal configuration:\n")
  cat("  Number of orgunits:", best_config$n_units, "\n")
  cat("  Number of months:", best_config$n_months, "\n")
  cat("  Time range:", best_config$start_time, "to", best_config$end_time, "\n")
  cat("  Time ID range:", best_config$start_timeid, "to", best_config$end_timeid, "\n")
  cat("  Total observations:", best_config$total_obs, "\n")
  
  # Load the recommended orgunits
  recommended_path <- "output/recommended_orgunits.csv"
  if (!file.exists(recommended_path)) {
    stop("Recommended orgunits file not found. Please run scripts/find_complete_timespan.py first.")
  }
  
  recommended <- read.csv(recommended_path, stringsAsFactors = FALSE)
  selected_units <- recommended$orgunitname
  
  cat("Loaded", length(selected_units), "recommended orgunits\n")
  
  # Filter data to selected orgunits and time span
  cat("Filtering data to selected orgunits and time span...\n")
  aug4 <- aug4[aug4$orgunitname %in% selected_units & 
               aug4$timeid >= best_config$start_timeid & 
               aug4$timeid <= best_config$end_timeid, ]
  
  cat("Rows after filtering:", nrow(aug4), "\n")
  
  # Save training metadata
  training_metadata <- list(
    n_units = best_config$n_units,
    n_months = best_config$n_months,
    start_time = best_config$start_time,
    end_time = best_config$end_time,
    start_timeid = best_config$start_timeid,
    end_timeid = best_config$end_timeid,
    total_obs = best_config$total_obs,
    selected_orgunits = selected_units
  )
  
  metadata_path <- "output/training_metadata.rds"
  saveRDS(training_metadata, file = metadata_path)
  cat("Saved training metadata to:", metadata_path, "\n")
  
  vars_all <- c(
    "disease_cases", "pop",
    "preci", "temp_max",
    "lag1_PRCP", "lag2_PRCP",
    "lag1_TEMPmax", "lag2_TEMPmax", "lag3_TEMPmax"
  )

  # Check data completeness without filtering
  cat("Checking data completeness...\n")
  cat("Rows before checks:", nrow(aug4), "\n")
  cat("Unique orgunits:", length(unique(aug4$orgunitname)), "\n")
  cat("Unique timeids:", length(unique(aug4$timeid)), "\n")
  
  # Check for balanced panel
  expected_rows <- length(unique(aug4$orgunitname)) * length(unique(aug4$timeid))
  cat("Expected rows (balanced panel):", expected_rows, "\n")
  
  if (nrow(aug4) != expected_rows) {
    cat("WARNING: Data is not a balanced panel!\n")
    # Show which units/times are missing
    all_combos <- expand.grid(
      orgunitname = unique(aug4$orgunitname),
      timeid = unique(aug4$timeid),
      stringsAsFactors = FALSE
    )
    aug4_check <- merge(all_combos, aug4, by = c("orgunitname", "timeid"), all.x = TRUE)
    missing_count <- sum(is.na(aug4_check$disease_cases))
    cat("Missing", missing_count, "orgunit-time combinations\n")
  }
  
  # Check for NAs in each variable
  cat("\nNA counts by variable:\n")
  for (v in vars_all) {
    na_count <- sum(is.na(aug4[[v]]))
    if (na_count > 0) {
      cat("  ", v, ":", na_count, "NAs\n")
    }
  }
  
  # Use data as-is without filtering
  dat <- aug4
  
  # Order data: first all K areas at time 1, then time 2, etc. (required by ST.CARar)
  dat <- dat[order(dat$timeid, dat$orgunitname), ]
  rownames(dat) <- NULL

  cat("\nRows in model dataset:", nrow(dat), "\n")
  cat("Distinct orgunits in model dataset:", length(unique(dat$orgunitname)), "\n")
  
  cat("Computing quartiles for spline knots...\n")
  quartiles <- compute_lag_quartiles(dat)

  # Load adjacency matrix W (indexed by orgunitname)
  cat("Loading adjacency matrix from:", W_path, "\n")
  W <- readRDS(W_path)

  # Restrict W to the orgunits actually present in the model data
  areas <- sort(unique(dat$orgunitname))

  dat <- dat[ , c("time", "orgunitname", vars_all), drop = FALSE]

  missing_in_W <- setdiff(areas, rownames(W))
  if (length(missing_in_W) > 0) {
    stop("The following orgunits are in the data but not in W:\n  ",
         paste(missing_in_W, collapse = ", "))
  }

  cat("Subsetting adjacency matrix to", length(areas), "selected orgunits\n")
  W <- W[areas, areas, drop = FALSE]
  cat("Adjacency matrix dimensions:", nrow(W), "x", ncol(W), "\n")

  # Add random neighbors to isolated areas (stopgap for zero-neighbor areas)
    cat("Checking for isolated areas (zero neighbors)...\n")
    neighbor_counts <- rowSums(W)
    isolated <- names(neighbor_counts)[neighbor_counts == 0]
    
    if (length(isolated) > 0) {
      cat("Found", length(isolated), "isolated areas with no neighbors\n")
      cat("Adding random neighbors as stopgap...\n")
      
      for (area in isolated) {
        # Get potential neighbors (all areas except itself)
        potential_neighbors <- setdiff(areas, area)
        
        # Randomly select 1-3 neighbors
        n_neighbors <- sample(1:3, 1)
        new_neighbors <- sample(potential_neighbors, min(n_neighbors, length(potential_neighbors)))
        
        # Add symmetric connections
        for (neighbor in new_neighbors) {
          W[area, neighbor] <- 1
          W[neighbor, area] <- 1
        }
        
        cat("  Added", length(new_neighbors), "random neighbor(s) to:", area, "\n")
      }
      
      # Verify all areas now have neighbors
      neighbor_counts_after <- rowSums(W)
      still_isolated <- names(neighbor_counts_after)[neighbor_counts_after == 0]
      if (length(still_isolated) == 0) {
        cat("All areas now have at least one neighbor\n")
      } else {
        warning("Some areas still have no neighbors: ", paste(still_isolated, collapse = ", "))
      }
    } else {
      cat("No isolated areas found - all areas have neighbors\n")
    }

  # (Optional) sanity check that W is symmetric and has zero diagonal
  if (!all(rownames(W) == colnames(W))) {
    stop("Row and column names of W do not match after subsetting.")
  }
  if (!all(W == t(W))) {
    warning("W is not perfectly symmetric after subsetting.")
  }
  diag(W) <- 0

  # Define spline-based formula using quartiles
  f <- disease_cases ~ offset(log(pop)) +
    ns(lag1_PRCP,    knots = c(quartiles$lag1_PRCP_q1,
                               quartiles$lag1_PRCP_q2,
                               quartiles$lag1_PRCP_q3)) +
    ns(lag2_PRCP,    knots = c(quartiles$lag2_PRCP_q1,
                               quartiles$lag2_PRCP_q2,
                               quartiles$lag2_PRCP_q3)) +
#    ns(temp_max,     knots = c(quartiles$TEMPmax_q1,
#                               quartiles$TEMPmax_q2,
#                               quartiles$TEMPmax_q3)) +
    temp_max + # TODO figure out why temp spline is causing collinearity error
    ns(lag1_TEMPmax, knots = c(quartiles$lag1_TEMPmax_q1,
                               quartiles$lag1_TEMPmax_q2,
                               quartiles$lag1_TEMPmax_q3)) +
    ns(lag2_TEMPmax, knots = c(quartiles$lag2_TEMPmax_q1,
                               quartiles$lag2_TEMPmax_q2,
                               quartiles$lag2_TEMPmax_q3)) +
    ns(lag3_TEMPmax, knots = c(quartiles$lag3_TEMPmax_q1,
                               quartiles$lag3_TEMPmax_q2,
                               quartiles$lag3_TEMPmax_q3))

  cat("Fitting ST.CARar model...\n")
  model <- CARBayesST::ST.CARar(
    formula  = f,
    family   = "poisson",
    data     = dat,
    W        = W,
    burnin   = config$burnin,   # feel free to increase for real runs
    n.sample = config$n_sample,  # feel free to increase for real runs
    thin     = config$thin,
    n.chains = 3,
    n.cores  = 3,
    AR       = 2,
    verbose  = TRUE
  )

  cat("Saving fitted model to:", model_fn, "\n")
  saveRDS(model, file = model_fn)

  quartiles_file <- paste0(model_fn, ".quartiles")
  saveRDS(quartiles, file = quartiles_file)

  cat("Done.\n")
}

# ---- CLI interface ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  csv_fn   <- args[1]
  model_fn <- args[2]

  train_chap(csv_fn, model_fn)
} else {
  cat("Usage: Rscript train.r <trainingData.csv> <output_model.rds>\n")
}
