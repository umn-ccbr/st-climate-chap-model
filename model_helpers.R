## augment_harmonized_data.R
## Load data/harmonized_data.csv and augment with lagged indicators.

augment_harmonized_data <- function(
  csv_path      = "data/harmonized_data.csv",
  max_lag_prcp  = 2,
  max_lag_temp  = 3
) {
  # ---- Load harmonized CSV ----
  df <- read.csv(csv_path, stringsAsFactors = FALSE)

  # ---- Normalise CHAP-written column names to the names used internally ----
  # CHAP writes: time_period (YYYY-MM), location, population
  # Internal names expected below: time (YYYYMM), orgunitname, pop
  if ("time_period" %in% names(df) && !"time" %in% names(df)) {
    # Convert "2017-01" → 201701
    df$time <- as.integer(gsub("-", "", df$time_period))
  }
  if ("location" %in% names(df) && !"orgunitname" %in% names(df)) {
    df$orgunitname <- df$location
  }
  if ("population" %in% names(df) && !"pop" %in% names(df)) {
    df$pop <- df$population
  }

  # CARBayesST Poisson requires integer response; CHAP's DataSet always writes float.
  if ("disease_cases" %in% names(df)) {
    df$disease_cases <- suppressWarnings(as.integer(round(as.numeric(df$disease_cases))))
  }

  required_cols <- c("time", "orgunitname", # "orgunitid",
                     "temp_max", "preci", "pop", "disease_cases")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing expected columns in harmonized data: ",
         paste(missing, collapse = ", "))
  }
  
  # ---- Basic indices for CARBayesST-style use later ----
  # spaceid based on orgunitname (since adjacency is keyed by orgunitname)
  df$spaceid <- as.integer(factor(df$orgunitname,
                                  levels = sort(unique(df$orgunitname))))
  
  # timeid based on sorted unique time values
  # (assuming time is like 201701, 201702, ...)
  df$timeid <- match(df$time, sort(unique(df$time)))
  
  # ---- Map to names used in your previous pipeline ----
  df$marlaria   <- df$disease_cases
  df$population <- df$pop
  
  # Center PRCP (preci) and TEMPmax (temp_max) across all rows (no scaling)
  df$PRCP    <- as.numeric(scale(df$preci,    scale = FALSE))
  df$TEMPmax <- as.numeric(scale(df$temp_max, scale = FALSE))
  
  # ---- Sort by orgunitname + time for consistent lagging ----
  df <- df[order(df$orgunitname, df$timeid), ]
  rownames(df) <- NULL
  
  # ---- Helper to add lags within each orgunitname ----
  add_lags_by_orgunit <- function(data, value_col, max_lag) {
    orgs <- unique(data$orgunitname)
    
    for (lag_k in seq_len(max_lag)) {
      new_col <- paste0("lag", lag_k, "_", value_col)
      data[[new_col]] <- NA_real_
    }
    
    for (org in orgs) {
      idx <- which(data$orgunitname == org)
      x   <- data[idx, value_col]
      n   <- length(x)
      
      for (lag_k in seq_len(max_lag)) {
        new_col <- paste0("lag", lag_k, "_", value_col)
        if (n > lag_k) {
          data[idx, new_col] <- c(rep(NA_real_, lag_k),
                                  x[1:(n - lag_k)])
        } else {
          data[idx, new_col] <- NA_real_
        }
      }
    }
    
    data
  }
  
  # ---- Add PRCP lags: lag1_PRCP, lag2_PRCP ----
  if (max_lag_prcp >= 1) {
    df <- add_lags_by_orgunit(df, value_col = "PRCP", max_lag = max_lag_prcp)
  }
  
  # ---- Add TEMPmax lags: lag1_TEMPmax, lag2_TEMPmax, lag3_TEMPmax ----
  if (max_lag_temp >= 1) {
    df <- add_lags_by_orgunit(df, value_col = "TEMPmax", max_lag = max_lag_temp)
  }
  
  # ---- Done ----
  df
}

compute_lag_quartiles <- function(df) {
  # Variables we care about based on your formula
  vars_prcp  <- c("lag1_PRCP", "lag2_PRCP")
  vars_temp  <- c("TEMPmax", "lag1_TEMPmax", "lag2_TEMPmax", "lag3_TEMPmax")
  vars_all   <- c(vars_prcp, vars_temp)
  
  # Check they exist
  missing <- setdiff(vars_all, names(df))
  if (length(missing) > 0) {
    stop("Missing expected columns in augmented data: ",
         paste(missing, collapse = ", "))
  }
  
  # Restrict to rows with complete data for all these variables
  complete_rows <- stats::complete.cases(df[ , vars_all])
  if (!any(complete_rows)) {
    stop("No rows with complete data for all lag variables; ",
         "check your lag construction / filtering.")
  }
  df_use <- df[complete_rows, , drop = FALSE]
  
  # Helper: 25th, 50th, 75th percentiles
  get_q <- function(x) {
    stats::quantile(x, probs = c(0.25, 0.50, 0.75),
                    na.rm = TRUE, names = FALSE)
  }
  
  quartiles <- list()
  
  for (v in vars_all) {
    q <- get_q(df_use[[v]])
    quartiles[[paste0(v, "_q1")]] <- q[1]
    quartiles[[paste0(v, "_q2")]] <- q[2]
    quartiles[[paste0(v, "_q3")]] <- q[3]
  }
  
  quartiles
}

default_model_config <- list(
  n_sample              = 50000,
  burnin                = 10000,
  thin                  = 40,
  prediction_time_point = 56,
  proposal_sd_phi       = 0.1,
  save_predictions      = TRUE,
  result_suffix         = "prediction"
)
