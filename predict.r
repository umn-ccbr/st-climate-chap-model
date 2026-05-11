options(warn=1)
## Not functional yet, still working on adjusting the model saving
## -Mitch

#' Predict Malaria Cases using Trained Model
#'
#' This function performs prediction/inference using a previously trained 
#' spatiotemporal malaria model. It implements the MCMC prediction approach
#' to generate forecasts for specified time periods.
#'
#' @param trained_model_path Path to saved trained model results, or trained model object
#' @param prediction_data Optional: custom prediction data. If NULL, uses training data with missing values
#' @param prediction_periods Vector of time periods to predict (default: c(60))
#' @param n_sample Number of MCMC samples for prediction (default: 20000)
#' @param burnin Number of burn-in samples (default: 4000) 
#' @param thin Thinning interval (default: 40)
#' @param prediction_time_point Time point to extract final predictions from (default: 56)
#' @param proposal_sd_phi Proposal standard deviation for spatial effects (default: 0.1)
#' @param helpers_source_file Path to carbayes helpers (default: 'carbayes_helpers.R')
#' @param save_predictions Whether to save prediction results (default: TRUE)
#' @param result_suffix Suffix for saved prediction files
#'
#' @return List containing:
#'   \item{fitted_predictions}{Vector of predicted cases for specified time point}
#'   \item{fitted_values_all}{Matrix of all fitted values (K x N)}
#'   \item{samples_fitted}{Array of MCMC samples}
#'   \item{prediction_diagnostics}{MCMC and model diagnostics}
#'
predict_malaria_trained <- function(trained_model_path = NULL,
                                   prediction_data = NULL,
                                   config = default_model_config,
                                   prediction_time_point = 56,
                                   proposal_sd_phi = 0.1,) {

  # Load required libraries
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(MASS)
  library(tidyr)
  library(CARBayesST)
  library(splines)
  
  # Source helper functions
  source("carbayes_helpers.R")
  
  # Load trained model components
  if (!is.null(trained_model_object)) {
    # Use provided model object
    summary_results <- trained_model_object$summary_results
    formula <- trained_model_object$formula
    training_data <- trained_model_object$processed_data
    adjacency_matrix <- trained_model_object$adjacency_matrix
    region_info <- trained_model_object$region_info
    region <- region_info$region
    
    cat("Using provided trained model object\n")
    
    # Debug formula
    cat("=== FORMULA DEBUG (from model object) ===\n")
    cat("Formula class:", class(formula), "\n")
    cat("Formula structure:\n")
    print(formula)
    cat("Formula environment variables:", ls(envir = environment(formula)), "\n")
    cat("==========================================\n")
  } else if (!is.null(trained_model_path)) {
    # Load from file path (assumes specific file naming convention)
    if (file.exists(paste0(trained_model_path, '_summary.RData'))) {
      load(paste0(trained_model_path, '_summary.RData'))
      load(paste0(trained_model_path, '_components.RData'))
      summary_results <- summary # or whatever the saved variable name is
      # Extract other components from loaded training_components
      formula <- training_components$formula
      training_data <- training_components$processed_data
      adjacency_matrix <- training_components$adjacency_matrix
      region <- training_components$region
      region_info <- list(region = region, quartiles = training_components$quartiles)
    } else {
      stop("Could not find trained model files at specified path")
    }
    cat("Loaded trained model from files\n")
    
    # Debug formula
    cat("=== FORMULA DEBUG (from saved files) ===\n")
    cat("Formula class:", class(formula), "\n")
    cat("Formula structure:\n")
    print(formula)
    cat("Formula environment variables:", ls(envir = environment(formula)), "\n")
    cat("Training data variables: ", ls(envir = environment(training_data)), "\n")
    cat("==========================================\n")
  } else {
    stop("Must provide either trained_model_path or trained_model_object")
  }
  
  # Prepare prediction data
  if (is.null(prediction_data)) {
    # Use training data and set specified periods to missing
    final_data_lagged4_aug <- training_data
    
    # Set prediction periods to NA
    for (period in prediction_periods) {
      final_data_lagged4_aug[final_data_lagged4_aug$timeid == period, 'marlaria'] <- NA
    }
    
    # Also set last few rows to NA as in original code
    n_rows <- nrow(final_data_lagged4_aug)
    final_data_lagged4_aug[(n_rows-1):n_rows, 'marlaria'] <- NA
    
    cat("Using training data with missing values for prediction periods:", prediction_periods, "\n")
  } else {
    final_data_lagged4_aug <- prediction_data
    cat("Using provided prediction data\n")
  }
  
  # Debug before using formula
  cat("=== PRE-FORMULA DEBUG ===\n")
  cat("Formula about to be used:\n")
  print(formula)
  cat("Prediction data column names:", colnames(final_data_lagged4_aug), "\n")
  cat("Prediction data dimensions:", dim(final_data_lagged4_aug), "\n")
  cat("Sample of prediction data:\n")
  print(head(final_data_lagged4_aug, 3))
  cat("========================\n")

  quartiles = trained_model_object$region_info$quartiles

  lag1_PRCP_q1 = quartiles$lag1_PRCP[1]
  lag1_PRCP_q2 = quartiles$lag1_PRCP[2]
  lag1_PRCP_q3 = quartiles$lag1_PRCP[3]
  lag2_PRCP_q1 = quartiles$lag2_PRCP[1]
  lag2_PRCP_q2 = quartiles$lag2_PRCP[2]
  lag2_PRCP_q3 = quartiles$lag2_PRCP[3]
  TEMPmax_q1 = quartiles$TEMPmax[1]
  TEMPmax_q2 = quartiles$TEMPmax[2]
  TEMPmax_q3 = quartiles$TEMPmax[3]
  lag1_TEMPmax_q1 = quartiles$lag1_TEMPmax[1]
  lag1_TEMPmax_q2 = quartiles$lag1_TEMPmax[2]
  lag1_TEMPmax_q3 = quartiles$lag1_TEMPmax[3]
  lag2_TEMPmax_q1 = quartiles$lag2_TEMPmax[1]
  lag2_TEMPmax_q2 = quartiles$lag2_TEMPmax[2]
  lag2_TEMPmax_q3 = quartiles$lag2_TEMPmax[3]
  lag3_TEMPmax_q1 = quartiles$lag3_TEMPmax[1]
  lag3_TEMPmax_q2 = quartiles$lag3_TEMPmax[2]
  lag3_TEMPmax_q3 = quartiles$lag3_TEMPmax[3]

# Set up model framework
tryCatch({
    frame.results <- common.frame(formula, final_data_lagged4_aug, "poisson")
}, error = function(e) {
    cat("ERROR in common.frame():\n")
    cat("Error message:", e$message, "\n")
    cat("Formula causing error:\n")
    print(formula)
    cat("\n=== DETAILED DIAGNOSTICS ===\n")

    df = final_data_lagged4_aug
    f = formula
    
    cat("Formula environment:\n")
    print(environment(f))
    training_componenti
    cat("\nAll variables in formula:\n")
    print(all.vars(f))
    
    cat("\nVariables in formula NOT in data:\n")
    missing_vars <- setdiff(all.vars(f), names(df))
    print(missing_vars)
    
    cat("\nVariables in data:\n")
    print(names(df))
    
    # See if the formula env has same-named objects that could be masking df columns
    fe <- environment(f)
    cat("\nColumn names that ALSO exist in formula environment:\n")
    overlapping <- intersect(names(df), ls(fe))
    print(overlapping)
    
    cat("\nFormula environment contents:\n")
    print(ls(fe))
    
    cat("============================\n")

    miss <- c("lag1_PRCP_q1","lag1_PRCP_q2","lag1_PRCP_q3",
            "lag2_PRCP_q1","lag2_PRCP_q2","lag2_PRCP_q3",
            "TEMPmax_q1","TEMPmax_q2","TEMPmax_q3",
            "lag1_TEMPmax_q1","lag1_TEMPmax_q2","lag1_TEMPmax_q3",
            "lag2_TEMPmax_q1","lag2_TEMPmax_q2","lag2_TEMPmax_q3",
            "lag3_TEMPmax_q1","lag3_TEMPmax_q2","lag3_TEMPmax_q3")

    fe <- environment(f)
    sapply(miss, function(s) if (exists(s, fe, inherits=FALSE)) length(get(s, fe)) else NA)
    paste(" nrow(final_data_lagged4_aug):", nrow(final_data_lagged4_aug))

    stop("Formula processing failed: ", e$message)
})

N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
  X.standardised <- frame.results$X.standardised
  X.sd <- frame.results$X.sd
  X.mean <- frame.results$X.mean
  X.indicator <- frame.results$X.indicator 
  offset <- frame.results$offset
  Y <- frame.results$Y
  which.miss <- frame.results$which.miss
  n.miss <- frame.results$n.miss  
  Y.DA <- Y
  
  cat("Model framework set up successfully. Missing observations:", n.miss, "\n")
  cat("Design matrix dimensions: X =", dim(X), "X.standardised =", dim(X.standardised), "\n")
  
  # Extract parameters from trained model
  rho <- summary_results['rho.S','Mean']
  alpha <- c(summary_results['rho1.T','Mean'], summary_results['rho2.T','Mean'])
  
  # Determine beta dimensions based on region
  beta_dimensions <- list(
    "north" = 25,
    "central" = 25,
    "south" = 25
  )
  
  beta_n <- beta_dimensions[[region]]
  beta <- summary_results[1:beta_n, 'Mean']
  
  cat("Parameters extracted. Beta dimension:", beta_n, ", rho:", round(rho, 3), "\n")
  
  # Setup CAR quantities
  W.quants <- common.Wcheckformat.leroux(adjacency_matrix)
  K <- W.quants$n
  N <- N.all / K
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  W.n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  # Compute beta blocking structure
  block.temp <- common.betablock(p)
  beta.beg  <- block.temp[[1]]
  beta.fin <- block.temp[[2]]
  n.beta.block <- block.temp[[3]]
  
  # Initialize spatial effects and other parameters
  log.Y <- log(Y)
  log.Y[Y==0] <- -0.1  
  res.temp <- log.Y - X.standardised %*% beta - offset
  res.sd <- sd(res.temp, na.rm=TRUE)/5
  phi <- rnorm(n=N.all, mean=0, sd = res.sd)
  tau2 <- summary_results['tau2','Mean']
  
  # Create matrix quantities
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
  phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   
  fitted <- exp(as.numeric(offset.mat + regression.mat + phi.mat))
  
  cat("Spatial structure set up. K =", K, ", N =", N, "\n")
  
  # Setup MCMC storage
  n.keep <- floor((n_sample - burnin)/thin)
  samples.fitted <- array(NA, c(n.keep, N.all))
  if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  # Initialize acceptance counters
  accept.all <- rep(0,6)
  accept <- accept.all
  
  # Handle missing data imputation
  if(n.miss > 0) {
    Y.DA[which.miss==0] <- tail(Y.DA[which.miss!=0], K)   
  }
  Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)
  
  cat("Starting MCMC prediction sampling...\n")
  cat("Samples:", n_sample, ", Burn-in:", burnin, ", Thin:", thin, ", Keep:", n.keep, "\n")
  
  # MCMC sampling loop
  start_time <- Sys.time()
  for(j in 1:n_sample) {
    
    # Sample from spatial random effects (phi)
    phi.offset <- offset.mat + regression.mat
    den.offset <- rho * W.triplet.sum + 1 - rho
    temp1 <- poissonar2carupdateRW(W.triplet, W.begfin, W.triplet.sum,  
                                   K, N, phi.mat, tau2, alpha[1], alpha[2], 
                                   rho, Y.DA.mat, proposal_sd_phi, 
                                   phi.offset, den.offset)      
    phi.temp <- temp1[[1]]
    phi <- as.numeric(phi.temp) - mean(as.numeric(phi.temp))
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    
    # Update fitted values
    lp <- as.numeric(offset.mat + regression.mat + phi.mat)
    fitted <- exp(lp)
    
    # Store results after burn-in and thinning
    if(j > burnin & (j-burnin)%%thin==0) {
      ele <- (j - burnin) / thin
      samples.fitted[ele, ] <- fitted
    }
    
    # Progress reporting
    if(j %% 2000 == 0) {
      elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units="mins")), 1)
      cat("Iteration:", j, "/", n_sample, "(", round(100*j/n_sample, 1), "%) -", 
          elapsed, "min elapsed\n")
    }
  }
  
  end_time <- Sys.time()
  total_time <- round(as.numeric(difftime(end_time, start_time, units="mins")), 1)
  cat("MCMC sampling completed in", total_time, "minutes\n")
  
  # Compute final predictions
  fitted.values <- apply(samples.fitted, 2, mean)
  fitted.values_mat <- matrix(fitted.values, nrow=K, ncol=N, byrow=FALSE)
  
  
  # Prepare diagnostics
  prediction_diagnostics <- list(
    n_sample = n_sample,
    burnin = burnin,
    thin = thin,
    n_keep = n.keep,
    n_missing = n.miss,
    K = K,
    N = N,
    region = region,
    prediction_periods = prediction_periods,
    prediction_time_point = prediction_time_point,
    total_time_minutes = total_time,
    rho = rho,
    alpha = alpha,
    tau2 = tau2
  )
  
  cat("Prediction completed successfully!\n")
  cat("Predicted", length(fitted_predictions), "values for time point", prediction_time_point, "\n")
  
  # Save results if requested
  if (save_predictions) {
    save(fitted_predictions, file = paste0("fitted_pred_", region, "_", result_suffix, ".RData"))
    save(fitted.values_mat, file = paste0("fitted_values_all_", region, "_", result_suffix, ".RData"))
    save(prediction_diagnostics, file = paste0("prediction_diagnostics_", region, "_", result_suffix, ".RData"))
    cat("Prediction results saved\n")
  }
  
  # Return results
  return(list(
    fitted_predictions = fitted_predictions,
    fitted_values_all = fitted.values_mat,
    samples_fitted = samples.fitted,
    prediction_diagnostics = prediction_diagnostics
  ))
}


#' Convenience function to run prediction with saved model files
#' 
#' @param region Region name ("north", "south", "central")
#' @param model_suffix Suffix used when saving model (default: "spline_full") 
#' @param ... Additional arguments passed to predict_malaria_trained
#'
predict_from_saved_model <- function(region, model_suffix = "spline_full", ...) {
  
  # Construct file paths based on naming convention from training
  summary_file <- paste0('fit_beta_lagged3_max_summary_', region, '_', model_suffix, '.RData')
  components_file <- paste0('training_components_', region, '_', model_suffix, '.RData')
  
  if (!file.exists(summary_file) || !file.exists(components_file)) {
    stop("Could not find required model files. Expected:\n", 
         summary_file, "\n", components_file)
  }
  
  # Load the files
  cat("Loading summary file:", summary_file, "\n")
  load(summary_file)  # loads 'summary_results'
  cat("Loading components file:", components_file, "\n")
  load(components_file)  # loads 'training_components'
  
  # Debug: Check what variables were loaded
  cat("All variables in environment:", ls(), "\n")
  cat("Variables matching 'summary':", ls(pattern="summary"), "\n")
  cat("Variables matching 'training':", ls(pattern="training"), "\n")
  
  # Check if expected variables exist
  if (!exists("summary_results")) {
    cat("WARNING: summary_results not found, available variables:", ls(), "\n")
  } else {
    cat("summary_results loaded successfully, class:", class(summary_results), "\n")
  }
  
  if (!exists("training_components")) {
    cat("WARNING: training_components not found, available variables:", ls(), "\n")
  } else {
    cat("training_components loaded successfully, structure:\n")
    cat("Names:", names(training_components), "\n")
    cat("Formula in training_components:", class(training_components$formula), "\n")
  }
  
  # Create a model object
  trained_model <- list(
    summary_results = summary_results,
    formula = training_components$formula,
    processed_data = training_components$processed_data,
    adjacency_matrix = training_components$adjacency_matrix,
    region_info = list(
      region = training_components$region,
      quartiles = training_components$quartiles
    )
  )
  
  # Debug the created model object
  cat("=== CREATED MODEL OBJECT DEBUG ===\n")
  cat("Model object names:", names(trained_model), "\n")
  cat("Formula from model object:\n")
  print(trained_model$formula)
  cat("=================================\n")
  
  # Run prediction
  return(predict_malaria_trained(trained_model_object = trained_model, ...))
}


## Base predict_chap template

predict_chap <- function(model_fn, historic_data_fn, future_climatedata_fn, predictions_fn) {
  df <- read.csv(future_climatedata_fn)
  X <- df[, c("rainfall", "mean_temperature"), drop = FALSE]
  model <- readRDS(model_fn)  # Assumes the model was saved using saveRDS

  y_pred <- predict(model, newdata = X)
  df$sample_0 <- y_pred
  write.csv(df, predictions_fn, row.names = FALSE)

  print(paste("Forecasted values:", paste(y_pred, collapse = ", ")))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 4) {
  model_fn <- args[1]
  historic_data_fn <- args[2]
  future_climatedata_fn <- args[3]
  predictions_fn <- args[4]
  
  predict_chap(model_fn, historic_data_fn, future_climatedata_fn, predictions_fn)
}# else {
#  stop("Usage: Rscript predict.R <model_fn> <historic_data_fn> <future_climatedata_fn> <predictions_fn>")
#}
