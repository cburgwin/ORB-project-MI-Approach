#---------------------------------------------------------
# FINAL script to send to cluster
#---------------------------------------------------------

# load functions and libraries
source("Code/00.Functions.R")
library(parallel)
library(dplyr)
library(readr)

# set seed for reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# ---------------------------------------------------------
# Global Parameters
# ---------------------------------------------------------
n_sim         <- 1900
M_imputations <- 200
n_cores       <- parallel::detectCores() - 32

# create folder to save each scenario
save_dir <- "Data/Sim_results"

if (dir.exists(save_dir)) {
  # Find all old scenario files inside the directory
  old_files <- list.files(save_dir, pattern = "scenario_.*\\.rds", full.names = TRUE)
  if (length(old_files) > 0) {
    # Delete them cleanly
    unlink(old_files, force = TRUE)
  }
} else {
  # Create the folder if it doesn't exist yet
  dir.create(save_dir, recursive = TRUE)
}
# create folder to save final csv output file
out_dir <- "Data_Sim"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---------------------------------------------------------
# Scenario grid 
# ---------------------------------------------------------
scenarios_grid <- expand.grid(
  K           = c(25),
  p1          = c(0.2),
  tau2_val    = c(0, 0.02, 0.06, 0.36),
  theta_1     = c(0.4),
  theta_2     = c(0.4),
  rho_b       = c(0.8),
  rho_w       = c(0.8),
  delta_sim   = seq(0, 1, by = 0.2),
  delta_est   = seq(0, 1, by = 0.2),
  select_type = c("zscore", "effect"),
  stringsAsFactors = FALSE
)  %>% unique() %>%
  dplyr::filter(`delta_sim` == `delta_est`)


total_scenarios <- nrow(scenarios_grid)
cat(sprintf("System verified. Starting evaluation of %d simulation scenarios across %d system cores...\n",
            total_scenarios, n_cores))


# ---------------------------------------------------------
# Wrappers to handle two cases of errors:
# 1) -> failure of convergence of rma inside imputation function (either for 
# numerical issues or more probably cause there are less than 2 reported studies)
# 2) -> failure of convergence of rma inside adjust function for numerical issues
# ---------------------------------------------------------

safe_adj_uni <- function(mi, delta, sel_type) {
  
  # it get passed a null df in the case of a failure in rma naive
  if (is.null(mi)) return(c(est = NA, ci_l = NA, ci_u = NA, ess = NA))
  
  tryCatch({
    res <- adj_univariate(mi,
                          delta = delta, 
                          select_type = sel_type, 
                          model_type = "REML",
                          track.ess = TRUE)
    
    return(c(est = res$Estimate[1],
             ci_l = res$CI_Lower[1],
             ci_u = res$CI_Upper[1],
             ess = res$ess[1]))
  }, error = function(e) c(est = NA, ci_l = NA, ci_u = NA, ess = NA))
}

safe_adj_biv <- function(mi, delta, sel_type) {
  
  # it get passed a null df in the case of a failure in rma naive
  if (is.null(mi)) return(c(est = NA, ci_l = NA, ci_u = NA, ess = NA, fail = NA))
  
  tryCatch({
    res <- adj_bivariate(mi, 
                         delta = delta,
                         select_type = sel_type)
    return(c(est = res$Estimate[1], 
             ci_l = res$CI_Lower[1],
             ci_u = res$CI_Upper[1],
             ess = res$ess[1], 
             fail = res$failed.proportion[1]))
  }, error = function(e) c(est = NA, ci_l = NA, ci_u = NA, ess = NA, fail = NA))
}


# ---------------------------------------------------------
# MAIN FUNCTION
# ---------------------------------------------------------

run_ORB <- function(scenario_idx) {
  
  output_file <- sprintf("%s/scenario_%05d.rds", save_dir, scenario_idx)
  if (file.exists(output_file)) return(NULL)
  
  # run through scenario grid
  s <- scenarios_grid[scenario_idx, ]
  
  K           <- s$K
  theta_1     <- s$theta_1
  theta_2     <- s$theta_2
  tau2_val    <- s$tau2_val
  rho_b       <- s$rho_b
  rho_w       <- s$rho_w
  delta_sim   <- s$delta_sim
  p1          <- s$p1
  select_type <- s$select_type
  delta_est   <- s$delta_est
  
  true_theta  <- theta_1
  
  sim_results <- vector("list", n_sim)
  
  # Number of admissible simulation repetitions
  n_eligible <- 0
  
  # Number of attempted repetitions, including failed attempts
  n_attempts <- 0
  n_unexpected_failures <- 0
  
  # Maximum number of outer attempts
  max_attempts <- 5 * n_sim
  
  # Redraw diagnostics
  redraws_total <- 0
  n_redraw_repetitions <- 0
  
  # Maximum number of redraws allowed for a single repetition
  max_redraws_per_rep <- 100
  
  # =========================================================
  # Compute Complete Estimate Bivariate
  # =========================================================
  compute_full <- function(full_data) {
    
    res_biv_long <- data.frame(
      Study_id = rep(full_data$Study_id, each = 2),
      outcome  = factor(rep(c("O1", "O2"), times = K)),
      yi       = as.numeric(t(as.matrix(full_data[, c("O1_yi", "O2_yi")]))),
      sei      = as.numeric(t(as.matrix(full_data[, c("O1_sei", "O2_sei")])))
    )
    
    V_list <- lapply(1:K, function(j) {
      v1 <- full_data$O1_sei[j]^2
      v2 <- full_data$O2_sei[j]^2
      cov12 <- rho_w * sqrt(v1) * sqrt(v2)
      matrix(c(v1, cov12, cov12, v2), 2, 2)
    })
    
    V <- as.matrix(Matrix::bdiag(V_list))
    
    res_biv <- rma.mv(yi,
                      V = V,
                      mods = ~ outcome - 1,
                      random = ~ outcome | Study_id,
                      struct = "UN",
                      data = res_biv_long,
                      method = "REML",
                      tau2 = NULL,
                      rho = NULL)
    
    list(
      biv = c(est = res_biv$beta[1],
              l   = res_biv$ci.lb[1],
              u   = res_biv$ci.ub[1])
    )
  }
  
  # =========================================================
  # Main loop
  # =========================================================
  while (n_eligible < n_sim && n_attempts < max_attempts) {
    
    n_attempts <- n_attempts + 1
    
    res_df <- tryCatch({
      
      n_redraws <- 0
      
      repeat {
        
        # Generate complete bivariate dataset
        full_data <- generate_bivariate_ma(K = K,
                                           theta = c(theta_1, theta_2),
                                           tau2 = c(tau2_val, tau2_val),
                                           rho_b = rho_b,
                                           rho_w = rho_w)
        
        # Introduce ORB
        obs_data <- impose_orb(full_data,
                               p1 = p1,
                               delta_sim = delta_sim,
                               select_type = select_type,
                               orb.se = TRUE,
                               theta_1 = theta_1,
                               tau2_val = tau2_val,
                               n_arm = 50)
        
        n_reported <- sum(!is.na(obs_data$O1_yi))
        
        # Admissible dataset
        if (n_reported >= 4) {
          break
        }
        
        n_redraws <- n_redraws + 1
        
        # Avoid an infinite loop in extremely difficult scenarios
        if (n_redraws >= max_redraws_per_rep) {
          stop("Maximum number of redraws exceeded.")
        }
      }
      
      redraws_total <- redraws_total + n_redraws
      
      if (n_redraws > 0) {
        n_redraw_repetitions <- n_redraw_repetitions + 1
      }
      
      # -----------------------------------------------------
      # 1. FULL (always from full_data)
      # -----------------------------------------------------
      full_res <- tryCatch(
        compute_full(full_data),
        error = function(e) {
          message("Full-data analysis failed: ", conditionMessage(e))
          NULL
        }
      )
      
      full_success <- !is.null(full_res) && all(is.finite(full_res$biv[c("est", "l", "u")]))
      
      # -----------------------------------------------------
      # 2. NAIVE (always from reported data)
      # -----------------------------------------------------
      naive_res <- tryCatch({
        
        observed_uni <- obs_data[!is.na(obs_data$O1_yi), ]
        K_observed <- nrow(observed_uni)
        
        complete_cases <- which(!is.na(observed_uni$O1_yi) & !is.na(observed_uni$O2_yi))
        n_complete_pairs <- length(complete_cases)
        
        # Correlation used for the naive analysis
        if (n_complete_pairs >= 4) {
          
          rho_hat <- cor(observed_uni$O1_yi[complete_cases],
                         observed_uni$O2_yi[complete_cases])
          
          rho_hat_fallback <- FALSE
        } else {
          rho_hat <- rho_w
          rho_hat_fallback <- TRUE
        }
        
        res_naive_biv_long <- data.frame(
          Study_id = rep(observed_uni$Study_id, each = 2),
          outcome = factor(rep(c("O1", "O2"), times = K_observed)),
          yi = as.numeric(t(as.matrix(observed_uni[, c("O1_yi", "O2_yi")]))),
          sei = as.numeric(t(as.matrix(observed_uni[, c("O1_sei", "O2_sei")]))))
        
        V_list <- lapply(1:K_observed, function(j) {
          
          v1 <- observed_uni$O1_sei[j]^2
          v2 <- observed_uni$O2_sei[j]^2
          cov12 <- rho_w * sqrt(v1) * sqrt(v2)
          matrix(c(v1, cov12,
                   cov12, v2), 2, 2)
        })
        
        V_naive <- as.matrix(Matrix::bdiag(V_list))
        
        res_naive_biv <- rma.mv(
          yi,
          V = V_naive,
          mods = ~ outcome - 1,
          random = ~ outcome | Study_id,
          struct = "UN",
          tau2 = NULL,
          rho = rho_hat,
          data = res_naive_biv_long,
          method = "REML",
          control = list(rel.tol = 1e-5,
                         maxiter = 200)
        )
        
        list(
          est = as.numeric(res_naive_biv$beta[1]),
          ci_l = as.numeric(res_naive_biv$ci.lb[1]),
          ci_u = as.numeric(res_naive_biv$ci.ub[1]),
          n_complete_pairs = n_complete_pairs,
          rho_hat = rho_hat,
          rho_hat_fallback = rho_hat_fallback
        )
        
      }, error = function(e) {
        
        message("Naive analysis failed: ", conditionMessage(e))
        
        list(est = NA_real_,
             ci_l = NA_real_,
             ci_u = NA_real_,
             n_complete_pairs = NA_integer_,
             rho_hat = NA_real_,
             rho_hat_fallback = NA)
      })
      
      naive_success <- is.finite(naive_res$est) &&
        is.finite(naive_res$ci_l) &&
        is.finite(naive_res$ci_u)
      
      # -----------------------------------------------------
      # 3. ADJUSTED
      # -----------------------------------------------------
      # impute unreported standard errors
      obs_data_imp <- tryCatch(
        impute_missing_se(obs_data,
                          "O1_yi",
                          "O1_sei",
                          "n_total"),
        error = function(e) {
          message("SE imputation failed: ", conditionMessage(e))
          NULL
        }
      )
      
      # univariate imputation
      mi_uni <- NULL
      
      if (!is.null(obs_data_imp)) {
        
        mi_uni <- tryCatch(
          run_univariate_imputation(obs_data_imp,
                                    theta_col = "O1_yi",
                                    se_col = "O1_sei",
                                    m = M_imputations),
          error = function(e) {
            message("Univariate imputation failed: ", conditionMessage(e))
            NULL
          }
        )
      }
      
      adj_uni <- safe_adj_uni(mi_uni,
                              delta_est,
                              select_type)
      
      uni_success <- all(is.finite(as.numeric(adj_uni[c("est", "ci_l", "ci_u")])))
      
      
      # bivariate imputation
      mi_biv <- NULL
      
      if (!is.null(obs_data_imp)) {
        
        mi_biv <- tryCatch(
          run_bivariate_imputation(obs_data_imp,
                                   theta_cols = c("O1_yi", "O2_yi"),
                                   se_cols = c("O1_sei", "O2_sei"),
                                   rho_w = "pearson",
                                   rho_b = NULL, # set them equal for kirkham global correlation
                                   tau2_val = NULL,
                                   m = M_imputations),
          error = function(e) {
            message("Bivariate imputation failed: ", conditionMessage(e))
            NULL
          }
        )
      }
      
      adj_biv <- safe_adj_biv(mi_biv,
                              delta_est,
                              select_type)
      
      biv_success <- all(is.finite(as.numeric(adj_biv[c("est", "ci_l", "ci_u")])))
      
      # -----------------------------------------------------
      # OUTPUT
      # -----------------------------------------------------
      output_row <- data.frame(
        
        # -----------------------------
        # Simulation diagnostics
        # -----------------------------
        n_redraws = n_redraws,
        full_success = full_success,
        naive_success = naive_success,
        uni_success = uni_success,
        biv_success = biv_success,
        
        n_complete_pairs = naive_res$n_complete_pairs,
        
        rho_hat_fallback = naive_res$rho_hat_fallback,
        
        # -----------------------------
        # Full
        # -----------------------------
        full = if (full_success)
          full_res$biv["est"] else NA_real_,
        
        full_ci_l = if (full_success)
          full_res$biv["l"] else NA_real_,
        
        full_ci_u = if (full_success)
          full_res$biv["u"] else NA_real_,
        
        # -----------------------------
        # Naive
        # -----------------------------
        naive_biv      = naive_res$est,
        naive_biv_ci_l = naive_res$ci_l,
        naive_biv_ci_u = naive_res$ci_u,
        
        # -----------------------------
        # Adjusted univariate
        # -----------------------------
        uni      = as.numeric(adj_uni["est"]),
        uni_ci_l = as.numeric(adj_uni["ci_l"]),
        uni_ci_u = as.numeric(adj_uni["ci_u"]),
        
        # -----------------------------
        # Adjusted bivariate
        # -----------------------------
        biv      = as.numeric(adj_biv["est"]),
        biv_ci_l = as.numeric(adj_biv["ci_l"]),
        biv_ci_u = as.numeric(adj_biv["ci_u"]),
        
        # -----------------------------
        # Diagnostics
        # -----------------------------
        u_ess = as.numeric(adj_uni["ess"]),
        b_ess = as.numeric(adj_biv["ess"]),
        b_f   = as.numeric(adj_biv["fail"])
      )
      output_row
    }, error = function(e) {
      
      message("Unexpected repetition-level error: ", conditionMessage(e))
      
      n_unexpected_failures <<- n_unexpected_failures + 1
      
      NULL
    })
    
    if (!is.null(res_df)) {
      n_eligible <- n_eligible + 1
      sim_results[[n_eligible]] <- res_df
    }
    
  }  
  
  # Combine all successful simulation repetitions
  res_df <- do.call(rbind,
                    sim_results[seq_len(n_eligible)])
  
  metric_fun <- function(est,
                         ci_l,
                         ci_u,
                         true_theta
  ) {
    
    valid <- complete.cases(est, ci_l, ci_u)
    
    list(N = sum(valid),
         Bias = mean(est[valid]) - true_theta,
         MSE = mean((est[valid] - true_theta)^2),
         Coverage = mean(ci_l[valid] <= true_theta &
                           ci_u[valid] >= true_theta),
         CI_Width = mean(ci_u[valid] - ci_l[valid]),
         Failure_Rate = mean(!valid))
  }
  
  # method-wise analysis
  full_m <- metric_fun(res_df$full,
                       res_df$full_ci_l,
                       res_df$full_ci_u,
                       true_theta)
  
  naive_m <- metric_fun(res_df$naive_biv,
                        res_df$naive_biv_ci_l,
                        res_df$naive_biv_ci_u,
                        true_theta)
  
  uni_m <- metric_fun(res_df$uni,
                      res_df$uni_ci_l,
                      res_df$uni_ci_u,
                      true_theta)
  
  biv_m <- metric_fun(res_df$biv,
                      res_df$biv_ci_l,
                      res_df$biv_ci_u,
                      true_theta)
  
  common_success <- with(
    res_df,
    full_success &
      naive_success &
      uni_success &
      biv_success
  )
  
  res_common <- res_df[common_success, ]
  
  # repetition-wise analysis
  full_common <- metric_fun(res_common$full,
                            res_common$full_ci_l,
                            res_common$full_ci_u,
                            true_theta)
  
  naive_common <- metric_fun(res_common$naive_biv,
                             res_common$naive_biv_ci_l,
                             res_common$naive_biv_ci_u,
                             true_theta)
  
  uni_common <- metric_fun(res_common$uni,
                           res_common$uni_ci_l,
                           res_common$uni_ci_u,
                           true_theta)
  
  biv_common <- metric_fun(res_common$biv,
                           res_common$biv_ci_l,
                           res_common$biv_ci_u,
                           true_theta)
  
  # failure_rate <- (n_attempts - n_success) / n_attempts
  
  avg_redraws <- if (n_eligible > 0) {
    redraws_total / n_eligible
  } else {
    NA_real_
  }
  
  summary_row <- data.frame(
    
    # Scenario parameters
    scenario_idx = scenario_idx,
    K = K,
    theta_1 = theta_1,
    theta_2 = theta_2,
    tau2_val = tau2_val,
    rho_b = rho_b,
    rho_w = rho_w,
    delta_sim = delta_sim,
    p1 = p1,
    select_type = select_type,
    delta_est = delta_est,
    
    # -----------------------------
    # Assessing missingness
    # -----------------------------
    N_Eligible = n_eligible,
    N_Attempts = n_attempts,
    N_Unexpected_Failures = n_unexpected_failures,
    N_Redraws = sum(res_df$n_redraws),
    N_Repeated_Redraw = sum(res_df$n_redraws > 0),
    Prop_Redrawn = mean(res_df$n_redraws > 0),
    Avg_Redraws = mean(res_df$n_redraws),
    Max_Redraws = max(res_df$n_redraws),
    
    # -----------------------------
    # Method-specific failures
    # -----------------------------
    N_Full_Valid = full_m$N,
    Failure_Rate_Full = full_m$Failure_Rate,
    
    N_Naive_Valid = naive_m$N,
    Failure_Rate_Naive = naive_m$Failure_Rate,
    
    N_Uni_Valid = uni_m$N,
    Failure_Rate_Uni = uni_m$Failure_Rate,
    
    N_Biv_Valid = biv_m$N,
    Failure_Rate_Biv = biv_m$Failure_Rate,
    
    # -----------------------------
    # Complete-data benchmark
    # -----------------------------
    
    Bias_Full = full_m$Bias,
    MSE_Full = full_m$MSE,
    Coverage_Full = full_m$Coverage,
    CI_Width_Full = full_m$CI_Width,
    
    # -----------------------------
    # Naive
    # -----------------------------
    
    Bias_Naive_Biv = naive_m$Bias,
    MSE_Naive_Biv = naive_m$MSE,
    Coverage_Naive_Biv = naive_m$Coverage,
    CI_Width_Naive_Biv = naive_m$CI_Width,
    
    # -----------------------------
    # Adjusted univariate
    # -----------------------------
    
    Bias_Adj_Uni = uni_m$Bias,
    MSE_Adj_Uni = uni_m$MSE,
    Coverage_Adj_Uni = uni_m$Coverage,
    CI_Width_Adj_Uni = uni_m$CI_Width,
    
    # -----------------------------
    # Adjusted bivariate
    # -----------------------------
    
    Bias_Adj_Biv = biv_m$Bias,
    MSE_Adj_Biv = biv_m$MSE,
    Coverage_Adj_Biv = biv_m$Coverage,
    CI_Width_Adj_Biv = biv_m$CI_Width,
    
    # -----------------------------
    # Common successful repetitions
    # -----------------------------
    
    N_Common = sum(common_success),
    Prop_Common = mean(common_success),
    
    Bias_Full_Common = full_common$Bias,
    MSE_Full_Common = full_common$MSE,
    Coverage_Full_Common = full_common$Coverage,
    CI_Width_Full_Common = full_common$CI_Width,
    
    Bias_Naive_Common = naive_common$Bias,
    MSE_Naive_Common = naive_common$MSE,
    Coverage_Naive_Common = naive_common$Coverage,
    CI_Width_Naive_Common = naive_common$CI_Width,
    
    Bias_Adj_Uni_Common = uni_common$Bias,
    MSE_Adj_Uni_Common = uni_common$MSE,
    Coverage_Adj_Uni_Common = uni_common$Coverage,
    CI_Width_Adj_Uni_Common = uni_common$CI_Width,
    
    Bias_Adj_Biv_Common = biv_common$Bias,
    MSE_Adj_Biv_Common = biv_common$MSE,
    Coverage_Adj_Biv_Common = biv_common$Coverage,
    CI_Width_Adj_Biv_Common = biv_common$CI_Width,
    
    # -----------------------------
    # Additional diagnostics
    # -----------------------------
    
    Prop_Rho_Fallback = mean(res_df$rho_hat_fallback, na.rm = TRUE),
    Mean_Complete_Pairs = mean(res_df$n_complete_pairs,na.rm = TRUE),
    Mean_ESS_Uni = mean(res_df$u_ess, na.rm = TRUE),
    Mean_ESS_Biv = mean(res_df$b_ess, na.rm = TRUE),
    Mean_Failed_Imputations_Biv = mean(res_df$b_f, na.rm = TRUE)
  )
  
  saveRDS(summary_row, file = output_file)
  return(NULL)
}


# Using native mclapply at scenario level to eliminate data transfer friction
mclapply(
  X = 1:total_scenarios,
  FUN = run_ORB,
  mc.cores = n_cores,
  mc.preschedule = FALSE # CRITICAL: Dynamic balancing so slow scenarios don't stall cores
)

# save data
cat("\nAll scenario files calculated. Merging to master data frame... ")
all_files <- list.files(save_dir, pattern = "scenario_.*\\.rds", full.names = TRUE)
final_metrics_df <- do.call(rbind, lapply(all_files, readRDS))
write_csv(final_metrics_df,  file.path(out_dir, "data_simulation_test.csv"))
cat("Complete!\n")


