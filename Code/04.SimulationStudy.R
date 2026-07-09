#---------------------------------------------------------
# Final script to send to cluster
#---------------------------------------------------------

# load functions and libraries
source("00.Functions.R")
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
save_dir <- "sim_results"

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
out_dir <- "data_Sim"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---------------------------------------------------------
# Scenario grid 
# ---------------------------------------------------------
scenarios_grid <- expand.grid(
  K           = c(6, 12, 25),
  p1          = c(0.2, 0.4),
  tau2_val    = c(0, 0.02, 0.06, 0.36),
  theta_1     = c(0, 0.4),
  theta_2     = c(0, 0.4),
  rho_b       = c(0, 0.4),
  rho_w       = c(0, 0.4),
  delta_sim   = seq(0, 1, by = 0.2),
  delta_est   = seq(0, 1, by = 0.2),
  select_type = c("zscore", "effect"),
  stringsAsFactors = FALSE
)  %>% unique() 


total_scenarios <- nrow(scenarios_grid)
cat(sprintf("System verified. Starting evaluation of %d simulation scenarios across %d system cores...\n",
            total_scenarios, n_cores))


# ---------------------------------------------------------
# Wrappers to handle two cases of errors:
# 1) -> failure of convergence of rma inside imputation function (either for 
# numerical issues or more probably cause there are less than 4 reported studies)
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
  
  true_theta <- theta_1
  
  sim_results <- vector("list", n_sim)
  # initialize vectors to assess numerical stability
  n_success <- 0
  n_attempts <- 0
  max_attempts <- 5 * n_sim
  redraws_total <- 0
  
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
    
    res_biv <- rma.mv(
      yi,
      V = V,
      mods = ~ outcome - 1,
      random = ~ outcome | Study_id,
      struct = "UN",
      data = res_biv_long,
      method = "REML",
      tau2 = tau2_val,
      rho = rho_b
    )
    
    list(
      biv = c(est = res_biv$beta[1],
              l   = res_biv$ci.lb[1],
              u   = res_biv$ci.ub[1])
    )
  }
  
  # =========================================================
  # Main loop
  # =========================================================
  while (n_success < n_sim && n_attempts < max_attempts) {
    
    n_attempts <- n_attempts + 1
    
    res_df <- tryCatch({
      
      n_redraws <- 0
      
      repeat {
        # generate bivariate meta-analysis for each scenario
        full_data <- generate_bivariate_ma(
          K = K,
          theta = c(theta_1, theta_2),
          tau2 = c(tau2_val, tau2_val),
          rho_b = rho_b,
          rho_w = rho_w
        )
        # introduce ORB on outcome 1
        obs_data <- impose_orb(
          full_data,
          p1 = p1,
          delta_sim = delta_sim,
          select_type = select_type,
          orb.se = TRUE,
          theta_1 = theta_1,
          tau2_val = tau2_val,
          n_arm = 50
        )
        
        if (sum(!is.na(obs_data$O1_yi)) >= 4) break
        n_redraws <- n_redraws + 1
      }
      
      redraws_total <- redraws_total + n_redraws
      
      # -----------------------------------------------------
      # 1. FULL (always from full_data)
      # -----------------------------------------------------
      full_res <- compute_full(full_data)
      
      # -----------------------------------------------------
      # 2. NAIVE (always from reported data)
      # -----------------------------------------------------
      # Isolate only the observed rows
      observed_uni <- obs_data[!is.na(obs_data$O1_yi), ]
      K_observed   <- nrow(observed_uni)
      
      # Match the complete cases logic from your imputation function:
      # Find rows where both outcomes are reported to calculate a clean Pearson correlation
      complete_cases <- which(!is.na(observed_uni$O1_yi) & !is.na(observed_uni$O2_yi))
      
      if (length(complete_cases) >= 4) {
        rho_hat <- cor(observed_uni$O1_yi[complete_cases],
                       observed_uni$O2_yi[complete_cases])
      } else {
        # Fallback default if selection cuts too many complete pairs
        rho_hat <- rho_b 
      }
      
      # Build long framework for the bivariate model
      res_naive_biv_long <- data.frame(
        Study_id = rep(observed_uni$Study_id, each = 2),
        outcome  = factor(rep(c("O1","O2"), times = K_observed)),
        yi       = as.numeric(t(as.matrix(observed_uni[, c("O1_yi","O2_yi")]))),
        sei      = as.numeric(t(as.matrix(observed_uni[, c("O1_sei","O2_sei")])))
      )
      
      V_list <- lapply(1:K_observed, function(j) {
        v1 <- observed_uni$O1_sei[j]^2
        v2 <- observed_uni$O2_sei[j]^2
        cov12 <- rho_w * sqrt(v1) * sqrt(v2)
        matrix(c(v1, cov12, cov12, v2), 2, 2)
      })
      V_naive <- as.matrix(Matrix::bdiag(V_list))
      
      # Fit bivariate model using the calculated rho_hat variable
      # and fixed tau2
      res_naive_biv <- rma.mv(
        yi,
        V = V_naive,
        mods = ~ outcome - 1,
        random = ~ outcome | Study_id,
        struct = "UN",
        tau2 = tau2_val, 
        rho = rho_hat,                
        data = res_naive_biv_long,
        method = "REML",
        control = list(rel.tol = 1e-5,
                       maxiter = 200)
      )
      
      # -----------------------------------------------------
      # 3. ADJUSTED
      # -----------------------------------------------------
      # impute unreported standard errors
      obs_data <- impute_missing_se(obs_data,
                                    "O1_yi",
                                    "O1_sei",
                                    "n_total")
      # univariate imputation
      mi_uni <- run_univariate_imputation(
        obs_data,
        theta_col = "O1_yi",
        se_col = "O1_sei",
        m = M_imputations
      )
      
      adj_uni <- safe_adj_uni(mi_uni, delta_est, select_type)
      
      # bivariate imputation
      mi_biv <- run_bivariate_imputation(
        obs_data,
        theta_cols = c("O1_yi","O2_yi"),
        se_cols = c("O1_sei","O2_sei"),
        rho_w = rho_w,
        tau2_val = tau2_val,
        m = M_imputations
      )
      
      adj_biv <- safe_adj_biv(mi_biv, delta_est, select_type)
      
      # -----------------------------------------------------
      # OUTPUT
      # -----------------------------------------------------
      output_row <- data.frame(
        
        full       = full_res$biv["est"],
        full_ci_l  = full_res$biv["l"],
        full_ci_u  = full_res$biv["u"],
        
        naive_biv       = as.numeric(res_naive_biv$beta[1]),
        naive_biv_ci_l  = as.numeric(res_naive_biv$ci.lb[1]),
        naive_biv_ci_u  = as.numeric(res_naive_biv$ci.ub[1]),
        
        uni       = as.numeric(adj_uni["est"]),
        uni_ci_l  = as.numeric(adj_uni["ci_l"]),
        uni_ci_u  = as.numeric(adj_uni["ci_u"]),
        
        biv       = as.numeric(adj_biv["est"]),
        biv_ci_l  = as.numeric(adj_biv["ci_l"]),
        biv_ci_u  = as.numeric(adj_biv["ci_u"]),
        
        u_ess = as.numeric(adj_uni["ess"]),
        b_ess = as.numeric(adj_biv["ess"]),
        b_f   = as.numeric(adj_biv["fail"])
      )
      
      output_row
      
  },error = function(e) {
    print(conditionMessage(e))
    NULL
  }
  )
    
    if (!is.null(res_df)) {
      n_success <- n_success + 1
      sim_results[[n_success]] <- res_df
    }
  }
  
  if (n_success == 0) return(NULL)
  
  res_df <- do.call(rbind, sim_results[1:n_success])
  
  failure_rate <- (n_attempts - n_success) / n_attempts
  
  avg_redraws <- if (n_success > 0) {
    redraws_total / n_success
  } else {
    NA_real_
  }
  
  summary_row <- data.frame(
    # scenario parameters
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
    
    # numerical stability parameters
    N_Successful = n_success,
    N_Attempts = n_attempts,
    Failure_Rate = failure_rate,
    Avg_Redraws = avg_redraws,
    
    # -------------------------------------------------------
    # Complete Data Estimate - benchmark
    # -------------------------------------------------------
    Bias_Full = mean(res_df$full, na.rm = TRUE) - true_theta,
    
    MSE_Full = mean((res_df$full - true_theta)^2, na.rm = TRUE),
    
    Coverage_Full = mean(res_df$full_ci_l <= true_theta &
                         res_df$full_ci_u >= true_theta, na.rm = TRUE),
    
    CI_Width_Full = mean(res_df$full_ci_u - res_df$full_ci_l, na.rm = TRUE),
    
    # -------------------------------------------------------
    # Naive estimate (biased under ORB)
    # -------------------------------------------------------
    
    Bias_Naive_Biv = mean(res_df$naive_biv, na.rm = TRUE) - true_theta,
    
    MSE_Naive_Biv = mean((res_df$naive_biv - true_theta)^2, na.rm = TRUE),
    
    Coverage_Naive_Biv = mean(res_df$naive_biv_ci_l <= true_theta &
                              res_df$naive_biv_ci_u >= true_theta, na.rm = TRUE),
    
    CI_Width_Naive_Biv = mean(res_df$naive_biv_ci_u - res_df$naive_biv_ci_l, na.rm = TRUE),
    
    # -------------------------------------------------------
    # Adjusted univariate and bivariate estimate (ORB method)
    # -------------------------------------------------------
    Bias_Adj_Uni = mean(res_df$uni, na.rm = TRUE) - true_theta,
    
    MSE_Adj_Uni = mean((res_df$uni - true_theta)^2, na.rm = TRUE),
    
    Coverage_Adj_Uni = mean(res_df$uni_ci_l <= true_theta &
                            res_df$uni_ci_u >= true_theta, na.rm = TRUE),
    
    CI_Width_Adj_Uni = mean(res_df$uni_ci_u - res_df$uni_ci_l, na.rm = TRUE),
    
    Bias_Adj_Biv = mean(res_df$biv, na.rm = TRUE) - true_theta,
    
    MSE_Adj_Biv = mean((res_df$biv - true_theta)^2, na.rm = TRUE),
    
    Coverage_Adj_Biv = mean(res_df$biv_ci_l <= true_theta &
                            res_df$biv_ci_u >= true_theta, na.rm = TRUE),
    
    CI_Width_Adj_Biv = mean(res_df$biv_ci_u - res_df$biv_ci_l, na.rm = TRUE),
    
    # -------------------------------------------------------
    # Diagnostics
    # -------------------------------------------------------
    Fail_Rate_Biv = mean(res_df$b_f, na.rm = TRUE)
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
write_csv(final_metrics_df,  file.path(out_dir, "data_simulation.csv"))
cat("Complete!\n")


