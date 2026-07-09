
library(MASS)
library(Matrix)
library(metafor)


#---------------------------------------
# Imputing missing se
# --------------------------------


impute_missing_se <- function(data,
                              target_theta_col,
                              target_se_col,
                              n_col) {
  
  # Identify reported and unreported indices
  rep_idx <- which(!is.na(data[[target_theta_col]]))
  unrep_idx <- which(is.na(data[[target_theta_col]]))
  
  # If there are no missing values for this outcome, return data 
  if(length(unrep_idx) == 0) {
    return(data)
  }
  
  # Calculate k_hat using only the reported studies (Eq. 6)
  precisions <- 1 / (data[[target_se_col]][rep_idx]^2)
  num <- sum(precisions)
  den <- sum(data[[n_col]][rep_idx])
  
  k_hat <- num / den
  
  # Impute the missing standard errors
  data[[target_se_col]][unrep_idx] <- sqrt(1 / (k_hat * data[[n_col]][unrep_idx]))
  
  return(data)
}



# ------------------------------------------
# Univariate
# -----------------------------------------


# Univariate Imputation under MAR (remember it assumes SE already imputed)
run_univariate_imputation <- function(data, 
                                      theta_col = "log_RR2", 
                                      se_col = "se_RR2",
                                      m = 1000, 
                                      model_type = "REML") {
  
  # Get reported and unreported indexes
  rep_idx_start <- which(!is.na(data[[theta_col]]))
  unrep_idx_start <- which(is.na(data[[theta_col]]))
  
  # Order correctly (first  rep, then unrep)
  data_rep <- data[rep_idx_start, ]
  data_unrep <- data[unrep_idx_start, ]
  data_ordered <- rbind(data_rep, 
                        data_unrep)
  
  # determine number of reported studies and unreported studies
  n_rep <- nrow(data_rep)
  n_unrep <- nrow(data_unrep)
  K <- nrow(data_ordered)

  
  # Indexes after reordering
  rep_idx <- 1:n_rep
  unrep_idx <- (n_rep + 1):K
  
  # if nothing missing, return null
  if (length(unrep_idx) == 0) return(NULL)
  
  # but also we need at least 4 reported studies
  if (length(rep_idx) < 4) {
    warning("Fewer than 4 reported studies. Cannot estimate heterogeneity. Skipping.")
    return(NULL) 
  }
  
  
  # Fit REML on the reported data (Naive Estimate), save estimate & SE
  res_naive <- metafor::rma(yi = data_rep[[theta_col]],
                            sei = data_rep[[se_col]],
                            method = model_type)
  theta_MA <- as.numeric(res_naive$beta)
  tau2_hat <- as.numeric(res_naive$tau2)
  SE2_mean <- res_naive$se^2
  
  # determine variance matrix
  J_K <- matrix(1, K, K) 
  
  #Total cov matrix
  
  Sigma <- diag(data_ordered[[se_col]]^2 + tau2_hat) + (SE2_mean * J_K) 
  
  
  # Partitioning (Eq. 10)
  Sigma_RR <- Sigma[rep_idx, rep_idx]
  Sigma_UU <- Sigma[unrep_idx, unrep_idx]
  Sigma_UR <- Sigma[unrep_idx, rep_idx]
  Sigma_RU <- Sigma[rep_idx, unrep_idx]
  
  # Conditionals (Eq. 11) 
  theta_R <- data_rep[[theta_col]]
  inv_Sigma_RR <- solve(Sigma_RR) 
  
  mu_cond <- as.numeric(theta_MA + Sigma_UR %*% inv_Sigma_RR %*% (theta_R - theta_MA))
  Sigma_cond <- Sigma_UU - Sigma_UR %*% inv_Sigma_RR %*% Sigma_RU
  
  # Generate M Imputations (we get a matrix m x n_unrep)
  imputed_draws <- MASS::mvrnorm(n = m, 
                                 mu = mu_cond, 
                                 Sigma = Sigma_cond)
  
  
  return(list(
    imp_draws  = imputed_draws,
    unrep_idx    = unrep_idx,
    res_naive    = res_naive,
    data_ordered = data_ordered,
    theta_col    = theta_col,  
    se_col       = se_col     
  ))
}



# --------------------------
#  Importance Sampling 
# --------------------------

adj_univariate <- function(mi_results, 
                           delta = 0.5, 
                           select_type = "zscore",
                           model_type = "REML", 
                           track.ess = TRUE) {
  
  imp_draws  <- mi_results$imp_draws
  data_ordered <- mi_results$data_ordered
  unrep_idx    <- mi_results$unrep_idx
  
  # ensure imp_draws is a matrix
  if (!is.matrix(imp_draws)) imp_draws <- matrix(imp_draws, ncol = length(unrep_idx))
  M <- nrow(imp_draws)
  
  
  # if nothing missing, return null, will handle better in the outer loop of the simulation
  if (length(unrep_idx) == 0) return(NULL)
  
  
  theta_col    <- mi_results$theta_col  #  from mi_results
  se_col       <- mi_results$se_col     #  from mi_results
  tau2_val <- mi_results$res_naive$tau2 # from mi_results
  
  # extract the se useful for the z scores
  sei_unrep <- data_ordered[[se_col]][unrep_idx]   
  
  # small check to check that dimensions match
  if (length(sei_unrep) != ncol(imp_draws) ) stop ("Something is wrong in the dimensions of the unreported studies")
  
  if (select_type == "zscore") {
    # divides every row in 'imp_draws' by the 'sei_unrep' vector to create z score
    z_matrix <- sweep(imp_draws, 
                      2, 
                      sei_unrep,
                      FUN = "/")
    log_weights <- - delta * rowSums(z_matrix)
  } else {
    log_weights <- - delta * rowSums(imp_draws)
  }
  
  # rma on each complete dataset 
  # vector of thetas that is like [c(reported), NA, NA,NA ...]
  
  yi_base <- data_ordered[[theta_col]]
  
  fits <- lapply(seq_len(M), function(m) {
    yi_complete <- yi_base
    yi_complete[unrep_idx] <- imp_draws[m, ] # complete the vector of thetas with the simulated draws from here
    tryCatch(
      rma(yi = yi_complete, 
          sei = data_ordered[[se_col]], 
          method = model_type,
          tau2 = tau2_val,
          control = list(stepadj = 0.1, 
                         rel.tol = 1e-5, 
                         maxiter = 200)),
      error = function(e) NULL
    )
  })
  
  # only take not failed 
  valid       <- !sapply(fits, is.null)
  fits        <- fits[valid]
  log_weights <- log_weights[valid]
  
  # extract estimates
  theta_MA_m <- sapply(fits, function(r) as.numeric(r$beta))
  var_MA_m   <- sapply(fits, function(r) r$se^2)
  
  # Normalize weights (log-sum-exp)
  max_log <- max(log_weights)
  w_norm <- exp(log_weights - max_log) / sum(exp(log_weights - max_log))
  
  
  # Adjusted (Eq. 16)
  theta_adj  <- sum(w_norm * theta_MA_m)
  
  # Rubin's Rules (Eq. 18)
  var_within <- sum(w_norm * var_MA_m)
  var_between <- sum(w_norm * (theta_MA_m - theta_adj)^2)
  total_var <- var_within + var_between
  
  
  out <- data.frame(
    Approach = paste("Selection on", select_type),
    Estimate = theta_adj,
    SE       = sqrt(total_var),
    CI_Lower = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper = theta_adj + 1.96 * sqrt(total_var)
  )
  
  if (track.ess) {
    ess <- 1 / sum(w_norm^2)
    out$ess <- ess
  }
  
  return(out)
  
}



# ---------------------------
## Bivariate functions 
# ----------------------------

run_bivariate_imputation <- function(data,
                                     theta_cols, 
                                     se_cols, 
                                     rho_w = 0.4, 
                                     m = 1000, 
                                     tau2_val,
                                     model_type = "REML") {
  
  
  K <- nrow(data)
  
  # --- Compute or Assign the Correlation Vector (r_vec) ---
  if (is.character(rho_w) && rho_w %in% c("estimated", "studyspecific")) {
    
    # Identify complete cases using the dynamic theta_cols passed into the function
    complete_cases <- which(!is.na(data[[theta_cols[1]]]) & !is.na(data[[theta_cols[2]]]))
    
    if (length(complete_cases) < 4) {
      stop("Not enough complete pairs in the data to estimate a Pearson correlation.")
    }
    
    # Estimate correlation and Fisher CI directly from your input data columns
    cor_test <- cor.test(data[[theta_cols[1]]][complete_cases],
                         data[[theta_cols[2]]][complete_cases],
                         method = "pearson")
    
    rho_hat  <- as.numeric(cor_test$estimate)
    
    CI <- c(cor_test$conf.int[1], cor_test$conf.int[2])  
    CIz <- atanh(CI) # Fisher z-transformation
    
    # Mean and SE on Fisher z-scale
    mu_z <- mean(CIz)
    se_z <- diff(CIz) / (2 * 1.96)
    
    if (rho_w == "estimated") {
      # One global correlation with uncertainty drawn once
      r_global <- tanh(rnorm(1, 
                             mean = mu_z,
                             sd = se_z)) 
      r_vec <- rep(r_global, K)
    } else if (rho_w == "studyspecific") {
      # One unique correlation sampled per study
      r_vec <- tanh(rnorm(K,
                          mean = mu_z, 
                          sd = se_z))
    }
    
  } else {
    # Use fixed value provided numerically (e.g., 0.4)
    rho_hat <- as.numeric(rho_w)
    if (rho_hat <= -1 || rho_hat >= 1) stop("Correlation must be strictly between -1 and 1.")
    r_vec <- rep(rho_hat, K)
  }
  
  # study-major vector of length 2K (remember that as.numeric flatten by columns)
  # so we have (Out1study1, out2study1, ...... out1studyK, out2studyK)
  y_vec <- as.numeric(t(as.matrix (data[, theta_cols]))) 
  se_vec <- as.numeric(t( as.matrix(data[, se_cols])))
  
  # Identify reported and unreported studies
  rep_idx <- which(!is.na(y_vec))
  unrep_idx <- which(is.na(y_vec))
  
  # if no study outcomes are unreported, return null
  if (length(unrep_idx) == 0) return(NULL)
  
  # but also we need at least 4 reported studies
  if (length(rep_idx) < 4) {
    warning("Fewer than 4 reported studies. Cannot estimate heterogeneity. Skipping.")
    return(NULL) 
  }
  
  # build the within study covariance 
  V_list <- list()
  for (i in 1:K) { # one for each study 
    v1 <- data[[se_cols[1]]][i]^2
    v2 <- data[[se_cols[2]]][i]^2
    cov_12 <- r_vec[i] * sqrt(v1) * sqrt(v2) 
    V_list[[i]] <- matrix(c(v1, cov_12, cov_12, v2), 2, 2)
  }
  
  V_full <- as.matrix(Matrix::bdiag(V_list))
  
  # Create Long Data 
  # --> df that has row 1 = study 1 outcome 1, row 2 = study 1 outcome 2, row 3 = study 2 outcome 1 ....
  data_long <- data.frame(
    Study_id = rep(data$Study_id, each = 2),
    outcome  = factor(rep(c("O1", "O2"), times = K)),
    yi       = y_vec,
    sei      = se_vec,
    obs_id   = 1:(2*K)
  )
  
  # built within covariance matrix only based on reported studies
  reported_long <- data_long[rep_idx, ]
  V_rep <- V_full[rep_idx, rep_idx]
  
  # Fit Naive Bivariate REML  (Eq. 12 & 13)
  res_naive <- rma.mv(yi,
                      V = V_rep, 
                      mods = ~ outcome - 1, 
                      random = ~ outcome | Study_id, 
                      struct = "UN",  # true effects varies for each study id (heterogeneity)
                      data = reported_long, 
                      method = model_type, 
                      rho = rho_hat,
                      tau2 = tau2_val,
                      control=list(rel.tol=1e-8,
                                   maxiter=200))
  
  # Construct Total Covariance Matrix Sigma
  Cov_theta_MA <- vcov(res_naive)
  # Create a K x K matrix of ones
  J_K <- matrix(1, K, K)
  # save estimates from naive bivariate MA
  tau2_1 <- res_naive$tau2[1]
  tau2_2 <- res_naive$tau2[2]
  rho_b  <- res_naive$rho
  
  # Extract Psi (Between-study variance-covariance matrix)  
  Psi <- matrix(c(
    tau2_1, rho_b * sqrt(tau2_1) * sqrt(tau2_2),
    rho_b * sqrt(tau2_1) * sqrt(tau2_2), tau2_2),
    nrow = 2, ncol = 2
    )
  I_K <- diag(K)
  
  # Use the Kronecker product to tile the 2x2 Vb_matrix across every study pair
  # This results in the 2K x 2K matrix  
  Sigma <- V_full + kronecker(I_K, Psi) + kronecker(J_K, Cov_theta_MA)
  
  # partitioning 
  Sigma_RR <- Sigma[rep_idx, rep_idx, drop = FALSE]
  Sigma_UU <- Sigma[unrep_idx, unrep_idx, drop = FALSE]
  Sigma_UR <- Sigma[unrep_idx, rep_idx, drop = FALSE]
  Sigma_RU <- Sigma[rep_idx, unrep_idx, drop = FALSE]
  
  #  Conditionals (Eq. 11)
  theta_MA_vec <- rep(coef(res_naive), times = K)  # vector 2K
  theta_R_MA <- theta_MA_vec[rep_idx] # vector K_R
  theta_U_MA <- theta_MA_vec[unrep_idx] # vector K_U
  theta_R <- y_vec[rep_idx] # get reported summary estimates
  
  inv_Sigma_RR <- solve(Sigma_RR)
  mu_cond <- theta_U_MA + Sigma_UR %*% inv_Sigma_RR %*% (theta_R - theta_R_MA)
  Sigma_cond <- Sigma_UU - Sigma_UR %*% inv_Sigma_RR %*% Sigma_RU
  
  # Generate M Imputations
  imputed_draws <- MASS::mvrnorm(n = m,
                                 mu = as.numeric(mu_cond),
                                 Sigma = Sigma_cond)
  if (length(unrep_idx) == 1) {
    imputed_draws <- matrix(imputed_draws, ncol = 1)
  }
  
  return(list(
    imp_draws    = imputed_draws,
    unrep_idx    = unrep_idx,
    res_naive    = res_naive,
    data_long    = data_long,
    V_full       = V_full,
    theta_cols   = theta_cols
  ))
}



adj_bivariate <- function(mi_results,
                          delta = 0.5,
                          select_type = "zscore",
                          model_type = "REML", 
                          track.failed.proportion = TRUE, 
                          track.ess = TRUE) {
  
  imp_draws <- mi_results$imp_draws
  data_long <- mi_results$data_long
  unrep_idx <- mi_results$unrep_idx
  V_full <- mi_results$V_full
  
  tau2_fixed <- mi_results$res_naive$tau2
  rho_fixed  <- mi_results$res_naive$rho

  if (!is.matrix(imp_draws)) imp_draws <- matrix(imp_draws, ncol = length(unrep_idx))
  M <- nrow(imp_draws)
  
  # extract the se useful for the z scores
  sei_unrep <- data_long$sei[unrep_idx]  
  
  # check dimensions
  if (ncol(imp_draws) != length(sei_unrep)) stop ("Something is wrong with the number of unreported studies")
  
  
  if (select_type == "zscore") {
    # divides every row in 'imp_draws' by the 'sei_unrep' vector
    z_matrix <- sweep(imp_draws, 2, sei_unrep, FUN = "/")
    log_weights <- -delta * rowSums(z_matrix)
  } else {
    log_weights <- -delta * rowSums(imp_draws)
  }
  
  # loop for each draw
  fits <- lapply(seq_len(M), function(m) {
    d <- data_long
    
    # create complete dataset
    d$yi[unrep_idx] <- imp_draws[m, ]
    
    #fit the results
    tryCatch( 
      rma.mv(yi,
             V = V_full,
             mods = ~ outcome - 1,
             random = ~ outcome | Study_id, 
             struct = "UN",
             data = d, 
             method = model_type,
             tau2 = tau2_fixed,
             rho = rho_fixed,
             control = list(stepadj = 0.1,
                            rel.tol = 1e-5, 
                            maxiter = 200)),
      error = function(e) NULL
    )
  })
  
  # extract valid
  valid       <- sapply(fits, function(x) !is.null(x))
  fits        <- fits[valid]
  log_weights <- log_weights[valid]
  
  if (length(fits) == 0) stop("All rma.mv fits failed")
  
  # Extract (note that sapply bind vectors together by columns)
  theta_MA_m <- t(sapply(fits, function(r) as.numeric(r$beta)))  # M x 2
  var_MA_m   <- t(sapply(fits, function(r) diag(r$vb)))           # M x 2
  
  # Normalize weights
  max_log <- max(log_weights)
  w_norm <- exp(log_weights - max_log) / sum(exp(log_weights - max_log))
  
  # Adjusted estimates 
  # remember w_norm is a vector of 1000 weights
  # theta_MA_m is a 1000 x 2 matrix
  theta_adj <- colSums(w_norm * theta_MA_m)               # this returns a VECTOR of length 2
  
  # Rubin rule
  var_within  <- colSums(w_norm * var_MA_m)
  var_between <- colSums(w_norm * sweep(theta_MA_m, 2, theta_adj, "-")^2)
  total_var   <- var_within + var_between #length 2 vector
  
  if (length(theta_adj) != length(total_var) ) stop ("Something is off in the dimension of outcomes")
  out <- data.frame(
    Outcome  = c("O1", "O2"),
    Approach = paste("Bivariate Selection on", select_type),
    Estimate = theta_adj,
    SE       = sqrt(total_var),
    CI_Lower = theta_adj - 1.96 * sqrt(total_var),
    CI_Upper = theta_adj + 1.96 * sqrt(total_var)
  )
  
  rownames(out) <- NULL
  
  if (track.ess) {
    ess <- 1 / sum(w_norm^2)
    out$ess <- ess
  }
  
  if (track.failed.proportion) { 
    failed.proportion <- sum(!valid) / M
    out$failed.proportion <- failed.proportion
  }
  return(out)
}




# --------------------------------
# Simulate data
# -------------------------------

# here we generate the full data 

# K = number of studies
# theta = vector containing the true values for the two outcomes
# tau2 = vector containing the true tau
# rho_b is the between study correlation
# rho_w withing study correlation 
# n_arm is fixed to 50 in the whole simulation 


generate_bivariate_ma <- function(K = 12,
                                  theta = c(0.4, 0.4),
                                  tau2 = c(0.06, 0.06), 
                                  rho_b = 0.4,
                                  rho_w = 0.4, 
                                  n_arm = 50) {
  
  # Between-study covariance matrix (Psi)
  cov_b <- rho_b * sqrt(tau2[1]) * sqrt(tau2[2])
  Psi <- matrix(c(tau2[1], cov_b, cov_b, tau2[2]), 2, 2)
  
  # True study-specific effects (bivariate normal)
  theta_i <- matrix(MASS::mvrnorm(n = K,
                                  mu = theta,
                                  Sigma = Psi), 
                    nrow = K, ncol = 2)  
  
  #  matrix V for Wishart distribution
  V_scale <- (1 / ((n_arm - 1) * n_arm)) * matrix(c(1, rho_w, rho_w, 1), 2, 2)
  df_wishart <- 2 * (n_arm - 1)
  
  # Initialize vectors for obs effects and standard errors
  y_obs <- matrix(NA, nrow = K, ncol = 2)
  se_obs <- matrix(NA, nrow = K, ncol = 2)
  
  for (i in 1:K) {
    #  within-study covariance matrix for study i using Wishart (Sigma_i)
    Sigma_i <- stats::rWishart(n = 1, df = df_wishart, Sigma = V_scale)[,,1]
    
    # Draw observed effects (bivariate)
    # sample the observed value from a Normal centerd in the true value (sampled above)
    # and with variance from Widhart
    y_obs[i, ] <- mvrnorm(n = 1, mu = theta_i[i, ], Sigma = Sigma_i)
    
    # Extract standard errors
    se_obs[i, 1] <- sqrt(Sigma_i[1, 1])
    se_obs[i, 2] <- sqrt(Sigma_i[2, 2])
  }
  
  data <- data.frame(
    Study_id = 1:K,
    n_total = 2 * n_arm,
    O1_yi = y_obs[, 1], O1_sei = se_obs[, 1],
    O2_yi = y_obs[, 2], O2_sei = se_obs[, 2]
  )
  
  attr(data, "theta1") <- theta[1]
  attr(data, "tau2_1") <- tau2[1]
  attr(data, "n_arm") <- n_arm
  
  return(data)
}





expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(p) log(p / (1 - p))

# Now we need to choose:
# p1 = prob of unreporting
# delta value used in the simulation 
# selection type 



impose_orb <- function(data,
                       p1 = 0.4,
                       delta_sim = 0.5, 
                       select_type = "zscore", 
                       orb.se = TRUE,
                       theta_1 = 0.4,       
                       tau2_val = 0.06,     
                       n_arm = 50) {
 
  
  n_total <- n_arm * 2 # we can safely assume that both arms have the same n
  
  
  # Determine the selection variable s_1 and its expected value E(s_1)
  if (select_type == "zscore") {
    s1 <- data$O1_yi / data$O1_sei #vector of real value 
    E_s1 <- theta_1 / sqrt((2/n_arm) + tau2_val) #expected score to define theoretical alpha
  } else if (select_type == "effect") {
    s1 <- data$O1_yi 
    E_s1 <- theta_1 # if we do selection on effect, the exp value is directly the true value
  } else {
    stop("select_type must be 'zscore' or 'effect'")
  }
  
  target_reported <- 1 - p1
  
  # Avoid logit of 1 or 0 if p1 is extreme (this should never happen in the code)
  if(target_reported >= 1) target_reported <- 0.999 
  if(target_reported <= 0) target_reported <- 0.001
  
  alpha_1 <- logit(target_reported) - (delta_sim * E_s1)
  
  # Calculate reporting probabilities for each outocme
  prob_report <- expit(alpha_1 + delta_sim * s1)
  
  # Simulate reporting or not (1 = reported, 0 = missing)
  reported_flag <- rbinom(n = nrow(data),
                          size = 1, 
                          prob = prob_report)
  
  # Apply missingness to Outcome 1
  data_orb <- data
  data_orb$O1_yi[reported_flag == 0] <- NA
  
  if (orb.se) data_orb$O1_sei[reported_flag == 0] <- NA 
  
  # add the sample size column
  data_orb$n_total <- n_total
  
  attr(data_orb, "actual_missing_rate") <- sum(reported_flag == 0) / nrow(data)
  
  return(data_orb)
}