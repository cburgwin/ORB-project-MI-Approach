# load packages
library(metafor)
library(dplyr)
library(MASS)



# read data in 
df_topiramate <- data.frame(
  Study_id = 1:12,
  Author = c("Ben-Menachem", "Elterman", "Faught", "Guberman", "Korean", 
             "Privitera", "Rosenfeld", "Sharief", "Tassinari", "Yen", 
             "Zhang", "Coles"),
  
  # Sample Sizes (T = Treated, C = Control)
  N_t = c(28, 41, 136, 171, 91, 143, 167, 23, 30, 23, 46, 52),
  N_c = c(28, 45, 45, 92, 86,  47, 42, 24, 30, 23, 40, 51),
  
  # Outcome 1: 50% Seizure Reduction
  T1 = c(12, 16, 54, 77, 45, 58, 86, 8, 14, 11, 22, NA),
  C1 = c(0, 9, 8, 22, 11, 4, 8, 2, 3, 3, 3, NA),
  
  # Outcome 2: Seizure Freedom
  T2 = c (NA, 4, NA, 10, 7, NA, NA, 2, 0, NA, 0, NA),
  C2 = c (NA, 2, NA,  2, 1, NA, NA, 0, 0, NA, 0, NA)
)

compute_effects <- function(data) {
  
  # ----------------------------------
  # OUTCOME 1: 50% Seizure Reduction
  # ----------------------------------
  #  Log Odds Ratio
  data <- escalc(measure = "OR", 
                 ai = T1, 
                 n1i = N_t, 
                 ci = C1,
                 n2i = N_c, 
                 data = data,
                 var.names = c("log_OR1", "var_OR1"), 
                 append = TRUE, 
                 add = 0.5, to = "only0")
  # Convert variance to SE
  data$se_OR1 <- sqrt(data$var_OR1) 
  
  # Log Relative Risk
  data <- escalc(measure = "RR", 
                 ai = T1, 
                 n1i = N_t, 
                 ci = C1, 
                 n2i = N_c, 
                 data = data, 
                 var.names = c("log_RR1", "var_RR1"), 
                 append = TRUE,
                 add = 0.5, 
                 to = "only0")
  # Convert variance to SE
  data$se_RR1 <- sqrt(data$var_RR1)
  
  # --------------------------------------
  # OUTCOME 2: Seizure Freedom
  # --------------------------------------
  # Log Odds Ratio
  data <- escalc(measure = "OR", 
                 ai = T2, 
                 n1i = N_t,
                 ci = C2, 
                 n2i = N_c, 
                 data = data, 
                 var.names = c("log_OR2", "var_OR2"),
                 append = TRUE,
                 add = 0.5, 
                 to = "only0")
  # Convert variance to SE
  data$se_OR2 <- sqrt(data$var_OR2)
  
  # Log Relative Risk
  data <- escalc(measure = "RR", 
                 ai = T2, 
                 n1i = N_t, 
                 ci = C2, 
                 n2i = N_c, 
                 data = data, 
                 var.names = c("log_RR2", "var_RR2"),
                 append = TRUE,
                 add = 0.5, 
                 to = "only0")
  # Convert variance to SE
  data$se_RR2 <- sqrt(data$var_RR2)
  
  data_out <- data |> 
    dplyr::mutate(n_total = N_t + N_c) |> 
    dplyr::select(Study_id,
                  n_total, 
                  log_OR1, 
                  se_OR1, 
                  log_RR1, 
                  se_RR1, 
                  log_OR2,
                  se_OR2, 
                  log_RR2, 
                  se_RR2)  
  
  return(data_out)
}

df_topiramate <- compute_effects(df_topiramate)






# slightly different function than the one in 00.functions.R
impute_missing_se <- function(data, 
                              target_theta_col, 
                              target_se_col) {
  
  # Identify reported and unreported indices
  rep_idx <- which(!is.na(data[[target_theta_col]]))
  unrep_idx <- which(is.na(data[[target_theta_col]]))
  
  # If there are no missing values for this outcome, return data 
  if(length(unrep_idx) == 0) {
    return(data)
  }
  
  # Calculate k_hat using only the reported studies (Eq. 6)
  # Precision is the inverse of the variance (1 / SE^2)
  precisions <- 1 / (data[[target_se_col]][rep_idx]^2)
  k_hat <- sum(precisions) / sum(data$n_total[rep_idx])
  
  # Impute the missing standard errors
  data[[target_se_col]][unrep_idx] <- sqrt(1 / (k_hat * data$n_total[unrep_idx]))
  
  return(data)
}


# Impute missing SEs for each of the outcomes and effect measures 

df_topiramate <- impute_missing_se(
  data = df_topiramate, 
  target_theta_col = "log_RR1", 
  target_se_col = "se_RR1"
)

df_topiramate <- impute_missing_se(
  data = df_topiramate, 
  target_theta_col = "log_RR2", 
  target_se_col = "se_RR2"
)

df_topiramate <- impute_missing_se(
  data = df_topiramate, 
  target_theta_col = "log_OR1", 
  target_se_col = "se_OR1"
)

df_topiramate <- impute_missing_se(
  data = df_topiramate, 
  target_theta_col = "log_OR2", 
  target_se_col = "se_OR2"
)

# save dataset 
saveRDS(df_topiramate, file = "df_topiramate.RDS")

