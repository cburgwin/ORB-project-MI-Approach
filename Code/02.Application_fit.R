# load source code and data
source("00.Functions.R")
df_topiramate <- readRDS("Data/df_topiramate.RDS")


# load libraries
library(dplyr)
library(tidyr)
library(tidyverse)


# set seed
set.seed(1)

# set parameters
measures <- c("OR", "RR")
rhos <- c(0, 0.7) # that needs to be used for multivariate
m_imputations <- 1000


# -------------------------------------------------------------------------
# Univariate Table Generation 
# -------------------------------------------------------------------------

results_list_uni <- list()

#loop the two measures
for (meas in measures) {

  #loop the two outcomes
  for (out_num in 1:2) {

    cat(sprintf("Running Univariate Imputation: Measure = %s, Outcome = %s\n", meas, out_num))

    theta_name <- paste0("log_", meas, out_num)
    se_name <- paste0("se_", meas, out_num)

    # imputation
    mi <- run_univariate_imputation(
      data = df_topiramate,
      theta_col = theta_name,
      se_col = se_name,
      m = m_imputations
      )

    #  Extract Naive into a df matching adj_univariate output
    naive_df <- data.frame(
      Approach = "Naive estimate",
      Estimate = as.numeric(mi$res_naive$beta),
      SE       = mi$res_naive$se,
      CI_Lower = mi$res_naive$ci.lb,
      CI_Upper = mi$res_naive$ci.ub
    )

    #  Adjusted estimate for - Effect measure and Z-score
    adj_eff <- adj_univariate(mi,
                              delta = 0.5,
                              select_type = "effect",
                              track.ess = FALSE)

    adj_z   <- adj_univariate(mi,
                              delta = 0.5,
                              select_type = "zscore",
                              track.ess = FALSE)

    # Stack
    tmp <- rbind(naive_df,
                 adj_eff,
                 adj_z)

    # add info columns (correlation is empty)
    tmp$Correlation <- ""
    tmp$Measure     <- paste0("log ", meas)
    tmp$Outcome     <- paste0("O", out_num) # e.g., "O1" or "O2"

    results_list_uni[[length(results_list_uni) + 1]] <- tmp
  }
}

final_uni_df <- do.call(rbind,
                        results_list_uni)


# -------------------------------------------------------------------------
# Formatting 
# -------------------------------------------------------------------------

table_2_uni_repro <- final_uni_df  |>
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),

    Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  )  |>
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95)  |>
  tidyr::pivot_wider(
    names_from = Outcome,
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  )  |>
  dplyr::select(Correlation, Measure, Approach,
                `Outcome 1_Est_SE`, `Outcome 1_CI_95`,
                `Outcome 2_Est_SE`, `Outcome 2_CI_95`)

# Change the Correlation column to say "Univariate" so it looks nice in the final table
table_2_uni_repro$Correlation <- "Univariate"




# -------------------------------------------------------------------------
# Bivariate Table Generation (r = 0 and r = 0.7)
# -------------------------------------------------------------------------

results_list_biv <- list()


# loop over the rhos
for (rho in rhos) {

  #loop over the two measures (here we dont loop for the outcomes since they are processed together)
  for (m in measures) {

    cat(sprintf("Running Bivariate Imputation: Measure = %s, rho_w = %s\n", m, rho))

    theta_cols <- paste0("log_", m, c(1, 2))
    se_cols <- paste0("se_", m, c(1, 2))

    # Run imputation
    mi_biv <- run_bivariate_imputation(
      data = df_topiramate,
      theta_cols = theta_cols,
      se_cols = se_cols,
      rho_w = rho,
      m = m_imputations
    )

    # create naive df similar to the one produced by adj_bivariate
    naive_df <- data.frame(
      Outcome   = c("O1", "O2"),
      Approach  = "Naive estimate",
      Estimate  = as.numeric(mi_biv$res_naive$beta),
      SE        = sqrt(diag(mi_biv$res_naive$vb)),
      CI_Lower  = mi_biv$res_naive$ci.lb,
      CI_Upper  = mi_biv$res_naive$ci.ub
    )

    # Calculate both Adjustments (each return a 2-row dataframe for O1 and O2)
    adj_eff <- adj_bivariate(mi_biv,
                             delta = 0.5,
                             select_type = "effect",
                             track.ess = FALSE,
                             track.failed.proportion = FALSE)

    adj_z   <- adj_bivariate(mi_biv,
                             delta = 0.5,
                             select_type = "zscore",
                             track.ess = FALSE,
                             track.failed.proportion = FALSE)

    # stack them
    tmp <- rbind(naive_df,
                 adj_eff,
                 adj_z)

    #add explaining columns
    tmp$Correlation <- paste0("r = ", rho)
    tmp$Measure     <- paste0("log ", m)

    results_list_biv[[length(results_list_biv) + 1]] <- tmp
  }
}

final_biv_df <- do.call(rbind, results_list_biv)

# -------------------------------------------------------------------------
# Formatting
# -------------------------------------------------------------------------

table_2_biv_repro <- final_biv_df  |>
  dplyr::mutate(
    Outcome = ifelse(Outcome == "O1", "Outcome 1", "Outcome 2"),

    Est_SE = paste0(sprintf("%.2f", Estimate), " (", sprintf("%.2f", SE), ")"),
    CI_95  = paste0(sprintf("%.2f", CI_Lower), " to ", sprintf("%.2f", CI_Upper))
  )  |>
  dplyr::select(Correlation, Measure, Approach, Outcome, Est_SE, CI_95)  |>
  tidyr::pivot_wider(
    names_from = Outcome,
    values_from = c(Est_SE, CI_95),
    names_glue = "{Outcome}_{.value}"
  )  |>
  dplyr::select(Correlation, Measure, Approach,
                `Outcome 1_Est_SE`, `Outcome 1_CI_95`,
                `Outcome 2_Est_SE`, `Outcome 2_CI_95`)


# ----------------------------------
# FINAL TABLE 2
# ----------------------------------
FULL_TABLE_2 <- rbind(table_2_uni_repro, table_2_biv_repro)

saveRDS(FULL_TABLE_2, file = "Table_2.rds")
cat("Table 2 created")



# -------------------------------------------------------------------------
# Run Code for Plots
# -------------------------------------------------------------------------

####### Univariate: only for seizure freedom


res_uni_plot <- list()
deltas <- seq(from = 0, to = 1.3, by = 0.1)
measures <- c("OR", "RR")
m_imputations <- 1000

# Loop 1: Iterate over measures first and impute once per measure
for (m in measures) {

  theta_cols <- paste0("log_", m, "1")
  se_cols <- paste0("se_", m, "1")

  cat(sprintf("Generating imputed data for Measure = %s\n", m))

  # Run imputation ONCE here. The random pool is now locked in.
  mi <- run_univariate_imputation(
    data = df_topiramate,
    theta_col = theta_cols,
    se_col = se_cols,
    m = m_imputations
  )

  # Extract the static Naive estimate for this measure
  naive_df <- data.frame(
    Approach = "Naive estimate",
    Estimate = as.numeric(mi$res_naive$beta),
    SE       = mi$res_naive$se,
    CI_Lower = mi$res_naive$ci.lb,
    CI_Upper = mi$res_naive$ci.ub
  )

  # Loop 2: Now test different deltas using the same imputed data
  for (delta in deltas) {

    cat(sprintf("  -> Calculating adjustments for delta = %s\n", delta))

    # Importance sampling calculations over the static 'mi' object
    adj_eff <- adj_univariate(mi,
                              delta = delta,
                              select_type = "effect",
                              track.ess = FALSE)

    adj_z   <- adj_univariate(mi,
                              delta = delta,
                              select_type = "zscore",
                              track.ess = FALSE)

    # Stack the results
    tmp <- rbind(naive_df, adj_eff, adj_z)

    # Add metadata
    tmp$Measure   <- paste0("log ", m)
    tmp$Selection <- delta  # Saved as numeric directly to avoid your earlier character bug!

    # Append to list
    res_uni_plot[[length(res_uni_plot) + 1]] <- tmp
  }
}

plot_uni_df <- do.call(rbind,
                       res_uni_plot)

# bring data into wide format
naive_df <- plot_uni_df %>%
  filter(Approach == "Naive estimate") %>%
  dplyr::select(Measure, Selection,
         Estimate_unadjusted = Estimate,
         CI_unadjusted_lower = CI_Lower,
         CI_unadjusted_upper = CI_Upper)

adjusted_df <- plot_uni_df %>%
  filter(Approach != "Naive estimate") %>%
  # Create the Z_score column your plot expects (TRUE if zscore, FALSE if effect)
  mutate(Z_score = ifelse(Approach == "Selection on zscore", TRUE, FALSE)) %>%
  dplyr::select(Measure, Selection, Z_score,
         Estimate_adjusted = Estimate,
         CI_adjusted_lower = CI_Lower,
         CI_adjusted_upper = CI_Upper)

final_uni_df <- adjusted_df %>%
  left_join(naive_df, by = c("Selection", "Measure"))

# save data
write_csv(final_uni_df, "data/data_uni_df_O1.csv")



############## Bivariate #####################

res_biv_df <- list()
rhos <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.67, "estimated", "studyspecific")
deltas <- seq(from = 0, to = 1.3, by = 0.1) 
m_imputations <- 1000


# loop over the rhos 
for (rho in rhos) {
  
  #loop over the two measures (here we dont loop for the outcomes since they are processed together)
  for (m in measures) {
    
    cat(sprintf("\n--- Running Bivariate Imputation Pool: Measure = %s, correlation setting = %s ---\n", m, rho))
    
    theta_cols <- paste0("log_", m, c(1, 2))
    se_cols <- paste0("se_", m, c(1, 2))
    
    # Run imputation
    mi_biv <- run_bivariate_imputation(
      data = df_topiramate, 
      theta_cols = theta_cols, 
      se_cols = se_cols, 
      rho_w = rho, 
      m = m_imputations
    )
    
    # create naive df similar to the one produced by adj_bivariate
    naive_df <- data.frame(
      Outcome   = c("O1", "O2"),
      Approach  = "Naive estimate",
      Estimate  = as.numeric(mi_biv$res_naive$beta),
      SE        = sqrt(diag(mi_biv$res_naive$vb)),
      CI_Lower  = mi_biv$res_naive$ci.lb,
      CI_Upper  = mi_biv$res_naive$ci.ub
    )
    
    for(delta in deltas) {
      
      cat(sprintf("  -> Calculating Bivariate Adjustments for delta = %s\n", delta))
    
      # Calculate both Adjustments (each return a 2-row dataframe for O1 and O2)
      adj_eff <- adj_bivariate(mi_biv, 
                               delta = delta, 
                               select_type = "effect",
                               track.ess = FALSE, 
                               track.failed.proportion = FALSE)
      
      adj_z   <- adj_bivariate(mi_biv, 
                               delta = delta, 
                               select_type = "zscore", 
                               track.ess = FALSE, 
                               track.failed.proportion = FALSE)
      
      # stack them
      tmp <- rbind(naive_df,
                   adj_eff,
                   adj_z)
      
      #add explaining columns
      tmp$Correlation <- paste0(rho)
      tmp$Measure     <- paste0("log ", m)
      tmp$Selection   <- delta
      
      res_biv_df[[length(res_biv_df) + 1]] <- tmp
    }
  }
}

plot_biv_df <- do.call(rbind,
                       res_biv_df)

# save data
write_csv(plot_biv_df, "data_biv_df.csv")


