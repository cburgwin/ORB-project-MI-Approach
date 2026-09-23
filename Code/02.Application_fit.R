# load source code and data
source("Code/00.Functions.R")
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

# proportion of failed refits 
failed_refits <- c()


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
    
    # refit the imputed datasets once, reused for both selection types
    fits_uni <- fit_imputations_uni(mi)
    failed_refits[paste0("Table2 uni ", meas, " O", out_num)] <- attr(fits_uni, "failed.proportion")

    #  Adjusted estimate for - Effect measure and Z-score
    adj_eff <- adj_univariate(mi,
                              delta = 0.5,
                              select_type = "effect",
                              track.ess = FALSE,
                              fits = fits_uni)

    adj_z   <- adj_univariate(mi,
                              delta = 0.5,
                              select_type = "zscore",
                              track.ess = FALSE,
                              fits = fits_uni)

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
    
    # refit the imputed datasets once, reused for both selection types
    fits_biv <- fit_imputations_biv(mi_biv)
    failed_refits[paste0("Table2 biv ", m, " rho_w=", rho)] <- mean(sapply(fits_biv, is.null))

    # Calculate both Adjustments (each return a 2-row dataframe for O1 and O2)
    adj_eff <- adj_bivariate(mi_biv,
                             delta = 0.5,
                             select_type = "effect",
                             track.ess = FALSE,
                             track.failed.proportion = FALSE,
                             fits = fits_biv)

    adj_z   <- adj_bivariate(mi_biv,
                             delta = 0.5,
                             select_type = "zscore",
                             track.ess = FALSE,
                             track.failed.proportion = FALSE,
                             fits = fits_biv)

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

saveRDS(FULL_TABLE_2, file = "Data/Table_2.rds")
cat("Table 2 created")



# -------------------------------------------------------------------------
# Run Code for Plots
# -------------------------------------------------------------------------

####### Univariate: both outcomes


deltas <- seq(from = 0, to = 1.3, by = 0.1)
measures <- c("OR", "RR")
m_imputations <- 1000
outcomes <- c(1,2)

#Loop 0: iterate over the two outcomes
for (o in outcomes) {
res_uni_plot <- list()   # reset per outcome, otherwise the O2 file also contains the O1 rows
# Loop 1: Iterate over measures first and impute once per measure
for (m in measures) {

  theta_cols <- paste0("log_", m, o)
  se_cols <- paste0("se_", m, o)

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
  
  # refit the imputed datasets once, reused for every delta and selection type
  fits_uni <- fit_imputations_uni(mi)
  failed_refits[paste0("Plot uni ", m, " O", o)] <- attr(fits_uni, "failed.proportion")

  # Loop 2: Now test different deltas using the same imputed data
  for (delta in deltas) {

    cat(sprintf("  -> Calculating adjustments for delta = %s\n", delta))

    # Importance sampling calculations over the static 'mi' object
    adj_eff <- adj_univariate(mi,
                              delta = delta,
                              select_type = "effect",
                              track.ess = FALSE,
                              fits = fits_uni)

    adj_z   <- adj_univariate(mi,
                              delta = delta,
                              select_type = "zscore",
                              track.ess = FALSE,
                              fits = fits_uni)

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

# file to save 
file_name <- paste0("Data/data_uni_df_O", o, ".csv")

# save data
write_csv(final_uni_df, file_name)
}


############## Bivariate #####################

res_biv_df <- list()
rhos <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.67, "estimated", "studyspecific")
deltas <- seq(from = 0, to = 1.3, by = 0.1) 
m_imputations <- 1000



set.seed(1)

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
    
    # refit the imputed datasets once, reused for every delta and selection type
    fits_biv <- fit_imputations_biv(mi_biv)
    failed_refits[paste0("Plot biv ", m, " rho_w=", rho)] <- mean(sapply(fits_biv, is.null))
    
    for(delta in deltas) {
      
      cat(sprintf("  -> Calculating Bivariate Adjustments for delta = %s\n", delta))
    
      # Calculate both Adjustments (each return a 2-row dataframe for O1 and O2)
      adj_eff <- adj_bivariate(mi_biv, 
                               delta = delta, 
                               select_type = "effect",
                               track.ess = FALSE, 
                               track.failed.proportion = FALSE,
                               fits = fits_biv)
      
      adj_z   <- adj_bivariate(mi_biv, 
                               delta = delta, 
                               select_type = "zscore", 
                               track.ess = FALSE, 
                               track.failed.proportion = FALSE,
                               fits = fits_biv)
      
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
write_csv(plot_biv_df, "Data/data_biv_df.csv")

# proportion of failed refits per imputation run (0 = every refit succeeded)
print(failed_refits)
saveRDS(failed_refits, file = "Data/failed_refits.rds")











############## Number of imputations M (Appendix) #####################

# Bivariate method, log OR and log RR, selection on the z-score and on the
# effect estimate. Within-study correlation fixed at rho_W = -0.3. tau2_1,
# tau2_2 and rho_B estimated.
# One analysis with M = 1000 imputations per measure. The M = 200 curve uses
# only the first 200 of these imputations.

library(ggplot2)

types_M   <- c("zscore", "effect")
deltas_ov <- seq(from = 0, to = 1.3, by = 0.1)
res_ov    <- list()
measures <- c("OR", "RR")

set.seed(1)

for (m in measures) {
  mi_ov   <- run_bivariate_imputation(df_topiramate, paste0("log_", m, c(1, 2)), paste0("se_", m, c(1, 2)),
                                      rho_w = -0.3, tau2_val = NULL, rho_b = NULL, m = 1000)
  fits_ov <- fit_imputations_biv(mi_ov)

  for (mm in c(200, 1000)) {
    mi_sub <- mi_ov
    mi_sub$imp_draws <- mi_ov$imp_draws[1:mm, , drop = FALSE]   # first mm imputations
    for (type in types_M) {
      for (delta in deltas_ov) {
        tmp <- adj_bivariate(mi_sub, delta = delta, select_type = type,
                             track.failed.proportion = FALSE, track.ess = FALSE,
                             fits = fits_ov[1:mm])
        res_ov[[length(res_ov) + 1]] <- data.frame(Measure = paste0("log ", m), M = mm,
                                                   Select_type = type, Selection = delta,
                                                   Outcome = tmp$Outcome, Estimate = tmp$Estimate,
                                                   CI_Lower = tmp$CI_Lower, CI_Upper = tmp$CI_Upper)
      }
    }
  }
}

dat_ov <- dplyr::bind_rows(res_ov) %>%
  mutate(M = factor(M, levels = c(200, 1000)),
         Measure = factor(Measure, levels = c("log OR", "log RR")),
         Outcome = factor(ifelse(Outcome == "O1", "50% Seizure Reduction", "Seizure Freedom"),
                          levels = c("50% Seizure Reduction", "Seizure Freedom")),
         Type = factor(ifelse(Select_type == "zscore", "Selection on z-score", "Selection on estimate"),
                       levels = c("Selection on z-score", "Selection on estimate")))
write_csv(dat_ov, "Data/data_M_overlay.csv")

# rows = measure x outcome, columns = selection type
ggplot(dat_ov, aes(x = Selection, y = Estimate, colour = M, fill = M)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(Measure + Outcome ~ Type) +
  scale_colour_manual(values = c("200" = "#D55E00", "1000" = "#0072B2")) +
  scale_fill_manual(values = c("200" = "#D55E00", "1000" = "#0072B2")) +
  labs(x = expression("Selection weight"~delta), y = "Adjusted estimate",
       colour = "Imputations M", fill = "Imputations M") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("Paper/figures/App_M_convergence.png", width = 9, height = 11)

# largest difference between the two curves, for the text
dat_ov |>
  select(M, Measure, Type, Selection, Outcome, Estimate) |>
  tidyr::pivot_wider(names_from = M, values_from = Estimate, names_prefix = "M_") |>
  group_by(Measure, Outcome, Type) |>
  summarise(max_abs_diff = round(max(abs(M_200 - M_1000)), 3), .groups = "drop")



# -------------------------------------------------------------------------
# Justification paragraph
# -------------------------------------------------------------------------
# For computational reasons, the simulation uses M = 200 imputations instead of
# the M = 1000 used in the application. See Appendix to check that this reduction does not
# change the results.


# --- In the appendix add this:
# To check that this reduction does not change the results, the bivariate
# adjustment was applied to the topiramate data with M = 1000 imputations and
# the adjusted estimate recomputed from the first 200 of them, on the log OR
# and on the log RR scale. The within-study correlation was fixed at
# rho_W = -0.3, while the between-study variances and correlation were
# estimated. Both outcomes were considered, 50% seizure reduction (1 study
# unreported) and seizure freedom (6 of 12 unreported), with selection on the
# z-score and on the effect estimate, for delta from 0 to 1.3.
#
# Figure caption: Adjusted estimate (line) and 95% confidence interval (band) as
# a function of the selection weight delta, computed from the first 200 (orange)
# and from all 1000 (blue) imputations of the same analysis. 
#
# For 50% seizure reduction the two curves coincide on both scales (largest
# difference 0.06, against confidence intervals about 0.6 to 0.7 wide). For
# seizure freedom they overlap up to delta = 0.5 and differ by at most 0.12 at
# larger delta, well within confidence intervals 1.3 to 1.8 wide.
# M = 200 therefore leads to the same conclusions as M = 1000.
