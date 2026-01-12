# ==============================================================================
# Robust Power-Moment Scale Estimation Simulation (R Version)
# ==============================================================================
# This script performs a Monte Carlo simulation to evaluate the performance
# of the power-moment based robust scale estimator T_n(alpha).
#
# It generates:
# 1. /robust_power_moment_simulation/tables/ : CSV summaries of MSE, Bias, Variance.
# 2. /robust_power_moment_simulation/figs/   : Plots of metrics vs Alpha.
# ==============================================================================

# -----------------------------
# 1. Setup and Libraries
# -----------------------------
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Define Output Directories
base_dir <- "robust_power_moment_simulation"
fig_dir <- file.path(base_dir, "figs")
tab_dir <- file.path(base_dir, "tables")

# Create directories (suppress warnings if they exist)
dir.create(base_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

cat("Directories created at:", normalizePath(base_dir, mustWork = FALSE), "\n")

# -----------------------------
# 2. Core Estimator Functions
# -----------------------------

# Compute theoretical moment m(alpha) = E(|U|^(2*alpha)) for U ~ N(0,1)
# Closed form: 2^alpha * Gamma(alpha + 0.5) / sqrt(pi)
get_m_alpha <- function(alpha) {
  # The moment is defined for 2*alpha > -1 => alpha > -0.5
  if (alpha <= -0.5) {
    return(NA_real_)
  }
  return((2^alpha) * gamma(alpha + 0.5) / sqrt(pi))
}

# Compute the estimator T_n(alpha) for a vector x
# T_n(alpha) = ( mean(|x - mean(x)|^(2*alpha)) / m(alpha) )^(1/alpha)
calc_Tn <- function(x, alpha, m_val) {
  # Input validation
  if (is.na(m_val) || m_val <= 0) return(NA_real_)
  if (abs(alpha) < 1e-9) return(NA_real_) # Undefined at alpha = 0
  
  z <- x - mean(x)
  absz <- abs(z)
  
  # If alpha is negative, any zero deviation causes infinite moment
  if (alpha < 0 && any(absz < 1e-12)) {
    return(NA_real_)
  }
  
  # Compute S_n(alpha)
  # Using absz^(2*alpha)
  terms <- absz^(2 * alpha)
  Sn <- mean(terms)
  
  if (!is.finite(Sn) || Sn <= 0) return(NA_real_)
  
  # Final estimator
  # Result can be NaN if base is negative (shouldn't be) or overflow
  Tn <- (Sn / m_val)^(1/alpha)
  return(Tn)
}

# -----------------------------
# 3. Data Samplers
# -----------------------------

# Scenario 1: Normal Mixture
# F ~ N(10, 2), G ~ N(20, 1)
sample_mixture_nn <- function(n, eps) {
  # Generate mixture mask
  n_contam <- rbinom(1, n, eps)
  n_main <- n - n_contam
  
  x_main <- rnorm(n_main, mean = 10, sd = 2)
  x_contam <- rnorm(n_contam, mean = 20, sd = 1)
  
  sample(c(x_main, x_contam)) # shuffle
}

# Scenario 2: Lognormal Mixture
# F is LogNormal matching mean=10, sd=2
# G is LogNormal shifted (heavier)
sample_mixture_lnl <- function(n, eps) {
  # Calculate params for F
  muF_target <- 10
  sigmaF_target <- 2
  
  varF <- sigmaF_target^2
  s2_F <- log(1 + varF / (muF_target^2))
  s_F <- sqrt(s2_F)
  m_F <- log(muF_target) - s2_F / 2
  
  # Params for G (Contaminant)
  # mG = mF + 0.7, sG = sF * 1.25
  m_G <- m_F + 0.7
  s_G <- s_F * 1.25
  
  n_contam <- rbinom(1, n, eps)
  n_main <- n - n_contam
  
  x_main <- rlnorm(n_main, meanlog = m_F, sdlog = s_F)
  x_contam <- rlnorm(n_contam, meanlog = m_G, sdlog = s_G)
  
  sample(c(x_main, x_contam))
}

# Scenario 3: Normal + Uniform
# F ~ N(10, 2), G ~ Unif(25, 35)
sample_mixture_nu <- function(n, eps) {
  n_contam <- rbinom(1, n, eps)
  n_main <- n - n_contam
  
  x_main <- rnorm(n_main, mean = 10, sd = 2)
  x_contam <- runif(n_contam, min = 25, max = 35)
  
  sample(c(x_main, x_contam))
}

# -----------------------------
# 4. Simulation Driver
# -----------------------------

run_simulation <- function(scenario_name, sampler_func, true_sigma2, 
                           eps_list, alpha_grid, n = 200, R = 2000, seed = 123) {
  
  set.seed(seed)
  cat(sprintf("Running Scenario: %s ...\n", scenario_name))
  
  # Precompute m(alpha)
  m_map <- setNames(sapply(alpha_grid, get_m_alpha), alpha_grid)
  
  results_list <- list()
  counter <- 1
  
  for (eps in eps_list) {
    # We will run R replications. 
    # To optimize R speed, we can generate a matrix (n x R) 
    # but calculating column-wise means/stats is often cleaner in a loop 
    # given the logic complexity of Tn.
    
    # Store estimates for all alphas across all R runs
    # Structure: list of vectors, key = alpha
    T_storage <- vector("list", length(alpha_grid))
    names(T_storage) <- as.character(alpha_grid)
    for(a_str in names(T_storage)) T_storage[[a_str]] <- numeric(R)
    
    for (r in 1:R) {
      x <- sampler_func(n, eps)
      
      for (i in seq_along(alpha_grid)) {
        a <- alpha_grid[i]
        a_str <- as.character(a)
        val <- calc_Tn(x, a, m_map[[a_str]])
        T_storage[[a_str]][r] <- val
      }
    }
    
    # Summarize results for this epsilon
    for (i in seq_along(alpha_grid)) {
      a <- alpha_grid[i]
      a_str <- as.character(a)
      
      estimates <- T_storage[[a_str]]
      
      # Filter finite values
      valid_est <- estimates[is.finite(estimates)]
      n_valid <- length(valid_est)
      valid_rate <- n_valid / R
      
      if (n_valid > 0) {
        mean_T <- mean(valid_est)
        var_T <- var(valid_est) # Uses denominator n-1
        bias <- mean_T - true_sigma2
        mse <- var_T + bias^2 # Bias-variance decomposition
      } else {
        mean_T <- NA
        var_T <- NA
        bias <- NA
        mse <- NA
      }
      
      results_list[[counter]] <- data.frame(
        scenario = scenario_name,
        n = n,
        R = R,
        eps = eps,
        alpha = a,
        m_alpha = m_map[[a_str]],
        valid_rate = valid_rate,
        mean_T = mean_T,
        bias = bias,
        var = var_T,
        mse = mse
      )
      counter <- counter + 1
    }
  }
  
  return(do.call(rbind, results_list))
}

# -----------------------------
# 5. Execute Simulation
# -----------------------------

# Parameters
n <- 200
R <- 2000
eps_list <- c(0.1, 0.2, 0.3)
alpha_grid <- seq(-0.9, 2.0, by = 0.1)
alpha_grid <- round(alpha_grid, 2) # avoid floating point naming issues

# Target variances
# Main component is always roughly N(10, 2^2) or Lognorm with sd=2
true_sigma2 <- 4.0 

# Run Scenarios
df_nn <- run_simulation(
  "Normal(main) + Normal(contam)",
  sample_mixture_nn,
  true_sigma2, eps_list, alpha_grid, n, R, seed = 2026
)

df_lnl <- run_simulation(
  "Lognormal(main) + Lognormal(contam)",
  sample_mixture_lnl,
  true_sigma2, eps_list, alpha_grid, n, R, seed = 2027
)

df_nu <- run_simulation(
  "Normal(main) + Uniform(contam)",
  sample_mixture_nu,
  true_sigma2, eps_list, alpha_grid, n, R, seed = 2028
)

# Combine Results
df_all <- rbind(df_nn, df_lnl, df_nu)

# -----------------------------
# 6. Save Outputs
# -----------------------------

# 6.1 Save Full Table
write_csv(df_all, file.path(tab_dir, "simulation_summary_all.csv"))

# 6.2 Save Best Alpha Table (Min MSE for positive alpha)
df_best <- df_all %>%
  filter(alpha > 0) %>%
  group_by(scenario, eps) %>%
  slice(which.min(mse)) %>%
  select(scenario, eps, alpha, mse, bias, var, valid_rate)

write_csv(df_best, file.path(tab_dir, "best_alpha_by_mse_positive_alpha.csv"))

# -----------------------------
# 7. Plotting
# -----------------------------

plot_metric <- function(df, scen_name, metric, y_lab) {
  # Filter data
  data_plot <- df %>% filter(scenario == scen_name)
  
  # File name sanitization
  safe_name <- gsub("[^A-Za-z0-9]", "_", scen_name)
  safe_name <- gsub("__+", "_", safe_name)
  filename <- paste0(metric, "_", safe_name, ".png")
  
  p <- ggplot(data_plot, aes(x = alpha, y = .data[[metric]], color = as.factor(eps), group = as.factor(eps))) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    labs(
      title = paste0(scen_name, ": ", y_lab, " vs Alpha"),
      x = "Alpha",
      y = y_lab,
      color = "Epsilon"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(fig_dir, filename), plot = p, width = 7, height = 5, bg = "white")
}

unique_scenarios <- unique(df_all$scenario)

cat("Generating plots...\n")
for (sc in unique_scenarios) {
  plot_metric(df_all, sc, "mse", "MSE")
  plot_metric(df_all, sc, "bias", "Bias")
  plot_metric(df_all, sc, "var", "Variance")
  plot_metric(df_all, sc, "valid_rate", "Valid Rate")
}

# -----------------------------
# 8. Create README
# -----------------------------
readme_text <- "Robust Power-Moment Scale Estimation Simulation (R Output)
=========================================================

Description
-----------
This folder contains results from an R simulation of the estimator:
T_n(alpha) = ( S_n(alpha) / m(alpha) )^(1/alpha)

Target: sigma^2 = 4 (Standard Deviation = 2) for the main component.

Scenarios
---------
1. Normal + Normal: F=N(10,2), G=N(20,1)
2. Lognormal + Lognormal: F params match mean=10,sd=2. G is heavier.
3. Normal + Uniform: F=N(10,2), G=U(25,35)

Files
-----
tables/simulation_summary_all.csv : Complete simulation results.
tables/best_alpha_by_mse...csv    : Optimal alpha minimizing MSE (alpha > 0).
figs/*.png                        : Visualization of performance metrics.

Notes
-----
- alpha <= -0.5 results in undefined theoretical moments for Normal baseline.
- alpha = 0 is undefined.
- Results generated using R with seed fixed for reproducibility.
"

write_file(readme_text, file.path(base_dir, "README.txt"))

cat("Simulation complete. All files saved in:", normalizePath(base_dir, mustWork = FALSE), "\n")