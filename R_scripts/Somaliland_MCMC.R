# Somaliland_MCMC.R
#
# Author: Bob Verity
# Date created: 2026-05-07
#
# Description:
# This script estimates deletion prevalence at individual sampling sites using a
# hierarchical beta-binomial model implemented in drjacoby. Site-specific
# prevalences are partially pooled around a shared global prevalence
# distribution. This mirrors exactly the model used in DRpower, the only
# difference being that here we explicitly model the site-level prevalences,
# whereas in DRpower we integrate over them analytically.
#
# Model structure:
#   n_deletions_i ~ Binomial(n_tested_i, p_site_i)
#   p_site_i      ~ Beta(alpha, beta)
#
# where:
#   alpha = p_global * (1/r - 1)
#   beta  = (1 - p_global) * (1/r - 1)
#
# Parameters:
#   p_global   Global mean prevalence
#   r          Overdispersion / heterogeneity parameter
#   p_site_i   Site-specific prevalence parameters
#
# Input files:
#   data_raw/Somaliland_2026.05.07.csv
#
# Output files:
#   outputs/mcmc.rds

library(here)
library(tidyverse)
#remotes::install_github("mrc-ide/drjacoby@v1.5.4")
library(drjacoby)

set.seed(1)

# ------------------------------------------------------------------

# import data
dat <- read.csv(here("data_raw", "Somaliland_2026.05.07.csv")) |>
  mutate(prev = round(n_deletions / n_tested * 100, digits = 2))

dat

# define parameters dataframe
n_sites <- nrow(dat)
df_params <- data.frame(name = c("p_global", "r", sprintf("p_site_%s", 1:n_sites)),
                        min = 0, max = 1)


# define log-likelihood function
r_loglike <- function(params, data, misc) {
  
  n_sites <- misc$n_sites
  p_site <- params[sprintf("p_site_%s", 1:n_sites)]
  alpha <- params["p_global"]*(1 / params["r"] - 1)
  beta <- (1 - params["p_global"])*(1 / params["r"] - 1)
  
  # beta-binomial model
  sum(dbinom(x = data$n_deletions,
         size = data$n_tested,
         prob = p_site,
         log = TRUE)) +
    sum(dbeta(p_site, shape1 = alpha, shape2 = beta, log = TRUE))
}

# define log-prior function
r_logprior <- function(params, misc) {
  dbeta(params["r"], shape1 = 1, shape2 = 9, log = TRUE)
}

# run MCMC
mcmc <- run_mcmc(data = dat,
                 df_params = df_params,
                 misc = list(n_sites = n_sites),
                 loglike = r_loglike,
                 logprior = r_logprior,
                 burnin = 500,
                 samples = 1e4,
                 chains = 10,
                 pb_markdown = FALSE)

# MCMC diagnostics
mcmc$diagnostics$rhat
mcmc$diagnostics$ess

# exploratory plots
plot_trace(mcmc, phase = "sampling")

plot_density(mcmc)

# save MCMC to file
saveRDS(mcmc, file = here("outputs", "mcmc.rds"))

