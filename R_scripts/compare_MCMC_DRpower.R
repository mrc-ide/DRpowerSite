# compare_MCMC_DRpower.R
#
# Author: Bob Verity
# Date created: 2026-05-07
#
# Description:
# This script compares posterior estimates of global deletion prevalence
# obtained using:
#   1. The analytical DRpower approach
#   2. Explicit Bayesian MCMC sampling of site-level prevalences
#
# Both approaches are based on the same hierarchical beta-binomial
# model structure. The key distinction is computational:
# - DRpower analytically integrates over the latent site-level
#   prevalences
# - The MCMC approach explicitly samples the site-level prevalences
#   using `drjacoby`
#
# The purpose of this script is to demonstrate that both methods
# produce equivalent posterior distributions for the global prevalence
# parameter (`p_global`).
#
# Input files:
#   outputs/mcmc.rds
#   data_raw/Somaliland_2026.05.07.csv
#
# Output:
#   outputs/posterior_compare.pdf
#
# ------------------------------------------------------------------

library(here)
library(tidyverse)
#remotes::install_github("mrc-ide/drjacoby@v1.5.4")
library(drjacoby)
#remotes::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)

# ------------------------------------------------------------------

# import MCMC results
mcmc <- readRDS(here("outputs", "mcmc.rds"))

# get posterior density from MCMC
mcmc_density <- mcmc$output |>
  filter(phase == "sampling") |>
  pull(p_global) |>
  density(from = 0, to = 0.3, n = 1000)

# import data
dat <- read.csv(here("data_raw", "Somaliland_2026.05.07.csv"))

# run standard DRpower method and extract full posterior
drp <- DRpower::get_prevalence(n = dat$n_deletions,
                               N = dat$n_tested,
                               post_full_on = TRUE,
                               post_full_breaks = mcmc_density$x)


# plot both distributions
data.frame(p = mcmc_density$x * 100,
           MCMC = mcmc_density$y,
           DRpower = drp$post_full[[1]]) |>
  pivot_longer(cols = 2:3, names_to = "Method") |>
  ggplot() + theme_bw() +
  geom_line(aes(x = p, y = value, col = Method)) +
  xlab("Global prevalence (%)") + ylab("Posterior density") +
  ggtitle("Comparison of DRpower vs. MCMC posterior distributions")

# write image to file
ggsave(filename = here("outputs", "posterior_compare.pdf"))
