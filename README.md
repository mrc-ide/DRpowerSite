# DRpowerSite

Bayesian estimation of site-level deletion prevalence using hierarchical beta-binomial models.

### Overview

DRpowerSite extends the methodology implemented in [DRpower](https://mrc-ide.github.io/DRpower/index.html) to estimate prevalence at individual sampling sites.

The original DRpower framework models site-level prevalences as latent random effects drawn from a shared beta distribution, but integrates over these analytically in order to efficiently estimate the global mean prevalence and between-site heterogeneity.

In contrast, DRpowerSite explicitly estimates the site-level prevalences using Bayesian MCMC implemented in drjacoby.

This allows:

- estimation of posterior prevalence distributions for each site
- uncertainty quantification at the site level
- shrinkage / partial pooling across sites
- borrowing of information for sparsely sampled locations

while retaining exactly the same underlying hierarchical model structure used in DRpower.
