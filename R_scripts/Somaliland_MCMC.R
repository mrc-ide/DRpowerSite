
library(here)
library(tidyverse)


dat <- read.csv(here("data_raw", "Somaliland_2026.05.07.csv")) |>
  mutate(prev = round(n_deletions / n_tested * 100, digits = 2))

dat
