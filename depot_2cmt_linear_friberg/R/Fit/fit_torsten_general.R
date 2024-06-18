rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

# set_cmdstan_path("~/Torsten/cmdstan")
set_cmdstan_path("/data/Random/mpi_example1/Torsten/v0.91.0/cmdstan/")

nonmem_data <- read_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_prop_friberg_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

n_subjects <- nonmem_data %>%  # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data)   # total number of records

i_obs <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  filter(evid == 0) %>%
  select(row_num) %>%
  deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(row_num) %>%
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  rate = nonmem_data$rate,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  location_tvcl = 0.75,
                  location_tvvc = 18,
                  location_tvq = 3,
                  location_tvvp = 35,
                  location_tvka = 0.8,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  location_tvmtt = 100,
                  location_tvcirc0 = 5,
                  location_tvgamma = 0.2,
                  location_tvalpha = 3e-4,
                  scale_tvmtt = 1,
                  scale_tvcirc0 = 1,
                  scale_tvgamma = 1,
                  scale_tvalpha = 1,
                  scale_omega_mtt = 0.4,
                  scale_omega_circ0 = 0.4,
                  scale_omega_gamma = 0.4,
                  scale_omega_alpha = 0.4,
                  lkj_df_omega_pd = 2,
                  scale_sigma_p_pd = 0.5,
                  prior_only = 0,
                  no_gq_predictions = 0)

write_stan_json(stan_data, file = "depot_2cmt_linear_friberg/Data/stan_data.json")

system("cd /data/Random/within-chain-parallelization-paper/Torsten/v0.91.0/cmdstan && make /data/Random/within-chain-parallelization-paper/depot_2cmt_linear_friberg/Stan/Fit/depot_2cmt_prop_friberg_prop_torsten_general")
system("cd /data/Random/within-chain-parallelization-paper/depot_2cmt_linear_friberg/Stan/Fit && ./fit_torsten_general.sh")
