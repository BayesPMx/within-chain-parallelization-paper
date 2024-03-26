rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_1cmt_linear/Data/depot_1cmt_prop.csv",
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
                  location_tvcl = 0.6,
                  location_tvvc = 15,
                  location_tvka = 1,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0)

write_stan_json(stan_data, file = "depot_1cmt_linear/Data/stan_data.json")

system("mpiexec -n 2 ./depot_1cmt_prop_torsten_general sample num_warmup=300 num_samples=200 algorithm=hmc engine=nuts \
       data file=../../Data/stan_data.json output file=torsten_general.csv \
       init=../../Data/Inits/inits_1_1.json")

#################################

mpi_cmd <- "/data/home/cadavis/mpich-install/bin/mpiexec"

n_jobs <- 24
run_number <- 1

init_files <- str_c("depot_1cmt_linear/Data/Inits/inits_", run_number, "_",
                    1:4, ".json")

model <- cmdstan_model(
  "depot_1cmt_linear/Stan/Fit/depot_1cmt_prop_torsten_general.stan",
  cpp_options = list(STAN_MPI = TRUE, 
                     CXX = "mpicxx", 
                     TBB_CXX_TYPE = "gcc"))

model <- cmdstan_model(
  "depot_1cmt_linear/Stan/Fit/depot_1cmt_prop_torsten_general.stan",
  cpp_options = list(STAN_MPI = TRUE, 
                     CXX = "mpicxx", 
                     TBB_CXX_TYPE = "gcc"))

sample_and_save_torsten_general <- function(n_jobs, run_number){
  
  fit <- model$sample_mpi(data = stan_data,
                          mpi_cmd = mpi_cmd,
                          mpi_args = list("bind-to" = "core", 
                                          "n" = n_jobs),
                          seed = 112356,
                          chains = 1,
                          iter_warmup = 500,
                          iter_sampling = 500,
                          adapt_delta = 0.8,
                          refresh = 50,
                          max_treedepth = 10,
                          output_dir = "depot_1cmt_linear/Stan/Fits",
                          output_basename = str_c("torsten_general_test_", n_jobs),
                          save_warmup = TRUE,
                          init = init_files[1])
  
  fit$save_object(str_c("depot_1cmt_linear/Stan/Fits/torstengeneral_",
                        n_jobs, "_jobs_run_", run_number, ".rds"))
  
}

expand_grid(n_jobs = c(24, 8, 4),
            run_number = 3) %>% 
  pwalk(.f = sample_and_save_torsten_general)

