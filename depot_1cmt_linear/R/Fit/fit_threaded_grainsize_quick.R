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

sample_and_save_all_quick <- function(grainsize, solver, run_number){
  
  print(str_c("solver = ", solver, ", ", 
              "grainsize = ", grainsize, ", ",  
              "run # ", run_number))
  
  init_files <- str_c("depot_1cmt_linear/Data/Inits/inits_", 
                      run_number, "_", 1:4, ".json")
  
  stan_data$solver <- solver
  stan_data$grainsize <- grainsize
  
  model <- cmdstan_model(
    "depot_1cmt_linear_grainsize/Stan/Fit/depot_1cmt_prop_all_solvers.stan",
    cpp_options = list(stan_threads = TRUE))
  
  fit <- model$sample(data = stan_data,
                      seed = 112356,
                      chains = 4,
                      parallel_chains = 4,
                      threads_per_chain = 4,
                      iter_warmup = 500,
                      iter_sampling = 1000,
                      adapt_delta = 0.8,
                      refresh = 500,
                      max_treedepth = 10,
                      init = init_files)
  
  solver_string <- case_when(solver == 1 ~ "analytical",
                             solver == 2 ~ "matexp",
                             solver == 3 ~ "rk45",
                             TRUE ~ NA_character_)
  
  fit$save_object(str_c("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/",
                        solver_string, "_grainsize_", grainsize, "_run_", 
                        run_number, "_threads_4.rds"))
  
}

expand_grid(solver = c(1, 2, 3),
            grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96), 
            run_number = 2) %>% 
  pwalk(.f = sample_and_save_all_quick)

fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_1_run_2_threads_4.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_2_run_2_threads_4.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_4_run_2_threads_4.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_8_run_2_threads_4.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_10_run_2_threads_4.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_12_run_2_threads_4.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_16_run_2_threads_4.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_24_run_2_threads_4.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_48_run_2_threads_4.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_96_run_2_threads_4.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))

fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_1_run_2_threads_4.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_2_run_2_threads_4.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_4_run_2_threads_4.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_8_run_2_threads_4.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_10_run_2_threads_4.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_12_run_2_threads_4.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_16_run_2_threads_4.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_24_run_2_threads_4.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_48_run_2_threads_4.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_96_run_2_threads_4.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))

fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_1_run_2_threads_4.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_2_run_2_threads_4.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_4_run_2_threads_4.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_8_run_2_threads_4.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_10_run_2_threads_4.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_12_run_2_threads_4.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_16_run_2_threads_4.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_24_run_2_threads_4.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_48_run_2_threads_4.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_96_run_2_threads_4.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))






fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_1_run_2.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_2_run_2.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_4_run_2.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_8_run_2.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_10_run_2.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_12_run_2.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_16_run_2.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_24_run_2.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_48_run_2.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/analytical_grainsize_96_run_2.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))

fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_1_run_2.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_2_run_2.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_4_run_2.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_8_run_2.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_10_run_2.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_12_run_2.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_16_run_2.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_24_run_2.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_48_run_2.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/matexp_grainsize_96_run_2.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))

fit_1 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_1_run_2.rds")
fit_2 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_2_run_2.rds")
fit_4 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_4_run_2.rds")
fit_8 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_8_run_2.rds")
fit_10 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_10_run_2.rds")
fit_12 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_12_run_2.rds")
fit_16 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_16_run_2.rds")
fit_24 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_24_run_2.rds")
fit_48 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_48_run_2.rds")
fit_96 <- read_rds("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/rk45_grainsize_96_run_2.rds")

tibble(grainsize = c(1, 2, 4, 8, 10, 12, 16, 24, 48, 96),
       time = c(fit_1$time()$total,
                fit_2$time()$total,
                fit_4$time()$total,
                fit_8$time()$total,
                fit_10$time()$total,
                fit_12$time()$total,
                fit_16$time()$total,
                fit_24$time()$total,
                fit_48$time()$total,
                fit_96$time()$total))
