rm(list = ls())
cat("\014")

library(trelliscopejs)
library(mrgsolve)
library(tidybayes)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

TVCL <- 0.5
TVVC <- 12
TVKA <- 1

omega_cl <- 0.3
omega_vc <- 0.3
omega_ka <- 0.3

R <- diag(rep(1, times = 3))
R[1, 2] <- R[2, 1] <- 0.4 # Put in some correlation between CL and VC

sigma_p <- 0.2
sigma_a <- 0

cor_p_a <- 0

n_subjects_per_dose <- 24

dosing_data_1 <- expand.ev(ID = 1:n_subjects_per_dose, addl = 6, ii = 24, 
                           cmt = 1, amt = 50, ss = 0, tinf = 0, 
                           evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  select(ID, TIME, everything()) 

dosing_data_2 <- expand.ev(ID = 1:n_subjects_per_dose, 
                           addl = 13, ii = 24, 
                           cmt = 1, amt = 100, ss = 0, tinf = 0, 
                           evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(ID = (n_subjects_per_dose + 1):(2*n_subjects_per_dose)) %>% 
  select(ID, TIME, everything()) 

dosing_data_3 <- expand.ev(ID = 1:n_subjects_per_dose, 
                           addl = 20, ii = 24, 
                           cmt = 1, amt = 200, ss = 0, tinf = 0, 
                           evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(ID = (2*n_subjects_per_dose + 1):(3*n_subjects_per_dose)) %>% 
  select(ID, TIME, everything()) 

dosing_data_4 <- expand.ev(ID = 1:n_subjects_per_dose, 
                           addl = 27, ii = 24, 
                           cmt = 1, amt = 400, ss = 0, tinf = 0, 
                           evid = 1) %>%
  as_tibble() %>% 
  rename_all(toupper) %>%
  mutate(ID = (3*n_subjects_per_dose + 1):(4*n_subjects_per_dose)) %>% 
  select(ID, TIME, everything()) 

sampling_times <- c(0.25, 0.5, 1, 2, 4, 8, 12, 24)
realistic_times_1 <- c(sampling_times, 72, 144, 144 + sampling_times)
realistic_times_2 <- c(sampling_times, 72, 144, 144 + sampling_times, 
                       312, 312 + sampling_times)
realistic_times_3 <- c(sampling_times, 72, 144, 144 + sampling_times, 
                       312, 312 + sampling_times,
                       480, 480 + sampling_times)
realistic_times_4 <- c(sampling_times, 72, 144, 144 + sampling_times, 
                       312, 312 + sampling_times,
                       480, 480 + sampling_times,
                       648, 648 + sampling_times)

times_to_simulate_1 <- realistic_times_1
times_to_simulate_2 <- realistic_times_2
times_to_simulate_3 <- realistic_times_3
times_to_simulate_4 <- realistic_times_4

nonmem_data_simulate_1 <- dosing_data_1 %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_1))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_1) %>% 
  ungroup() %>%
  bind_rows(dosing_data_1) %>% 
  arrange(ID, TIME, AMT)

nonmem_data_simulate_2 <- dosing_data_2 %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_2))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_2) %>% 
  ungroup() %>%
  bind_rows(dosing_data_2) %>% 
  arrange(ID, TIME, AMT)

nonmem_data_simulate_3 <- dosing_data_3 %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_3))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_3) %>% 
  ungroup() %>%
  bind_rows(dosing_data_3) %>% 
  arrange(ID, TIME, AMT)

nonmem_data_simulate_4 <- dosing_data_4 %>% 
  group_by(ID) %>% 
  slice(rep(1, times = length(times_to_simulate_4))) %>% 
  mutate(AMT = 0,
         ADDL = 0,
         II = 0,
         CMT = 2,
         EVID = 0,
         # RATE = 0,
         TIME = times_to_simulate_4) %>% 
  ungroup() %>%
  bind_rows(dosing_data_4) %>% 
  arrange(ID, TIME, AMT)

nonmem_data_simulate <- mget(str_subset(ls(), 
                                        "^nonmem_data_simulate_[0-9]")) %>%
  bind_rows() %>% 
  rowwise() %>% 
  mutate(drop = if_else(EVID == 1, 0L, rbinom(1, 1, 0.1))) %>% 
  ungroup() %>% 
  filter(drop == 0)

n_subjects <- nonmem_data_simulate %>%  # number of individuals to simulate
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_total <- nrow(nonmem_data_simulate) # total number of time points at which to predict

subj_start <- nonmem_data_simulate %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total) 

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  time = nonmem_data_simulate$TIME,
                  amt = nonmem_data_simulate$AMT,
                  cmt = nonmem_data_simulate$CMT,
                  evid = nonmem_data_simulate$EVID,
                  rate = nonmem_data_simulate$RATE,
                  ii = nonmem_data_simulate$II,
                  addl = nonmem_data_simulate$ADDL,
                  ss = nonmem_data_simulate$SS,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  TVCL = TVCL,
                  TVVC = TVVC,
                  TVKA = TVKA,
                  omega_cl = omega_cl,
                  omega_vc = omega_vc,
                  omega_ka = omega_ka,
                  R = R,
                  sigma_p = sigma_p,
                  sigma_a = sigma_a,
                  cor_p_a = cor_p_a,
                  solver = 2) # analytical = 1, mat exp = 2, rk45 = 3, bdf = 4

model <- cmdstan_model("depot_1cmt_linear_ragged/Stan/Simulate/depot_1cmt_ppa.stan") 

simulated_data <- model$sample(data = stan_data,
                               fixed_param = TRUE,
                               seed = 11235,
                               iter_warmup = 0,
                               iter_sampling = 1,
                               chains = 1,
                               parallel_chains = 1)

params_ind <- simulated_data$draws(c("CL", "VC", "KA")) %>% 
  spread_draws(c(CL, VC, KA)[ID]) %>%
  ungroup() %>%
  select(ID, CL, VC, KA)

data <- simulated_data$draws(c("dv", "ipred")) %>% 
  spread_draws(c(dv, ipred)[i]) %>% 
  ungroup() %>% 
  inner_join(nonmem_data_simulate %>% 
               mutate(i = 1:n()),
             by = "i") %>%
  select(ID, AMT, II, ADDL, RATE, CMT, EVID, SS, TIME, 
         DV = "dv", IPRED = "ipred") %>% 
  mutate(LLOQ = 1, 
         BLOQ = case_when(EVID == 1 ~ NA_real_,
                          DV <= LLOQ ~ 1,
                          DV > LLOQ ~ 0,
                          TRUE ~ NA_real_),
         DV = if_else(EVID == 1 | BLOQ == 1, NA_real_, DV),
         MDV = if_else(is.na(DV), 1, 0)) %>% 
  relocate(DV, .after = last_col()) %>% 
  relocate(TIME, .before = DV)

(p_1 <- ggplot(data %>% 
                 group_by(ID) %>% 
                 mutate(Dose = factor(max(AMT))) %>% 
                 ungroup() %>% 
                 filter(!is.na(DV))) +
    geom_point(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = TIME, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(data$TIME), by = 24*7),
                       labels = seq(0, max(data$TIME)/(24*7), by = 1),
                       limits = c(0, max(data$TIME))) +
    labs(color = "Dose (mg)") +
    theme(legend.position = "bottom"))

data %>%
  select(-IPRED) %>% 
  # write_csv("depot_1cmt_linear_ragged/Data/depot_1cmt_ppa.csv", na = ".")
  write_csv("depot_1cmt_linear_ragged/Data/depot_1cmt_prop.csv", na = ".")

params_ind %>%
  # write_csv("depot_1cmt_linear_ragged/Data/depot_1cmt_ppa_params_ind.csv")
  write_csv("depot_1cmt_linear_ragged/Data/depot_1cmt_prop_params_ind.csv")


