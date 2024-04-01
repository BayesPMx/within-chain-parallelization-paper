rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(furrr)
library(tidyverse)


convert_csvs_to_cmdstanmcmc <- function(threads_per_chain, solver, run_number){
  
  solver_string <- case_when(solver == 1 ~ "analytical",
                             solver == 2 ~ "matexp",
                             solver == 3 ~ "rk45",
                             solver == 4 ~ "torsten_general",
                             solver == 5 ~ "torsten_general_mpi",
                             TRUE ~ NA_character_)
  
  files <- str_c("depot_1cmt_linear/Stan/Fits/Torsten_General/", 
                 solver_string, "_", threads_per_chain, "_jobs_run_", 
                 run_number, "_", 1:4, ".csv")
  
  fit <- as_cmdstan_fit(files = files)
  
  # fit$save_object(str_c("depot_1cmt_linear/Stan/Fits/", solver_string, "_",
  #                       threads_per_chain, "_jobs_run_", run_number, ".rds"))
  
  return(fit)
  
}
  
get_time_and_main_variables <- function(threads_per_chain, solver, run_number){
  
  solver_string <- case_when(solver == 1 ~ "analytical",
                             solver == 2 ~ "matexp",
                             solver == 3 ~ "rk45",
                             solver == 4 ~ "torsten_general",
                             solver == 5 ~ "torsten_general_mpi",
                             TRUE ~ NA_character_)
  
  if(solver <= 3){
    
    file_name <- str_c("depot_1cmt_linear/Stan/Fits/", solver_string, "_",
                       threads_per_chain, "_threads_run_", run_number, ".rds")
    
    fit <- read_rds(file_name)
    
  }else{
    
    fit <- convert_csvs_to_cmdstanmcmc(threads_per_chain, solver, run_number)
    
  }
  
  # The ones fit explicitly with cmdstanr and saved directly as CmdStanMCMC
  # objects (the multi-threaded ones) have fit$time()$total along with the 
  # total time by chain, but the ones fit and saved in CSVs and then converted
  # into CmdStanMCMC objects (Torsten's group solver) have the total time by 
  # chain but fit$time()$total = NA. Therefore, to do the same thing for all 
  # fits, I'll take the maximum time of the 4 chains as the elapsed time.
  return(fit$draws(c("TVCL", "TVVC", "TVKA", "omega", "sigma_p"), 
                   format = "draws_df") %>% 
           mutate(solver = solver_string, 
                  run_number = run_number, 
                  threads_per_chain = threads_per_chain,
                  time = max(fit$time()$chains$total)))
  
}

plan(multisession, workers = parallel::detectCores())

zoom <- expand_grid(threads_per_chain = c(0, 1, 2, 4, 8, 12, 24, 48), 
                    solver = c(1:4), # c(1:5) if I get the cluster going
                    run_number = 1:10) %>%
  filter(!(threads_per_chain == 0 & solver %in% c(4, 5))) %>% 
  future_pmap_dfr(.f = get_time_and_main_variables)

plan(sequential)

zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time)) %>% 
  mutate(across(ends_with("time"), ~ ./.[1]))

