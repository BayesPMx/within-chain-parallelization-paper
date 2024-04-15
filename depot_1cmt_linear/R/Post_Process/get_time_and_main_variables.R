rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(furrr)
library(tidyverse)
library(patchwork)


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

# plan(multisession, workers = parallel::detectCores())
# 
# zoom <- expand_grid(threads_per_chain = c(0, 1, 2, 4, 8, 12, 24, 48), 
#                     solver = c(1:4), 
#                     run_number = 1:10) %>%
#   filter(!(threads_per_chain == 0 & solver %in% c(4, 5))) %>% 
#   future_pmap_dfr(.f = get_time_and_main_variables)
# 
# plan(sequential)
# 
# zoom %>% 
#   write_rds("depot_1cmt_linear/Results/time_and_main_variables.rds")

zoom <- read_rds("depot_1cmt_linear/Results/time_and_main_variables.rds")

p_1 <- zoom %>% 
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel")),
         solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>% 
  ggplot() +
  geom_density(aes(x = TVCL, color = parallel), alpha = 0.1) +
  theme_bw() +
  scale_color_manual(name = "Parallel",
                     values = c("Parallel" = "blue",
                                "Not Parallel" = "red")) +
  facet_wrap(~solver, nrow = 1) +
  theme(legend.position = "top")

p_2 <- zoom %>% 
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel")),
         solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>% 
  ggplot() +
  geom_density(aes(x = TVCL, color = solver), alpha = 0.1) +
  theme_bw() +
  scale_color_manual(name = "Solver",
                     values = c("Analytical Threaded" = "blue",
                                "Linear ODE Threaded" = "red",
                                "General ODE Threaded" = "orange",
                                "Torsten Group ODE MPI" = "magenta")) +
  facet_wrap(~parallel, nrow = 1) +
  theme(legend.position = "bottom")
    
p_1 /
  p_2

zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time)) %>% 
  ungroup() %>%
  bind_rows({
    data_1 <- .
    data_1 %>% 
      filter(solver == "rk45", threads_per_chain == 0) %>% 
      mutate(solver = factor("torsten_general"))
  }) %>% 
  group_by(solver) %>% 
  arrange(threads_per_chain, .by_group = TRUE) %>% 
  mutate(across(ends_with("time"), ~ ./.[1]),
         speedup_mean = 1/mean_time,
         speedup_min = 1/min_time,
         speedup_max = 1/max_time)

p_time_vs_processes <- zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time)) %>% 
  ungroup() %>% 
  mutate(type = if_else(solver %in% c("rk45", "torsten_general"), "ode", 
                        solver) %>% 
           factor(),
         processes = if_else(threads_per_chain == 0, "Not Parallel", 
                             as.character(threads_per_chain)) %>% 
           factor(levels = c("Not Parallel", as.character(unique(threads_per_chain)))),
         across(ends_with("time"), ~./60),
         solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>%
  ggplot(aes(x = processes, y = mean_time, 
             group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = min_time, ymax = max_time), width = 0.2) +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     values = c("Analytical Threaded" = "blue",
                                "Linear ODE Threaded" = "red",
                                "General ODE Threaded" = "orange",
                                "Torsten Group ODE MPI" = "magenta")) +
  scale_y_continuous(name = "Elapsed Time (m)",
                     trans = "log10") +
  annotation_logticks(sides = "l") +
  scale_x_discrete(name = "Processes per Chain") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom")
               

p_speedup_vs_processes <- zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time)) %>% 
  ungroup() %>%
  bind_rows({
    data_1 <- .
    data_1 %>% 
      filter(solver == "rk45", threads_per_chain == 0) %>% 
      mutate(solver = factor("torsten_general"))
  }) %>% 
  group_by(solver) %>% 
  arrange(threads_per_chain, .by_group = TRUE) %>% 
  mutate(across(ends_with("time"), ~ ./.[1]),
         speedup_mean = 1/mean_time,
         speedup_min = 1/min_time,
         speedup_max = 1/max_time) %>% 
  ungroup() %>% 
  mutate(type = if_else(solver %in% c("rk45", "torsten_general"), "ode", 
                        solver) %>% 
           factor(),
         processes = if_else(threads_per_chain == 0, "Not Parallel", 
                             as.character(threads_per_chain)) %>% 
           factor(levels = c("Not Parallel", as.character(unique(threads_per_chain)))),
         across(ends_with("time"), ~./60),
         solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>% 
  ggplot(aes(x = processes, y = speedup_mean, group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = speedup_min, ymax = speedup_max), width = 0.2) +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     values = c("Analytical Threaded" = "blue",
                                "Linear ODE Threaded" = "red",
                                "General ODE Threaded" = "orange",
                                "Torsten Group ODE MPI" = "magenta")) +
  scale_y_continuous(name = "Speedup",
                     trans = "identity",
                     breaks = seq(1, 9, by = 2),
                     labels = seq(1, 9, by = 2)) +
  scale_x_discrete(name = "Processes per Chain") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom")

p_time_vs_processes /
  p_speedup_vs_processes +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## Do the grainsize stuff

read_quick_grainsize <- function(solver, threads_per_chain, grainsize){

  solver_string <- case_when(solver == 1 ~ "analytical",
                             solver == 2 ~ "matexp",
                             solver == 3 ~ "rk45",
                             TRUE ~ NA_character_)
  
  file_name_base <- str_c("depot_1cmt_linear/Stan/Fits/Grainsize/Quick_Test/", 
                          solver_string, "_grainsize_", grainsize, "_run_2")
  
  file_name <- if_else(threads_per_chain == 4,
                       str_c(file_name_base, "_threads_4.rds"),
                       str_c(file_name_base, ".rds"))
  
  fit <- read_rds(file_name)

  # The ones fit explicitly with cmdstanr and saved directly as CmdStanMCMC
  # objects (the multi-threaded ones) have fit$time()$total along with the 
  # total time by chain, but the ones fit and saved in CSVs and then converted
  # into CmdStanMCMC objects (Torsten's group solver) have the total time by 
  # chain but fit$time()$total = NA. Therefore, to do the same thing for all 
  # fits, I'll take the maximum time of the 4 chains as the elapsed time.
  return(fit$draws(c("TVCL", "TVVC", "TVKA", "omega", "sigma_p"), 
                   format = "draws_df") %>% 
           mutate(solver = solver_string, 
                  threads_per_chain = threads_per_chain,
                  grainsize = grainsize,
                  time = max(fit$time()$chains$total)))
  
}

plan(multisession, workers = min(parallel::detectCores(), 8))

quick_grainsize <- expand_grid(threads_per_chain = c(4, 24), 
                               solver = c(1:3),
                               grainsize = c(1, 2, 4, 8, 10, 12, 
                                             16, 24, 48, 96)) %>%
  future_pmap_dfr(.f = read_quick_grainsize)

plan(sequential)


summary_of_not_parallel_threaded <- zoom %>% 
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel")),
         solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI")),
         time = time/60) %>% 
  ungroup() %>% 
  filter(parallel == "Not Parallel") %>% 
  group_by(solver, threads_per_chain, parallel, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time))

quick_grainsize %>% 
  distinct(solver, grainsize, threads_per_chain, time) %>% 
  mutate(solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI")),
         processes = factor(threads_per_chain),
         grainsize = factor(grainsize),
         time = time/60) %>% 
  ggplot(aes(x = grainsize, y = time, 
             group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     values = c("Analytical Threaded" = "blue",
                                "Linear ODE Threaded" = "red",
                                "General ODE Threaded" = "orange")) +
  scale_y_continuous(name = "Elapsed Time (m)",
                     trans = "log10") +
  annotation_logticks(sides = "l") +
  scale_x_discrete(name = "Grainsize") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom") +
  geom_hline(data = summary_of_not_parallel_threaded,
             mapping = aes(yintercept = mean_time, color = solver),
             linetype = "dashed") +
  facet_wrap(~ processes, 
             labeller = as_labeller(function(x) str_c("Processes per Chain: ", 
                                                      x)))

quick_grainsize %>% 
  distinct(solver, grainsize, threads_per_chain, time) %>% 
  mutate(solver = case_when(solver == "analytical" ~ "Analytical Threaded",
                            solver == "matexp" ~ "Linear ODE Threaded",
                            solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("Analytical Threaded",
                             "Linear ODE Threaded",
                             "General ODE Threaded",
                             "Torsten Group ODE MPI")),
         processes = factor(threads_per_chain),
         grainsize = factor(grainsize),
         time = time/60) %>% 
  group_by(solver, threads_per_chain) %>% 
  filter(grainsize == 1 | time == min(time)) %>% 
  mutate(speedup = time[1]/time)




# zoom %>% 
#   group_by(solver, threads_per_chain, run_number, .chain) %>% 
#   mutate(.chain = cur_group_id()) %>% 
#   ungroup() %>% 
#   select(-solver, -run_number, -threads_per_chain, -time) %>% 
#   summarize_draws()
