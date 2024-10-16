rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(furrr)
library(tidyverse)
library(patchwork)

## Do the grainsize stuff

read_grainsize <- function(solver, threads_per_chain, grainsize){
  
  solver_string <- case_when(solver == 1 ~ "analytical",
                             solver == 2 ~ "matexp",
                             solver == 3 ~ "rk45",
                             TRUE ~ NA_character_)
  
  file_name <- str_c("depot_1cmt_linear_ragged/Stan/Fits/Grainsize/", 
                     solver_string, "_threads_", threads_per_chain,
                     "_grainsize_", grainsize, "_run_2.rds")
  
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

quick_grainsize_2 <- tibble(threads_per_chain = 0, grainsize = 0, 
                          solver = 1:3) %>% 
  bind_rows(expand_grid(threads_per_chain = c(4, 24), 
                        solver = c(1:3),
                        grainsize = c(1, 2, 4, 8, 10, 12, 
                                      16, 24, 48, 96))) %>%
  future_pmap_dfr(.f = read_grainsize)

plan(sequential)


data_for_plot_2 <- quick_grainsize_2 %>% 
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
         time = time/60)

p_grainsize_dataset_2 <- data_for_plot_2 %>%
  filter(threads_per_chain != 0) %>% 
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
        legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  geom_hline(data = data_for_plot_2 %>%
               filter(threads_per_chain == 0) %>% 
               select(solver, time),
             mapping = aes(yintercept = time, color = solver),
             linetype = "dashed") +
  facet_wrap(~ processes, 
             labeller = as_labeller(function(x) str_c("Processes per Chain: ", 
                                                      x)))

data_for_plot_2 %>% 
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel"))) %>%
  group_by(solver) %>% 
  mutate(baseline_no_parallel = time[parallel == "Not Parallel"]) %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(baseline_grainsize_1 = case_when(all(grainsize == 0) ~ NA_real_,
                                          TRUE ~ time[1])) %>% 
  summarize(grainsize_1_vs_no_parallel = baseline_no_parallel/baseline_grainsize_1) %>% 
  ungroup() %>% 
  distinct() %>% 
  drop_na() %>% 
  pivot_wider(names_from = threads_per_chain, 
              values_from = grainsize_1_vs_no_parallel)


data_for_plot_2 %>% 
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel"))) %>%
  group_by(solver) %>% 
  mutate(baseline_no_parallel = time[parallel == "Not Parallel"]) %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(baseline_grainsize_1 = case_when(all(grainsize == 0) ~ NA_real_,
                                          TRUE ~ time[1]),
         optimum_grainsize_time = min(time),
         optimum_grainsize = grainsize[time == optimum_grainsize_time]) %>% 
  ungroup() %>% 
  select(solver, threads_per_chain, starts_with("baseline"),
         starts_with("optimum")) %>% 
  distinct() %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(speedup_vs_grainsize_1 = baseline_grainsize_1/optimum_grainsize_time,
         speedup_vs_no_threading = baseline_no_parallel/optimum_grainsize_time,
         speedup_grainsize_1_vs_no_threading = baseline_no_parallel/baseline_grainsize_1) %>% 
  filter(threads_per_chain > 0) %>% 
  arrange(solver) %>% 
  select(solver, optimum_grainsize, starts_with("speedup"))




blah <- data_for_plot_2 %>%  
  mutate(parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel"))) %>% 
  # ungroup() %>% 
  group_by(solver, threads_per_chain, parallel) %>% 
  distinct(time) %>% 
  # ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  summarize(mean_time = mean(time), min_time = min(time), 
            max_time = max(time)) %>% 
  filter(str_detect(solver, "Threaded")) %>% 
  mutate(grainsize = if_else(parallel == "Parallel", 1, 0)) %>% 
  select(solver, threads_per_chain, parallel, grainsize, time = mean_time) %>% 
  filter(threads_per_chain %in% c(0, 4, 24)) %>% 
  bind_rows(quick_grainsize %>% 
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
                     parallel = "Parallel",
                     time = time/60) %>% 
              group_by(solver, threads_per_chain) %>% 
              filter(time == min(time)) %>% 
              ungroup()) %>% 
  arrange(solver, threads_per_chain, grainsize)




t_1 <- blah %>% 
  group_by(solver) %>% 
  mutate(baseline_no_parallel = time[parallel == "Not Parallel"]) %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(baseline_grainsize_1 = case_when(all(grainsize == 0) ~ NA_real_,
                                          TRUE ~ time[1])) %>% 
  summarize(grainsize_1_vs_no_parallel = baseline_no_parallel/baseline_grainsize_1) %>% 
  ungroup() %>% 
  distinct() %>% 
  drop_na() %>% 
  pivot_wider(names_from = threads_per_chain, 
              values_from = grainsize_1_vs_no_parallel)

t_2 <- blah %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(baseline_grainsize_1 = case_when(all(grainsize == 0) ~ NA_real_,
                                          TRUE ~ time[1]),
         optimum_grainsize = case_when(all(grainsize == 0) ~ NA_real_,
                                       TRUE ~ time[2])) %>% 
  summarize(optimum_grainsize_vs_grainsize_1 = baseline_grainsize_1/optimum_grainsize) %>% 
  ungroup() %>% 
  distinct() %>% 
  drop_na() %>% 
  pivot_wider(names_from = threads_per_chain, 
              values_from = optimum_grainsize_vs_grainsize_1)  


t_3 <- blah %>% 
  group_by(solver) %>% 
  mutate(baseline_no_parallel = time[parallel == "Not Parallel"]) %>% 
  group_by(solver, threads_per_chain) %>% 
  mutate(optimum_grainsize = case_when(all(grainsize == 0) ~ NA_real_,
                                       TRUE ~ time[2])) %>% 
  summarize(optimum_grainsize_vs_no_parallel = 
              baseline_no_parallel/optimum_grainsize) %>% 
  ungroup() %>% 
  distinct() %>% 
  drop_na() %>% 
  pivot_wider(names_from = threads_per_chain, 
              values_from = optimum_grainsize_vs_no_parallel)

t_1
t_2
t_3



