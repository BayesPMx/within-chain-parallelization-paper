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

tmp <- zoom %>% 
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
                             "Torsten Group ODE MPI")))

(p_1_tvcl <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVCL") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = parallel), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Parallel",
                       values = c("Parallel" = "blue",
                                  "Not Parallel" = "red")) +
    facet_grid(Variable ~ solver, scales = "free") +
    # facet_wrap(Variable ~ solver, scales = "free") +
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()))

(p_1_tvvc <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVVC") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = parallel), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Parallel",
                       values = c("Parallel" = "blue",
                                  "Not Parallel" = "red")) +
    facet_grid(Variable ~ solver, scales = "free") +
    # facet_wrap(Variable ~ solver, scales = "free") +
    theme(legend.position = "none",
          axis.title.x = element_blank()))

(p_1_tvka <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVKA") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = parallel), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Parallel",
                       values = c("Parallel" = "blue",
                                  "Not Parallel" = "red")) +
    facet_grid(Variable ~ solver, scales = "free") +
    # facet_wrap(Variable ~ solver, scales = "free") +
    theme(legend.position = "none",
          axis.title.y = element_blank()))

(p_2_tvcl <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVCL") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = solver), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Solver",
                       values = c("Analytical Threaded" = "purple",
                                  "Linear ODE Threaded" = "darkgreen",
                                  "General ODE Threaded" = "orange",
                                  "Torsten Group ODE MPI" = "magenta")) +
    facet_grid(Variable ~ parallel, scales = "free") +
    # facet_wrap(Variable ~ parallel, scales = "free") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()))

(p_2_tvvc <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVVC") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = solver), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Solver",
                       values = c("Analytical Threaded" = "purple",
                                  "Linear ODE Threaded" = "darkgreen",
                                  "General ODE Threaded" = "orange",
                                  "Torsten Group ODE MPI" = "magenta")) +
    facet_grid(Variable ~ parallel, scales = "free") +
    # facet_wrap(Variable ~ parallel, scales = "free") +
    theme(legend.position = "none",
          axis.title.x = element_blank()))

(p_2_tvka <- tmp %>% 
    bind_rows(tmp %>% 
                filter(solver == "General ODE Threaded", 
                       parallel == "Not Parallel") %>% 
                mutate(solver = "Torsten Group ODE MPI")) %>%
    mutate(solver = factor(solver, levels = c("Analytical Threaded",
                                              "Linear ODE Threaded",
                                              "General ODE Threaded",
                                              "Torsten Group ODE MPI")),
           parallel = factor(parallel, levels = c("Not Parallel", "Parallel"))) %>%
    select(starts_with("TV"), solver, parallel) %>% 
    pivot_longer(cols = starts_with("TV"), 
                 names_to = "Variable", 
                 values_to = "Value") %>% 
    mutate(Variable = factor(Variable, levels = c("TVCL", "TVVC", "TVKA"))) %>% 
    filter(Variable == "TVKA") %>% 
    ggplot() +
    geom_density(aes(x = Value, color = solver), alpha = 0.1) +
    theme_bw() +
    scale_color_manual(name = "Solver",
                       values = c("Analytical Threaded" = "purple",
                                  "Linear ODE Threaded" = "darkgreen",
                                  "General ODE Threaded" = "orange",
                                  "Torsten Group ODE MPI" = "magenta")) +
    facet_grid(Variable ~ parallel, scales = "free") +
    # facet_wrap(Variable ~ parallel, scales = "free") +
    theme(legend.position = "bottom",
          axis.title.y = element_blank()))


p_1_tvcl /
  p_1_tvvc /
  p_1_tvka /
  p_2_tvcl /
  p_2_tvvc /
  p_2_tvka

ggsave(filename = "Figures/figure_3.pdf", device = "pdf", dpi = 1200, 
       width = 5, height = 6.5, units = "in", scale = 1.5)

boo <- zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  inner_join(zoom, 
             by = c("solver", "run_number", "threads_per_chain", "time"))

all_summary_tables <- zoom %>% 
  group_by(solver, threads_per_chain, run_number) %>% 
  distinct(time) %>% 
  ungroup(run_number) %>%
  arrange(time, .by_group = TRUE) %>% 
  slice(-c(1, n())) %>% 
  inner_join(zoom, 
             by = c("solver", "run_number", "threads_per_chain", "time")) %>% 
  group_by(solver, threads_per_chain, run_number, time) %>% 
  group_modify(~ { 
    .x %>%
      select(-contains("omega"), -contains("sigma")) %>%
      as_draws_array() %>%
      summarize_draws(mean, median, sd, mcse_mean,
                      ~quantile2(.x, probs = c(0.025, 0.975)), rhat,
                      ess_bulk, ess_tail)
  }) %>% 
  ungroup()

# all_summary_tables %>%
#   mutate(Variable = factor(variable, levels = c("TVCL", "TVVC", "TVKA"))) %>%
#   pivot_longer(cols = c(mean, ess_bulk, ess_tail), values_to = "Value",
#                names_to = "Estimate") %>%
#   ggplot(aes(x = Variable,
#              y = Value)) +
#   geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
#   scale_y_continuous(name = "Posterior Mean",
#                      trans = "identity") +
#   theme_bw() +
#   facet_wrap(variable ~ Estimate, scales = "free")

data_for_mcse_ess_plot <- all_summary_tables %>%
  select(solver, threads_per_chain, run_number, time, variable, mean, mcse_mean, 
         starts_with("ess")) %>% 
  pivot_longer(cols = c(mean, mcse_mean, ess_bulk, ess_tail), values_to = "Value",
               names_to = "Estimate") %>%
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
         parallel = threads_per_chain > 0,
         parallel = if_else(parallel, "Parallel", "Not Parallel") %>% 
           factor(levels = c("Not Parallel", "Parallel")),
         Variable = factor(variable, levels = c("TVCL", "TVVC", "TVKA")))

data_for_mcse_ess_plot <- data_for_mcse_ess_plot %>% 
  bind_rows(data_for_mcse_ess_plot %>% 
              filter(solver == "General ODE Threaded", 
                     threads_per_chain == 0) %>% 
              mutate(solver = "Torsten Group ODE MPI")) %>% 
  mutate(solver = factor(solver, 
                         levels = c("Analytical Threaded",
                                    "Linear ODE Threaded",
                                    "General ODE Threaded",
                                    "Torsten Group ODE MPI")))

p_mean <- data_for_mcse_ess_plot %>% 
  filter(Estimate == "mean") %>% 
  ggplot(aes(x = solver,
             y = Value)) +
  # geom_boxplot(color = "black", size = 0.4, alpha = 0.9) +
  geom_boxplot(aes(color = parallel), size = 1, alpha = 0.9) +
  scale_y_continuous(name = "Posterior Mean",
                     trans = "identity") +
  theme_bw() +
  theme(legend.position = "top",
        # strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_manual(name = "Parallel",
                     values = c("Parallel" = "blue",
                                "Not Parallel" = "red")) +
  # scale_color_manual(name = "Solver",
  #                    guide = "none",
  #                    values = c("Analytical Threaded" = "purple",
  #                               "Linear ODE Threaded" = "darkgreen",
  #                               "General ODE Threaded" = "orange",
  #                               "Torsten Group ODE MPI" = "magenta")) +
  facet_wrap(~ Variable, nrow = 1, scales = "free")

p_mcse <- data_for_mcse_ess_plot %>% 
  filter(Estimate == "mcse_mean") %>% 
  ggplot(aes(x = solver,
             y = Value)) +
  # geom_boxplot(color = "black", size = 0.4, alpha = 0.9) +
  geom_boxplot(aes(color = parallel), size = 1, alpha = 0.9) +
  scale_y_continuous(name = "MCSE for the Mean",
                     trans = "identity") +
  theme_bw() +
  theme(# legend.position = "top",
    # strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_color_manual(name = "Parallel",
                     guide = "none",
                     values = c("Parallel" = "blue",
                                "Not Parallel" = "red")) +
  # scale_color_manual(name = "Solver",
  #                    guide = "none",
  #                    values = c("Analytical Threaded" = "purple",
  #                               "Linear ODE Threaded" = "darkgreen",
  #                               "General ODE Threaded" = "orange",
  #                               "Torsten Group ODE MPI" = "magenta")) +
  facet_wrap(~ Variable, nrow = 1, scales = "free")

p_ess_bulk <- data_for_mcse_ess_plot %>% 
  filter(Estimate == "ess_bulk") %>% 
  ggplot(aes(x = solver,
             y = Value)) +
  # geom_boxplot(color = "black", size = 0.4, alpha = 0.9) +
  geom_boxplot(aes(color = parallel), size = 1, alpha = 0.9) +
  scale_y_continuous(name = latex2exp::TeX("$ESS_{bulk}"),
                     trans = "identity") +
  theme_bw() +
  theme(# legend.position = "top",
    # strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_color_manual(name = "Parallel",
                     guide = "none",
                     values = c("Parallel" = "blue",
                                "Not Parallel" = "red")) +
  # scale_color_manual(name = "Solver",
  #                    guide = "none",
  #                    values = c("Analytical Threaded" = "purple",
  #                               "Linear ODE Threaded" = "darkgreen",
  #                               "General ODE Threaded" = "orange",
  #                               "Torsten Group ODE MPI" = "magenta")) +
  facet_wrap(~ Variable, nrow = 1, scales = "free")

p_ess_tail <- data_for_mcse_ess_plot %>% 
  filter(Estimate == "ess_tail") %>% 
  ggplot(aes(x = solver,
             y = Value)) +
  # geom_boxplot(color = "black", size = 0.4, alpha = 0.9) +
  geom_boxplot(aes(color = parallel), size = 1, alpha = 0.9) +
  scale_y_continuous(name = latex2exp::TeX("$ESS_{tail}"),
                     trans = "identity") +
  theme_bw() +
  theme(# legend.position = "bottom",
    # strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = -60, hjust = 0.25, vjust = 0)) +
  scale_color_manual(name = "Parallel",
                     guide = "none",
                     values = c("Parallel" = "blue",
                                "Not Parallel" = "red")) +
  # scale_color_manual(name = "Solver",
  #                    guide = "none",
  #                    values = c("Analytical Threaded" = "purple",
  #                               "Linear ODE Threaded" = "darkgreen",
  #                               "General ODE Threaded" = "orange",
  #                               "Torsten Group ODE MPI" = "magenta")) +
  facet_wrap(~ Variable, nrow = 1, scales = "free")

p_mean /
  p_mcse /
  p_ess_bulk /
  p_ess_tail +
  theme(plot.margin = unit(c(0, 0.35, 0, 0), 
                           "inches"))

ggsave(filename = "Figures/figure_4.pdf", device = "pdf", dpi = 1200,
       width = 5, height = 6.5, units = "in", scale = 3)

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
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, 
                                  size = 20),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  ggtitle("Dataset 1") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


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
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

p_dataset_1 <- p_time_vs_processes /
  p_speedup_vs_processes +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

read_threaded_fits <- function(threads_per_chain){
  
  file_name <- str_c("depot_2cmt_linear_friberg/Stan/Fits/", 
                     threads_per_chain, "_threads.rds")
  
  fit <- read_rds(file_name)
  
  return(fit$draws(c("TVCL", "TVVC", "TVQ", "TVVP", "TVKA",
                     "TVMTT", "TVCIRC0", "TVGAMMA", "TVALPHA",
                     "omega", "omega_pd",
                     "sigma_p", "sigma_p_pd"), 
                   format = "draws_df") %>% 
           mutate(threads_per_chain = threads_per_chain,
                  solver = "rk45",
                  time = max(fit$time()$chains$total)))
  
}

read_mpi_fits <- function(threads_per_chain){
  
  files <- str_c("depot_2cmt_linear_friberg/Stan/Fits/Torsten_General/torsten_general_", 
                 threads_per_chain, "_jobs_", 1:4, ".csv")
  
  fit <- as_cmdstan_fit(files = files)
  
  return(fit$draws(c("TVCL", "TVVC", "TVQ", "TVVP", "TVKA",
                     "TVMTT", "TVCIRC0", "TVGAMMA", "TVALPHA",
                     "omega", "omega_pd",
                     "sigma_p", "sigma_p_pd"), 
                   format = "draws_df") %>% 
           mutate(threads_per_chain = threads_per_chain,
                  solver = "torsten_general",
                  time = max(fit$time()$chains$total)))
  
}

threaded_fits <- c(0, 1, 4, 12, 24) %>% 
  map_dfr(.f = read_threaded_fits)

mpi_fits <- c(1, 4, 12, 24) %>%
  # mpi_fits <- c(1, 4, 12) %>% 
  map_dfr(.f = read_mpi_fits) %>% 
  bind_rows(threaded_fits %>% 
              filter(threads_per_chain == 0) %>% 
              mutate(solver = "torsten_general")) %>% 
  arrange(threads_per_chain)

fits <- threaded_fits %>% 
  bind_rows(mpi_fits)

tmp <- fits %>% 
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
                             "Torsten Group ODE MPI")))

fits %>% 
  group_by(solver, threads_per_chain) %>% 
  distinct(time) %>% 
  ungroup(threads_per_chain) %>% 
  mutate(timeh = time/3600,
         across(ends_with("time"), ~ ./.[1]),
         speedup = 1/time)

p_time_vs_processes <- fits %>% 
  group_by(solver, threads_per_chain) %>% 
  distinct(time) %>% 
  ungroup() %>% 
  mutate(processes = if_else(threads_per_chain == 0, "Not Parallel", 
                             as.character(threads_per_chain)) %>% 
           factor(levels = c("Not Parallel", as.character(unique(threads_per_chain)))),
         across(ends_with("time"), ~./3600),
         solver = case_when(solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>%
  ggplot(aes(x = processes, y = time, 
             group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     values = c("General ODE Threaded" = "orange",
                                "Torsten Group ODE MPI" = "magenta")) +
  scale_y_continuous(name = "Elapsed Time (h)",
                     trans = "log10") +
  annotation_logticks(sides = "l") +
  scale_x_discrete(name = "Processes per Chain") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, 
                                  size = 20),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  ggtitle("Dataset 3")


p_speedup_vs_processes <- fits %>% 
  group_by(solver, threads_per_chain) %>% 
  distinct(time) %>% 
  ungroup() %>%
  group_by(solver) %>% 
  arrange(threads_per_chain, .by_group = TRUE) %>% 
  mutate(across(ends_with("time"), ~ ./.[1]),
         speedup = 1/time) %>% 
  ungroup() %>% 
  mutate(processes = if_else(threads_per_chain == 0, "Not Parallel", 
                             as.character(threads_per_chain)) %>% 
           factor(levels = c("Not Parallel", as.character(unique(threads_per_chain)))),
         solver = case_when(solver == "rk45" ~ "General ODE Threaded",
                            solver == "torsten_general" ~ 
                              "Torsten Group ODE MPI",
                            TRUE ~ NA_character_) %>% 
           factor(levels = c("General ODE Threaded",
                             "Torsten Group ODE MPI"))) %>% 
  ggplot(aes(x = processes, y = speedup, group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     values = c("General ODE Threaded" = "orange",
                                "Torsten Group ODE MPI" = "magenta")) +
  scale_y_continuous(name = "Speedup",
                     trans = "identity",
                     breaks = seq(1, 9, by = 2),
                     labels = seq(1, 9, by = 2)) +
  scale_x_discrete(name = "Processes per Chain") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))

p_dataset_3 <- p_time_vs_processes /
  p_speedup_vs_processes +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_dataset_1 | p_dataset_3

ggsave(filename = "Figures/figure_5.pdf", device = "pdf", dpi = 1200,
       width = 6.5, height = 4.5, units = "in", scale = 1.85)


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

p_grainsize_dataset_1 <- quick_grainsize %>% 
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
        legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  geom_hline(data = summary_of_not_parallel_threaded,
             mapping = aes(yintercept = mean_time, color = solver),
             linetype = "dashed") +
  facet_wrap(~ processes, 
             labeller = as_labeller(function(x) str_c("Processes per Chain: ", 
                                                      x)))

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

read_grainsize <- function(threads_per_chain, grainsize){
  
  if(threads_per_chain > 0){
    if(grainsize > 1){
      
      file_name <- str_c("depot_2cmt_linear_friberg/Stan/Fits/", 
                         threads_per_chain, "_threads_", 
                         grainsize, "_grainsize.rds")
      
    }else{ # grainsize == 1
      
      file_name <- str_c("depot_2cmt_linear_friberg/Stan/Fits/", 
                         threads_per_chain, "_threads.rds")
      
    }
  }else{
    
    file_name <- str_c("depot_2cmt_linear_friberg/Stan/Fits/0_threads.rds")
    
  }
  
  fit <- read_rds(file_name)
  
  # The ones fit explicitly with cmdstanr and saved directly as CmdStanMCMC
  # objects (the multi-threaded ones) have fit$time()$total along with the 
  # total time by chain, but the ones fit and saved in CSVs and then converted
  # into CmdStanMCMC objects (Torsten's group solver) have the total time by 
  # chain but fit$time()$total = NA. Therefore, to do the same thing for all 
  # fits, I'll take the maximum time of the 4 chains as the elapsed time.
  return(fit$draws(c("TVCL", "TVVC", "TVQ", "TVVP", "TVKA", 
                     "TVMTT", "TVCIRC0", "TVGAMMA", "TVALPHA", 
                     "omega", "omega_pd", "sigma_p", "sigma_p_pd"), 
                   format = "draws_df") %>% 
           mutate(threads_per_chain = threads_per_chain,
                  grainsize = grainsize,
                  time = max(fit$time()$chains$total)))
  
}

plan(multisession, workers = min(parallel::detectCores(), 8))

quick_grainsize_3 <- tibble(threads_per_chain = 0, grainsize = 0) %>% 
  bind_rows(expand_grid(threads_per_chain = c(4, 24), 
                        grainsize = c(1, 2, 8, 12))) %>%
  future_pmap_dfr(.f = read_grainsize)

plan(sequential)


data_for_plot_3 <- quick_grainsize_3 %>% 
  distinct(grainsize, threads_per_chain, time) %>% 
  mutate(solver = "General ODE Threaded" %>% 
           factor(levels = c("General ODE Threaded")),
         processes = factor(threads_per_chain),
         grainsize = factor(grainsize),
         time = time/3600)

p_grainsize_dataset_3 <- data_for_plot_3 %>%
  filter(threads_per_chain != 0) %>% 
  ggplot(aes(x = grainsize, y = time, 
             group = solver, color = solver)) +
  geom_point() +
  geom_line() +
  theme_bw() + 
  scale_color_manual(name = "Solver",
                     guide = "none",
                     values = c("General ODE Threaded" = "orange")) +
  scale_y_continuous(name = "Elapsed Time (h)",
                     trans = "log10") +
  annotation_logticks(sides = "l") +
  scale_x_discrete(name = "Grainsize") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.position = "bottom") +
  geom_hline(data = data_for_plot_3 %>%
               filter(threads_per_chain == 0) %>% 
               select(solver, time),
             mapping = aes(yintercept = time, color = solver),
             linetype = "dashed") +
  facet_wrap(~ processes, 
             labeller = as_labeller(function(x) str_c("Processes per Chain: ", 
                                                      x)))

p_grainsize_dataset_1 /
  p_grainsize_dataset_2 /
  p_grainsize_dataset_3 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = "Figures/figure_6.pdf", device = "pdf", dpi = 1200,
       width = 6.5, height = 7.5, units = "in", scale = 1.25)


## Figures in Supplementary Info

nonmem_data_1 <- read_csv("depot_1cmt_linear/Data/depot_1cmt_prop.csv",
                          na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(dataset = "Dataset 1",
         DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

nonmem_data_2 <- read_csv("depot_1cmt_linear_ragged/Data/depot_1cmt_prop.csv",
                          na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(dataset = "Dataset 2",
         DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

nonmem_data <- nonmem_data_1 %>% 
  bind_rows(nonmem_data_2)

ggplot(nonmem_data %>%
         group_by(ID) %>%
         mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
         ungroup() %>%
         filter(mdv == 0)) +
  geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
  geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
  scale_color_discrete(name = "Dose (mg)") +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     limits = c(NA, NA),
                     trans = "log10") +
  scale_x_continuous(name = "Time (d)",
                     breaks = seq(0, 672, by = 168),
                     labels = seq(0, 28, by = 7),
                     limits = c(0, NA)) +
  theme_bw(18) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(size = 2),
        legend.position = "bottom") +
  facet_wrap(~ dataset, scales = "free_x")

nonmem_data <- read_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_prop_friberg_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

(p_pk <- ggplot(nonmem_data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(amt)/1000)) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), cmt == 2, bloq == 0) %>% 
                  mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
                                                cmt == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(ng/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 216, by = 24),
                       labels = seq(0, 216/24, by = 24/24),
                       limits = c(0, 216)) +
    scale_color_manual(name = "Dose (mg)",
                       values = c("10" = "red", "20" = "blue")) +
    # facet_wrap(~ cmt, scales = "free") +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom"))

(p_pd <- ggplot(nonmem_data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(amt)/1000)) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), cmt == 4, bloq == 0) %>% 
                  mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
                                                cmt == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Neutrophils $(\\times 10^9/L)"),
                       trans = "identity") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 672, by = 168),
                       labels = seq(0, 28, by = 7),
                       limits = c(0, 672)) +
    scale_color_manual(name = "Dose (mg)",
                       values = c("10" = "red", "20" = "blue")) +
    # facet_wrap(~ cmt, scales = "free") +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 2),
          legend.position = "bottom"))

p_pk + 
  p_pd +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")