rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(furrr)
library(tidyverse)
library(patchwork)

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

## Work from here down

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
