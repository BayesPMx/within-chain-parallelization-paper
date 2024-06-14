rm(list = ls())
cat("\014")

library(tidyjson)
library(tidyverse)

n_subjects <- read_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_prop_friberg_prop.csv",
                       na = ".") %>% 
  distinct(ID) %>%
  count() %>%
  deframe()

write_inits <- function(chain){
  list(TVCL = rlnorm(1, log(0.6), 0.3),
       TVVC = rlnorm(1, log(18), 0.3),
       TVQ = rlnorm(1, log(2), 0.3),
       TVVP = rlnorm(1, log(40), 0.3),
       TVKA = rlnorm(1, log(1.3), 0.3),
       omega = rlnorm(5, log(0.3), 0.3),
       L = diag(5),
       sigma_p = rlnorm(1, log(0.2), 0.3),
       Z = matrix(rnorm(n_subjects*5), ncol = n_subjects, nrow = 5),
       TVMTT = rlnorm(1, log(120), 0.3),
       TVCIRC0 = rlnorm(1, log(5), 0.3),
       TVGAMMA = rlnorm(1, log(0.2), 0.3),
       TVALPHA = rlnorm(1, log(3e-4), 0.3),
       omega_pd = rlnorm(4, log(0.35), 0.3),
       L_pd = diag(4),
       sigma_p_pd = rlnorm(1, log(0.2), 0.3),
       Z_pd = matrix(rnorm(n_subjects*4), ncol = n_subjects, nrow = 4)) %>% 
    map(round, 5) %>% 
    write_stan_json(str_c("depot_2cmt_linear_friberg/Data/Inits/inits_", 
                          chain, ".json"))

}

expand_grid(chain = 1:4) %>% 
  {map(.$chain, write_inits)}
