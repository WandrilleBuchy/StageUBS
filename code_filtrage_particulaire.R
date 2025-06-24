rm(list)
library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(RVCompare)


source("utils_simulation_functions.R")
source("utils_kalman_functions.R")

pars_X0_1d <- list(dim_x = 1,
             dim_y = 1,
             m_x = 70,
             S_x = 1,
             A = 1,
             b = 0,
             matF = 1,
             S_y = 1)
pars_1d <- list(dim_x = 1, #pas de valeur moyenne 
                dim_y = 1,
                S_x = 2,
                A = 1,
                b = 0,
                matF = 0.995,
                S_y = 0.1)

data_1d <- get_data(n_sim = 2, pars_X0_1d, pars_1d) #Données 1d
filter_1d <- get_filtering_kalman(pars_X0_1d, pars_1d, data_1d)





## ----SNIS-------------------------------------------------------------------------------------------------------------------------------------------------

get_MCMC_SNIS <- function(f, u_target, q_proposal, observation, pars, m = 1000){
  samples <- q_proposal$get_samples(m, pars)
  W <- u_target(samples, observation, pars)/ q_proposal$get_density(samples, pars)
  W_tilde <- W / sum(W)
  I_hat <- sum(W_tilde * f(samples))
  return(list(W_tilde = W_tilde, X = samples, I_hat = I_hat))
}

q_proposal <- list(get_density = function(x, pars) return(mixtools::dmvnorm(x,
                                                                      pars$m_x,
                                                                      pars$S_x)),
               get_samples = function(n, pars) 
                 return(mixtools::rmvnorm(n, pars$m_x, as.matrix(pars$S_x))))

u_target <-  function(x, y, pars){
  sapply(x, function(x_){
    mixtools::dmvnorm(x_,
                      pars$m_x,
                      pars$S_x) * 
      mixtools::dmvnorm(y, as.numeric(pars$A %*% x_), pars$S_y)
  })
}
# Densite non normalisee

SNIS_X0 <- get_MCMC_SNIS(f = function(x) x, u_target, q_proposal, 
                         observation = data_1d$y[1, ], pars = pars_X0_1d,  m = 1000)
