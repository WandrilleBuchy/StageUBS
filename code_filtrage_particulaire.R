
library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(RVCompare)


source("D:/Users/buchy/Documents/Stage/utils_simulation_functions.R")
source("D:/Users/buchy/Documents/Stage/utils_kalman_functions.R")

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
                S_y = 1)

n_obs <- 2
data_1d <- get_data(n_sim = n_obs, pars_X0_1d, pars_1d) #Données 1d
filter_1d <- get_filtering_kalman(pars_X0_1d, pars_1d, data_1d)





## ----SNIS-------------------------------------------------------------------------------------------------------------------------------------------------

get_SNIS <- function(f, u_target, q_proposal, observation, pars, m = 1000){
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

SNIS_X0 <- get_SNIS(f = function(x) x, u_target, q_proposal, 
                    observation = data_1d$y[1, ], pars = pars_X0_1d,  m = 1000)


# Affichage de la vraie loi

plot(objectif, filter_1d$mu_hat[1,] - 10, filter_1d$mu_hat[1,] + 10)
points(SNIS_X0$X, rep(0, 1000), cex = 100 * SNIS_X0$W_tilde)

df = data.frame(x = SNIS_X0$X, w =  SNIS_X0$W_tilde)

ggplot(df) + 
  geom_point(aes(x = x, size = w), y = 0) +
  stat_function(fun = function(x) dnorm(x, 70, sqrt(0.5))) +
  stat_function(fun = function(x) dnorm(x, filter_1d$mu_hat[2], 
                                        sqrt(filter_1d$V_hat[,,2])), color = "red")


#fonction SIS ------ brouillon

get_SIS <- function(SNIS_X0, u_target, q_proposal, data, m = 1000){
  
  W_tilde <- matrix(data = rep(0,m*n_obs), nrow = n_obs)
  W_tilde[1,] <- SNIS_X0$W_tilde
  W_tilde_out <- array(data = rep(0, m))
  
  for (i in 2:n_obs){ 
    
    q_proposal <- 
    SNIS_i <- get_SNIS(f = function(x) x,
                       u_proposal = dnorm(x, filter_1d$mu_hat[i], sqrt(filter_1d$V_hat[,,i])),
                       q_proposal, 
                       observation = data$y[i, ], m)
    alpha[i,] <- 
    W_tilde[i,] <- W_tilde[i-1,] * alpha[i]
    W_tilde[]
    
    
  }
  
}



#Neff_hat <- 1/(sum(SNIS_X0$W_tilde^2))
