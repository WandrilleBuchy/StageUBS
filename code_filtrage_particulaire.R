
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

n_obs <- 4
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

get_SIS <- function(f, u_target, q_proposal, data, pars_X0, pars, m = 1000){
  
  SNIS_X0 <- get_SNIS(f, u_target, q_proposal, 
                      observation = data$y[1, ], pars_X0, m)
  traj_mat <- W_tilde <- matrix(data = rep(0, m * n_obs), nrow = n_obs)
  W_tilde[1,] <- SNIS_X0$W_tilde
  traj_mat[1,] <- SNIS_X0$X

  
  for (i in 2:n_obs){
  
    for (j in 1:m){ 
      
      pars$m_X <- traj_mat[i-1,j]
      traj_mat[i,j] <- q_proposal$get_samples(1,pars)
      
      if(q_proposal$get_density(traj_mat[i,j], pars) == 0) print(i,j)
      
      #normalisation
      W_tilde[i,j] <- W_tilde[i-1,j] * u_target(traj_mat[i,j], data$y[i,], pars) / q_proposal$get_density(traj_mat[i,j], pars)
      
      print(u_target(traj_mat[i,j], data$y[i,], pars))
      
    }
    W_tilde[i,] <- W_tilde[i,]/sum(W_tilde[i,])
  }
  return(list(traj_mat = traj_mat, W_tilde = W_tilde))
}


test <- get_SIS(f = function(x) x, u_target, q_proposal, data_1d, pars_X0_1d, pars_1d, m = 1000)

