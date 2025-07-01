
library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
#library(RVCompare)
#library(survey)


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

n_obs <- 30
data_1d <- get_data(n_sim = n_obs, pars_X0_1d, pars_1d) #Données 1d
filter_1d <- get_filtering_kalman(pars_X0_1d, pars_1d, data_1d)

objectif <- function(x) return(dnorm(x, filter_1d$mu_hat[1,], sqrt(filter_1d$V_hat[1])))





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

#plot(objectif, filter_1d$mu_hat[1,] - 10, filter_1d$mu_hat[1,] + 10)
#points(SNIS_X0$X, rep(0, 1000), cex = 100 * SNIS_X0$W_tilde)

df = data.frame(x = SNIS_X0$X, w =  SNIS_X0$W_tilde)

ggplot(df) + 
  geom_point(aes(x = x, size = w), y = 0) +
  stat_function(fun = function(x) dnorm(x, 70, sqrt(0.5))) +
  stat_function(fun = function(x) dnorm(x, filter_1d$mu_hat[2], 
                                        sqrt(filter_1d$V_hat[,,2])), color = "red")



#fonction SIS-------------------------------------------------------------------------------------------------------------------------------------------------

get_SIS <- function(f, u_target, q_proposal, data, pars_X0, pars, n_obs, m = 1000){
  
  SNIS_X0 <- get_SNIS(f, u_target, q_proposal, 
                      observation = data$y[1, ], pars_X0, m)
  traj_mat <- W_tilde <- matrix(data = rep(0, m * n_obs), nrow = n_obs)
  I_hat <- rep(0, n_obs)
  W_tilde[1,] <- SNIS_X0$W_tilde
  traj_mat[1,] <- SNIS_X0$X
  
  
  for (i in 2:n_obs){
    
    for (j in 1:m){ 
      
      pars_j <- pars
      pars_j$m_x <- traj_mat[i-1, j]
      traj_mat[i,j] <- q_proposal$get_samples(1, pars_j)
      
      #normalisation
      W_tilde[i,j] <- W_tilde[i-1,j] * u_target(traj_mat[i,j], data$y[i,], pars_j) / q_proposal$get_density(traj_mat[i,j], pars_j)
    }
    W_tilde[i,] <- W_tilde[i,]/sum(W_tilde[i,])
  }
  for (k in 1:n_obs){ 
    
    I_hat[k] <- sum(W_tilde[k,] * f(traj_mat[k,]))
  }
  return(list(traj_mat = traj_mat, W_tilde = W_tilde, I_hat = I_hat))
}

SIS <- get_SIS(f = function(x) x, u_target, q_proposal, data_1d, pars_X0_1d, pars_1d, n_obs, m = 1000)


mean_W_tilde <- rep(0, dim(SIS$W_tilde)[1])
for (i in 1:dim(SIS$W_tilde)[1]) mean_W_tilde[i] <- mean(SIS$W_tilde[,i])

plot(mean_W_tilde)


# SIS Resampling -------------------------------------------------------------------------------------------------------------------------------------------------

get_SIR <- function(f, u_target, q_proposal, data, pars_X0, pars, N_thres, n_obs, m = 1000){

  SIS <- get_SIS(f, u_target, q_proposal, data, pars_X0, pars, n_obs, m)
  traj_mat <- SIS$traj_mat
  W_tilde <- SIS$W_tilde
  I_hat <- rep(0, n_obs)
  
  for (i in 1:n_obs){
    
    #Création de la condition 
    N_eff_hat <- 1/sum(SIS$W_tilde[,i]^2)
    if (N_eff_hat < N_thres){
      
      #On sample
      sample_i <- sample(x = c(1:m), size = m, prob = SIS$W_tilde[i,], replace = TRUE)
      
      for (j in 1:m){
        
        traj_mat[,j] <- traj_mat[,sample_i[j]]
        W_tilde[,j] <- rep(1/m, n_obs)
        
      }
    }
  }
  for (k in 1:n_obs){ 
    
    I_hat[k] <- sum(W_tilde[k,] * f(traj_mat[k,]))
  }
  
  return(list(traj_mat = traj_mat, W_tilde = W_tilde, I_hat = I_hat))
}


N_thres <- 10000
SIR <- get_SIR(f = function(x) x, u_target, q_proposal, data_1d, pars_X0_1d, pars_1d, N_thres, n_obs, m = 1000)

mean_W_tilde_SIR <- rep(0, dim(SIR$W_tilde)[1])
for (i in 1:dim(SIR$W_tilde)[1]) mean_W_tilde_SIR[i] <- mean(SIR$W_tilde[,i])

plot(mean_W_tilde_SIR)
abline(h = 1/1000)

