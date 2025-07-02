
library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(mixtools)



source("D:/Users/buchy/Documents/Stage/utils_simulation_functions.R")
source("D:/Users/buchy/Documents/Stage/utils_kalman_functions.R")

# Modele
# X_0 ~ N(m0_x, S_x)
# X_{t} ~ N(F_x, X_{t-1}, S_x) t >= 1
# Y_t ~ N(A_y X_{t}, S_y) t >=0

full_pars_list <- list(dim_x = 1,
                       dim_y = 1,
                       m0_x = as.numeric(70),
                       S_x = as.matrix(2), # Must be a matrix
                       F_x = as.matrix(0.995), 
                       A_y = as.matrix(1), # must be a matrix...
                       S_y = as.matrix(1))

n_obs <- 500
data_1d <- get_data(n_steps = n_obs, full_pars_list) #Données 1d
filter_1d <- get_filtering_kalman(pars = full_pars_list, y_obs = data_1d$y)

objectif <- function(x) return(dnorm(x, filter_1d$mu_hat[1,], sqrt(filter_1d$V_hat[1])))



true_model <- list(get_initial_density = function(x, pars){
  dmvnorm(as.matrix(x), pars$m0_x, pars$S_x)
}, 
get_transition_density = function(x, ancestors, pars){
  # ancestors is either a N vector (if d_X = 1)
  # or  N * d_X matrix
  x <- as.matrix(x)
  ancestors <- as.matrix(ancestors)
  sapply(1:nrow(x),
         function(i){
           dmvnorm(x[i, ], as.numeric(pars$F_x %*% ancestors[i, ]), pars$S_x)
         })
},
get_emission_density = function(x, y, pars){
  x <- as.matrix(x) # N * d_X
  # y is single observation
  sapply(1:nrow(x),
         function(i){
           dmvnorm(y, as.numeric(pars$A_y %*% x[i, ]), pars$S_y)
         })
})

target_model <- list(get_x0_knowing_y0 = function(x, y, pars){
  x <- as.matrix(x)
  true_model$get_initial_density(x, pars) * 
    sapply(1:nrow(x),
           function(i){
             true_model$get_emission_density(x = x[i, ], y = y, pars = pars)
           })
  }, 
  get_xt_knowing_y_xtm1 = function(x, y, ancestors, pars){
  x <- as.matrix(x) # N * d_X
  ancestors <- as.matrix(ancestors)
  # y is single observation
  sapply(1:nrow(x),
         function(i){
           true_model$get_emission_density(x = x[i, ], y = y, pars = pars)
         }) *
   true_model$get_transition_density(x = x, ancestors = ancestors, pars = pars)
})



q_proposal <- list(get_initial_density = true_model$get_initial_density,
                   get_initial_samples = function(n, pars)
                     return(rmvnorm(n, as.numeric(pars$m0_x), as.matrix(pars$S_x))),
                   get_transition_density = function(x, ancestors, y, pars){
                     # Ne dépend pas de y pour le moment, mais pourrait...
                     true_model$get_transition_density(x = x, 
                                                       ancestors = ancestors, 
                                                       pars = pars)
                     },
                   get_transition_samples = function(ancestors, y, pars){
                     ancestors <- as.matrix(ancestors)
                     apply(ancestors, 1, function(anc_){
                       rmvnorm(1, as.numeric(pars$F_x %*% anc_), pars$S_x)
                     }, simplify = FALSE) |> 
                       do.call(what = rbind)
                     })


#fonction SIS-------------------------------------------------------------------------------------------------------------------------------------------------

get_SIS <- function(target_model, # Alist with elements ..... 
                    q_proposal,  # A list with elements ....
                    y_obs, # A matrix of size n_obs * dim_y
                    pars, # A list with all model elements
                    n_particles = 100){
  # target_model is a list having at least two elements which are functions
  # - get_x0_knowing_y0
  # - get_xt_knowing_y_xtm1
  
  # Initialisation des objets
  n_obs <- nrow(y_obs)
  particles <-  array(NA, dim = c(n_particles, pars$dim_x, n_obs))
  filtering_weights <- matrix(data = rep(0, n_particles * n_obs), nrow = n_obs)
  
  # Initialisation du filtre particulaire
  
  particles[,,1] <- q_proposal$get_initial_samples(n_particles, pars)
  W <- target_model$get_x0_knowing_y0(particles[,,1], 
                                      y_obs[1, ], pars) / 
    q_proposal$get_initial_density(particles[,,1], pars)
  filtering_weights[1,] <- W / sum(W)
  for (i in 2:n_obs){
    #normalisation
    particles[,,i] <- q_proposal$get_transition_samples(ancestors = particles[,, i - 1], 
                                                        y = y_obs[i, ], pars)
    W <- target_model$get_xt_knowing_y_xtm1(x = particles[,,i], y = y_obs[i, ], 
                                            ancestors = particles[,, i - 1], pars = pars) / 
        q_proposal$get_transition_density(x = particles[,,i], y = y_obs[i, ], 
                                          ancestors = particles[,, i - 1],
                                          pars = pars)
    filtering_weights[i,] <- W / sum(W)
  }
  return(list(particles = particles, 
              filtering_weights= filtering_weights))
}

SIS <- get_SIS(target_model = target_model,
               q_proposal = q_proposal,
               y_obs =  data_1d$y,
               pars = full_pars_list,
               n_particles = 100)

get_df_from_SIS_1d <- function(res_SIS){
  particle_df <- res_SIS$particles[,1,] |> 
    t() |> 
    as.data.frame() |> 
    rowid_to_column(var = "t") |> 
    pivot_longer(cols = -c("t"),
                 names_to = "Particle_index", 
                 values_to = "Particle_value",
                 names_prefix = "V")
  particle_weight <- res_SIS$filtering_weights |> 
    as.data.frame() |> 
    rowid_to_column(var = "t") |> 
    pivot_longer(cols = -c("t"),
                 names_to = "Particle_index", 
                 values_to = "Particle_weight",
                 names_prefix = "V")
  return(inner_join(particle_df, particle_weight, by = c("t", "Particle_index")))
}

SIS_df <- get_df_from_SIS_1d(SIS)

kalman_mean <- filter_1d$mu_hat

ggplot(SIS_df) + 
  aes(x = t, y = Particle_value) + 
  geom_line(aes(group = Particle_index)) +
  geom_line(data = data_1d$x |> 
              as.data.frame() |> 
              rename(x = V1) |> 
              rowid_to_column(var = "t"), 
            aes(y = x, color = "Truth")) + 
  geom_line(data = data.frame(t = 1:nrow(data_1d$x),  y = kalman_mean[,1]),
            aes(y = y, color = "Kalman"))


mean_filtering_weights <- colMeans(SIS$filtering_weights)

plot(mean_filtering_weights)

get_ESS <- function(W){
  # W is a matrix with T rows (number of time steps) and N columns (number of MC samples)
  apply(W, 1, function(w) 1 / sum(w^2))
}

apply(SIS$filtering_weights, 1, function(w) 1 / sum(w^2))


# SIS Resampling -------------------------------------------------------------------------------------------------------------------------------------------------


get_SIR <- function(target_model, # A list with elements ..... 
                    q_proposal,  # A list with elements ....
                    y_obs, # A matrix of size n_obs * dim_y
                    pars, # A list with all model elements
                    n_particles = 100){
  # target_model is a list having at least two elements which are functions
  # - get_x0_knowing_y0
  # - get_xt_knowing_y_xtm1
  
  SIS <- get_SIS(target_model, 
                 q_proposal, 
                 y_obs, 
                 pars,
                 n_particles = 100)
  
  particles <- SIS$particles
  filtering_weights <- SIS$filtering_weights

  for (i in 1:n_obs){

      #On sample
      sample_i <- sample(x = c(1:n_particles), size = n_particles, prob = SIS$filtering_weights[i,], replace = TRUE)

      for (j in 1:n_particles){

        particles[,,j] <- particles[,,sample_i[j]]
        filtering_weights[,j] <- rep(1/n_particles, n_obs)

      }
    }
  
  return(list(particles = particles, filtering_weights = filtering_weights))
}


N_thres <- 10
SIR <- get_SIR(target_model = target_model,
               q_proposal = q_proposal,
               y_obs =  data_1d$y,
               pars = full_pars_list,
               n_particles = 100)

# mean_W_tilde_SIR <- rep(0, dim(SIR$W_tilde)[1])
# for (i in 1:dim(SIR$W_tilde)[1]) mean_W_tilde_SIR[i] <- mean(SIR$W_tilde[,i])
# 
# plot(mean_W_tilde_SIR)
# abline(h = 1/100)
