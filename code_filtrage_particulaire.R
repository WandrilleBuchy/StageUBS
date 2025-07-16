
library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(mixtools)


source("utils_get_toy_dataset_linear_gaussian_HMM.R")
source("utils_kalman_functions.R")
source("utils_particle_filter_functions.R")


# Definition of models for particle filter --------------------------------


true_model <- list(get_initial_density = function(x, pars){
  mixtools::dmvnorm(as.matrix(x), pars$m0_x, pars$S_x)
}, get_transition_density = function(x, ancestors, pars){
  # ancestors is either a N vector (if d_X = 1)
  # or  N * d_X matrix
  x <- matrix(x, ncol = pars$dim_x)
  ancestors <- matrix(ancestors, ncol = pars$dim_x)
  
  sapply(1:nrow(x),
         function(i){
           mixtools::dmvnorm(x[i, ], as.numeric(pars$F_x %*% ancestors[i, ]), pars$S_x)
         })
},
get_emission_density = function(x, y, pars){
  x <- as.matrix(x) # N * d_X
  # y is single observation
  sapply(1:nrow(x),
         function(i){
           mixtools::dmvnorm(y, as.numeric(pars$A_y %*% x[i, ]), pars$S_y)
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
                     return(mixtools::rmvnorm(n, as.numeric(pars$m0_x), as.matrix(pars$S_x))),
                   get_transition_density = function(x, ancestors, y, pars){
                     # Ne dépend pas de y pour le moment, mais pourrait...
                     true_model$get_transition_density(x = x, 
                                                       ancestors = ancestors, 
                                                       pars = pars)
                     },
                   get_transition_samples = function(ancestors, y, pars){
                     ancestors <- as.matrix(ancestors)
                     apply(ancestors, 1, function(anc_){
                       mixtools::rmvnorm(1, as.numeric(pars$F_x %*% anc_), pars$S_x)
                     }, simplify = FALSE) |> 
                       do.call(what = rbind)
                     })



# Performing particle filter ----------------------------------------------

n_particles = 500
n_replicates = 100 # Nombre de réplicats des estimateurs
methods <- list(SIS = list(threshold = 0),
                SIR = list(threshold = n_particles + 1),
                SIR_bis = list(threshold = n_particles / 2 ))
all_results <- lapply(methods,
                      function(my_method){
                        replicate(n_replicates,
                                  {
                                    get_SIR(target_model = target_model,
                                            q_proposal = q_proposal,
                                            y_obs =  sim_data$y,
                                            pars = full_pars_list,
                                            n_particles = n_particles,
                                            threshold = my_method$threshold)
                                  }, simplify = FALSE)
                      })

saveRDS(all_results, file = "results_particle_filter.rds")



