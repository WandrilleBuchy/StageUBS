source("utils_kalman_functions.R")
source("utils_particle_filter_functions.R")

library(matrixStats)
  
  


estimate_loglikelihood_y0n <- function(particle, particle_weight, y_obs){
  #particle and particle_weight are n_obs * n_particles matrixes
  #y_obs is a n_obs vector
  
  estimate_yk_knowing_y0k_prev <- rowSums(true_model$get_emission_density(x = particle,
                                                                         y = y_obs,
                                                                         pars = full_pars_list) * particle_weight)
  loglikelihood <- logSumExp(estimate_yk_knowing_y0k_prev)
  return(loglikelihood)
}

SIR_loglikelihood <- estimate_loglikelihood_y0n(particle = all_results$SIR[[1]]$particles[,1,],
                                                particle_weight = all_results$SIR[[1]]$filtering_weights, 
                                                y_obs = sim_data$y)


true_likelihood <- function(){
  
}