



get_ESS <- function(w){
  # W is a matrix with T rows (number of time steps) and N columns (number of MC samples)
  apply(w, 1, function(w) 1 / sum(w^2))
}


# One code for different method,
# SIS corresponds to a threshopld of 0 leading to no resampling
# threshold to n_particles + 1 leads to always resample.
get_SIR <- function(target_model, # A list with elements ..... 
                    q_proposal,  # A list with elements ....
                    y_obs, # A matrix of size n_obs * dim_y
                    pars, # A list with all model elements
                    n_particles = 100,
                    threshold = NULL){
  if(is.null(threshold)){ # Cas où on réechantillonne tout le temps
    threshold = n_particles + 1
  }
  # target_model is a list having at least two elements which are functions
  # - get_x0_knowing_y0
  # - get_xt_knowing_y_xtm1
  
  n_obs <- nrow(y_obs)
  dim_x <- pars$dim_x
  
  # On initialise
  particles <- array(NA, dim = c(n_particles, dim_x, n_obs))
  unnormed_log_filtering_weights <- filtering_weights <- matrix(0, nrow = n_obs, ncol = n_particles)
  
  particles[,,1] <- q_proposal$get_initial_samples(n_particles, pars)
  log_w0 <- log(target_model$get_x0_knowing_y0(particles[,,1], y_obs[1,], pars)) -
    log(q_proposal$get_initial_density(particles[,,1], pars))
  log_sum_unnormed_w <- max(log_w0) + log(sum(exp(log_w0 - max(log_w0)))) # log_sum_exp_trick
  unnormed_log_filtering_weights[1,] <- log_w0
  filtering_weights[1,] <- exp(unnormed_log_filtering_weights[1, ] - log_sum_unnormed_w)
  log_likelihood <- log_sum_unnormed_w - log(n_particles)
  old_log_sum_unnormed_w <- log_sum_unnormed_w
  for (i in 2:n_obs) {
    log_old_weights <- unnormed_log_filtering_weights[i - 1, ]
    current_ESS <-  1 / sum(filtering_weights[i - 1, ]^2)
    ancestor_indexes <- 1:n_particles
    if(current_ESS <= threshold){
      log_old_weights <- rep(0, n_particles) 
      ancestor_indexes <- sample(1:n_particles,
                                 size = n_particles,
                                 replace = TRUE,
                                 prob = filtering_weights[i - 1, ])
    }
    particles[,,i] <- q_proposal$get_transition_samples(ancestors = particles[ancestor_indexes, ,i - 1], 
                                                        y_obs[i,], 
                                                        pars)
    
    log_w <- log_old_weights + 
      log(target_model$get_xt_knowing_y_xtm1(x = particles[,,i], y = y_obs[i,],
                                             ancestors = particles[ancestor_indexes,,i - 1], pars)) -
      log(q_proposal$get_transition_density(x = particles[,,i], 
                                            y = y_obs[i,],
                                            ancestors = particles[ancestor_indexes,,i - 1], pars))
    log_sum_unnormed_w <- max(log_w) + log(sum(exp(log_w - max(log_w)))) # Log sum exp trick
    unnormed_log_filtering_weights[i, ] <- log_w
    filtering_weights[i,] <- exp(unnormed_log_filtering_weights[i, ] - log_sum_unnormed_w)
    log_likelihood <- log_likelihood + log_sum_unnormed_w - 
      ifelse(current_ESS <= threshold, 
             log(n_particles),
             old_log_sum_unnormed_w)
    print(log_sum_unnormed_w - 
            ifelse(current_ESS <= threshold, 
                   log(n_particles),
                   old_log_sum_unnormed_w))
    old_log_sum_unnormed_w <- log_sum_unnormed_w
  }
  
  return(list(particles = particles, filtering_weights =  filtering_weights,
              log_likelihood = log_likelihood))
}

