



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
  if(is.null(threshold)){ # Always resample
    threshold = n_particles + 1
  }
  # target_model is a list having at least two elements which are functions
  # - get_x0_knowing_y0
  # - get_xt_knowing_y_xtm1
  
  # initialization
  n_obs <- nrow(y_obs)
  dim_x <- pars$dim_x

  particles <- array(NA, dim = c(n_particles, dim_x, n_obs)) # 3 dim matrix
  unnormed_log_filtering_weights <- filtering_weights <- matrix(0, nrow = n_obs, ncol = n_particles) # n_obs * n_particles matrices
  
  particles[,,1] <- q_proposal$get_initial_samples(n_particles, pars) # get n_particles initials samples cf n_particles x_0
  
  log_w0 <- log(target_model$get_x0_knowing_y0(particles[,,1], y_obs[1,], pars)) -
    log(q_proposal$get_initial_density(particles[,,1], pars)) # get the initials log weights for the n_particles x_0 ------- the log is needed to process the log likelihood
  
  log_sum_unnormed_w <- max(log_w0) + log(sum(exp(log_w0 - max(log_w0)))) # log_sum_exp_trick to prevent numerical explosions cf board pictures
  unnormed_log_filtering_weights[1,] <- log_w0 # saving the log initials weights
  filtering_weights[1,] <- exp(unnormed_log_filtering_weights[1, ] - log_sum_unnormed_w) # get by the exp the real normed weights from the log weights
  log_likelihood <- log_sum_unnormed_w - log(n_particles) # get the log likelihood
  old_log_sum_unnormed_w <- log_sum_unnormed_w # save the old log weights sums
  
  for (i in 2:n_obs) { # start of the loop
    log_old_weights <- unnormed_log_filtering_weights[i - 1, ] # keep the old unnormed weights
    current_ESS <-  1 / sum(filtering_weights[i - 1, ]^2) # get ESS in case of SIR or SIR_bis
    ancestor_indexes <- 1:n_particles # creation of indexes in case of SIR or SIR_bis
    
    if(current_ESS <= threshold){ # compare threshold and ESS -------- threshold = 0 -> SIS | threshold = n_particles + 1 or NULL -> SIR | 0 < threshold < n_particles + 1 -> SIR_bis
      log_old_weights <- rep(0, n_particles) # reinitialization of a n_particles vector
      ancestor_indexes <- sample(1:n_particles,
                                 size = n_particles,
                                 replace = TRUE,
                                 prob = filtering_weights[i - 1, ]) # sample n_particles from the index ancestor_indexes with probability filtering_weights[i-1,] for each number from the index
    }
    particles[,,i] <- q_proposal$get_transition_samples(ancestors = particles[ancestor_indexes, ,i - 1], # no change in the case of SIS or SIR_bis with unrealized condition
                                                        y_obs[i,], 
                                                        pars) #creation of the next generation of points after the resampling or not  
    
    log_w <- log_old_weights + 
      log(target_model$get_xt_knowing_y_xtm1(x = particles[,,i], y = y_obs[i,],
                                             ancestors = particles[ancestor_indexes,,i - 1], pars)) -
      log(q_proposal$get_transition_density(x = particles[,,i], 
                                            y = y_obs[i,],
                                            ancestors = particles[ancestor_indexes,,i - 1], pars)) # updating the log weights
    
    log_sum_unnormed_w <- max(log_w) + 
      log(sum(exp(log_w - max(log_w)))) # Log sum exp trick to prevent numerical explosions cf board pictures
    unnormed_log_filtering_weights[i, ] <- log_w
    filtering_weights[i,] <- exp(unnormed_log_filtering_weights[i, ] - log_sum_unnormed_w) # get by the exp the normed real weights from the log weights
    
    # get the log likelihood
    log_likelihood <- log_likelihood + log_sum_unnormed_w - 
      ifelse(current_ESS <= threshold, 
             log(n_particles), # if condition not realized ie no resampling, the weights are the same -> classic MC estimator
             old_log_sum_unnormed_w) # if condition realized ie resampling, the weights are not the same -> estimator SIR with weights' sum
    
    old_log_sum_unnormed_w <- log_sum_unnormed_w # saving the old weights for the next log likelihood
  }
  
  return(list(particles = particles, filtering_weights =  filtering_weights,
              log_likelihood = log_likelihood))
}

