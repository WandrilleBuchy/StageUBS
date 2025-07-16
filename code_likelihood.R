source("utils_kalman_functions.R")
source("utils_particle_filter_functions.R")

#all_result <- read_rds("result_particle_filter.rds")

library(matrixStats)


estimate_kalman_loglikelihood_y0n <- function(kalman_result, full_pars_list){
  #ancestors and Pred_kalm are n_obs vectors
  
  log_p_yk_knowing_ykm1 <- sapply(2:n_obs, function(i){
    mixtools::logdmvnorm(full_pars_list$A_y %*% full_pars_list$F_x %*% kalman_result$mu_hat[i-1],
                         full_pars_list$A_y %*% (full_pars_list$F_x %*% kalman_result$Pred_kalm[i] 
                                %*% t(full_pars_list$F_x) %*% t(full_pars_list$A_y) + full_pars_list$S_x) + full_pars_list$S_y)})
  
  log_p_y1 <- mixtools::logdmvnorm(full_pars_list$A_y %*% full_pars_list$m0_x, 
                                   full_pars_list$S_y + 
                                     full_pars_list$A_y %*% full_pars_list$S_x %*% t(full_pars_list$A_y))
  
  log_likelihood <- sum(log_p_yk_knowing_ykm1) + log_p_y1  
  
  return(list(log_likelihood = log_likelihood,
              log_p_yk_knowing_ykm1 = c(log_p_y1,log_p_yk_knowing_ykm1)))
}


kalman_likelihood <- estimate_kalman_loglikelihood_y0n(kalman_result = kalman_result,
                                                       full_pars_list = full_pars_list)

plot(kalman_likelihood$log_p_yk_knowing_ykm1, type = "l")

SIS_log_likelihoods <- sapply(1:n_replicates, function(i) {
  all_results$SIS[[i]]$log_likelihood
})
SIR_log_likelihoods <- sapply(1:n_replicates, function(i) {
  all_results$SIR[[i]]$log_likelihood
})
SIR_bis_log_likelihoods <- sapply(1:n_replicates, function(i) {
  all_results$SIR_bis[[i]]$log_likelihood
})


methods_mean_log_likelihood <- list(SIS_mean_log_likelihood = mean(SIS_log_likelihoods),
                                    SIR_mean_log_likelihood = mean(SIR_log_likelihoods),
                                    SIR_bis_mean_log_likelihood = mean(SIR_bis_log_likelihoods))