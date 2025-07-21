source("utils_kalman_functions.R")
source("utils_particle_filter_functions.R")
source("utils_get_toy_dataset_linear_gaussian_HMM.R")



library(tidyverse)
library(readr)
all_results <- read_rds("results_particle_filter.rds")


estimate_kalman_loglikelihood_y0n <- function(sim_data, kalman_result, full_pars_list){
  #ancestors and Pred_kalm are n_obs vectors
  
  log_p_yk_knowing_ykm1 <- sapply(2:n_obs, function(i){
    mixtools::logdmvnorm(y = sim_data$y[i, ],
      mu = full_pars_list$A_y %*% full_pars_list$F_x %*% kalman_result$mu_hat[i-1],
      sigma = full_pars_list$A_y %*% (full_pars_list$F_x %*% 
                                        kalman_result$V_hat[,,i-1] %*% t(full_pars_list$F_x) +
                                        full_pars_list$S_x) %*% t(full_pars_list$A_y) + full_pars_list$S_y)})
  
  log_p_y1 <- mixtools::logdmvnorm(y = sim_data$y[1,],
                                   mu = full_pars_list$A_y %*% full_pars_list$m0_x, 
                                   sigma = full_pars_list$S_y + 
                                     full_pars_list$A_y %*% full_pars_list$S_x %*% t(full_pars_list$A_y))
  
  log_likelihood <- sum(log_p_yk_knowing_ykm1) + log_p_y1  
  
  return(list(log_likelihood = log_likelihood,
              log_p_yk_knowing_ykm1 = c(log_p_y1,log_p_yk_knowing_ykm1)))
}
kalman_result <- get_filtering_kalman(pars = full_pars_list, y_obs = sim_data$y)


kalman_likelihood <- estimate_kalman_loglikelihood_y0n(sim_data = sim_data,
                                                       kalman_result = kalman_result,
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

bind_rows(data.frame(method = "SIS", ll = SIS_log_likelihoods),
          data.frame(method = "SIR", ll = SIR_log_likelihoods),
          data.frame(method = "SIR_adaptative", ll = SIR_bis_log_likelihoods)) |> 
  mutate(likelihood = exp(ll)) |> 
  group_by(method) |> 
  mutate(mean_likelihood = mean(likelihood)) |> 
  ungroup() |> 
  # filter(method != "SIS") |> 
  ggplot(aes(x = method, y = likelihood)) +
  geom_boxplot() + 
  geom_point(aes(y = mean_likelihood), color = "red") +
  geom_hline(yintercept = exp(kalman_likelihood$log_likelihood), color = "blue") +
  labs(x = "Méthodes utilisées", y = "Vraisemblance")

methods_mean_log_likelihood <- list(SIS_mean_log_likelihood = mean(SIS_log_likelihoods),
                                    SIR_mean_log_likelihood = mean(SIR_log_likelihoods),
                                    SIR_bis_mean_log_likelihood = mean(SIR_bis_log_likelihoods))

