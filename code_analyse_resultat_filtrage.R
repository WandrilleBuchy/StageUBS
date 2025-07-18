
library(tidyverse)

source("utils_get_toy_dataset_linear_gaussian_HMM.R")
source("utils_kalman_functions.R")

# Loading particle filter results
if(!file.exists("results_particle_filter.rds")){
  source("code_filtrage_particulaire.R")
}

all_results <- readRDS("results_particle_filter.rds")
kalman_result <- get_filtering_kalman(pars = full_pars_list, y_obs = sim_data$y)
kalman_df <- data.frame(t = 1:nrow(sim_data$y),
                        x = kalman_result$mu_hat[,1],
                        IC_low = kalman_result$mu_hat[,1] - 1.96 * sqrt(kalman_result$V_hat[1, 1, ]),
                        IC_sup = kalman_result$mu_hat[,1] + 1.96 * sqrt(kalman_result$V_hat[1, 1, ]))

# Get estimates of E[X_t | Y_{0:t}] ---------------------------------------

all_estimates_E_Xt <- imap_dfr(all_results,
                               function(x, nm){
                                 map2_dfr(.x = x, .y = 1:n_replicates,
                                          .f = function(res_PF, idx){
                                            # On chope une matrice de taille n_obs * n_particles
                                            X <- t(res_PF$particles[, 1, ]) # Code valable en dimension 1, sinon utiliser apply
                                            est_X <- rowSums(res_PF$filtering_weights * X)
                                            data.frame(t = 0:(nrow(X) - 1),
                                                       x = est_X,
                                                       replicate = idx,
                                                       method = nm)
                                          })
                               })

# Cheching variance of estimators -----------------------------------------


all_estimates_E_Xt |> 
  group_by(method, t) |> 
  summarise(V_est = var(x)) |>
  ggplot(aes(x = t, y = V_est, color = method)) + 
  lims(x = c(0, 2), y = c(0, .1)) +
  geom_line() 

# Comparing trajectories

all_estimates_E_Xt |> 
  ggplot(aes(x = t, y = x)) +
  geom_line(aes(color = method, group = interaction(method, replicate))) +
  geom_line(data = sim_data$x |> 
              as.data.frame() |> 
              rename(x = V1) |> 
              mutate(t = 0:(nrow(sim_data$x) - 1)),
            aes(color = "True X")) +
  geom_line(data = kalman_df, aes(color = "Kalman"))
  
