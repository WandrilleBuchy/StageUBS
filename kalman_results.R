source("get_libraries.R")
source("utils_simulation_functions.R")
source("utils_kalman_functions.R")

filter_1d <- get_filtering_kalman(pars = full_pars_list,
                                  y_obs = sim_data$y)
filter_2d <- get_filtering_kalman(pars = full_pars_list_2d,
                                  y_obs = sim_data_2d$y)

smoother_1d <- get_smoothering_RTS(ini = filter_1d,
                                   pars = full_pars_list)
smoother_2d <- get_smoothering_RTS(ini = filter_2d,
                                   pars = full_pars_list_2d)

#différence des estimations
estimation_variation <- abs(smoother_1d$mu_hat - filter_1d$mu_hat)
print(mean(estimation_variation))