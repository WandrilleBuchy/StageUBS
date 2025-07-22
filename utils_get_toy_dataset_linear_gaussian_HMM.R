# Modele
# X_0 ~ N(m0_x, S_x)
# X_{t} ~ N(F_x, X_{t-1}, S_x) t >= 1
# Y_t ~ N(A_y X_{t}, S_y) t >=0

source("get_libraries.R")
source("utils_simulation_functions.R")

full_pars_list <- list(dim_x = 1,
                       dim_y = 1,
                       m0_x = as.numeric(70),
                       S_x = as.matrix(2), # Must be a matrix
                       F_x = as.matrix(0.9), 
                       A_y = as.matrix(1), # must be a matrix...
                       S_y = as.matrix(1))

full_pars_list_2d <- list(dim_x = 2,
                          dim_y = 2,
                          m0_x = as.vector(rep(70, 2)),
                          S_x = diag(2, 2), # Must be a matrix
                          F_x = diag(0.9,2), 
                          A_y = diag(1,2), # must be a matrix...
                          S_y = diag(1,2))


n_obs <- 2000
set.seed(1234)
sim_data <- get_data(n_steps = n_obs, full_pars_list) #Données 1d
sim_data_2d <- get_data(n_steps = n_obs, full_pars_list_2d) #Données 2d
n_particles = 50
n_replicates = 100 # Nombre de réplicats des estimateurs