source("D:/Users/buchy/Documents/Stage/utils_simulation_functions.R")
source("D:/Users/buchy/Documents/Stage/utils_kalman_functions.R")
source("D:/Users/buchy/Documents/Stage/code_filtrage_particulaire.R")

## ----SNIS-------------------------------------------------------------------------------------------------------------------------------------------------

get_SNIS <- function(f, u_target, q_proposal, observation, pars, m = 100){
  samples <- q_proposal$get_samples(m, pars)
  W <- u_target(samples, observation, pars)/ q_proposal$get_density(samples, pars)
  W_tilde <- W / sum(W)
  I_hat <- sum(W_tilde * f(samples))
  return(list(W_tilde = W_tilde, X = samples, I_hat = I_hat))
}


#Densite non normalisee

SNIS_X0 <- get_SNIS(f = function(x) x,
                    u_target,
                    observation = data_1d$y[1, ],
                    pars = full_pars_list,
                    m = 100)


# Affichage de la vraie loi

plot(objectif, filter_1d$mu_hat[1,] - 10, filter_1d$mu_hat[1,] + 10)
points(SNIS_X0$X, rep(0, 1000), cex = 100 * SNIS_X0$W_tilde)

df = data.frame(x = SNIS_X0$X, w =  SNIS_X0$W_tilde)

ggplot(df) +
  geom_point(aes(x = x, size = w), y = 0) +
  stat_function(fun = function(x) dnorm(x, 70, sqrt(0.5))) +
  stat_function(fun = function(x) dnorm(x, filter_1d$mu_hat[2],
                                        sqrt(filter_1d$V_hat[,,2])), color = "red")