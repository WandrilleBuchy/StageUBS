source("utils_simulation_functions.R")
source("utils_kalman_functions.R")
source("code_filtrage_particulaire.R")

# Loi cible, vraie loi de X_0 | Y_0
objectif <- function(x) return(dnorm(x, filter_1d$mu_hat[1,], sqrt(filter_1d$V_hat[1])))

# Approximations données par des particules (stockées dans SNIS_X0$X), 
# pondérées par leur poids SNIS_X0$W_tilde

# Affichage de la vraie loi
plot(objectif, filter_1d$mu_hat[1,] - 10, filter_1d$mu_hat[1,] + 10)

# Affichage de l'estimation de la densité 
lines(density(SNIS_X0$X, weights = SNIS_X0$W_tilde), col = "red")

# Vraie probabilité que  X_0 < 69 | Y_0
pnorm(69, filter_1d$mu_hat[1,], sqrt(filter_1d$V_hat[1]))

# Probabilité estimée que  X_0 < 69 | Y_0
sum(SNIS_X0$W_tilde * (SNIS_X0$X < 69))

# Estimation de E[X_0 | Y_0]
Esp_X0_Y0 <- sum(SNIS_X0$W_tilde * SNIS_X0$X)

# Estimation de V[X_0 | Y_0]
Var_X0_Y0 <- sum(SNIS_X0$W_tilde * SNIS_X0$X^2) - Esp_X0_Y0^2

# Vraie valeur:
filter_1d$V_hat[,,1]

# Estimation de l'intervalle de confiance à 95% pour X_0| Y_0 = y_0
Esp_X0_Y0 + c(-1, 1) * 1.96 * sqrt(Var_X0_Y0)

# verite
qnorm(c(0.025, 0.975),  filter_1d$mu_hat[1,], sqrt(filter_1d$V_hat[1]))


