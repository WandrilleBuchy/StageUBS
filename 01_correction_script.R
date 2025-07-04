
# Exo 3.2 -----------------------------------------------------------------


## Approche directe --------------------------------------------------------

m <- 1000
echantillons <- rnorm(m, mean = 0, sd = sqrt(2))
F_chapeau <- mean(2 * sqrt(pi) * echantillons^2 * cos(echantillons))
sigma2_F_chapeau <- var(2 * sqrt(pi) * echantillons^2 * cos(echantillons))
IC <- F_chapeau + c(-1, 1) * 1.96 * sqrt(sigma2_F_chapeau / m)

# Réponse à 3.3

phi_echantillons <- 2 * sqrt(pi) * echantillons^2 * cos(echantillons)
estimation_sequentielle <- cumsum(phi_echantillons) / (1:m)
plot(estimation_sequentielle)

# Cette manière de faire est bonne ponctuellement mais il vaut mieux
# coder grâce à des fonctions

## Approche programmation --------------------------------------------------

rm(list = ls()) # On supprime 

# On crée des fonctions

get_phi_F <- function(x){ # Fonction phi speciale pour F
  2 * sqrt(pi) * x^2 * cos(x)
} 

get_samples_F <- function(m){ # Fonction d'échantillonnage pour F
  rnorm(m, 0, sqrt(2))
}

# Fonction pour obtenir estimation Monte Carlo (MC)

get_MC_estimate <- function(m, # Effort Monte Carlo 
                            phi, # Fonction phi 
                            sampling_method){ # Manière d'obtenir échantillons
  samples <- sampling_method(m)
  phi_samples <- phi(samples)
  estimation <- mean(phi_samples)
  sigma2_est <- var(phi_samples)
  IC <- estimation + c(-1, 1) * 1.96 * sqrt(sigma2_est / m)
  # Estimations séquentielles
  sequence_estimates <- cumsum(phi_samples) / (1:m)
  sequence_IC_inf <- sequence_estimates - 1.96 * sqrt(sigma2_est / (1:m))
  sequence_IC_sup <- sequence_estimates + 1.96 * sqrt(sigma2_est / (1:m))
  results_table <- data.frame(Effort = 1:m,
                              Estimation = sequence_estimates,
                              IC_inf = sequence_IC_inf,
                              IC_sup = sequence_IC_sup)
  return(list(Estimation = estimation,
              Sigma2_chapeau = sigma2_est,
              IC = IC,
              Tableau = results_table))
}

resultats_F <- get_MC_estimate(m = 1000, phi = get_phi_F, 
                sampling_method = get_samples_F)

library(ggplot2)
ggplot(resultats_F$Tableau) + 
  aes(x = Effort, y = Estimation) + 
  geom_line() + # On dessine la séquence
  geom_ribbon(aes(ymin = IC_inf, ymax = IC_sup), # Rajout bande de confiance
              fill = "lightblue", alpha = .5) # alpha: transparence


# Exo 3.4 -----------------------------------------------------------------

# On définit les fonctions phi et la manière d'échantillonner
# et on procède avec la fonction de tout à l'heure

get_phi_pi <- function(X){ # X matrice m * 2 avec des échantillons
  a <- X[, 1]
  b <- X[, 2]
  4 * ((a^2 + b^2) <= 1)
}

get_samples_pi <- function(m){
  a <- runif(m, -1, 1) # Abscisse
  b <- runif(m, -1, 1) # Ordonnée
  X <- cbind(a, b) # Renvoie une matrice m * 2
}

# On utilise get-MC-estimate qu'on a défini plus haut

results_pi <- get_MC_estimate(m = 5e4, get_phi_pi, get_samples_pi)

ggplot(results_pi$Tableau) + 
  aes(x = Effort, y = Estimation) + 
  geom_line() + # On dessine la séquence
  geom_ribbon(aes(ymin = IC_inf, ymax = IC_sup), # Rajout bande de confiance
              fill = "lightblue", alpha = .5) +
  geom_hline(yintercept = pi, color = "red")# Ajout vraie valeur


