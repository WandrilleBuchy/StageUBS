

## ----création de la fonction de génération des données pour tout choix de dimensions----------------------------------------------------------------------

#Fonction racine carrée de matrice
get_mat_sqrt <- function(A){
  if(is.null(dim(A))){
    A <- matrix(A)
  }
  y = eigen(A)
  MCV <- matrix(y$vectors,byrow=F,ncol=dim(A)[1]) #matrice orthogonale de vecteurs propres
  MCVinv <- t(MCV) #transposée de matrice de vecteur propre
  Diag <- diag(y$values, nrow = dim(A)[1]) #matrice diagonale des valeurs propres
  DiagCalc <- sqrt(Diag) # matrice diagonale racine carrée des valeurs propres
  MCV %*% DiagCalc
}


#Fonction génération d'une donnée de n dimensions

get_initial_x_y <- function(pars){
  x <- pars$m0_x + get_mat_sqrt(pars$S_x) %*% rnorm(pars$dim_x)
  y <- pars$A_y %*% x +  get_mat_sqrt(pars$S_y) %*% rnorm(pars$dim_y)
  return(list(x = x, y = y))
}
get_next_x_y <- function(old_x, pars){
  x <- pars$F_x %*% old_x + get_mat_sqrt(pars$S_x) %*% rnorm(pars$dim_x)
  y <- pars$A_y %*% x +  get_mat_sqrt(pars$S_y) %*%  rnorm(pars$dim_y)
  return(list(x = x, y = y))
}


#Fonction de création des N données
get_data <- function(n_steps, pars){
  # Initialisation des objets
  x_states <- matrix(NA, nrow = n_steps, ncol = pars$dim_x)
  y_obs <- matrix(NA, nrow = n_steps, ncol = pars$dim_y)
  # Initialisation de la série temporelle
  xy_init <- get_initial_x_y(pars)
  x_states[1, ] <- xy_init$x
  y_obs[1, ] <- xy_init$y
  # Propagation
  for (i in 1:(n_steps-1)){ #Boucle de la création des données suivantes
    xy_next <- get_next_x_y(old_x = x_states[i], pars)
    x_states[i  + 1, ] <- xy_next$x
    y_obs[i + 1, ] <- xy_next$y
  }
  return(list(x = x_states, y = y_obs))
}

