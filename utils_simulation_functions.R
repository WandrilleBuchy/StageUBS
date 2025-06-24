

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
get_simulation <- function(n_sim, pars){
  x_white_noise <- matrix(rnorm(n_sim * pars$dim_x),
                          nrow = pars$dim_x,
                          ncol = n_sim)
  x <- pars$matF %*% pars$m_x + get_mat_sqrt(pars$S_x) %*% x_white_noise
  y_white_noise <- matrix(rnorm(n_sim * pars$dim_y),
                          nrow = pars$dim_y,
                          ncol = n_sim)
  y <- pars$A %*% x + pars$b +  get_mat_sqrt(pars$S_y) %*% y_white_noise
  return(list(x = t(x), y = t(y)))
}


#Fonction de création des N données
get_data <- function(n_sim, pars_X0, pars){
  
  Data <- get_simulation(1, pars_X0) #Création de la première donnée
  X_loop <- Data$x #Récupération de la donnée pour l'utiliser en argument de la suivante
  pars$m_x <- t(X_loop) #Ajout dans les paramètre de la donnée récupérée
  
  
  for (i in 1:(n_sim-1)){ #Boucle de la création des données suivantes
    
    append <- get_simulation(1, pars) #Création de la donnée en fonction de la précédente
    X_loop <- append$x #Récupération de la donnée pour l'utiliser en argument de la suivante
    pars$m_x <- as.matrix(t(X_loop)) #Ajout dans les paramètre de la donnée récupérée
    Data <- bind_rows(Data,append) #Enregistrement des valeurs de X et Y dans le Data
  }
  return(Data)
}

