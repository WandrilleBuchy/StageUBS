
## kalman filtering function
get_filtering_kalman <- function(pars, y_obs){
  
  #initialisation du tableau de sortie
  n <- nrow(y_obs)
  
  mu_hat = matrix(nrow = n, ncol = pars$dim_x)
  V_hat = P_hat = array(dim = c(pars$dim_x, pars$dim_x, n))
  
  #initialisation de l'algorithme pour le premier terme
  K <- pars$S_x %*% t(pars$A_y) %*% solve(pars$A_y %*% pars$S_x %*% t(pars$A_y) + pars$S_y)
  mu <- pars$m0_x + K %*% (y_obs[1,] - pars$A_y %*% pars$m0_x)
  V <- (diag(pars$dim_x) - K %*% pars$A_y) %*% pars$S_x
  Pred_kalm <- pars$F_x %*% V %*% t(pars$F_x) + pars$S_x
  
  mu_hat[1,] <- as.numeric(mu)
  V_hat[,,1] <- V
  P_hat[,,1] <- Pred_kalm
  
  #boucle de l'algorithme
  for (i in 2:n){
    K <- Pred_kalm %*% t(pars$A_y) %*% solve(pars$A_y %*% Pred_kalm %*% t(pars$A_y) + pars$S_y)
    mu <- pars$F_x %*% mu + K %*% (y_obs[i,] - pars$A_y %*% pars$F_x %*% mu)
    V <- (diag(pars$dim_x) - K %*% pars$A_y) %*% Pred_kalm
    Pred_kalm <- pars$F_x %*% V %*% t(pars$F_x) + pars$S_x
    
    mu_hat[i,] <- as.numeric(mu)
    V_hat[,,i] <- V
    P_hat[,,i] <- Pred_kalm
    
  }
  
  colnames(mu_hat) <- paste0("M", 1:pars$dim_x)
  means <- as.data.frame(mu_hat)
  
  marginal_variances <- matrix(NA_real_, nrow = n, ncol = pars$dim_x)
  for (i in 1:n) {
    marginal_variances[i, ] <- diag(matrix(V_hat[,,i], nrow = pars$dim_x))
  }
  
  
  colnames(marginal_variances) <- paste0("V", 1:pars$dim_x)
  marginal_variances <- as.data.frame(marginal_variances)
  
  formatted_data <- bind_cols(means, marginal_variances)  %>%
    mutate(index = 1:nrow(.))  %>% 
    pivot_longer(
      cols = -index,
      names_to = c(".value", "dimension"),
      names_pattern = "([MV])(\\d+)"
    )  %>% 
    rename(mean = M, Var = V)  %>%
    mutate(lower = mean - 1.96 * sqrt(Var),
           upper = mean + 1.96 * sqrt(Var))
  
  return(list(mu_hat = mu_hat, V_hat = V_hat, result_df = formatted_data, Pred_kalm = P_hat, index = 1:n))
}



#RTS smoothering function
get_smoothering_RTS <- function(ini, pars){
  
  #initialisation du tableau de sortie
  n <- nrow(ini$mu_hat)
  mu_hat = matrix(nrow = n,
                  ncol = pars$dim_x)
  V_hat = array(dim = c(pars$dim_x,
                        pars$dim_x,
                        n))
  
  mu_hat[n,] <- t(ini$mu_hat[n,])
  V_hat[,,n] <- ini$V_hat[,,n]
  
  #boucle de l'algorithme
  for (i in (n-1):1){
    
    C <- ini$V_hat[,,i] %*% t(pars$F_x) %*%
      solve(ini$Pred_kalm[,,i])
    mu_hat_loop <- ini$mu_hat[i,] + C %*%
      (mu_hat[i+1,] - pars$F_x %*% ini$mu_hat[i,])
    V_hat_loop <- ini$V_hat[,,i] + C %*%
      (V_hat[,,i+1] - ini$Pred_kalm[,,i]) %*% t(C)
    
    mu_hat[i,] <- as.matrix(mu_hat_loop)
    V_hat[,,i] <- V_hat_loop
    
  }
  marginal_variances <- matrix(NA_real_,
                               nrow = n,
                               ncol = pars$dim_x)
  for (i in 1:n) {
    marginal_variances[i, ] <- diag(matrix(V_hat[,,i],
                                           nrow = pars$dim_x))
  }
  
  colnames(mu_hat) <- paste0("M", 1:pars$dim_x)
  colnames(marginal_variances) <- paste0("V", 1:pars$dim_x)
  
  means <- as.data.frame(mu_hat)
  marginal_variances <- as.data.frame(marginal_variances)
  
  formatted_data <- bind_cols(means, marginal_variances)%>% 
    mutate(index = 1:nrow(.))%>% 
    pivot_longer(
      cols = -index,
      names_to = c(".value", "dimension"),
      names_pattern = "([MV])(\\d+)"
    ) %>%
    rename(mean = M, Var = V) %>%
    mutate(lower = mean - 1.96 * sqrt(Var),
           upper = mean + 1.96 * sqrt(Var))
  return(list(mu_hat = mu_hat,
              V_hat = V_hat,
              result_df = formatted_data))
}
