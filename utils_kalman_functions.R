get_filtering_kalman <- function(pars_X0, pars, data){
  
  #initialisation du tableau de sortie
  n <- nrow(data)
  
  mu_hat = matrix(nrow = n, ncol = pars$dim_x)
  V_hat = P_hat = array(dim = c(pars$dim_x, pars$dim_x, n))
  
  #initialisation de l'algorithme pour le premier terme
  K <- pars_X0$S_x %*% t(pars$A) %*% solve(pars$A %*% pars_X0$S_x %*% t(pars$A) + pars$S_y)
  mu <- pars_X0$m_x + K %*% (data$y[1,] - pars$A %*% pars_X0$m_x)
  V <- (diag(pars$dim_x) - K %*% pars$A) %*% pars_X0$S_x
  Pred_kalm <- pars$matF %*% V %*% t(pars$matF) + pars$S_x
  
  mu_hat[1,] <- as.numeric(mu)
  V_hat[,,1] <- V
  P_hat[,,1] <- Pred_kalm
  
  #boucle de l'algorithme
  for (i in 2:n){
    K <- Pred_kalm %*% t(pars$A) %*% solve(pars$A %*% Pred_kalm %*% t(pars$A) + pars$S_y)
    mu <- pars$matF %*% mu + K %*% (data$y[i,] - pars$A %*% pars$matF %*% mu)
    V <- (diag(pars$dim_x) - K %*% pars$A) %*% Pred_kalm
    Pred_kalm <- pars$matF %*% V %*% t(pars$matF) + pars$S_x
    
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
  
  formatted_data <- bind_cols(means, marginal_variances) %>%
    mutate(index = 1:nrow(.)) %>%
    pivot_longer(
      cols = -index,
      names_to = c(".value", "dimension"),
      names_pattern = "([MV])(\\d+)"
    ) %>%
    rename(mean = M, Var = V) %>%
    mutate(lower = mean - 1.96 * sqrt(Var),
           upper = mean + 1.96 * sqrt(Var))
  
  return(list(mu_hat = mu_hat, V_hat = V_hat, result_df = formatted_data, Pred_kalm = P_hat, index = 1:n))
}