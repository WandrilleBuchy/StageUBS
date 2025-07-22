## fonction representation graphs

plot_dimension <- function(result_df, data, filter_result_df = NULL, dim_to_plot = 1) {
  #Filtrer sur la dimension choisie
  df_dim <- result_df %>% filter(dimension == as.character(dim_to_plot))
  
  #Extraction l'état réel et l'observation correspondants
  true_state <- data.frame(index = 1:nrow(data$x), 
                           value = data$x[, dim_to_plot], 
                           type = paste0("x (état réel ", dim_to_plot, ")"))
  
  observations <- data.frame(index = 1:nrow(data$y), 
                             value = data$y[, dim_to_plot], 
                             type = paste0("y (observation ", dim_to_plot, ")"))
  
  #Récupération du filtrage
  filter_df <- NULL
  if (!is.null(filter_result_df)) {
    filter_df <- filter_result_df %>% filter(dimension == as.character(dim_to_plot)) %>%
      mutate(type = paste0("mu_hat (filtrage ", dim_to_plot, ")")) %>%
      select(index, mean, type)
  }
  
  #Récupération du lissage
  if (!is.null(filter_result_df)) {
    smoother_df <- df_dim %>% 
      select(index, mean, lower, upper) %>%
      mutate(type = paste0("mu_hat (lissage ", dim_to_plot, ")"))
  }
  
  if (is.null(filter_result_df)) {
    smoother_df <- df_dim %>% 
      select(index, mean, lower, upper) %>%
      mutate(type = paste0("mu_hat (filtrage ", dim_to_plot, ")"))
  }
  
  #Données pour ribbon IC
  ribbon_df <- df_dim %>% 
    select(index, lower, upper) %>%
    mutate(type = "Intervalle de confiance")
  
  #Préparation des données pour les points
  points_df <- bind_rows(
    true_state,
    observations,
    if (!is.null(filter_df)) filter_df %>% rename(value = mean),
    smoother_df %>% rename(value = mean)
  ) %>%
    mutate(type = factor(type, levels = c(
      paste0("x (état réel ", dim_to_plot, ")"),
      paste0("y (observation ", dim_to_plot, ")"),
      if (!is.null(filter_df)) paste0("mu_hat (filtrage ", dim_to_plot, ")"),
      if (!is.null(filter_df)) paste0("mu_hat (lissage ", dim_to_plot, ")"),
      if (is.null(filter_df)) paste0("mu_hat (filtrage ", dim_to_plot, ")")
    )))
  
  
  # Plot
  colors <- c("black", "blue", "red", "darkgreen", "cyan4")
  names(colors) <- c(
    paste0("x (état réel ", dim_to_plot, ")"),
    paste0("y (observation ", dim_to_plot, ")"),
    paste0("mu_hat (filtrage ", dim_to_plot, ")"),
    paste0("mu_hat (lissage ", dim_to_plot, ")"),
    "Intervalle de confiance"
  )
  
  shapes <- c(16, 4, 1, 1, NA)
  names(shapes) <- c(
    paste0("x (état réel ", dim_to_plot, ")"),
    paste0("y (observation ", dim_to_plot, ")"),
    paste0("mu_hat (filtrage ", dim_to_plot, ")"),
    paste0("mu_hat (lissage ", dim_to_plot, ")"),
    "Intervalle de confiance"
  )
  
  fills <- c("cyan4")
  names(fills) <- c("Intervalle de confiance")
  
  plot <- ggplot() +
    geom_ribbon(data = ribbon_df, aes(x = index, ymin = lower, ymax = upper, fill = type), alpha = 0.2) +
    geom_line(data = smoother_df, aes(x = index, y = mean, color = type), linewidth = 1, alpha = 0.6) +
    geom_point(data = points_df, aes(x = index, y = value, color = type, shape = type), size = 2, na.rm = TRUE) +
    scale_color_manual(name = "Légende", values = colors) +
    scale_fill_manual(name = "Légende", values = fills) +
    scale_shape_manual(name = "Légende", values = shapes) +
    labs(x = "Numéro de l'observation", y = "Valeur")
  
  if (!is.null(filter_result_df)) {
    plot <- plot + geom_line(data = filter_df, aes(x = index, y = mean, color = type), linewidth = 1) + 
      labs(title = paste0("Estimation et lissage de la dim ", dim_to_plot, " de x avec intervalle de confiance à 95%"))
  }
  
  if (is.null(filter_result_df)) {
    plot <- plot + labs(title = paste0("Estimation de la dim ", dim_to_plot, " de x avec intervalle de confiance à 95%"))
  }
  
  print(plot)
}
