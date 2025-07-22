## -----------------------------------------------1d representations

#représentation 1d
df <- data.frame(
  x = as.vector(sim_data$x),
  y = as.vector(sim_data$y)
)

df$index <- 1:nrow(df)

ggplot(df, aes(x = index, y = y)) +
  geom_line() +
  geom_point(color = 'red', size = .5, alpha = .8) +
  labs(x = "Instants temporels", y = "Données observées", title = "Courbe des données observées en fonction des différents instants temporels t")




data_histo <- rep(1,1000)
for (i in 1:1000){
  data_histo[i] <- get_data(2, full_pars_list)$y[1,]
}

#fonction gaussienne de la répartition des tirages Y_0
Y_densite <- function(x, pars){
  return(dnorm(x, mean = pars$A_y * pars$m0_x,
               sd = sqrt(pars$A_y^2 * pars$S_x + pars$S_y)))
}


#représentation graphique de l'échantillon créé
ggplot(data.frame(data_histo = data_histo),
       aes(x = data_histo)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = .5,
                 fill = "lightblue3",
                 col = 'black') +
  ggtitle("Histogramme des observations initiales effectivement observées") +
  
  #densité réelle
  geom_density(aes(color = "Densité réelle du tirage de Y initial"),
               linewidth = 1) +
  
  #densité théorique
  stat_function(fun = Y_densite,
                args = list(pars = full_pars_list),
                aes(color = "Densité théorique du tirage de Y initial"), 
                linewidth = 1) +
  
  #création de la légende
  scale_color_manual(name = "Légende",
                     values = c("Densité théorique du tirage de Y initial" = "red", 
                                "Densité réelle du tirage de Y initial" = "blue")) +
  
  labs(x = "Valeur des observations Y initiales",
       y = "Densité du comptage des observations Y initiales")


## ----------------------------------------------- 2d representations

#Séparation des données en intervalles
x_c <- cut(sim_data_2d$y[,1], 70)
y_c <- cut(sim_data_2d$y[,2], 70)

#Jonction des données
z <- table(x_c, y_c)

#Représentation de l'histogramme 3D
hist3D(z=z, border="black",
       main = "Histogramme 3D des données effectivement observées",
       xlab = "Valeur dim 1 de Y",
       ylab = "Valeur dim 2 de Y",
       zlab = "Comptage")

#Représentation de la heatmap
image2D(z=z, border="black",main = "Heatmap des données effectivement observées",
        xlab = "Valeur dim 1 de Y",
        ylab = "Valeur dim 2 de Y")


#représentation du graphe fleché
df_x <- tibble(x = sim_data_2d$x[,1], y = sim_data_2d$x[,2])
df_y <- tibble(x = sim_data_2d$y[,1], y = sim_data_2d$y[,2])

df_arrows_x <- tibble(
  x_start = df_x$x[-nrow(df_x)],
  y_start = df_x$y[-nrow(df_x)],
  x_end   = df_x$x[-1],
  y_end   = df_x$y[-1]
)

df_arrows_y <- tibble(
  x_start = df_y$x[-nrow(df_y)],
  y_start = df_y$y[-nrow(df_y)],
  x_end   = df_y$x[-1],
  y_end   = df_y$y[-1]
)

df_points <- bind_rows(
  mutate(df_x, Legende = "Points de X"),
  mutate(df_y, Legende = "Points de Y")
)

ggplot() +
  geom_segment(data = df_arrows_y,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               color = "lightblue4", show.legend = FALSE) +
  geom_segment(data = df_arrows_x,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               color = "pink3", show.legend = FALSE) +
  geom_point(data = df_points,
             aes(x = x, y = y, color = Legende),
             size = 0.4) +
  scale_color_manual(values = c("Points de X" = "red", "Points de Y" = "blue")) +
  labs(x = "Coordonnée abscisse", y = "Coordonnée ordonnée", title = "Trajectoire des points Y₁ → Y₂ → ... → Y₂₀₀₀ et X₁ → X₂ → ... → X₂₀₀₀") 
