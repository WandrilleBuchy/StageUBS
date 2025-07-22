source("get_libraries.R")
source("kalman_results.R")
source("utils_representation_graph.R")

#representation du filtre de kalman en 1d
plot_dimension(result_df = filter_1d$result_df,
               data = sim_data)

#representation du filtre de kalman en 2d
plot_dimension(result_df = filter_2d$result_df,
               data = sim_data_2d,
               dim_to_plot = 1)

plot_dimension(result_df = filter_2d$result_df,
               data = sim_data_2d,
               dim_to_plot = 2)

#representation du lisseur RTS en 1d
plot_dimension(result_df = smoother_1d$result_df,
               data = sim_data,
               filter_result_df = filter_1d$result_df)

#representation du lisseur RTS en 2d
plot_dimension(result_df = smoother_2d$result_df,
               data = sim_data_2d,
               filter_result_df = filter_2d$result_df,
               dim_to_plot = 1)

plot_dimension(result_df = smoother_2d$result_df,
               data = sim_data_2d,
               filter_result_df = filter_2d$result_df,
               dim_to_plot = 2)


# representaion ellipses de confiance

q <- qchisq(.95, full_pars_list_2d$dim_x)


point_filter_16 <- as.data.frame(t(filter_2d$mu_hat[16, ]))
colnames(point_filter_16) <- c("x", "y")

point_filter_500 <- as.data.frame(t(filter_2d$mu_hat[500, ]))
colnames(point_filter_500) <- c("x", "y")

point_smoother_16 <- as.data.frame(t(smoother_2d$mu_hat[16, ]))
colnames(point_smoother_16) <- c("x", "y")

point_smoother_500 <- as.data.frame(t(smoother_2d$mu_hat[500, ]))
colnames(point_smoother_500) <- c("x", "y")

mu_filter_16 <- as.numeric(filter_2d$mu_hat[16, ])  
sigma_filter_16 <- filter_2d$V_hat[,,16]

mu_filter_500 <- as.numeric(filter_2d$mu_hat[500, ])  
sigma_filter_500 <- filter_2d$V_hat[,,500]

mu_smoother_16 <- as.numeric(smoother_2d$mu_hat[16, ])  
sigma_smoother_16 <- smoother_2d$V_hat[,,16]

mu_smoother_500 <- as.numeric(smoother_2d$mu_hat[500, ])  
sigma_smoother_500 <- smoother_2d$V_hat[,,500]

x_seq_16 <- seq(filter_2d$mu_hat[16,1] - 10, filter_2d$mu_hat[16,1] + 10, by = .1)
y_seq_16 <- seq(filter_2d$mu_hat[16,2] - 10, filter_2d$mu_hat[16,2] + 10, by = .1)

x_seq_500 <- seq(filter_2d$mu_hat[500,1] - 10, filter_2d$mu_hat[500,1] + 10, by = .1)
y_seq_500 <- seq(filter_2d$mu_hat[500,2] - 10, filter_2d$mu_hat[500,2] + 10, by = .1)

grille_points1 <- expand.grid(x = x_seq_16, y = y_seq_16)
grille_points2 <- expand.grid(x = x_seq_500, y = y_seq_500)

dist_maha_filter_16 <- mapply(get_mahalanobis, 
                              x = grille_points1$x, 
                              y = grille_points1$y,
                              MoreArgs = list(mu = mu_filter_16, sigma = sigma_filter_16))

dist_maha_filter_500 <- mapply(get_mahalanobis, 
                               x = grille_points2$x, 
                               y = grille_points2$y, 
                               MoreArgs = list(mu = mu_filter_500, sigma = sigma_filter_500))

dist_maha_smoother_16 <- mapply(get_mahalanobis, 
                                x = grille_points1$x, 
                                y = grille_points1$y,
                                MoreArgs = list(mu = mu_smoother_16, sigma = sigma_smoother_16))

dist_maha_smoother_500 <- mapply(get_mahalanobis, 
                                 x = grille_points2$x, 
                                 y = grille_points2$y, 
                                 MoreArgs = list(mu = mu_smoother_500, sigma = sigma_smoother_500))

area_filter_16 <- grille_points1[dist_maha_filter_16 <= q, ]

area_filter_500 <- grille_points2[dist_maha_filter_500 <= q, ]

area_smoother_16 <- grille_points1[dist_maha_smoother_16 <= q, ]

area_smoother_500 <- grille_points2[dist_maha_smoother_500 <= q, ]


ggplot()+
  geom_point(data = filter_2d$mu_hat,
             aes(x = filter_2d$mu_hat[,1], y = filter_2d$mu_hat[,2]),
             size = 0.4)+
  geom_point(data = point_filter_16,
             aes(x = x, y = y),
             size = 2, col = 'red') +
  geom_point(data = point_filter_500,
             aes(x = x, y = y),
             size = 2, col = 'blue') +
  geom_tile(data = area_filter_16, aes(x = x, y = y), fill = "red", alpha = .4)+
  geom_tile(data = area_filter_500, aes(x = x, y = y), fill = "blue", alpha = .4)+
  geom_tile(data = grille_points1, aes(x = x, y = y), fill = "lightblue", alpha = 0.5)+
  geom_tile(data = grille_points2, aes(x = x, y = y), fill = "lightblue", alpha = 0.5)+
  labs(x = "Dimension 1 du filtrage", y = "Dimension 2 du filtrage")+
  ggtitle("Représentation des points estimés par la méthode du filtrage")



ggplot()+
  geom_point(data = smoother_2d$mu_hat,
             aes(x = smoother_2d$mu_hat[,1], y = smoother_2d$mu_hat[,2]),
             size = 0.4)+
  geom_point(data = point_smoother_16,
             aes(x = x, y = y),
             size = 2, col = 'red') +
  geom_point(data = point_smoother_500,
             aes(x = x, y = y),
             size = 2, col = 'blue') +
  geom_tile(data = area_smoother_16, aes(x = x, y = y), fill = "red", alpha = .4)+
  geom_tile(data = area_smoother_500, aes(x = x, y = y), fill = "blue", alpha = .4)+
  geom_tile(data = grille_points1, aes(x = x, y = y), fill = "lightblue", alpha = 0.5)+
  geom_tile(data = grille_points2, aes(x = x, y = y), fill = "lightblue", alpha = 0.5)+
  labs(x = "Dimension 1 du lissage", y = "Dimension 2 du lissage")+
  ggtitle("Représentation des points estimés par la méthode du lissage")




#zoom ellispes

x_min_16 <- point_filter_16$x - 5  
x_max_16 <- point_filter_16$x + 5  
y_min_16 <- point_filter_16$y - 5  
y_max_16 <- point_filter_16$y + 5  

ggplot() +
  geom_point(data = filter_2d$mu_hat,
             aes(x = filter_2d$mu_hat[,1], y = filter_2d$mu_hat[,2]),
             size = 0.4) +
  geom_point(data = point_filter_16,
             aes(x = x, y = y),
             size = 2, col = 'red') +
  geom_tile(data = area_filter_16,
            aes(x = x, y = y),
            fill = "red", alpha = 0.3) +
  geom_tile(data = area_smoother_16,
            aes(x = x, y = y),
            fill = "green3", alpha = 0.3) +
  coord_cartesian(xlim = c(x_min_16, x_max_16), ylim = c(y_min_16, y_max_16))+
  labs(x = "Dimension 1 du filtrage", y = "Dimension 2 du filtrage")+
  ggtitle("Ellipses des confiances des méthodes du lissage et du filtrage pour le point 16")


x_min_500 <- point_filter_500$x - 5  
x_max_500 <- point_filter_500$x + 5  
y_min_500 <- point_filter_500$y - 5  
y_max_500 <- point_filter_500$y + 5  

ggplot() +
  geom_point(data = filter_2d$mu_hat,
             aes(x = filter_2d$mu_hat[,1], y = filter_2d$mu_hat[,2]),
             size = 0.4) +
  geom_point(data = point_filter_500,
             aes(x = x, y = y),
             size = 2, col = 'blue') +
  geom_tile(data = area_filter_500,
            aes(x = x, y = y),
            fill = "blue", alpha = 0.4) +
  geom_tile(data = area_smoother_500,
            aes(x = x, y = y),
            fill = "grey3", alpha = 0.3) +
  coord_cartesian(xlim = c(x_min_500, x_max_500), ylim = c(y_min_500, y_max_500))+
  labs(x = "Dimension 1 du filtrage", y = "Dimension 2 du filtrage")+
  ggtitle("Ellipses des confiances des méthodes du lissage et du filtrage pour le point 500")
