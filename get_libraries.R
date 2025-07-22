library(ggplot2)
library(plot3D)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(RVCompare)
library(mixtools)


get_mahalanobis <- function(x, y, mu, sigma) {
  vect <- c(x, y) - mu
  return(t(vect) %*% solve(sigma) %*% vect)
}