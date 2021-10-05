source("functions.R")

y <- simulate_data(N = 200)
output <- mcmc(y$y)

colMeans(output$mu_save)

