source("functions.R")

y <- simulate_data()
output <- mcmc(y$y)

colSums(output$mu_save)

