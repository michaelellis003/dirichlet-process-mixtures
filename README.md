# Dirichlet Process Mixture Models

Bayesian nonparametric clustering using Dirichlet process mixture models with Gibbs sampling in R. The model automatically learns the number of clusters from data using the stick-breaking construction.

## Model

Implements a truncated Dirichlet process Gaussian mixture model:

```
y_i ~ sum_{k=1}^{K} pi_k * N(mu_k, sigma_k^2)

pi_k = V_k * prod_{j=1}^{k-1} (1 - V_j)    (stick-breaking)
V_k  ~ Beta(1, lambda)
mu_k ~ N(mu_0, sigma_0^2)
sigma_k^2 ~ InvGamma(alpha_0, beta_0)
```

The stick-breaking prior (Sethuraman, 1994) allows the model to use fewer components than the truncation level K, effectively discovering the number of clusters.

## Features

- Gibbs sampling with conjugate updates for all parameters
- Stick-breaking construction for mixture weights
- Automatic cluster discovery (truncated at K = 20)
- MCMC diagnostic trace plots for all parameters

## Dependencies

```r
install.packages(c("ggplot2", "reshape2"))
```

## Usage

```r
source("functions.R")

# Simulate 3-component mixture data
y <- simulate_data(N = 200, K = 3,
                   mu = c(-20, 5, 25),
                   sigmasq = c(0.5, 0.5, 0.5))

# Run Gibbs sampler (5000 iterations, 2500 burn-in)
output <- mcmc(y$y)

# Posterior means of component locations
colMeans(output$mu_save)

# Trace plots for diagnostics
plot_mcmc(output)
```

See [`main.R`](main.R) for a complete example.

## References

- Sethuraman, J. (1994). A constructive definition of Dirichlet priors. *Statistica Sinica*, 4, 639-650.
- Neal, R. M. (2000). Markov chain sampling methods for Dirichlet process mixture models. *Journal of Computational and Graphical Statistics*, 9(2), 249-265.
- Ferguson, T. S. (1973). A Bayesian analysis of some nonparametric problems. *The Annals of Statistics*, 1(2), 209-230.
