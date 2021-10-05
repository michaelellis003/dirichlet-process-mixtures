# simulate data ----------------------------------------------------------------
simulate_data <- function(
    N = 100,    # sample size
    K = 3,      # Number of mixtures
    mu = c(-10, 0, 10),
    sigmasq = c(0.5, 0.5, 0.5),
    probs = c(1/3, 1/3, 1/3),
    plot = TRUE
) {
    library(ggplot2)
    
    sample_clusters <- sample(1:K, size = N, replace = TRUE)
    y <- rnorm(N, mean = mu[sample_clusters], sd = sqrt(sigmasq[sample_clusters]))
    
    if(plot) {
        df <- data.frame(
            Sample = 1:N,
            y = y
        )
        
        ggplot2::ggplot(df, aes(x = y)) + 
            geom_density(kernel = "gaussian")
    }
    
    return(
        list(
            y = y,
            N = N, 
            K = K,
            mu = mu,
            sigmasq = sigmasq,
            probs = probs
        )
    )
}

# Stick breaking function ------------------------------------------------------
stick_breaking <- function(V, K) { 
    
    pi <- rep(NA, K)
    
    for(k in 1:K) {
        if(k == 1){
            pi[k] <- V[k]
        } else {
            stick_prob <- V[k]
            for(j in 1:(k-1)){
                stick_prob <- stick_prob*(1 - V[j])
            }
            pi[k] <- stick_prob
        }
    }
    
    return(pi)
}

# mcmc function ----------------------------------------------------------------
mcmc <- function(
    y,
    K = 20,
    lambda0 = 0.1,
    mu0 = 0,
    sigmasq0 = 10^2,
    alpha0 = 0.001,
    beta0 = 0.001,
    n_mcmc = 5000,
    burnin = 2500,
    n_message = 500
) { 
    
    N <- length(y) # sample size
    n_save <- n_mcmc-burnin
    
    ## initial parameters
    mu <- rep(mu0, K) 
    sigmasq <- rep(sigmasq0, K)
    V <- c(rbeta(K-1, 1, lambda0), 1)
    pi <- stick_breaking(V, K)
    z <- sample(c(1:2, 5), size = N, replace = TRUE)
    
    ## save parameters
    mu_save <- matrix(NA, nrow = n_save, ncol = K)
    sigmasq_save <- matrix(NA, nrow = n_save, ncol = K)
    pi_save <- matrix(NA, nrow = n_save, ncol = K)
    z_save <- matrix(NA, nrow = n_save, ncol = N)
    
    for(i in 1:n_mcmc) {
        if(i %% n_message == 0) {
            message("Iteration ", i , " out of ", n_mcmc)
        }
        
        clusters <- unique(z)
        counts <- table(z)
        y_sum_k <- data.frame(y=y, z=z)
        y_sum_k <- aggregate(x = y_sum_k$y,
                             by = list(z = y_sum_k$z), 
                             FUN = sum)
        
        for(k in 1:K) {
            if(k %in% clusters) {
                n_k <- as.vector(counts[names(counts) == k])
                
                ## update parameters for mu
                a <- 1/(1/sigmasq0 + n_k/sigmasq[k])
                b <- mu0/sigmasq0 + y_sum_k[y_sum_k$z == k, 2]/sigmasq[k]
                
                ## sample mu
                mu[k] <- rnorm(1, a*b, sqrt(a))
                
                ## updated parameters for sigma squared
                SS <- data.frame(z=z, res=(y-mu[k])^2)
                SS <- aggregate(x = SS$res,
                                by = list(z = SS$z), 
                                FUN = sum)
                SS <- SS[SS$z == k, 2]
                updated_alpha <- alpha0 + n_k/2
                updated_beta <- beta0 + SS/2
                
                ## sample sigma squared
                sigmasq[k] <- 1/rgamma(1, updated_alpha, updated_beta)
                
            } else { ## sample from prior
                mu[k] <- rnorm(1, mu0, sqrt(sigmasq0))
                sigmasq[k] <- 1/rgamma(1, alpha0, beta0)
            }
        }
        
        max_cluster <- max(clusters)
        for(k in 1:(K-1)) {
            if(k %in% clusters) {
                n_k <- as.integer(counts[names(counts) == k])
            } else {
                n_k <- 0
            }
            
            if(k <= max_cluster) {
                sum_counts <- sum(unname(counts[as.integer(names(counts)) > k]))
            }
            
            ## updated parameters for pi tilde
            updated_shape <- 1 + n_k
            updated_lambda <- lambda0 + sum_counts
            
            V[k] <- rbeta(1, updated_shape, updated_lambda)
        }
        V[K] <- 1
        
        ## update pi
        pi <- stick_breaking(V, K)
        
        ## sample z
        for(j in 1:N) {
            ## probabilities for z = k
            z_k_dens <- rep(NA, K)
            z_k_probs <- rep(NA, K)
            for(k in 1:K) {
                z_k_dens[k] <- pi[k]*dnorm(y[j], mu[k], sqrt(sigmasq[k]))
            }
            for(k in 1:K){
                z_k_probs[k] <- z_k_dens[k]/sum(z_k_dens)
            }
            
            z[j] <- sample(1:K, size = 1, replace = TRUE, prob = z_k_probs)
        }
    
        
        if(i > burnin) {
            index <- i-burnin
            mu_save[index, ] <- mu
            sigmasq_save[index, ] <- sigmasq
            pi_save[index, ] <- pi
            z_save[index, ] <- z
        }
    }
    
    return(
        list(
            mu_save = mu_save,
            sigmasq_save = sigmasq_save,
            pi_save = pi_save,
            z_save = z_save
        )
    )
    
}