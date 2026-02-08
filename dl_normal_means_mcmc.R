# Dirichlet-Laplace prior MCMC for normal means model (redundancy-free Gibbs)

rinvgauss_msh <- function(mu, lambda) {
  if (any(mu <= 0) || lambda <= 0) {
    stop("mu and lambda must be positive")
  }
  v <- rnorm(length(mu))^2
  mu2 <- mu^2
  term <- mu + (mu2 * v) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * v + mu2 * v^2)
  u <- runif(length(mu))
  x <- ifelse(u <= mu / (mu + term), term, mu2 / term)
  x
}

rgig_wrapper <- function(n, lambda, chi, psi) {
  if (!requireNamespace("GIGrvg", quietly = TRUE)) {
    stop("Package 'GIGrvg' is required for GIG sampling. Please install it.")
  }
  GIGrvg::rgig(n = n, lambda = lambda, chi = chi, psi = psi)
}

dl_gibbs_sampler <- function(
  y,
  a,
  n_iter = 3000,
  burnin = 1000,
  thin = 2,
  seed = 123
) {
  n <- length(y)
  set.seed(seed)

  theta <- rep(0, n)
  delta <- rgamma(n, shape = a, rate = 1 / 2)
  psi <- rexp(n, rate = 1 / 2)

  n_keep <- floor((n_iter - burnin) / thin)
  theta_draws <- matrix(NA_real_, nrow = n_keep, ncol = n)
  tau_draws <- numeric(n_keep)

  keep_idx <- 0
  for (iter in seq_len(n_iter)) {
    sigma2 <- 1 / (1 + 1 / (psi * delta^2))
    mu <- sigma2 * y
    theta <- rnorm(n, mean = mu, sd = sqrt(sigma2))

    abs_theta <- pmax(abs(theta), 1e-12)
    chi_val <- pmax(2 * abs_theta, 1e-12)
    delta <- rgig_wrapper(n, lambda = a - 1, chi = chi_val, psi = 1)

    u <- rinvgauss_msh(mu = delta / abs_theta, lambda = 1)
    psi <- 1 / u

    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      keep_idx <- keep_idx + 1
      theta_draws[keep_idx, ] <- theta
      tau_draws[keep_idx] <- sum(delta)
    }
  }

  list(
    theta_draws = theta_draws,
    tau_draws = tau_draws
  )
}

rinvgamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

hs_gibbs_sampler <- function(
  y,
  n_iter = 3000,
  burnin = 1000,
  thin = 2,
  seed = 123
) {
  n <- length(y)
  set.seed(seed)

  theta <- rep(0, n)
  lambda2 <- rep(1, n)
  nu <- rep(1, n)
  tau2 <- 1
  xi <- 1

  n_keep <- floor((n_iter - burnin) / thin)
  theta_draws <- matrix(NA_real_, nrow = n_keep, ncol = n)

  keep_idx <- 0
  for (iter in seq_len(n_iter)) {
    sigma2 <- 1 / (1 + 1 / (lambda2 * tau2))
    mu <- sigma2 * y
    theta <- rnorm(n, mean = mu, sd = sqrt(sigma2))

    rate_lambda <- (1 / nu) + (theta^2 / (2 * tau2))
    lambda2 <- rinvgamma(n, shape = 1, rate = rate_lambda)

    rate_nu <- 1 + 1 / lambda2
    nu <- rinvgamma(n, shape = 0.5, rate = rate_nu)

    rate_tau <- (1 / xi) + sum(theta^2 / (2 * lambda2))
    tau2 <- rinvgamma(1, shape = (n + 1) / 2, rate = rate_tau)

    rate_xi <- 1 + 1 / tau2
    xi <- rinvgamma(1, shape = 0.5, rate = rate_xi)

    if (iter > burnin && ((iter - burnin) %% thin == 0)) {
      keep_idx <- keep_idx + 1
      theta_draws[keep_idx, ] <- theta
    }
  }

  list(theta_draws = theta_draws)
}

simulate_dl_normal_means <- function(
  n,
  q,
  A,
  a,
  n_iter = 3000,
  burnin = 1000,
  thin = 2,
  seed = 123
) {
  if (q > n) {
    stop("q must be <= n")
  }
  set.seed(seed)

  theta_true <- c(rep(A, q), rep(0, n - q))
  y <- rnorm(n, mean = theta_true, sd = 1)

  dl_draws <- dl_gibbs_sampler(
    y = y,
    a = a,
    n_iter = n_iter,
    burnin = burnin,
    thin = thin,
    seed = seed
  )

  theta_mean <- colMeans(dl_draws$theta_draws)
  theta_median <- apply(dl_draws$theta_draws, 2, median)
  theta_ci <- apply(dl_draws$theta_draws, 2, quantile, probs = c(0.025, 0.975))
  mse_mean <- sum((theta_mean - theta_true)^2)
  mse_median <- sum((theta_median - theta_true)^2)

  list(
    theta_true = theta_true,
    y = y,
    theta_mean = theta_mean,
    theta_median = theta_median,
    theta_ci = theta_ci,
    mse_mean = mse_mean,
    mse_median = mse_median,
    tau_draws = dl_draws$tau_draws,
    theta_draws = dl_draws$theta_draws
  )
}

run_se_comparison <- function(
  n_set,
  q_set,
  A_set,
  R = 30,
  n_iter = 3000,
  burnin = 1000,
  thin = 2,
  seed = 1
) {
  for (n in n_set) {
    for (q in q_set) {
      if (q > n) {
        next
      }
      results <- matrix(NA_real_, nrow = 3, ncol = length(A_set))
      rownames(results) <- c("HS", "DL_1n", "DL_01")
      colnames(results) <- paste0("A=", A_set)

      for (a_idx in seq_along(A_set)) {
        A <- A_set[a_idx]
        se_hs <- numeric(R)
        se_dl_1n <- numeric(R)
        se_dl_01 <- numeric(R)

        for (r in seq_len(R)) {
          seed_r <- seed + 10000 * n + 100 * q + r
          set.seed(seed_r)
          theta_true <- c(rep(A, q), rep(0, n - q))
          y <- rnorm(n, mean = theta_true, sd = 1)

          dl_1n_draws <- dl_gibbs_sampler(
            y = y,
            a = 1 / n,
            n_iter = n_iter,
            burnin = burnin,
            thin = thin,
            seed = seed_r
          )
          dl_01_draws <- dl_gibbs_sampler(
            y = y,
            a = 0.1,
            n_iter = n_iter,
            burnin = burnin,
            thin = thin,
            seed = seed_r + 1
          )
          hs_draws <- hs_gibbs_sampler(
            y = y,
            n_iter = n_iter,
            burnin = burnin,
            thin = thin,
            seed = seed_r + 2
          )

          theta_dl_1n <- apply(dl_1n_draws$theta_draws, 2, median)
          theta_dl_01 <- apply(dl_01_draws$theta_draws, 2, median)
          theta_hs <- apply(hs_draws$theta_draws, 2, median)

          se_dl_1n[r] <- sum((theta_dl_1n - theta_true)^2)
          se_dl_01[r] <- sum((theta_dl_01 - theta_true)^2)
          se_hs[r] <- sum((theta_hs - theta_true)^2)
        }

        results["HS", a_idx] <- mean(se_hs)
        results["DL_1n", a_idx] <- mean(se_dl_1n)
        results["DL_01", a_idx] <- mean(se_dl_01)
      }

      cat(sprintf("\nn = %d, q = %d\n", n, q))
      print(round(results, 2))
    }
  }
}

if (sys.nframe() == 0) {
  demo <- simulate_dl_normal_means(
    n = 200,
    q = 20,
    A = 8,
    a = 1 / 200,
    n_iter = 3000,
    burnin = 1000,
    thin = 2,
    seed = 42
  )

  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
  plot(demo$theta_true, demo$theta_median,
       main = "True theta vs Posterior Median",
       xlab = "True theta", ylab = "Posterior median",
       pch = 16, col = "steelblue")
  abline(0, 1, col = "firebrick", lty = 2)

  plot(demo$y, demo$theta_median,
       main = "y vs Posterior Median",
       xlab = "y", ylab = "Posterior median",
       pch = 16, col = "darkgreen")
  abline(0, 1, col = "firebrick", lty = 2)

  hist(demo$tau_draws, breaks = 30,
       main = "Histogram of tau",
       xlab = "tau", col = "gray")

  message(sprintf("MSE (mean): %.3f", demo$mse_mean))
  message(sprintf("MSE (median): %.3f", demo$mse_median))
}
