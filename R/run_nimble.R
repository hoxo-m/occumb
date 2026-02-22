run_nimble <- function(data, inits = NULL, parameters.to.save, model.file,
                       const, list_covs_phi, list_covs_theta, list_covs_psi,
                       n.chains, n.adapt = NULL, n.iter, n.burnin = 0, n.thin = 1,
                       modules = c("glm"), factories = NULL, parallel = FALSE,
                       n.cores = NULL, DIC = TRUE, store.data = FALSE, codaOnly = FALSE,
                       seed = NULL, bugs.format = FALSE, verbose = TRUE, ...) {
  dat <- data
  len_m_phi   <- length(dat$m_phi)
  len_m_theta <- length(dat$m_theta)
  len_m_psi   <- length(dat$m_psi)
  pad2 <- function(x) if (length(x) == 1L) c(x, x) else x
  dat$m_phi   <- pad2(dat$m_phi)
  dat$m_theta <- pad2(dat$m_theta)
  dat$m_psi   <- pad2(dat$m_psi)
  params <- parameters.to.save
  model_code <- model.file

  body <- body(inits)
  n_rho <- as.integer(dat$M * (dat$M - 1) / 2)
  body[[2]][names(body[[2]]) == "rho"][[1]] <- substitute(double(n), list(n = n_rho))
  body(inits) <- body

  # Run MCMC in NIMBLE
  const_nimble <- c(
    const[c("I", "J", "K", "N")],
    dat[c("M", "M_phi_shared", "M_theta_shared", "M_psi_shared",
          "m_phi", "m_theta", "m_psi", "prior_prec", "prior_ulim")],
    M_phi = len_m_phi, M_theta = len_m_theta, M_psi = len_m_psi)
  const_nimble$rho_index <- make_rho_index(const_nimble$M)
  const_nimble <- const_nimble[!is.na(names(const_nimble))]
  dat_nimble <- dat[c("y", "cov_phi", "cov_theta", "cov_psi",
                      "cov_phi_shared", "cov_theta_shared", "cov_psi_shared")]
  dat_nimble <- dat_nimble[!is.na(names(dat_nimble))]
  start_time <- Sys.time()
  if (parallel) {
    fit <- run_nimble_parallel(code = model_code, const = const_nimble,
                               data = dat_nimble, inits = inits,
                               monitors = params,
                               n.iter = n.iter, n.burnin = n.burnin,
                               n.thin = n.thin, n.chains = n.chains,
                               n.cores = n.cores)
  } else {
    fit <- nimble::nimbleMCMC(code = model_code, constants = const_nimble,
                              data = dat_nimble, inits = inits,
                              monitors = params,
                              thin = n.thin, niter = n.iter, nburnin = n.burnin,
                              nchains = n.chains, WAIC = FALSE, ...)
  }
  elapsed_mins <- round(as.numeric(Sys.time() - start_time, units = "mins"),
                        digits = 3)
  message("Summarizing MCMC samples...")
  fit <- nimbleSummary(fit)
  message("Finished")
  rownames(fit$summary) <- gsub(pattern = "\\s", replacement = "",
                                rownames(fit$summary))
  fit$parallel <- parallel
  fit$parameters <- params
  fit$model <- model_code
  # fit$modfile
  fit$run.date <- start_time
  fit$mcmc.info$n.burnin <- n.burnin
  fit$mcmc.info$n.thin <- n.thin
  fit$mcmc.info$elapsed.mins <- elapsed_mins

  fit 
}

run_nimble_parallel <- function(code, const, data, inits, dimensions, monitors,
                                n.iter, n.burnin, n.thin, n.chains, n.cores) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for this function. Please install it.", call. = FALSE)
  }
  if (is.null(n.cores)) n.cores <- parallel::detectCores()
  n.cores <- min(n.cores, n.chains)
  cluster <- parallel::makeCluster(n.cores)
  results <- parallel::parLapply(cl = cluster, X = seq_len(n.chains),
                                 fun = run_nimble_MCMC, code = code,
                                 const = const,
                                 data = data, inits = inits,
                                 dimensions = dimensions, monitors = monitors,
                                 n.iter = n.iter, n.burnin = n.burnin,
                                 n.thin = n.thin, n.chains = 1L)
  parallel::stopCluster(cl = cluster)

  results
}

run_nimble_MCMC <- function(seed, code, const, data, inits, dimensions, monitors,
                            n.iter, n.burnin, n.thin, n.chains) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required for this function. Please install it.", call. = FALSE)
  }
  
  model <- nimble::nimbleModel(code = code, constants = const, data = data,
                               inits = inits(), dimensions = dimensions)
  Cmodel <- nimble::compileNimble(model)
  MCMC <- nimble::buildMCMC(Cmodel, monitors = monitors)
  CMCMC <- nimble::compileNimble(MCMC)
  result <- nimble::runMCMC(CMCMC, niter = n.iter, nburnin = n.burnin,
                            thin = n.thin, nchains = n.chains, inits = inits,
                            setSeed = seed, WAIC = FALSE)
  result
}

# Auto-generate JAGS model code
write_nimble_model <- function(phi, theta, psi,
                               phi_shared, theta_shared, psi_shared,
                               M_cov_phi, M_cov_phi_shared,
                               M_cov_theta, M_cov_theta_shared,
                               M_cov_psi, M_cov_psi_shared) {

  model <- readLines(system.file("nimble",
                                 "occumb_template1.nimble",
                                 package = "occumb"))

  if (phi == "i")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i], 1)")
  if (phi == "ij")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i, j], 1)")
  if (phi == "ijk")
    model <- c(model,
               "                r[i, j, k] ~ dgamma(phi[i, j, k], 1)")

  model <- c(model,
             readLines(system.file("nimble",
                                   "occumb_template2.nimble",
                                   package = "occumb")))

  if (theta == "i")
    model <- c(model,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i])")
  if (theta == "ij")
    model <- c(model,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j])")
  if (theta == "ijk")
    model <- c(model,
               "                u[i, j, k] ~ dbern(z[i, j] * theta[i, j, k])")

  model <- c(model,
             readLines(system.file("jags",
                                   "occumb_template3.jags",
                                   package = "occumb")))

  if (psi == "i")
    model <- c(model,
               "            z[i, j] ~ dbern(psi[i])")
  if (psi == "ij")
    model <- c(model,
               "            z[i, j] ~ dbern(psi[i, j])")

  model <- c(model,
             readLines(system.file("jags",
                                   "occumb_template4.jags",
                                   package = "occumb")))

  if (phi_shared) {
    if (phi == "i") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, 1:M_phi_shared])"
      }
      model <- c(model, paste0(
                 "        log(phi[i]) <- ", term1, " + ", term2))
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, 1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, 1:M_phi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                 "            log(phi[i, j]) <- ", term1, " + ", term2),
                 "        }")
    } else if (phi == "ijk") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, k, 1]"
      } else {
        term1 <- "inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(alpha_shared[1:M_phi_shared], cov_phi_shared[i, j, k, 1:M_phi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                 "                log(phi[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
    }
  } else {
    if (phi == "i") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        log(phi[i]) <- alpha[i, 1] * cov_phi[1]")
      } else {
        model <- c(model,
                   "        log(phi[i]) <- inprod(alpha[i, 1:M_phi], cov_phi[1:M_phi])")
      }
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- alpha[i, 1] * cov_phi[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, 1:M_phi])",
                   "        }")
      }
    } else if (phi == "ijk") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- alpha[i, 1] * cov_phi[j, k, 1]",
                   "            }",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                log(phi[i, j, k]) <- inprod(alpha[i, 1:M_phi], cov_phi[j, k, 1:M_phi])",
                   "            }",
                   "        }")
      }
    }
  }

  if (theta_shared) {
    if (theta == "i") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, 1:M_theta_shared])"
      }
      model <- c(model, paste0(
                 "        logit(theta[i]) <- ", term1, " + ", term2))
    } else if (theta == "ij") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, 1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, 1:M_theta_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                 "            logit(theta[i, j]) <- ", term1, " + ", term2),
                 "        }")
    } else if (theta == "ijk") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, k, 1]"
      } else {
        term1 <- "inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(beta_shared[1:M_theta_shared], cov_theta_shared[i, j, k, 1:M_theta_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {",
                 "            for (k in 1:K) {", paste0(
                 "                logit(theta[i, j, k]) <- ", term1, " + ", term2),
                 "            }",
                 "        }")
    }
  } else {
    if (theta == "i")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        logit(theta[i]) <- beta[i, 1] * cov_theta[1]")
      } else {
        model <- c(model,
                   "        logit(theta[i]) <- inprod(beta[i, 1:M_theta], cov_theta[1:M_theta])")
      }
    if (theta == "ij")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- beta[i, 1] * cov_theta[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(theta[i, j]) <- inprod(beta[i, 1:M_theta], cov_theta[j, 1:M_theta])",
                   "        }")
      }
    if (theta == "ijk")
      if (M_cov_theta == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- beta[i, 1] * cov_theta[j, k, 1]",
                   "            }",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            for (k in 1:K) {",
                   "                logit(theta[i, j, k]) <- inprod(beta[i, 1:M_theta], cov_theta[j, k, 1:M_theta])",
                   "            }",
                   "        }")
      }
  }

  if (psi_shared) {
    if (psi == "i") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[1]"
      } else {
        term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, 1]"
      } else {
        term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, 1:M_psi_shared])"
      }
      model <- c(model,
                 paste0("        logit(psi[i]) <- ", term1, " + ", term2))
  
    } else if (psi == "ij") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[j, 1]"
      } else {
        term1 <- "inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(gamma_shared[1:M_psi_shared], cov_psi_shared[i, j, 1:M_psi_shared])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                 "            logit(psi[i, j]) <- ", term1, " + ", term2),
                 "        }")
    }
  } else {
    if (psi == "i") {
      if (M_cov_psi == 1) {
        model <- c(model,
                   "        logit(psi[i]) <- gamma[i, 1] * cov_psi[1]")
      } else {
        model <- c(model,
                   "        logit(psi[i]) <- inprod(gamma[i, 1:M_psi], cov_psi[1:M_psi])")
      }
    } else if (psi == "ij") {
      if (M_cov_psi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- gamma[i, 1] * cov_psi[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            logit(psi[i, j]) <- inprod(gamma[i, 1:M_psi], cov_psi[j, 1:M_psi])",
                   "        }")
      }
    }
  }

  model <- c(model,
             readLines(system.file("nimble",
                                   "occumb_template5.nimble",
                                   package = "occumb")))

  if (phi_shared)
    model <- c(model,
               "    for (m in 1:M_phi_shared) {",
               "        alpha_shared[m] ~ dnorm(0, prior_prec)",
               "    }")
  if (theta_shared)
    model <- c(model,
               "    for (m in 1:M_theta_shared) {",
               "        beta_shared[m] ~ dnorm(0, prior_prec)",
               "    }")
  if (psi_shared)
    model <- c(model,
               "    for (m in 1:M_psi_shared) {",
               "        gamma_shared[m] ~ dnorm(0, prior_prec)",
               "    }")

  model <- c(model, "}", "")

  model
}

make_rho_index <- function(M) {
  index <- matrix(0L, M, M)
  ind <- 1L
  for (m1 in 1:(M - 1)) {
    for (m2 in (m1 + 1):M) {
      index[m1, m2] <- ind
      ind <- ind + 1L
    }
  }
  index
}
