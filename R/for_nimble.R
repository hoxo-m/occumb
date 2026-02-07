# Auto-generate JAGS model code
write_nimble_model <- function(phi, theta, psi,
                               phi_shared, theta_shared, psi_shared,
                               M_cov_phi, M_cov_phi_shared,
                               M_cov_theta, M_cov_theta_shared,
                               M_cov_psi, M_cov_psi_shared) {
  
  model <- readLines(system.file("jags",
                                 "occumb_template1.jags",
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
             readLines(system.file("jags",
                                   "occumb_template2.jags",
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
        term1 <- "inprod(alpha[i, ], cov_phi[])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, 1]"
      } else {
        term2 <- "inprod(alpha_shared[], cov_phi_shared[i, ])"
      }
      model <- c(model, paste0(
                 "        log(phi[i]) <- ", term1, " + ", term2))
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, 1]"
      } else {
        term1 <- "inprod(alpha[i, ], cov_phi[j, ])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(alpha_shared[], cov_phi_shared[i, j, ])"
      }
      model <- c(model, 
                 "        for (j in 1:J) {", paste0(
                 "            log(phi[i, j]) <- ", term1, " + ", term2),
                 "        }")
    } else if (phi == "ijk") {
      if (M_cov_phi == 1) {
        term1 <- "alpha[i, 1] * cov_phi[j, k, 1]"
      } else {
        term1 <- "inprod(alpha[i, ], cov_phi[j, k, ])"
      }
      if (M_cov_phi_shared == 1) {
        term2 <- "alpha_shared[1] * cov_phi_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(alpha_shared[], cov_phi_shared[i, j, k, ])"
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
                   "        log(phi[i]) <- inprod(alpha[i, ], cov_phi[])")
      }
    } else if (phi == "ij") {
      if (M_cov_phi == 1) {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- alpha[i, ] * cov_phi[j, 1]",
                   "        }")
      } else {
        model <- c(model,
                   "        for (j in 1:J) {",
                   "            log(phi[i, j]) <- inprod(alpha[i, ], cov_phi[j, ])",
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
                   "                log(phi[i, j, k]) <- inprod(alpha[i, ], cov_phi[j, k, ])",
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
        term1 <- "inprod(beta[i, ], cov_theta[])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, 1]"
      } else {
        term2 <- "inprod(beta_shared[], cov_theta_shared[i, ])"
      }
      model <- c(model, paste0(
                 "        logit(theta[i]) <- ", term1, " + ", term2))
    } else if (theta == "ij") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, 1]"
      } else {
        term1 <- "inprod(beta[i, ], cov_theta[j, ])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, 1]"
      } else {
        term2 <- "inprod(beta_shared[], cov_theta_shared[i, j, ])"
      }
      model <- c(model,
                 "        for (j in 1:J) {", paste0(
                 "            logit(theta[i, j]) <- ", term1, " + " + term2),
                 "        }")
    } else if (theta == "ijk") {
      if (M_cov_theta == 1) {
        term1 <- "beta[i, 1] * cov_theta[j, k, 1]"
      } else {
        term1 <- "inprod(beta[i, ], cov_theta[j, k, ])"
      }
      if (M_cov_theta_shared == 1) {
        term2 <- "beta_shared[1] * cov_theta_shared[i, j, k, 1]"
      } else {
        term2 <- "inprod(beta_shared[], cov_theta_shared[i, j, k, ])"
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
                   "        logit(theta[i]) <- inprod(beta[i, ], cov_theta[])")
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
                   "            logit(theta[i, j]) <- inprod(beta[i, ], cov_theta[j, ])",
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
                   "                logit(theta[i, j, k]) <- inprod(beta[i, ], cov_theta[j, k, ])",
                   "            }",
                   "        }")
      }
  }
  
  if (psi_shared) {
    if (psi == "i") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[1]"
      } else {
        term1 <- "inprod(gamma[i, ], cov_psi[])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, 1]"
      } else {
        term2 <- "inprod(gamma_shared[], cov_psi_shared[i, ])"
      }
      model <- c(model,
                 paste0("        logit(psi[i]) <- ", term1, " + ", term2))
    
    } else if (psi == "ij") {
      if (M_cov_psi == 1) {
        term1 <- "gamma[i, 1] * cov_psi[j, 1]"
      } else {
        term1 <- "inprod(gamma[i, ], cov_psi[j, ])"
      }
      if (M_cov_psi_shared == 1) {
        term2 <- "gamma_shared[1] * cov_psi_shared[i, j, 1]"
      } else {
        term2 <- "inprod(gamma_shared[], cov_psi_shared[i, j, ])"
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
                   "        logit(psi[i]) <- inprod(gamma[i, ], cov_psi[])")
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
                   "            logit(psi[i, j]) <- inprod(gamma[i, ], cov_psi[j, ])",
                   "        }")
      }
    }
  }
  
  model <- c(model,
             readLines(system.file("jags",
                                   "occumb_template5.jags",
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
  
  to_nimble_model_code(model)
}

to_nimble_model_code <- function(model_code) {
  model_code[1] <- "{"
  model_code <- paste0(model_code, collapse = "\n")
  model_code <- str2lang(model_code)
  model_code <- replace_calls(model_code, list(`dmnorm.vcov(Mu[1:M], Sigma)` = quote(dmnorm(Mu[1:M], cov = Sigma[1:M, 1:M]))))
  model_code <- replace_calls(model_code, list(`length(m_phi)` = quote(len_m_phi)))
  model_code <- replace_calls(model_code, list(`length(m_theta)` = quote(len_m_theta)))
  model_code <- replace_calls(model_code, list(`length(m_psi)` = quote(len_m_psi)))
  model_code <- replace_calls(model_code, list(`ifelse` = quote(.nimble_ifelse)))
  model_code
}

replace_calls <- function(expr, rules) {
  if (is.call(expr)) {
    if (as.character(expr[[1]]) %in% names(rules)) {
      expr[[1]] <- rules[[as.character(expr[[1]])]]
      return(expr)
    }
    
    key <- paste(deparse(expr), collapse = "")
    
    if (key %in% names(rules)) {
      return(rules[[key]])
    }
    
    for (i in seq_along(expr)) {
      expr[[i]] <- replace_calls(expr[[i]], rules)
    }
  }
  expr
}

