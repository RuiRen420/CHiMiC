##' @noRd
##'
P_mcp <- function(zj, lambda, gamma = 3) {
  zj <- as.numeric(zj)[1]
  lambda <- as.numeric(lambda)[1]
  gamma <- as.numeric(gamma)[1]
  if (!is.finite(zj) || !is.finite(lambda) || !is.finite(gamma)) {
    return(0)
  }
  if (lambda < 0) lambda <- 0
  if (abs(zj) <= gamma * lambda) {
    s <- sign(zj) * pmax(0, abs(zj) - lambda)
    s / (1 - 1 / gamma)
  } else {
    zj
  }
}


##' @noRd
##'
one_hot_expand <- function(mat) {
  n <- nrow(mat)
  p <- ncol(mat)
  one_hot_list <- vector("list", p)
  for (j in seq_len(p)) {
    Kj <- max(mat[, j])
    m <- matrix(0, nrow = n, ncol = Kj)
    m[cbind(seq_len(n), mat[, j])] <- 1
    one_hot_list[[j]] <- m[, -1] # use column 1 as the control level
  }
  do.call(cbind, one_hot_list)
}


##' @noRd
##'
R2_CHiMiC <- function(Z, q1, q2, TE = NULL, w = 0, G, lambda1, lambda2, nfolds = 3,
                      max_iter = 100, inner_iter = 30, tol = 1e-2) {
  Y <- TE
  n <- nrow(Z)
  q <- ncol(Z)
  if (q1 > 0) {
    Z1 <- scale(Z[, 1:q1], center = TRUE, scale = TRUE)
  }
  if (q2 > 0) Z2 <- Z[, (q - q2 + 1):q]

  if (q1 == 0) Z <- Z2
  if (q2 == 0) Z <- Z1
  if (q1 > 0 & q2 > 0) Z <- cbind(Z1, Z2)

  if (!is.null(Y)) {
    Y <- scale(Y)
  }
  if (is.null(Y)) {
    w <- 0
  }

  folds <- createFolds(1:n, k = nfolds)
  cv_results <- lapply(folds, function(fold_index) {
    if (!is.null(Y)) {
      Ytrain <- Y[-fold_index]
    } else {
      Ytrain <- NULL
      Xtrain <- NULL
    }
    Ztrain <- Z[-fold_index, ]

    if (!is.null(Y)) {
      Ytest <- Y[fold_index]
    } else {
      Ytest <- NULL
      Xtest <- NULL
    }
    Ztest <- Z[fold_index, ]
    if (q1 > 0) {
      Z1train <- Z1[-fold_index, ]
      Z1test <- Z1[fold_index, ]
    }
    if (q2 > 0) {
      Z2train <- Z2[-fold_index, ]
      Z2test <- Z2[fold_index, ]
    }

    ntest <- nrow(Ztest)

    result_temp <- CHiMiC(
      Z = Ztrain, q1 = q1, q2 = q2, TE = Ytrain, w = w, G = G,
      lambda1 = lambda1, lambda2 = lambda2, max_iter = max_iter,
      inner_iter = inner_iter, tol = tol
    )
    if (q1 > 0) {
      mu_est <- result_temp$cons_mu
      Sigma_est <- matrix(result_temp$cons_Sigma, ncol = q1, nrow = G, byrow = TRUE)
      selected_features1 <- result_temp$selected_features1
    }
    if (q2 > 0) {
      pj_est <- result_temp$cate_prob
      selected_features2 <- result_temp$selected_features2
    }
    if (!is.null(Y)) {
      beta_est <- result_temp$te_beta
      sigma_est <- result_temp$te_sigma
    }
    pi_est <- result_temp$pi

    sum_logfij <- matrix(0, nrow = ntest, ncol = G)
    for (g in 1:G) {
      # Continuous part
      if (q1 > 0) {
        mu_mat <- matrix(mu_est[g, ], nrow = ntest, ncol = q1, byrow = TRUE)
        sigma_mat <- matrix(Sigma_est[g, ], nrow = ntest, ncol = q1, byrow = TRUE)
        log_cont <- rowSums(dnorm(Z1test, mean = mu_mat, sd = sigma_mat, log = TRUE))
      } else {
        log_cont <- 0
      }
      # Categorical part
      if (q2 > 0) {
        probg <- lapply(1:q2, function(l) pj_est[[l]][g, ])
        cate_mat <- sapply(1:q2, function(j) probg[[j]][Z2test[, j]])
        log_cat_prob <- rowSums(log(cate_mat))
      } else {
        log_cat_prob <- 0
      }
      # CATE part
      if (!is.null(Ytest)) {
        log_y <- dnorm(Ytest, mean = beta_est[g], sd = sigma_est[g], log = TRUE)
      } else {
        log_y <- 0
      }
      # Log likelihood
      sum_logfij[, g] <- (1 - w) * (log_cont + log_cat_prob) + w * log_y
    }
    logf_weighted <- sum_logfij + matrix(log(pi_est), nrow = ntest, ncol = G, byrow = TRUE)
    max_logf <- apply(logf_weighted, 1, max)
    logf_shifted <- sweep(logf_weighted, 1, max_logf, "-")
    log_denom <- log(rowSums(exp(logf_shifted))) + max_logf
    wi_test <- exp(logf_weighted - matrix(log_denom, nrow = ntest, ncol = G))

    clusters_test <- apply(wi_test, 1, which.max)

    # CATE part
    if (!is.null(Ytest)) {
      beta_test <- beta_est[clusters_test]
      sigma_test <- sigma_est[clusters_test]
      loglik_test <- sum(dnorm(Ytest, mean = beta_test, sd = sigma_test, log = TRUE))
      loglik_null <- sum(dnorm(Ytest, mean = mean(Ytest), sd = sd(Ytest), log = TRUE))
      R2y <- (loglik_null - loglik_test) / loglik_null
    } else {
      R2y <- 0
    }

    # Continuous part
    if (q1 > 0) {
      mu_mat_test <- mu_est[clusters_test, ]
      sigma_mat_test <- Sigma_est[clusters_test, ]
      log_cont_test <- colSums(dnorm(Z1test, mean = mu_mat_test, sd = sigma_mat_test, log = TRUE))

      mu_mat_null <- matrix(colMeans(Z1test), nrow = ntest, ncol = q1, byrow = TRUE)
      sigma_mat_null <- matrix(apply(Z1test, 2, sd), nrow = ntest, ncol = q1, byrow = TRUE)
      log_cont_null <- colSums(dnorm(Z1test, mean = mu_mat_null, sd = sigma_mat_null, log = TRUE))

      R2z1 <- (log_cont_null - log_cont_test) / log_cont_null
    } else {
      R2z1 <- 0
    }

    # Categorical part
    if (q2 > 0) {
      Kj <- apply(Z2, 2, max)
      pj_mean <- lapply(1:q2, function(l) {
        tab <- tabulate(factor(Z2test[, l], levels = 1:Kj[l]), nbins = Kj[l])
        tab / max(ntest, 1)
      })

      calc_multiclass_brier <- function(y_true, y_prob) {
        n <- length(y_true)
        K <- ncol(y_prob)
        o_mat <- matrix(0, nrow = n, ncol = K)
        o_mat[cbind(1:n, y_true)] <- 1
        mean(rowSums((y_prob - o_mat)^2))
      }

      R2z2_list <- lapply(1:q2, function(l) {
        pjk <- pj_est[[l]][clusters_test, ]
        o_mat <- matrix(0, nrow = nrow(pjk), ncol = ncol(pjk))
        o_mat[cbind(1:nrow(pjk), Z2test[, l])] <- 1
        nll_cat_test <- sum(log(pmax(rowSums(pjk * o_mat), 1e-6)))

        prob_mat <- matrix(pj_mean[[l]], ncol = Kj[l], nrow = ntest, byrow = TRUE)
        nll_cat_null <- sum(log(pmax(rowSums(prob_mat * o_mat), 1e-6)))

        (nll_cat_null - nll_cat_test) / nll_cat_null
      })
      R2z2 <- unlist(R2z2_list)
    } else {
      R2z2 <- 0
    }

    # Final R2
    if (q1 > 0 & q2 > 0) {
      R2_comp <- sum(R2z1[selected_features1]) + sum(R2z2[selected_features2])
    } else if (q1 > 0 & q2 == 0) {
      R2_comp <- sum(R2z1[selected_features1])
    } else if (q1 == 0 & q2 > 0) {
      R2_comp <- sum(R2z2[selected_features2])
    }

    R2 <- R2y + R2_comp

    return(list(R2 = R2))
  })

  R2_final <- mean(sapply(cv_results, function(x) x$R2))

  return(list(R2 = R2_final))
}


##' @noRd
##'
mBIC_CHiMiC <- function(Z, q1, q2, TE = NULL, ParameterSet, max_iter = 100,
                        inner_iter = 30, tol = 1e-2) {
  result_temp <- vector("list", nrow(ParameterSet))
  for (j in seq_len(nrow(ParameterSet))) {
    G_val <- ParameterSet$G[j]
    w_val <- ParameterSet$w[j]
    lambda1_val <- ParameterSet$lambda1[j]
    lambda2_val <- ParameterSet$lambda2[j]
    result_temp[[j]] <- tryCatch(
      {
        result_list <- CHiMiC(
          Z = Z, q1, q2, TE = TE, w = w_val, G = G_val, lambda1 = lambda1_val,
          lambda2 = lambda2_val, max_iter = max_iter, inner_iter = inner_iter, tol = tol
        )
        list(G = G_val, w = w_val, lambda1 = lambda1_val, lambda2 = lambda2_val, mBIC = result_list$mBIC)
      },
      error = function(e1) {
        list(G = G_val, w = w_val, lambda1 = lambda1_val, lambda2 = lambda2_val, mBIC = 1e10)
      }
    )
  }

  do.call(rbind, Filter(Negate(is.null), result_temp))
}


##' @noRd
##'
CV_CHiMiC <- function(Z, q1, q2, TE = NULL, ParameterSet, nfolds = 3,
                      max_iter = 100, inner_iter = 30, tol = 1e-2) {
  result_temp <- vector("list", nrow(ParameterSet))
  for (j in seq_len(nrow(ParameterSet))) {
    G_val <- ParameterSet$G[j]
    w_val <- ParameterSet$w[j]
    lambda1_val <- ParameterSet$lambda1[j]
    lambda2_val <- ParameterSet$lambda2[j]
    result_temp[[j]] <- tryCatch(
      {
        result_list <- R2_CHiMiC(
          Z = Z, q1, q2, TE = TE, w = w_val, G = G_val, lambda1 = lambda1_val,
          lambda2 = lambda2_val, max_iter = max_iter, inner_iter = inner_iter, tol = tol
        )
        list(G = G_val, w = w_val, lambda1 = lambda1_val, lambda2 = lambda2_val, R2 = result_list$R2)
      },
      error = function(e1) {
        list(G = G_val, w = w_val, lambda1 = lambda1_val, lambda2 = lambda2_val, R2 = -1e10)
      }
    )
  }

  R2 <- sapply(result_temp, function(x) x$R2)
  ind.final <- which.max(R2)[1]
  w_selected <- ParameterSet$w[ind.final]
  G_selected <- ParameterSet$G[ind.final]
  lambda1_selected <- ParameterSet$lambda1[ind.final]
  lambda2_selected <- ParameterSet$lambda2[ind.final]

  result <- CHiMiC(Z,
    q1 = q1, q2 = q2, TE = TE, w = w_selected, G = G_selected, lambda1 = lambda1_selected,
    lambda2 = lambda2_selected, max_iter = max_iter, inner_iter = inner_iter, tol = tol
  )
  return(result)
}


##' @noRd
##'
HiMiC <- function(Z, q1, q2, G, lambda1, lambda2, max_iter = 100,
                  inner_iter = 30, tol = 1e-2) {
  n <- nrow(Z)
  q <- ncol(Z)
  if (q1 > 0) {
    Z1 <- scale(Z[, 1:q1], center = TRUE, scale = TRUE)
  }
  if (q2 > 0) {
    Z2 <- Z[, (q - q2 + 1):q]
  }

  if (q1 == 0) {
    Z <- Z2
    Z0 <- Z2 - 1
  }
  if (q2 == 0) {
    Z <- Z1
    Z0 <- Z1
  }
  if (q1 > 0 & q2 > 0) {
    Z <- cbind(Z1, Z2)
    Z0 <- cbind(Z1, Z2 - 1)
  }

  km.out <- KMeansSparseCluster(x = Z0, K = G, maxiter = max_iter, wbounds = 1.2, silent = TRUE)
  cluster0 <- km.out[[1]][["Cs"]]

  # Initial values
  # Continuous part
  if (q1 > 0) {
    mu_g <- mu_gg <- matrix(rnorm(G * q1, mean = 0, sd = 0.001), G, q1)
    Sigma_g <- Sigma_gg <- Sigma_mat <- matrix(rnorm(G * q1, mean = 0, sd = 0.001), G, q1)
    for (g in 1:G) {
      idx <- which(cluster0 == g)
      Z1_g <- Z1[idx, , drop = FALSE]
      mu_g[g, ] <- mu_gg[g, ] <- colMeans(Z1_g)
      Sigma_g[g, ] <- Sigma_gg[g, ] <- apply(Z1_g, 2, sd)
    }
  }
  # Categorical part
  if (q2 > 0) {
    Kj <- apply(Z2, 2, max)
    nj_g <- nj_gg <- lapply(1:q2, function(l) {
      mat <- matrix(0, G, Kj[l])
      for (g in 1:G) {
        idx <- which(cluster0 == g)
        mat[g, ] <- tabulate(factor(Z2[idx, l], levels = 1:Kj[l]))
      }
      mat
    })
    pj_g <- pj_gg <- lapply(1:q2, function(l) {
      mat <- matrix(0, G, Kj[l])
      for (g in 1:G) {
        idx <- which(cluster0 == g)
        mat[g, ] <- tabulate(factor(Z2[idx, l], levels = 1:Kj[l]), nbins = Kj[l]) / max(length(idx), 1)
      }
      pmin(pmax(mat, 1e-6), 1 - 1e-6)
    })
    theta0_g <- theta0_gg <- lapply(1:q2, function(l) rep(0, Kj[l]))
    thetaG_g <- thetaG_gg <- lapply(1:q2, function(l) matrix(0, G, Kj[l]))
  }
  # Probabilities
  pi_g <- pi_gg <- as.numeric(table(cluster0) / n)
  wi_g <- wi_gg <- matrix(0, n, G)

  loglik_old <- 1e10
  last_diff <- 0
  for (iter in 1:max_iter) {
    # ---------- E step ----------
    sum_logfij <- matrix(0, nrow = n, ncol = G)
    for (g in 1:G) {
      # Continuous part
      if (q1 > 0) {
        mu_mat <- matrix(mu_g[g, ], nrow = n, ncol = q1, byrow = TRUE)
        sigma_mat <- matrix(Sigma_g[g, ], nrow = n, ncol = q1, byrow = TRUE)
        log_cont <- rowSums(dnorm(Z1, mean = mu_mat, sd = sigma_mat, log = TRUE))
      } else {
        log_cont <- 0
      }
      # Categorical part
      if (q2 > 0) {
        probg <- lapply(1:q2, function(l) pj_g[[l]][g, ])
        cate_mat <- sapply(1:q2, function(j) probg[[j]][Z2[, j]])
        log_cat_prob <- rowSums(log(cate_mat))
      } else {
        log_cat_prob <- 0
      }
      # Log likelihood
      sum_logfij[, g] <- log_cont + log_cat_prob
    }
    max_sumlog <- apply(sum_logfij, 1, max)
    logsum <- exp(sum_logfij - matrix(max_sumlog, nrow = n, ncol = G)) * matrix(pi_g, nrow = n, ncol = G, byrow = TRUE)
    wi_gg <- logsum / rowSums(logsum)

    logf_weighted <- sum_logfij + matrix(log(pi_g), nrow = n, ncol = G, byrow = TRUE)
    max_logf <- apply(logf_weighted, 1, max)
    logf_shifted <- sweep(logf_weighted, 1, max_logf, "-")
    prod_hat <- log(rowSums(exp(logf_shifted))) + max_logf
    loglik_current <- sum(prod_hat)

    # ---------- M step ----------
    clusters <- apply(wi_gg, 1, which.max)
    N_g <- colSums(wi_gg)
    pi_gg <- N_g / n

    # Continuous part
    if (q1 > 0) {
      Sigma_mat <- Sigma_gg
      for (g in 1:G) {
        mu_temp <- colSums(matrix(rep(wi_gg[, g], times = q1), ncol = q1) * Z1) / N_g[g]
        lambda_temp <- lambda1 * Sigma_g[g, ]^2 / N_g[g]
        mu_gg[g, ] <- sapply(1:length(mu_temp), function(k) P_mcp(zj = mu_temp[k], lambda = lambda_temp[k], gamma = 3))
        Sigma_mat[g, ] <- colSums(matrix(rep(wi_gg[, g], times = q1), ncol = q1) * (Z1 - matrix(rep(mu_gg[g, ], each = n), ncol = q1))^2) / n
      }
      Sigma_gg <- matrix(rep(sqrt(colSums(Sigma_mat)), each = G), ncol = q1)
    }

    # Categorical part
    if (q2 > 0) {
      nj_gg <- lapply(1:q2, function(l) {
        mat <- matrix(0, n, Kj[l])
        mat[cbind(1:n, Z2[, l])] <- 1
        t(wi_gg) %*% mat
      })

      pj_freq <- lapply(1:q2, function(l) {
        t(sapply(1:G, function(g) {
          idx <- which(clusters == g)
          tab <- tabulate(factor(Z2[idx, l], levels = 1:Kj[l]), nbins = Kj[l])
          pmin(pmax(tab / max(length(idx), 1), 1e-6), 1 - 1e-6)
        }))
      })

      for (j in 1:q2) {
        theta0_new <- theta0_old <- theta0_g[[j]]
        thetaG_new <- thetaG_old <- thetaG_g[[j]]
        for (sub_iter in 1:inner_iter) {
          if (iter == 1 & sub_iter == 1) {
            pj_old <- pj_freq[[j]]
          } else {
            theta_sum <- sweep(thetaG_old, 2, theta0_old, FUN = "+")
            theta_centered <- theta_sum - apply(theta_sum, 1, max)
            exp_theta <- exp(theta_centered)
            row_sums <- rowSums(exp_theta)
            safe_row_sums <- ifelse(row_sums == 0, 1e-12, row_sums)
            pj_old <- pmin(pmax(exp_theta / safe_row_sums, 1e-3), 1 - 1e-3)
          }
          for (g in 1:G) {
            for (k in 2:Kj[j]) {
              theta0_new[k] <- theta0_old[k] +
                sum(nj_gg[[j]][, k] - N_g * pj_old[, k]) / sum(N_g * pj_old[, k] * (1 - pj_old[, k]))
              tilde_thetaG <- thetaG_old[g, k] +
                (nj_gg[[j]][g, k] - N_g[g] * pj_old[g, k]) / (N_g[g] * pj_old[g, k] * (1 - pj_old[g, k]))
              lambda_temp <- (1 * lambda2) / (Kj[j] * N_g[g] * pj_old[g, k] * (1 - pj_old[g, k]))
              thetaG_new[g, k] <- P_mcp(zj = tilde_thetaG, lambda = lambda_temp, gamma = 3)
            }
          }
          max_diff_0 <- max(abs(theta0_new - theta0_old))
          max_diff_G <- max(abs(thetaG_new - thetaG_old))
          max_val <- suppressWarnings(max(c(max_diff_0, max_diff_G), na.rm = TRUE))
          if (!is.infinite(max_val) && max_val < tol) {
            break
          }
        }
        theta0_gg[[j]] <- theta0_new
        thetaG_gg[[j]] <- thetaG_new
      }
      pj_gg <- lapply(1:q2, function(l) {
        theta_sum <- sweep(thetaG_gg[[l]], 2, theta0_gg[[l]], FUN = "+")
        theta_centered <- theta_sum - apply(theta_sum, 1, max)
        exp_theta <- exp(theta_centered)
        row_sums <- rowSums(exp_theta)
        safe_row_sums <- ifelse(row_sums == 0, 1e-12, row_sums)
        pmin(pmax(exp_theta / safe_row_sums, 1e-3), 1 - 1e-3)
      })
    }


    # ---------- Judgment conditions ----------
    if (any(N_g < 1)) {
      break
    }

    if (q1 > 0) {
      par_cons_g <- c(mu_g, Sigma_g[1, ])
      par_cons_gg <- c(mu_gg, Sigma_gg[1, ])
    } else {
      par_cons_g <- 0
      par_cons_gg <- 0
    }
    if (q2 > 0) {
      par_cate_g <- c(unlist(theta0_g), unlist(thetaG_g))
      par_cate_gg <- c(unlist(theta0_gg), unlist(thetaG_gg))
    } else {
      par_cate_g <- 0
      par_cate_gg <- 0
    }
    par_g <- c(par_cons_g, par_cate_g)
    par_gg <- c(par_cons_gg, par_cate_gg)

    par_diff <- norm(par_gg - par_g, type = "2") / norm(par_g, type = "2")

    current_diff <- abs(loglik_current - loglik_old)
    diff_change <- abs(current_diff - last_diff)
    if (iter > 1 && (current_diff < tol || diff_change < tol || par_diff < tol)) {
      break
    }
    last_diff <- current_diff
    loglik_old <- loglik_current

    if (q1 > 0) {
      mu_g <- mu_gg
      Sigma_g <- Sigma_gg
    }
    if (q2 > 0) {
      pj_g <- pj_gg
      theta0_g <- theta0_gg
      thetaG_g <- thetaG_gg
    }
    pi_g <- pi_gg

    if (q1 > 0) {
      selected_features1 <- which(apply(abs(mu_gg), 2, sum) > 0)
    }
    if (q2 > 0) {
      features2_vec <- sapply(pj_gg, function(mat) {
        all_same <- all(mat == matrix(mat[1, ], nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE))
        as.integer(!all_same)
      })
      selected_features2 <- which(features2_vec == 1)
    }
  }

  # Formatting results
  result <- list(
    G = G, clusters = clusters, lambda1 = lambda1, lambda2 = lambda2,
    iterations = iter, pi_g = pi_gg, wi_g = wi_gg
  )
  if (q1 > 0) {
    result$cons_mu <- mu_gg
    result$cons_Sigma <- Sigma_gg
    result$features1_len <- length(selected_features1)
    result$selected_features1 <- selected_features1
    df_cons <- length(which(as.vector(mu_gg) != 0)) + length(which(as.vector(Sigma_gg[1, ]) != 0))
  } else {
    df_cons <- 0
  }
  if (q2 > 0) {
    result$cate_prob <- pj_gg
    result$theta0 <- theta0_gg
    result$thetaG <- thetaG_gg
    result$features2_len <- length(selected_features2)
    result$selected_features2 <- selected_features2
    df_cate <- length(which(unlist(theta0_gg) != 0)) + length(which(unlist(thetaG_gg) != 0))
  } else {
    df_cate <- 0
  }


  # mBIC calculation (used only w = 0)
  sum_logfij <- matrix(0, nrow = n, ncol = G)
  for (g in 1:G) {
    # Continuous part
    if (q1 > 0) {
      mu_mat <- matrix(mu_gg[g, ], nrow = n, ncol = q1, byrow = TRUE)
      sigma_mat <- matrix(Sigma_gg[g, ], nrow = n, ncol = q1, byrow = TRUE)
      log_cont <- rowSums(dnorm(Z1, mean = mu_mat, sd = sigma_mat, log = TRUE))
    } else {
      log_cont <- 0
    }
    # Categorical part
    if (q2 > 0) {
      probg <- lapply(1:q2, function(l) pj_gg[[l]][g, ])
      cate_mat <- sapply(1:q2, function(j) probg[[j]][Z2[, j]])
      log_cat_prob <- rowSums(log(cate_mat))
    } else {
      log_cat_prob <- 0
    }
    # Log likelihood
    sum_logfij[, g] <- log_cont + log_cat_prob
  }
  logf_weighted <- sum_logfij + matrix(log(pi_gg), nrow = n, ncol = G, byrow = TRUE)
  max_logf <- apply(logf_weighted, 1, max)
  logf_shifted <- sweep(logf_weighted, 1, max_logf, "-")
  prod_hat <- log(rowSums(exp(logf_shifted))) + max_logf

  logL_hat <- sum(prod_hat)
  mBIC <- log(n) * (df_cons + df_cate) - 2 * logL_hat
  result$mBIC <- mBIC

  return(result)
}
