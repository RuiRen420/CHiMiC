##' @title Causal High-dimensional Mixed-type Clustering
##'
##' @description This function performs clustering on mixed-type clinical data Z (including continuous and categorical variables) guided by treatment effects (TE). It also identifies relevant clinical variables using an MCP (minimax concave penalty)–based variable selection method.
##'
##' @param Z A matrix of clinical covariates with rows representing samples and columns representing variables (q1 + q2 in total).
##' @param q1 An integer specifying that the first q1 columns of Z are continuous clinical covariates.
##' @param q2 An integer specifying that the last q2 columns of Z are categorical clinical covariates. Each categorical variable should be coded as 1, 2, 3, ... (one integer per category).
##' @param TE A numeric vector of conditional average treatment effects used to guide clustering. Default is NULL; in this case, the method reduces to an unguided clustering approach that uses only the clinical covariates Z.
##' @param w A numeric value specifying the TE guidance weight. Default is 0; in this case, the method reduces to an unguided clustering approach that uses only the clinical covariates Z.
##' @param G An integer (>= 2) specifying the number of clusters.
##' @param lambda1 A numeric value specifying the sparsity penalty parameter for continuous variables.
##' @param lambda2 A numeric value specifying the sparsity penalty parameter for categorical variables.
##' @param max_iter An integer specifying the maximum number of iterations for clustering. Default is 100.
##' @param inner_iter An integer specifying the maximum number of inner iterations for estimating parameters in the categorical-data subroutine. Default is 30.
##' @param tol A numeric value specifying the convergence tolerance.
##'
##' @return A "CHiMiC" object containing the following components:
##'   \item{G}{Number of clusters.}
##'   \item{clusters}{A numeric vector of cluster assignments.}
##'   \item{lambda1}{Sparsity penalty parameter for continuous variables.}
##'   \item{lambda2}{Sparsity penalty parameter for categorical variables.}
##'   \item{w}{TE guidance weight.}
##'   \item{iterations}{Total number of iterations.}
##'   \item{pi}{A numeric vector of estimated mixing proportions.}
##'   \item{cons_mu}{Estimated means for continuous variables.}
##'   \item{cons_Sigma}{Estimated diagonal elements of the covariance matrix for continuous variables.}
##'   \item{features_len1}{Number of selected (relevant) continuous variables.}
##'   \item{selected_features1}{Indices of selected continuous variables.}
##'   \item{cate_prob}{A list of estimated categorical probabilities. Each element corresponds to one categorical variable and contains a matrix with rows representing clusters and columns representing category probabilities.}
##'   \item{theta0}{A list of estimated \eqn{\theta_0} values for categorical variables. Each element corresponds to one categorical variable and is a vector as defined in the algorithm.}
##'   \item{thetaG}{A list of estimated \eqn{\theta_G} values for categorical variables. Each element corresponds to one categorical variable and contains a matrix with rows representing clusters and columns representing the estimated \eqn{\theta_G} values as defined in the algorithm.}
##'   \item{features_len2}{Number of selected (relevant) categorical variables.}
##'   \item{selected_features2}{Indices of selected categorical variables.}
##'   \item{te_beta}{Estimated cluster-specific means for the treatment-effect guidance component.}
##'   \item{te_sigma}{Estimated cluster-specific standard deviations for the treatment-effect guidance component.}
##'   \item{mBIC}{Modified BIC computed using only the clinical covariates \eqn{Z}.}
##' @seealso See Also as \code{\link{simData}}.
##'
##' @import MASS
##' @import stats
##' @import caret
##' @import sparcl
##' @export
##' @examples
##' set.seed(123)
##' simdata <- simData(n = 500, s1 = 20, s2 = 5, ratio = 0.2)
##' Z <- simdata$Z
##' G <- simdata$G
##' q1 <- simdata$q1
##' q2 <- simdata$q2
##' cate <- simdata$CATE
##' w <- 0.5
##' lambda1 <- 20
##' lambda2 <- 15
##'
##' # For demonstration purposes, max_iter and inner_iter are set to small values.
##' # In practice, larger iterations (e.g., 100) are needed for precise estimation.
##' chimic_obj <- CHiMiC(Z, q1, q2, TE = cate, w, G, lambda1, lambda2,
##'                      max_iter = 5, inner_iter = 5, tol = 1e-2)
##'
CHiMiC <- function(Z, q1, q2, TE = NULL, w = 0, G, lambda1, lambda2,
                   max_iter = 100, inner_iter = 30, tol = 1e-2) {
  Y <- TE
  n <- nrow(Z)
  q <- ncol(Z)
  if (q1 > 0) {
    Z1 <- scale(Z[, 1:q1], center = TRUE, scale = TRUE)
  }
  if (q2 > 0) {
    Z2 <- Z[, (q - q2 + 1):q]
  }
  if (!is.null(Y)) {
    Y <- scale(Y)
  }

  if (is.null(Y)) {
    res.init <- HiMiC(
      Z = Z, q1 = q1, q2 = q2, G = G, lambda1 = lambda1, lambda2 = lambda2,
      max_iter = max_iter, inner_iter = inner_iter, tol = tol
    )
    result <- res.init
  } else {
    if (w == 0) {
      Zc <- Z
      q1c <- q1
    } else {
      Zc <- cbind(Y, Z)
      q1c <- q1 + 1
    }
    res.init <- HiMiC(
      Z = Zc, q1 = q1c, q2 = q2, G = G, lambda1 = lambda1, lambda2 = lambda2,
      max_iter = max_iter, inner_iter = inner_iter, tol = tol
    )
    cluster0 <- res.init$clusters

    # Initial values
    # Continuous part
    if (q1 > 0) {
      if (w == 0) {
        mu_g <- mu_gg <- res.init$cons_mu
        Sigma_g <- Sigma_gg <- Sigma_mat <- res.init$cons_Sigma
      } else {
        mu_g <- mu_gg <- res.init$cons_mu[, -1]
        Sigma_g <- Sigma_gg <- Sigma_mat <- res.init$cons_Sigma[, -1]
      }
    }
    # Categorical part
    if (q2 > 0) {
      Kj <- apply(Z2, 2, max)
      pj_g <- pj_gg <- res.init$cate_prob
      theta0_g <- theta0_gg <- res.init$theta0
      thetaG_g <- thetaG_gg <- res.init$thetaG
    }
    # Probabilities
    pi_g <- pi_gg <- res.init$pi_g
    wi_g <- wi_gg <- res.init$wi_g
    # CATE part
    beta_g <- beta_gg <- rep(0, G)
    sigma_g <- sigma_gg <- rep(1, G)
    for (g in 1:G) {
      beta_g[g] <- beta_gg[g] <- mean(Y[which(cluster0 == g)])
      sigma_g[g] <- sigma_gg[g] <- sd(Y[which(cluster0 == g)])
    }

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
        # CATE part
        log_y <- dnorm(Y, mean = beta_g[g], sd = sigma_g[g], log = TRUE)
        # Log likelihood
        sum_logfij[, g] <- (1 - w) * (log_cont + log_cat_prob) + w * log_y
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
      if (any(!is.finite(N_g)) || any(N_g < max(1, 0.01 * n))) break

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
            pmin(pmax(tab / max(length(idx), 1), 1e-3), 1 - 1e-3)
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

      # Outcome part
      for (g in 1:G) {
        beta_gg[g] <- sum(wi_gg[, g] * Y) / sum(wi_gg[, g])
        sigma_gg[g] <- sqrt(sum(wi_gg[, g] * (Y - beta_gg[g])^2) / sum(wi_gg[, g]))
      }

      # ---------- Judgment conditions ----------
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
      par_w0_g <- c(par_cons_g, par_cate_g)
      par_w0_gg <- c(par_cons_gg, par_cate_gg)
      par_w1_g <- c(beta_g, sigma_g)
      par_w1_gg <- c(beta_gg, sigma_gg)

      par_g <- if (w == 1) par_w1_g else if (w == 0) par_w0_g else c(par_w0_g, par_w1_g)
      par_gg <- if (w == 1) par_w1_gg else if (w == 0) par_w0_gg else c(par_w0_gg, par_w1_gg)
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
      beta_g <- beta_gg
      sigma_g <- sigma_gg
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
      w = w, iterations = iter, pi = pi_gg
    )
    if (q1 > 0) {
      result$cons_mu <- mu_gg
      result$cons_Sigma <- Sigma_gg[1, ]
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

    result$te_beta <- beta_gg
    result$te_sigma <- sigma_gg

    # mBIC calculation
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
      # CATE part
      log_y <- dnorm(Y, mean = beta_g[g], sd = sigma_g[g], log = TRUE)
      # Log likelihood
      sum_logfij[, g] <- (log_cont + log_cat_prob)
    }
    logf_weighted <- sum_logfij + matrix(log(pi_gg), nrow = n, ncol = G, byrow = TRUE)
    max_logf <- apply(logf_weighted, 1, max)
    logf_shifted <- sweep(logf_weighted, 1, max_logf, "-")
    prod_hat <- log(rowSums(exp(logf_shifted))) + max_logf

    logL_hat <- sum(prod_hat)
    mBIC <- log(n) * (df_cons + df_cate) - 2 * logL_hat
    result$mBIC <- mBIC
  }
  return(result)
}
