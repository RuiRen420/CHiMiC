##' @title Generate simulated data for \code{\link{CHiMiC}} and \code{\link{CV.CHiMiC}}
##'
##' @description Generate simulated data for applying \code{\link{CHiMiC}} and \code{\link{CV.CHiMiC}}.
##'
##' @param n An integer specifying the number of samples to generate.
##' @param s1 An integer specifying the number of relevant continuous variables.
##' @param s2 An integer specifying the number of relevant categorical variables.
##' @param ratio A numeric value between 0 and 1 specifying the proportion of relevant variables among the total continuous/categorical variables considered.
##'
##' @return A "simData" object containing the following components:
##' \item{Z}{ A simulated matrix of clinical covariates with rows representing samples and columns representing variables (q1 + q2 in total).}
##' \item{q1}{ An integer indicating that the first q1 columns of Z are continuous covariates.}
##' \item{q2}{ An integer indicating that the last q2 columns of Z are categorical covariates.}
##' \item{CATE}{ A simulated numeric vector used to guide clustering.}
##' \item{G}{ Number of clusters.}
##' \item{clusters}{ A numeric vector of cluster assignments.}
##' \item{features_ind}{ Indices of relevant variables in Z.}
##' @seealso See Also as \code{\link{CHiMiC}}, \code{\link{CV.CHiMiC}}.
##'
##' @importFrom MASS mvrnorm
##' @importFrom stats rnorm
##' @export
##' @examples
##' set.seed(123)
##' simdata <- simData(n = 1000, s1 = 20, s2 = 5, ratio = 0.2)
##' Z <- simdata$Z
##' G <- simdata$G
##' q1 <- simdata$q1
##' q2 <- simdata$q2
##' cate <- simdata$CATE
##' hist(cate, breaks = 200)
simData <- function(n, s1 = 20, s2 = 5, ratio = 0.2) {
  # Parameters
  n <- n
  G <- 3
  q1 <- s1 / ratio
  q2 <- s2 / ratio
  xi1 <- 1
  xi2 <- 0.8
  xi3 <- 0.4

  group <- sample(1:3, size = n, replace = TRUE, prob = c(1 / 3, 1 / 3, 1 / 3))

  Tmean <- c(xi1, 2 * xi1, 3 * xi1)
  tau <- numeric(n)
  for (g in 1:3) {
    idx <- which(group == g)
    tau[idx] <- rnorm(n = length(idx), mean = Tmean[g], sd = 0.4)
  }

  # Generate Z Covariates
  # Z1
  true_mu <- matrix(0, nrow = G, ncol = q1)
  true_mu[1, 1:10] <- rep(xi2, 10)
  true_mu[2, 6:15] <- rep(xi2, 10)
  true_mu[3, 11:20] <- rep(xi2, 10)
  Z1 <- matrix(0, nrow = n, ncol = q1)

  # Z2
  prob_vec <- c(0.5 - xi3, 0.5, 0.5 + xi3)
  Z2 <- matrix(0, nrow = n, ncol = s2)
  for (g in 1:3) {
    idx <- which(group == g)
    Z1[idx, ] <- mvrnorm(n = length(idx), mu = true_mu[g, ], Sigma = diag(1, q1))
    Z2[idx, ] <- matrix(rbinom(n = length(idx) * ncol(Z2), size = 1, prob = prob_vec[g]), nrow = length(idx)) + 1
  }
  Z3 <- matrix(rbinom(n = n * (q2 - s2), size = 1, prob = 0.5), nrow = n) + 1
  Z <- cbind(Z1, Z2, Z3)
  colnames(Z) <- 1:(q1 + q2)

  true_ind <- c(rep(1, s1), rep(0, q1 - s1), rep(1, s2), rep(0, q2 - s2))

  return(list(Z = Z, q1 = q1, q2 = q2, CATE = tau, G = 3, clusters = group, features_ind = true_ind))
}
