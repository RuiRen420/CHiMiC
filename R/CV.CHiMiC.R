##' @title Causal High-dimensional Mixed-type Clustering with tuning parameters automatically chosen
##'
##' @description This function performs automatically tuning parameters selection for CHiMiC algorithm.
##'
##' @param Z A matrix of clinical covariates with rows representing samples and columns representing variables (q1 + q2 in total).
##' @param q1 An integer specifying that the first q1 columns of Z are continuous clinical covariates.
##' @param q2 An integer specifying that the last q2 columns of Z are categorical clinical covariates. Each categorical variable should be coded as 1, 2, 3, ... (one integer per category).
##' @param TE A numeric vector of conditional average treatment effects used to guide clustering. Default is NULL; in this case, the method reduces to an unguided clustering approach that uses only the clinical covariates Z.
##' @param ParameterGrid A data.frame specifying the candidate tuning parameters, including:
##' \describe{
##' \item{G:}{ Integer candidates for the number of clusters.}
##' \item{w:}{ Numeric candidate values for the treatment-effect guidance weight.}
##' \item{lambda1:}{ Numeric candidate values for the sparsity penalty parameter on continuous variables.}
##' \item{lambda2:}{ Numeric candidate values for the sparsity penalty parameter on categorical variables.}
##' }
##' @param nfolds An integer (>= 2) specifying the number of folds. Default is 3.
##' @param max_iter An integer specifying the maximum number of iterations for clustering. Default is 100.
##' @param inner_iter An integer specifying the maximum number of inner iterations for estimating parameters in the categorical-data subroutine. Default is 30.
##' @param tol A numeric value specifying the convergence tolerance.
##'
##' @return A "CV.CHiMiC" object containing the following components:
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
##' @seealso See Also as \code{\link{simData}}, \code{\link{CHiMiC}}.
##'
##' @import MASS
##' @import stats
##' @import caret
##' @import sparcl
##' @import data.table
##' @export
##' @examples
##' set.seed(123)
##' simdata <- simData(n = 500, s1 = 20, s2 = 5, ratio = 0.2)
##' Z <- simdata$Z
##' G <- simdata$G
##' q1 <- simdata$q1
##' q2 <- simdata$q2
##' cate <- simdata$CATE
##' nfolds <- 3
##' ParameterGrid <- expand.grid(G = 3, w = c(0, 0.5), lambda1 = 20, lambda2 = 15)
##'
##' # For demonstration purposes, max_iter and inner_iter are set to small values.
##' # In practice, larger iterations (e.g., 100) are needed for precise estimation.
##' cv.chimic.obj <- CV.CHiMiC(Z, q1, q2, TE = cate, ParameterGrid, nfolds,
##'                            max_iter = 5, inner_iter = 5, tol = 1e-2)
##'

CV.CHiMiC <- function(Z, q1, q2, TE, ParameterGrid, nfolds = 3, max_iter = 100, inner_iter = 30, tol = 1e-2) {
  .get_col <- function(df, nm) {
    pos <- which(names(df) == nm)
    if (length(pos) == 0) stop("Column '", nm, "' not found.")
    x <- df[[pos[1]]] # if duplicated, use the first
    if (is.data.frame(x) || is.matrix(x)) x <- x[, 1]
    x <- as.vector(x)
    if (length(x) != nrow(df)) {
      stop(
        "Column '", nm, "' has length ", length(x),
        " but nrow(df) = ", nrow(df), ". Likely duplicated/malformed column."
      )
    }
    x
  }
  ## ---------------------------
  ## 0) Prepare ParameterGrid (data.frame + numeric conversion where needed)
  ## ---------------------------
  ParameterGrid0 <- as.data.frame(ParameterGrid)

  for (nm in intersect(c("w", "G", "lambda1", "lambda2"), names(ParameterGrid0))) {
    x <- ParameterGrid0[[nm]]
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(ParameterGrid0[[nm]] <- as.numeric(x))
  }

  w_vec <- .get_col(ParameterGrid0, "w")

  key_cols <- c("G", "lambda1", "lambda2")
  missing_keys <- setdiff(key_cols, names(ParameterGrid0))
  if (length(missing_keys) > 0) stop("Missing columns in ParameterGrid: ", paste(missing_keys, collapse = ", "))

  ## ---------------------------
  ## 1) Identify which w have >1 unique (G,lambda1,lambda2)
  ## ---------------------------
  idx_by_w <- split(seq_len(nrow(ParameterGrid0)), w_vec)

  n_combo_per_w <- vapply(idx_by_w, function(ii) {
    nrow(unique(ParameterGrid0[ii, key_cols, drop = FALSE]))
  }, integer(1))

  w_multi <- names(n_combo_per_w)[n_combo_per_w > 1]
  w_single <- names(n_combo_per_w)[n_combo_per_w == 1]

  sub2 <- ParameterGrid0[w_vec %in% as.numeric(w_single), , drop = FALSE]
  if (nrow(sub2) > 0) {
    w2 <- .get_col(sub2, "w")
    keep2 <- unlist(tapply(seq_len(nrow(sub2)), w2, function(ii) ii[1]), use.names = FALSE)
    sub2 <- sub2[keep2, , drop = FALSE]
  }

  ## ---------------------------
  ## 2) sub1: for w with multiple combos, run mBIC_CHiMiC on sub1 only
  ## ---------------------------
  if (length(w_multi) > 0) {
    sub1 <- ParameterGrid0[w_vec %in% as.numeric(w_multi), , drop = FALSE]

    ParameterSet0 <- mBIC_CHiMiC(
      Z = Z, q1 = q1, q2 = q2, TE = TE,
      ParameterSet = sub1,
      max_iter = max_iter, inner_iter = inner_iter, tol = tol
    )
    ParameterSet0 <- as.data.frame(ParameterSet0)

    ParameterSet0[] <- lapply(ParameterSet0, function(x) {
      if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
    })

    idx <- unlist(
      tapply(seq_len(nrow(ParameterSet0)), ParameterSet0$w, function(ii) {
        ii[which.min(ParameterSet0$mBIC[ii])]
      }),
      use.names = FALSE
    )

    sub1_best <- ParameterSet0[idx, , drop = FALSE]
    row.names(sub1_best) <- NULL
  } else {
    sub1_best <- ParameterGrid0[0, , drop = FALSE] # empty
  }

  ## ---------------------------
  ## 3) Combine: sub2 + selected sub1_best  => ParameterSet
  ## ---------------------------
  needed_cols <- c("w", "G", "lambda1", "lambda2")
  ParameterSet <- rbind(sub2[, needed_cols, drop = FALSE], sub1_best[, needed_cols, drop = FALSE])
  row.names(ParameterSet) <- NULL


  ## ---------------------------
  ## 4) Final CV tuning
  ## ---------------------------
  CV_CHiMiC(
    Z = Z, q1 = q1, q2 = q2, TE = TE, ParameterSet = ParameterSet, nfolds = nfolds,
    max_iter = max_iter, inner_iter = inner_iter, tol = tol
  )
}
