#' Build the design matrix and the penalty matrix based on a formula and a data set
#'
#' @param formula right side of a formula as used in mgcv
#' @param data data frame containing the variables in the formula
#' @param knots optional list containing user specified knot values to be used for basis construction. 
#' For most bases the user simply supplies the knots to be used, which must match up with the k value supplied (note that the number of knots is not always just k).
#' See mgcv documentation for more details.
#'
#' @return a list containing the design matrix Z, the penalty matrix S, the formula, the data and the knots
#' @export
#' 
#' @import mgcv
#' @importFrom stats update predict
#'
#' @examples
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
make_matrices = function(formula, data, knots = NULL){
  gam_setup = gam(formula = update(formula, dummy ~ .),
                  data = cbind(dummy = 1, data), 
                  knots = knots,
                  fit = FALSE)
  
  Z = gam_setup$X
  S = gam_setup$S
  formula = gam_setup$formula
  
  return(list(Z = Z, S = S, formula = formula, data = data, knots = knots))
}

#' Build the prediction design matrix based on new data and model_matrices object created by \code{make_matrices()}
#'
#' @param model_matrices model_matrices object as returned from \code{make_matrices()}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#'
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#' @export
#' 
#' @import mgcv
#' @importFrom stats update predict
#'
#' @examples
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
#' pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
pred_matrix = function(model_matrices, newdata) {
  gam_setup0 = gam(model_matrices$form, 
                   data = cbind(dummy = 1, model_matrices$data),
                   knots = model_matrices$knots)
  
  predict(gam_setup0, newdata = cbind(dummy = 1, newdata), type = "lpmatrix")
}


#' Build a standardized P-Spline design matrix and the associated P-Spline penalty matrix
#' 
#' This function builds the B-spline design matrix for a given data vector x. 
#' The B-spline basis functions are normalized such that the integral of each basis function is 1, hence this basis can be used for spline-based density estimation.
#'
#' @param x data vector
#' @param K number of basis functions
#' @param degree degree of the B-spline basis functions, defaults to cubic B-Splines
#' @param npoints number of points used in the numerical integration for normalizing the B-spline basis functions
#' @param diff_order order of the difference used for the penalty matrix. Defaults to second-order differences.
#'
#' @return a list containing the design matrix Z, the penalty matrix S, the knots, the normalization weights w, the degree of the B-spline basis functions and the basis positions
#' @export
#'
#' @examples
#' modmat = make_matrices_dens(x = 1:100, K = 20)
make_matrices_dens = function(x, K, degree = 3, npoints = 1e4, diff_order = 2){
  
  ## building the design matrix
  rangex = range(x, na.rm = TRUE)
  nObs = length(x)
  ord = degree + 1
  nrknots = K - (degree-1) 
  d = (rangex[2] - rangex[1]) / nrknots
  bm = c(rangex[1] - degree*d, rangex[2] + degree*d)
  knots = seq(bm[1], bm[2], length = nrknots + 2*degree)
  
  # numerical integration for normalizing the B-spline basis functions
  xseq =  seq(bm[1], bm[2], length = npoints)
  B0 = splines::spline.des(knots, xseq, degree+1, outer.ok=T)$design # unnormalized
  w = rep(NA, K)
  h = diff(c(knots[1], knots[length(knots)])) / npoints
  for (k in 1:K){
    w[k] = (h* sum(B0[,k]))^(-1) 
    # this computes the integrals of the B-spline basis functions (which are then standardized below)
  } 
  
  # actual data design matrix
  B = matrix(NA, nrow = nObs, ncol = K)
  ind = which(!is.na(x))
  B[ind,] = t(t(splines::spline.des(knots, x[ind], degree+1, outer.ok = TRUE)$design) * w) 
  
  # basis positions
  basis_pos = knots[(degree+1):(length(knots)-degree+1)]
  
  ## building the penalty matrix
  L = diff(diag(K), differences = diff_order) # second-order difference matrix
  S = t(L[,-1])%*%L[,-1] # leaving out first column
  cat("Leaving out first column of S, fix first column of parameter matrix at zero!")
  
  list(Z=B, S = S, knots=knots, w=w, degree = degree, basis_pos = basis_pos)
}

#' Build the prediction design matrix for a new data vector based on model_matrices object created by \code{make_matrices_dens()}
#'
#' @param model_matrices model_matrices object as returned from \code{make_matrices_dens()}
#' @param xnew new data vector for which to evaluate the basis
#'
#' @return prediction design matrix for \code{xnew} with the same basis as used for \code{model_matrices}
#' @export
#'
#' @examples
#' modmat = make_matrices_dens(x = 1:100, K = 20)
#' pred_matrix_dens(modmat, xnew = 1:100 - 0.5)
pred_matrix_dens = function(model_matrices, xnew){
  knots = model_matrices$knots
  degree = model_matrices$degree
  w = model_matrices$w
  
  B = splines::spline.des(knots, xnew, degree+1, outer.ok=T)$design
  sweep(B, 2, w, FUN = "*")
}