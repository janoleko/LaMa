
# Regression setting ------------------------------------------------------

# hidden helper function 
# -> just to turn cosinor(x, period) terms into sine/ cosine terms
process_cosinor <- function(formula){
  # Extract formula terms
  Terms <- stats::terms(formula, specials = "cosinor")
  term_names <- attr(Terms, "term.labels")
  var_names <- all.vars(formula)
  
  # Identify cosinor terms
  cosInd <- grep("cosinor\\(", term_names)
  
  # Separate regular terms
  mainpart <- if(length(cosInd) > 0) term_names[-cosInd] else term_names
  
  # Expand cosinor terms while preserving interactions
  expanded_terms <- c()
  for(name in term_names[cosInd]){
    # Extract the cosinor(...) call using regex
    match <- regmatches(name, regexpr("cosinor\\(.*\\)", name))
    
    # Check if ", eval = TRUE" is present
    if (grepl(", eval = TRUE", match)) {
      # Replace ", eval = TRUE" with ", eval = FALSE"
      new_match <- sub(", eval = TRUE", ", eval = FALSE", match)
    } else if (!grepl(", eval = FALSE", match)) {
      # If ", eval = FALSE" is not present, insert it
      new_match <- sub("(cosinor\\(.*)(\\))", "\\1, eval = FALSE\\2", match)
    } else {
      # If ", eval = FALSE" is already present, keep it as is
      new_match <- match
    }
    
    # Evaluate the expression using the data environment
    replace <- eval(parse(text = new_match), envir = list(cosinor = cosinor))
    
    # Preserve interactions (all are converted to ":" by terms())
    if (grepl(":", name)) {
      parts <- unlist(strsplit(name, ":"))  # Split interaction terms
      other_factors <- parts[parts != match]  # Extract non-cosinor parts
      
      # Rebuild interactions with expanded terms
      expanded_terms <- c(expanded_terms, sapply(replace, function(rt) paste(c(other_factors, rt), collapse = ":")))
    } else {
      expanded_terms <- c(expanded_terms, replace)
    }
  }
  
  # Final expanded formula
  final_terms <- c(mainpart, expanded_terms)
  expanded_formula <- stats::as.formula(paste("~", paste(final_terms, collapse = " + ")))
  
  expanded_formula
}

#' Evaluate trigonometric basis expansion
#' 
#' This function can be used to evaluate a trigonometric basis expansion for a given periodic variable and period. 
#' It can also be used in formulas passed to \code{\link{make_matrices}}.
#' 
#' The returned basis can be used for linear predictors of the form
#' \deqn{ 
#'  \eta^{(t)} = \beta_0 + \sum_{k} \bigl( \beta_{1k} \sin(\frac{2 \pi t}{period_k}) + \beta_{2k} \cos(\frac{2 \pi t}{period_k}) \bigr). 
#' }
#' This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing the length of \code{period}).
#'  
#' @param x vector of periodic variable values
#' @param period vector of period length. For example for time of day \code{period = 24}, or \code{period = c(24,12)} for more flexibility.
#' @param eval logical, should not be changed. If \code{TRUE} the function returns the evaluated cosinor terms, if \code{FALSE} the function returns the terms as strings which is used internally form formula evaluation.
#'
#' @return either a desing matrix with the evaluated cosinor terms (\code{eval = TRUE}) or a character vector with the terms as strings (\code{eval = FALSE}).
#' @export
#'
#' @examples
#' ## Evaluate cosinor terms
#' # builds design matrix
#' X = cosinor(1:24, period = 24)
#' X = cosinor(1:24, period = c(24, 12, 6))
#' 
#' ## Usage in model formulas
#' # e.g. frequencies of 24 and 12 hours + interaction with temperature
#' form = ~ x + temp * cosinor(hour, c(24, 12)) 
#' data = data.frame(x = runif(24), temp = rnorm(24,20), hour = 1:24)
#' modmat = make_matrices(form, data = data)
cosinor = function(x = 1:24, period = 24, eval = TRUE){
  # get the name of input varible
  xname = deparse(substitute(x))
  
  if(eval == FALSE){
    out = c()
    # Loop over periods and construct sine and cosine strings
    for(p in period){
      out = c(out,
              paste0("sin(2*pi*", xname, "/", p, ")"),
              paste0("cos(2*pi*", xname, "/", p, ")"))
    }
    return(out)
  } else{
    # Evaluate the cosinor terms
    # x might be a vector
    out = matrix(NA, nrow = length(x), ncol = 0)
    names = c()
    for(p in period){
      out = cbind(out,
                  sin(2*pi*x/p),
                  cos(2*pi*x/p))
      
      names = c(names,
                paste0("sin(2*pi*", xname, "/", p, ")"),
                paste0("cos(2*pi*", xname, "/", p, ")"))
    }
    colnames(out) = names
    return(out)
  }
}


#' Build the design and the penalty matrix for models involving penalised splines based on a formula and a data set
#'
#' @param formula right side of a formula as used in \code{mgcv}
#' @param data data frame containing the variables in the formula
#' @param knots optional list containing user specified knot values to be used for basis construction
#' 
#' For most bases the user simply supplies the \code{knots} to be used, which must match up with the \code{k} value supplied (note that the number of knots is not always just \code{k}).
#' See \code{mgcv} documentation for more details.
#'
#' @return a list containing the design matrix \code{Z}, a (potentially nested) list of penalty matrices \code{S}, the \code{formula}, the \code{data}, the \code{knots}, and the original \code{mod} object returned by \code{mgcv}.
#' Note that for tensorproduct smooths, the corresponding list entry is itself a list, containing the d marginal penalty matrices if d is the dimension of the tensor product.
#' @export
#' 
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom stats update
#'
#' @examples
#' data = data.frame(x = runif(100), 
#'                   y = runif(100),
#'                   g = factor(rep(1:10, each = 10)))
#'
#' # unvariate thin plate regression spline
#' modmat = make_matrices(~ s(x), data)
#' # univariate P-spline
#' modmat = make_matrices(~ s(x, bs = "ps"), data)
#' # adding random intercept
#' modmat = make_matrices(~ s(g, bs = "re") + s(x, bs = "ps"), data)
#' # tensorproduct of x and y
#' modmat = make_matrices(~ s(x) + s(y) + ti(x,y), data)
make_matrices = function(formula, 
                         data, 
                         knots = NULL
                         ){
  
  ## Potenially expand cosinor terms
  formula = process_cosinor(formula)
  
  ## setting up the model using mgcv
  gam_setup = mgcv::gam(formula = stats::update(formula, dummy ~ .),
                        data = cbind(dummy = 1, data), 
                        knots = knots,
                        fit = FALSE)
  
  
  ## assiging design matrix
  term_names = gam_setup$term.names
  Z = gam_setup$X
  colnames(Z) = term_names

  
  ## dealing with the penalty matrices
  term_labels = sapply(gam_setup$smooth, function(x) x$label)
  
  # first option: tensorproducts -> save marginal penalty matrices
  S = lapply(gam_setup$smooth, function(x){
    if(is.null(x$margin)){ # univariate smooth -> one penalty matrix
      return(x$S[[1]])
    } else{ # multivariate smooth -> several penalty matrices
      S_sublist = lapply(x$margin, function(y) y$S[[1]])
      margin_names = sapply(x$margin, function(y) y$term)
      names(S_sublist) = margin_names
      return(S_sublist)
    }
  })
  names(S) = term_labels
  
  # second option: tensorproduct -> save blown-up marginal penalty matrices (with constraints baked in)
  S2 = list()
  counter = 1
  for(i in seq_along(gam_setup$smooth)){
    sm = gam_setup$smooth[[i]]
    if(is.null(sm$margin)){
      S2[[i]] = gam_setup$S[[counter]]
      counter = counter + 1
    } else{
      nPenMat = length(sm$margin)
      S_sublist = gam_setup$S[counter:(counter + nPenMat - 1)]
      margin_names = sapply(sm$margin, function(y) y$term)
      names(S_sublist) = margin_names
      S2[[i]] = S_sublist
      counter = counter + nPenMat
    }
  }
  names(S2) = term_labels
  
  pardim <- list(fixed_eff = gam_setup$nsdf)
  
  pardim_smooth = sapply(S2, function(x){
    if(is.matrix(x)){
      return(nrow(x))
    } else{
      return(nrow(x[[1]]))
    }
  })
  
  pardim = c(pardim, pardim_smooth)
  
  return(list(Z = Z, 
              S = S, 
              S2 = S2,
              pardim = pardim,
              formula = gam_setup$formula, 
              data = data, 
              knots = knots,
              mod = gam_setup))
}


#' Build the prediction design matrix based on new data and model_matrices object created by \code{\link{make_matrices}}
#'
#' @param model_matrices model_matrices object as returned from \code{\link{make_matrices}}
#' @param newdata data frame containing the variables in the formula and new data for which to evaluate the basis
#' @param exclude optional vector of terms to set to zero in the predicted design matrix. Useful for predicting main effects only when e.g. \code{sd(..., bs = "re")} terms are present. See \code{mgcv::predict.gam} for more details.
#' @return prediction design matrix for \code{newdata} with the same basis as used for \code{model_matrices}
#' @export
#' 
#' @importFrom mgcv gam
#' @importFrom mgcv s
#' @importFrom mgcv predict.gam
#'
#' @examples
#' modmat = make_matrices(~ s(x), data.frame(x = 1:10))
#' Z_predict = pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
pred_matrix = function(model_matrices, 
                       newdata,
                       exclude = NULL) {
  gam_setup0 = mgcv::gam(model_matrices$formula, 
                         data = cbind(dummy = 1, model_matrices$data),
                         knots = model_matrices$knots)
  
  predict.gam(gam_setup0, 
              newdata = cbind(dummy = 1, newdata), 
              type = "lpmatrix",
              exclude = exclude)
}



# Density estimation setting ----------------------------------------------


#' Build a standardised P-Spline design matrix and the associated P-Spline penalty matrix
#' 
#' This function builds the B-spline design matrix for a given data vector. 
#' Importantly, the B-spline basis functions are normalised such that the integral of each basis function is 1, hence this basis can be used for spline-based density estimation, when the basis functions are weighted by non-negative weights summing to one.
#'
#' @param x data vector
#' @param k number of basis functions
#' @param type type of the data, either \code{"real"} for data on the reals, \code{"positive"} for data on the positive reals or \code{"circular"} for circular data like angles.
#' @param degree degree of the B-spline basis functions, defaults to cubic B-splines
#' @param npoints number of points used in the numerical integration for normalizing the B-spline basis functions
#' @param diff_order order of differencing used for the P-Spline penalty matrix for each data stream. Defaults to second-order differences.
#' @param pow power for polynomial knot spacing
#' 
#' Such non-equidistant knot spacing is only used for \code{type = "positive"}.
#'
#' @return list containing the design matrix \code{Z}, the penalty matrix \code{S}, the prediction design matrix \code{Z_predict}, the prediction grid \code{xseq}, and details for the basis expansion.
#' @export
#' @importFrom splines spline.des
#' @importFrom splines bs
#' @importFrom mgcv cSplineDes
#'
#' @examples
#' modmat = make_matrices_dens(x = (-50):50, k = 20)
#' modmat = make_matrices_dens(x = 1:100, k = 20, type = "positive")
#' modmat = make_matrices_dens(x = seq(-pi,pi), k = 20, type = "circular")
make_matrices_dens = function(x, # data vector
                              k, # number of basis functions
                              type = "real", # type of the data
                              degree = 3, # degree of the B-Spline basis
                              npoints = 1e4, # number of points for numerical integration
                              diff_order = 2, # order of the differences for the penalty matrix
                              pow = 0.5 # power for polynomial knot spacing for positive values
){
  nObs = length(x)
  ord = degree + 1
  
  if(type != "circular"){
    ## building the design matrix
    rangex = range(x, na.rm = TRUE)
    nrknots = k - (degree-1) 
    
    if(type == "real"){ # equidistant knots, range in reals
      d = (rangex[2] - rangex[1]) / nrknots
      bm = c(rangex[1] - degree*d, rangex[2] + degree*d)
      knots = seq(bm[1], bm[2], length = nrknots + 2*degree)
      
      # numerical integration for normalizing the B-spline basis functions
      xseq = seq(bm[1], bm[2], length = npoints)
      B0 = splines::spline.des(knots, xseq, degree+1, outer.ok = TRUE)$design # unnormalized
      w = rep(NA, k)
      h = diff(c(knots[1], knots[length(knots)])) / npoints
      for (i in 1:k){
        w[i] = (h* sum(B0[,i]))^(-1) 
        # this computes the integrals of the B-spline basis functions (which are then standardized below)
      } 
      
      # rescaling prediction matrix
      B0 = t(t(B0)*w)
      
      # actual data design matrix
      B = matrix(NA, nrow = nObs, ncol = k)
      ind = which(!is.na(x))
      B[ind,] = t(t(splines::spline.des(knots, x[ind], degree+1, outer.ok = TRUE)$design) * w) 
      
      # basis positions
      basis_pos = knots[(degree):(length(knots)-degree+1)]
      
      ## building the penalty matrix
      L = diff(diag(k), differences = diff_order) # second-order difference matrix
      
    } else if(type == "positive") { # non-equidistant knots, no mass on < 0
      if(min(rangex) <= 0) stop("When positive = TRUE, x can only contain positive values")
      
      # square-root spacing
      xseq = seq(0, max(x, na.rm = TRUE) * 1.05, length = npoints)
      knots = seq(min(x, na.rm = TRUE)^pow, max(x, na.rm = TRUE)^pow, length = k - degree)^(1/pow)
      
      B0 = splines::bs(x = xseq, 
                       knots = knots,
                       degree = degree) # spline design matrix
      
      boundary_knots = attr(B0, "Boundary.knots")
      allknots = c(boundary_knots[1], knots, boundary_knots[2])
      
      w = rep(NA, ncol(B0))
      h = (max(xseq) - min(xseq)) / npoints
      for (i in 1:(ncol(B0))){
        w[i] = (h * sum(B0[,i]))^(-1)
        # this computes the integrals of the B-spline basis functions (which are then standardized below)
      }
      B0 = t(t(B0)*w)
      B0[is.nan(B0)] = 0
      
      # numerically computing the average for each basis function
      basis_pos = colSums(xseq * t(t(B0)/rowSums(t(B0))))
      
      B = matrix(NA, nrow = nObs, ncol = k)
      ind = which(!is.na(x))
      B[ind,] = splines::bs(x = x[ind], 
                            knots = knots,
                            degree = degree) # spline design matrix
      B = t(t(B)*w)
      
      ## building the penalty matrix
      L = diff(diag(k), differences = diff_order) # second-order difference matrix
    } else {
      stop("type must be on of real, continuous or circular")
    }
    
  } else if(type == "circular"){
    
    knots = seq(-pi, pi, length = k+1)
    xseq = seq(-pi, pi, length = npoints)
    B0 = mgcv::cSplineDes(xseq, knots, ord = ord)
    # numerical integration
    w = rep(NA, k)
    h = 2* pi / npoints
    for (i in 1:k){
      w[i] = (h* sum(B0[,i]))^(-1) 
      # this computes the integrals of the B-spline basis functions (which are then standardized below)
    } 
    # rescaling prediction matrix
    B0 = t(t(B0)*w)
    
    # actual data design matrix
    B = matrix(NA, nrow = nObs, ncol = k)
    ind = which(!is.na(x))
    B[ind,] = t(t(mgcv::cSplineDes(x[ind], knots, ord = ord)) * w) 
    
    # circular P-Spline penalty matrix
    L = diff(rbind(diag(k), diag(k)[1:diff_order,]), differences = diff_order) # second-order difference matrix
    basis_pos = knots[c(k, 1:(k-1))]
  }
  
  # constructing penalty matrix
  S = t(L[,-k]) %*% L[,-k] # leaving out first column because parameter set to zero
  cat("Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!\n")
  
  basis = list(type = type, knots = knots, w = w, degree = degree, basis_pos = basis_pos)
  list(Z = B, 
       S = S, 
       Z_predict = B0[(1:500) * (npoints/500) -(npoints/500)/2,], 
       xseq = xseq[(1:500)* (npoints/500) -(npoints/500)/2], 
       basis = basis)
}

# helper function, not exported
make_splinecoef = function(model_matrices, 
                           type = "real", 
                           par){
  basis_pos = model_matrices$basis$basis_pos
  k = length(basis_pos)
  
  if(type == "real"){ # if density has support on the reals -> use normal distribution to initialize
    beta = sapply(basis_pos[-k], dnorm, mean = par$mean, sd = par$sd, log = TRUE)
  } else if(type == "positive") { # if density has support on positive continuous -> use gamma distribution
    # transformation to scale and shape
    shape = par$mean^2 / par$sd^2
    scale = par$sd^2 / par$mean
    beta = sapply(basis_pos[-k], dgamma, shape = shape, scale = scale, log = TRUE)
    # rescaling to account for non-equidistant knot spacing
    beta = beta - log(apply(model_matrices$Z, 2, max)[-k])
  } else if(type == "circular") {
    beta = sapply(basis_pos[-k], LaMa::dvm, mu = par$mean, kappa = par$concentration, log = TRUE)
  }
  
  if(is.vector(beta)){
    beta = matrix(beta, nrow = 1, ncol = length(beta))
  }
  beta = beta - beta[,k-1]
  cat("Parameter matrix excludes the last column. Fix this column at zero!\n")
  return(beta)
}

#' Build the design and penalty matrices for smooth density estimation
#' 
#' @description
#' This high-level function can be used to prepare objects needed to estimate mixture models of smooth densities using P-Splines.
#' 
#' @details
#' Under the hood, \code{\link{make_matrices_dens}} is used for the actual construction of the design and penalty matrices.
#'
#' You can provide one or multiple data streams of different types (real, positive, circular) and specify initial means and standard deviations/ concentrations for each data stream. This information is then converted into suitable spline coefficients.
#' \code{buildSmoothDens} then constructs the design and penalty matrices for standardised B-splines basis functions (integrating to one) for each data stream.
#' For types \code{"real"} and \code{"circular"} the knots are placed equidistant in the range of the data, for type \code{"positive"} the knots are placed using polynomial spacing.
#'
#' @param data named data frame of different data streams
#' @param type type of each data stream, either \code{"real"} for data on the reals, \code{"positive"} for data on the positive reals or \code{"circular"} for angular data. Needs to be a vector corresponding to the number of data streams in \code{data}.
#' @param par nested named list of initial means and sds/concentrations for each data stream
#' @param k number of basis functions for each data stream
#' @param degree degree of the B-spline basis functions for each data stream, defaults to cubic B-splines
#' @param diff_order order of differencing used for the P-Spline penalty matrix for each data stream. Defaults to second-order differences.
#'
#' @return a nested list containing the design matrices \code{Z}, the penalty matrices \code{S}, the initial coefficients \code{coef} the prediction design matrices \code{Z_predict}, the prediction grids \code{xseq}, and details for the basis expansion for each data stream.
#' @export
#'
#' @examples
#' ## 3 data streams, each with one distribution
#' # normal data with mean 0 and sd 1
#' x1 = rnorm(100, mean = 0, sd = 1)
#' # gamma data with mean 5 and sd 3
#' x2 = rgamma2(100, mean = 5, sd = 3)
#' # circular data
#' x3 = rvm(100, mu = 0, kappa = 2)
#' 
#' data = data.frame(x1 = x1, x2 = x2, x3 = x3)
#' 
#' par = list(x1 = list(mean = 0, sd = 1),
#'            x2 = list(mean = 5, sd = 3),
#'            x3 = list(mean = 0, concentration = 2))
#' 
#' SmoothDens = buildSmoothDens(data, 
#'                              type = c("real", "positive", "circular"),
#'                              par)
#'                              
#' # extracting objects for x1
#' Z1 = SmoothDens$Z$x1
#' S1 = SmoothDens$S$x1
#' coefs1 = SmoothDens$coef$x1
#' 
#' ## one data stream, but mixture of two distributions
#' # normal data with mean 0 and sd 1
#' x = rnorm(100, mean = 0, sd = 1)
#' data = data.frame(x = x)
#' 
#' # now parameters for mixture of two normals
#' par = list(x = list(mean = c(0, 5), sd = c(1,1)))
#' 
#' SmoothDens = buildSmoothDens(data, par = par)
#' 
#' # extracting objects 
#' Z = SmoothDens$Z$x
#' S = SmoothDens$S$x
#' coefs = SmoothDens$coef$x
buildSmoothDens = function(data,
                           type = "real",
                           par,
                           k = 20,
                           degree = 3,
                           diff_order = 2
){
  if(!is.data.frame(data)){
    stop("datastreams must be a data frame")
  }
  nStreams = ncol(data)
  nObs = nrow(data)
  varnames = colnames(data)
  if(length(k) == 1){
    k = rep(k, nStreams)
  } else if(length(k) != nStreams){
    stop("k must be a scalar or a vector of length equal to the number of datastreams")
  }
  if(length(type) == 1){
    type = rep(type, nStreams)
  } else if(length(type) != nStreams){
    stop("type must be of length 1 or equal to the number of datastreams")
  }
  if(length(degree) == 1){
    degree = rep(degree, nStreams)
  } else if(length(degree) != nStreams){
    stop("degree must be of length 1 or equal to the number of datastreams")
  }
  if(length(diff_order) == 1){
    diff_order = rep(diff_order, nStreams)
  } else if(length(diff_order) != nStreams){
    stop("diff_order must be of length 1 or equal to the number of datastreams")
  }
  
  listseed = vector("list", length = nStreams)
  names(listseed) = varnames
  Z = S = Z_predict = xseq = betastart = basis = listseed
  
  for(i in seq_len(nStreams)){
    thisname = varnames[i]
    
    cat(thisname, "\n")
    
    modmat = make_matrices_dens(x = data[[thisname]], type = type[i], k = k[i], 
                                degree = degree[i], diff_order = diff_order[i], npoints = 1e4)
    Z[[thisname]] = modmat$Z
    S[[thisname]] = modmat$S
    Z_predict[[thisname]] = modmat$Z_predict
    xseq[[thisname]] = modmat$xseq
    basis[[thisname]] = modmat$basis
    betastart[[thisname]] = make_splinecoef(modmat, type = type[i], par = par[[thisname]])
  }
  
  return(list(
    Z = Z,
    S = S,
    coef = betastart,
    Z_predict = Z_predict,
    xseq = xseq,
    basis = basis
  ))
}


#' Compute the design matrix for a trigonometric basis expansion
#'
#' Given a periodically varying variable such as time of day or day of year and the associated cycle length, this function performs a basis expansion to efficiently calculate a linear predictor of the form
#' \deqn{ 
#'  \eta^{(t)} = \beta_0 + \sum_{k=1}^K \bigl( \beta_{1k} \sin(\frac{2 \pi k t}{L}) + \beta_{2k} \cos(\frac{2 \pi k t}{L}) \bigr). 
#'  }
#'  This is relevant for modeling e.g. diurnal variation and the flexibility can be increased by adding smaller frequencies (i.e. increasing \eqn{K}).
#'  
#' @param tod equidistant sequence of a cyclic variable
#' 
#' For time of day and e.g. half-hourly data, this could be 1, ..., L and L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.
#' @param L length of one cycle on the scale of the time variable. For time of day, this would be 24.
#' @param degree degree K of the trigonometric link above. Increasing K increases the flexibility.
#'
#' @return design matrix (without intercept column), ordered as sin1, cos1, sin2, cos2, ...
#' @export
#'
#' @examples
#' ## hourly data
#' tod = rep(1:24, 10)
#' Z = trigBasisExp(tod, L = 24, degree = 2)
#' 
#' ## half-hourly data
#' tod = rep(1:48/2, 10) # in [0,24] -> L = 24
#' Z1 = trigBasisExp(tod, L = 24, degree = 3)
#' 
#' tod = rep(1:48, 10) # in [1,48] -> L = 48
#' Z2 = trigBasisExp(tod, L = 48, degree = 3)
#' 
#' all(Z1 == Z2)
#' # The latter two are equivalent specifications!
trigBasisExp = function(tod, L = 24, degree = 1){
  n = length(tod)
  Z = matrix(nrow = n, ncol = 2*degree)
  inner = 2*pi*tod/L
  for(k in seq_len(degree)){
    Z[,2*(k-1)+1:2] = cbind(sin(inner*k), cos(inner*k))
  }
  colnames(Z) = paste0(c("sin_", "cos_"), rep(1:degree, each = 2))
  Z
}
