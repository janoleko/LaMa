#' Computes penalty based on quadratic form
#'
#' @description
#' This function computes quadratic penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.
#' 
#' It is intended to be used inside the \strong{penalised negative log-likelihood function} when fitting models with penalised splines or simple random effects via \strong{quasi restricted maximum likelihood} (qREML) with the \code{\link{qreml}} function.
#' For \code{\link{qreml}} to work, the likelihood function needs to be compatible with the \code{RTMB} R package to enable automatic differentiation.
#' 
#' @seealso \code{\link{qreml}} for the \strong{qREML} algorithm
#' 
#' @details
#' \strong{Caution:} The formatting of \code{re_coef} needs to match the structure of the parameter list in your penalised negative log-likelihood function, 
#' i.e. you cannot have two random effect vectors of different names (different list elements in the parameter list), combine them into a matrix inside your likelihood and pass the matrix to \code{penalty}.
#' If these are seperate random effects, each with its own name, they need to be passed as a list to \code{penalty}. Moreover, the ordering of \code{re_coef} needs to match the character vector \code{random} specified in \code{\link{qreml}}.
#' 
#'
#' @param re_coef coefficient vector/ matrix or list of coefficient vectors/ matrices
#'
#' Each list entry corresponds to a different smooth/ random effect with its own associated penalty matrix in \code{S}.
#' When several smooths/ random effects of the same kind are present, it is convenient to pass them as a matrix, where each row corresponds to one smooth/ random effect. 
#' This way all rows can use the same penalty matrix.
#' @param S fixed penalty matrix or list of penalty matrices matching the structure of \code{re_coef} and also the dimension of the individuals smooths/ random effects
#' @param lambda penalty strength parameter vector that has a length corresponding to the \strong{total number} of random effects/ spline coefficients in \code{re_coef}
#'
#' E.g. if \code{re_coef} contains one vector and one matrix with 4 rows, then \code{lambda} needs to be of length 5.
#'
#' @return returns the penalty value and reports to \code{\link{qreml}}.
#' @export
#' 
#' @import RTMB
#'
#' @examples
#' # Example with a single random effect
#' re = rep(0, 5)
#' S = diag(5)
#' lambda = 1
#' penalty(re, S, lambda)
#'
#' # Example with two random effects, 
#' # where one element contains two random effects of similar structure
#' re = list(matrix(0, 2, 5), rep(0, 4))
#' S = list(diag(5), diag(4))
#' lambda = c(1,1,2) # length = total number of random effects
#' penalty(re, S, lambda)
#' 
#' # Full model-fitting example
#' data = trex[1:1000,] # subset
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 1)), # step mean
#'            logsigma = log(c(0.2, 0.7)), # step sd
#'            beta0 = c(-2,-2), # state process intercept
#'            betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
#'           
#' # data object with initial penalty strength lambda
#' dat = list(step = data$step, # step length
#'            tod = data$tod, # time of day covariate
#'            N = 2, # number of states
#'            lambda = rep(10,2)) # initial penalty strength
#'
#' # building model matrices
#' modmat = make_matrices(~ s(tod, bs = "cp"), 
#'                        data = data.frame(tod = 1:24), 
#'                        knots = list(tod = c(0,24))) # wrapping points
#' dat$Z = modmat$Z # spline design matrix
#' dat$S = modmat$S # penalty matrix
#'
#' # penalised negative log-likelihood function
#' pnll = function(par) {
#'   getAll(par, dat) # makes everything contained available without $
#'   Gamma = tpm_g(Z, cbind(beta0, betaspline)) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1) # initial distribution
#'   mu = exp(logmu) # step mean
#'   sigma = exp(logsigma) # step sd
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step)) # only for non-NA obs.
#'   for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
#'   -forward_g(delta, Gamma[,,tod], allprobs) +
#'       penalty(betaspline, S, lambda) # this does all the penalization work
#' }
#'
#' # model fitting
#' mod = qreml(pnll, par, dat, random = "betaspline")
penalty = function(re_coef, S, lambda) {
  # Capture the argument name used in the call to `penalty`
  # current_name <- as.character(substitute(lambda))
  # 
  # # Try to recover the original name by searching the parent frame
  # recover_original_name <- function(value, env) {
  #   # Search for all objects in the parent frame
  #   all_objects <- ls(env)
  #   for (obj in all_objects) {
  #     # Check if the object matches the value of `lambda`
  #     if (identical(get(obj, envir = env), value)) {
  #       return(obj)  # Return the name of the matching object
  #     }
  #   }
  #   NULL  # Return NULL if no match is found
  # }
  # 
  # # Attempt to recover the original name
  # original_name <- recover_original_name(lambda, parent.frame())
  # if (!is.null(original_name)) {
  #   argname_lambda <- original_name
  # } else {
  #   argname_lambda <- current_name
  # }
  # 
  # # Store the recovered name in the penalty metadata environment
  # assign("argname_lambda", argname_lambda, envir = penalty_metadata)
  
  # Convert re_coef to a list of matrices (even if originally a vector)
  if (!is.list(re_coef)) {
    re_coef = list(re_coef)
  }
  
  re_coef = lapply(re_coef, function(x) {
    if (is.null(dim(x))) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  # Get number of distinct random effects (of the same structure)
  n_re = length(re_coef)
  
  # Ensure S is a list of length n_re, replicating it if necessary
  if (!is.list(S)) {
    S = list(S)
  }
  if (length(S) == 1) {
    S = rep(S, n_re)
  }
  
  # transpose if necessary to match S
  re_coef <- lapply(seq_len(n_re), function(i) {
    if (ncol(re_coef[[i]]) != nrow(S[[i]])) {
      t(re_coef[[i]])
    } else if (nrow(re_coef[[i]]) != nrow(S[[i]])) {
      re_coef[[i]]
    } else{
      stop("The coefficient structure does not match the provided penalty matrices.")
    }
  })
  
  # Get the number of similar random effects for each distinct random effect
  re_lengths = sapply(re_coef, nrow)  # All elements are matrices now
  
  # Precompute start and end indices for lambda
  end = cumsum(re_lengths)
  start = c(1, end[-length(end)] + 1)
  
  RTMB::REPORT(S) # Report penalty matrix list
  
  # Initialize penalty variables
  Pen = vector("list", n_re)
  pen = 0
  
  # check if re_coef and S match
  # if(any(sapply(re_coef, ncol) != sapply(S, nrow))){
  #   stop("The coefficient structure does not match the provided penalty matrices.\n Are the coefficients arranged by row?")
  # }
  
  # Loop over distinct random effects - each now a matrix
  for (i in 1:n_re) {
    current_re = re_coef[[i]]  # current_re is always a matrix now
    
    # Vectorized calculation of penalty for each random effect
    quadform = rowSums(current_re %*% S[[i]] * current_re)
    Pen[[i]] = quadform
    
    # Apply lambda directly using precomputed indices
    pen = pen + sum(lambda[start[i]:end[i]] * quadform)
  }
  
  RTMB::REPORT(Pen) # Report the penalty list for qreml update
  
  0.5 * pen
}

#' Quasi restricted maximum likelihood (qREML) algorithm for models with penalised splines or simple i.i.d. random effects
#'
#' @description
#' This algorithm can be used very flexibly to fit statistical models that involve \strong{penalised splines} or simple \strong{i.i.d. random effects}, i.e. that have penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.
#'
#' The \strong{qREML} algorithm is typically much faster than REML or marginal ML using the full Laplace approximation method, but may be slightly less accurate regarding the estimation of the penalty strength parameters.
#'
#' Under the hood, \code{qreml} uses the R package \code{RTMB} for automatic differentiation in the inner optimisation.
#' The user has to specify the \strong{penalised negative log-likelihood function} \code{pnll} structured as dictated by \code{RTMB} and use the \code{\link{penalty}} function to compute the quadratic-form penalty inside the likelihood.
#' 
#' @seealso \code{\link{penalty}} to compute the penalty inside the likelihood function
#'
#' @param pnll penalised negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{\link{penalty}} function from \code{LaMa} to compute the penalty
#'
#' Needs to be a function of the named list of initial parameters \code{par} only.
#' @param par named list of initial parameters
#'
#' The random effects/ spline coefficients can be vectors or matrices, the latter summarising several random effects of the same structure, each one being a row in the matrix.
#' @param dat initial data list that contains the data used in the likelihood function, hyperparameters, and the \strong{initial penalty strength} vector
#'
#' If the initial penalty strength vector is \strong{not} called \code{lambda}, the name it has in \code{dat} needs to be specified using the \code{psname} argument below.
#' Its length needs to match the to the total number of random effects.
#' @param random vector of names of the random effects/ penalised parameters in \code{par}
#' 
#' \strong{Caution:} The ordering of \code{random} needs to match the order of the random effects passed to \code{\link{penalty}} inside the likelihood function.
#' @param map optional map argument, containing factor vectors to indicate parameter sharing or fixing.
#' 
#' Needs to be a named list for a subset of fixed effect parameters or penalty strength parameters. 
#' For example, if the model has four penalty strength parameters, \code{map[[psname]]} could be \code{factor(c(NA, 1, 1, 2))} to fix the first penalty strength parameter, estimate the second and third jointly, and estimate the fourth separately.
#' @param psname optional name given to the penalty strength parameter in \code{dat}. Defaults to \code{"lambda"}.
#' @param alpha optional hyperparamater for exponential smoothing of the penalty strengths.
#'
#' For larger values smoother convergence is to be expected but the algorithm may need more iterations.
#' @param smoothing optional scaling factor for the final penalty strength parameters
#' 
#' Increasing this beyond one will lead to a smoother final model. Can be an integer or a vector of length equal to the length of the penalty strength parameter.
#' @param maxiter maximum number of iterations in the outer optimisation over the penalty strength parameters.
#' @param tol Convergence tolerance for the penalty strength parameters.
#' @param control list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimisation. Here, \code{optim} uses the \code{BFGS} method which cannot be changed.
#' 
#' We advise against changing the default values of \code{reltol} and \code{maxit} as this can decrease the accuracy of the Laplace approximation.
#' @param silent integer silencing level: 0 corresponds to full printing of inner and outer iterations, 1 to printing of outer iterations only, and 2 to no printing.
#' @param joint_unc logical, if \code{TRUE}, joint \code{RTMB} object is returned allowing for joint uncertainty quantification
#' @param saveall logical, if \code{TRUE}, then all model objects from each iteration are saved in the final model object.
#' # @param epsilon vector of two values specifying the cycling detection parameters. If the relative change of the new penalty strength to the previous one is larger than \code{epsilon[1]} but the change to the one before is smaller than \code{epsilon[2]}, the algorithm will average the two last values to prevent cycling.
#'
#' @return returns a model list influenced by the users report statements in \code{pnll}
#' @export
#'
#' @import RTMB
#'
#' @examples
#' data = trex[1:1000,] # subset
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 1)), # step mean
#'            logsigma = log(c(0.2, 0.7)), # step sd
#'            beta0 = c(-2,-2), # state process intercept
#'            betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
#'           
#' # data object with initial penalty strength lambda
#' dat = list(step = data$step, # step length
#'            tod = data$tod, # time of day covariate
#'            N = 2, # number of states
#'            lambda = rep(10,2)) # initial penalty strength
#'
#' # building model matrices
#' modmat = make_matrices(~ s(tod, bs = "cp"), 
#'                        data = data.frame(tod = 1:24), 
#'                        knots = list(tod = c(0,24))) # wrapping points
#' dat$Z = modmat$Z # spline design matrix
#' dat$S = modmat$S # penalty matrix
#'
#' # penalised negative log-likelihood function
#' pnll = function(par) {
#'   getAll(par, dat) # makes everything contained available without $
#'   Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = TRUE) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1, ad = TRUE) # initial distribution
#'   mu = exp(logmu) # step mean
#'   sigma = exp(logsigma) # step sd
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step)) # only for non-NA obs.
#'   for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
#'   -forward_g(delta, Gamma[,,tod], allprobs) +
#'       penalty(betaspline, S, lambda) # this does all the penalization work
#' }
#'
#' # model fitting
#' mod = qreml(pnll, par, dat, random = "betaspline")
qreml = function(pnll, # penalized negative log-likelihood function
                 par, # initial parameter list
                 dat, # initial dat object, currently needs to be called dat!
                 random, # names of parameters in par that are random effects/ penalized
                 map = NULL, # map for fixed effects
                 psname = "lambda", # name given to the psname parameter in dat
                 alpha = 0.25, # exponential smoothing parameter
                 smoothing = 1,
                 maxiter = 100, # maximum number of iterations
                 tol = 1e-4, # tolerance for convergence
                 control = list(reltol = 1e-10, maxit = 1000), # control list for inner optimization
                 silent = 1, # print level
                 joint_unc = TRUE, # should joint object be returned?
                 saveall = FALSE)# , # save all intermediate models?
                 #epsilon = c(1e-2, 1e-1)) # cycling detection parameters 
{
  
  # setting the argument name for par because later updated par is returned
  argname_par = as.character(substitute(par))
  argname_dat = as.character(substitute(dat))
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re = length(random) 
  
  # list to save all model objects
  allmods = list() 
  
  # initial lambda locally
  # Define a global environment to store the captured names
  # penalty_metadata <- new.env(parent = emptyenv())
  # if(is.null(psname)){
  #   pnll(par) # call once to get the name of the lambda parameter
  #   psname = get("argname_lambda", envir = penalty_metadata)
  # }
  
  lambda = dat[[psname]]
  lambda0 = lambda # saving initial lambda so that fixed pars can always be refilled
  
  # creating the objective function as wrapper around pnll to pull lambda from local
  f = function(par){
    environment(pnll) = environment()
    
    # overloading assignment operators, currently necessary
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    getLambda = function(x) lambda
    
    dat[[psname]] = DataEval(getLambda, rep(advector(1), 0))
    
    # assigning dat to whatever it is called in pnll() (hopefully)
    assign(argname_dat, dat, envir = environment())
    
    pnll(par)
  }
  
  # creating the RTMB objective function
  if(silent %in% 0:1){
    cat("Creating AD function\n")
  } 
  
  ## mapping
  # map can contain fixed effects -> just passed to MakeADFun
  # and it can contain penalty strength parameters
  if(!is.null(map)){
    # check that no random effects are fixed
    if(any(names(map) %in% random)){
      msg <- "'map' cannot contain random effects or spline parameters"
      stop(msg)
    }
    # make factor
    if(!all(sapply(map, is.factor))){
      message("Converting map to factor")
    }
    map = lapply(map, factor)
  } else {
    map[[psname]] = factor(seq_along(lambda))
  }
  lambda_map = map[[psname]]
  if(length(lambda_map) != length(lambda)){
    msg = paste0("Length of map argument for ", psname, " has wrong length.")
    stop(msg)
  }
  
  # pop lambda_map from map list
  map = map[names(map) != psname]
  if(length(map) == 0) map = NULL
  
  # deal with mapping of penalty strength parameters
  lambda_mapped = map_lambda(lambda, lambda_map)
  if(length(lambda_mapped) < 1){
    stop("No penalty paramters estimated as all are fixed.")
  }
  
  obj = MakeADFun(func = f, 
                  parameters = par, 
                  silent = TRUE,
                  map = map) # silent and replacing with own prints
  
  newpar = obj$par # saving initial parameter value as vector to initialize optimization in loop
  
  # own printing of maximum gradient component if silent = 0
  if(silent == 0){
    newgrad = function(par){
      gr = obj$gr(par)
      cat(" inner mgc:", max(abs(gr)), "\n")
      gr
    }
  } else{
    newgrad = obj$gr
  }
  
  # prepwork
  mod0 = obj$report() # getting all necessary information from penalty report
  S = mod0$S # penalty matrix/ matrices in list format
  # S_dims = sapply(S, nrow)
  
  # finding the indices of the random effects to later index Hessian
  re_inds = list() 
  for(i in seq_len(n_re)){
    if(is.vector(par[[random[i]]])){
      re_dim = c(1, length(par[[random[i]]]))
    } else if(is.matrix(par[[random[i]]])){
      re_dim = dim(par[[random[i]]])
    } else{
      stop(paste0(random[i], " must be a vector or matrix"))
    }
    
    byrow = FALSE
    if(re_dim[1] == nrow(S[[i]])){
      byrow = TRUE
    }
    #re_dim = dim(as.matrix(par[[random[i]]]))
    
    # which(re_dim == nrow(S[[i]])) - 1
    
    # if(re_dim[2] == S_dims[i]){
    #   byrow = FALSE
    # } else{
    #   byrow = TRUE
    # }
    re_inds[[i]] = matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])#, byrow = byrow)
    if(byrow) {
      re_inds[[i]] = t(re_inds[[i]]) # if byrow, then transpose
    }
    # if(dim(re_inds[[i]])[2] == 1) re_inds[[i]] = t(re_inds[[i]]) # if only one column, then transpose
  }
  
  # get number of similar random effects for each distinct random effect (of same structure)
  re_lengths = sapply(re_inds, function(x) if (is.vector(x)) 1 else nrow(x))
  
  # initialize list of penalty strength parameters
  Lambdas = list()
  Lambdas[[1]] = reshape_lambda(re_lengths, lambda) # reshaping to match structure of random effects
  
  if(silent < 2){
    cat("Initialising with", paste0(psname, ":"), round(lambda, 3), "\n")
  }
  
  # computing rank deficiency for each penalty matrix to use in correction term
  m = numeric(length(S)) 
  for(i in seq_len(length(m))) {
    m[i] = nrow(S[[i]]) - Matrix::rankMatrix(S[[i]])
  } 
  
  # initialising convergence check index (initially for all lambdas)
  convInd <- seq_along(unlist(Lambdas[[1]]))
  
  ### updating algorithm
  # loop over outer iterations until convergence or maxiter
  for(k in seq_len(maxiter)){
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    opt = stats::optim(newpar, obj$fn, newgrad, 
                       method = "BFGS", hessian = TRUE, # return hessian in the end
                       control = control)
    
    # setting new optimum par for next iteration
    newpar = opt$par 
    
    # reporting to extract penalties
    mod = obj$report() 
    
    # evaluating current Hessian
    # J = obj$he()
    J = opt$hessian
    
    # computing inverse Hessian
    J_inv = MASS::ginv(J) 
    
    # saving entire model object
    if(saveall){
      allmods[[k]] = mod
    }
    
    ## updating all lambdas
    lambdas_k = list() # temporary lambda list
    
    # looping over distinct random effects (matrices)
    edoF = rep(NA, length(lambda0)) # initialise edoF vector
    pens = rep(NA, length(lambda0)) # initialise penalty vector
    l = 1 # counter for lambda
    
    for(i in 1:n_re){
      # initializing lambda vector for i-th random effect
      lambdas_k[[i]] = numeric(nrow(re_inds[[i]]))
      
      # looping over similar random effects
      for(j in 1:nrow(re_inds[[i]])){
        idx = re_inds[[i]][j,] # indices of this random effect
        
        # effective degrees of freedom for this random effect: J^-1_p J
        edoF[l] = nrow(S[[i]]) - m[i] - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
        
        # penalty
        pens[l] = mod$Pen[[i]][j]
        
        l = l+1
        
      #   # effective degrees of freedom for this random effect: J^-1_p J
      #   edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
      #   
      #   # calculating new lambda based on updating rule
      #   lambda_new = as.numeric((edoF - m[i]) / mod$Pen[[i]][j]) # m is correction if S_i does not have full rank
      #   
      #   # potentially smoothing new lambda
      #   lambdas_k[[i]][j] = (1-alpha) * lambda_new + alpha * Lambdas[[k]][[i]][j]
      #   
      #   # check for cycling behaviour
      #   if(k > 2){
      #     if(abs((lambdas_k[[i]][j] - Lambdas[[k-1]][[i]][j]) / Lambdas[[k-1]][[i]][j]) < epsilon[1] & # change to lambda_t-2 is small
      #        abs((lambdas_k[[i]][j] - Lambdas[[k]][[i]][j]) / Lambdas[[k]][[i]][j]) > epsilon[2]) # but change to lambda_t-1 is large
      #     {
      #       cat("Cycling detected - averaging for faster convergence\n")
      #       # replacing with mean to prevent cycling
      #       lambdas_k[[i]][j] = (lambdas_k[[i]][j] + Lambdas[[k]][[i]][j]) / 2 
      #     }
      #   }
      }
    }
    
    # now loop over actual lambda_mapped to update
    for(i in seq_along(lambda_mapped)){
      this_level = levels(lambda_map)[i]
      this_ind = which(lambda_map == this_level)
      
      this_edoF = sum(edoF[this_ind])
      this_pen = sum(pens[this_ind])
      
      lambda_new = this_edoF / this_pen
      
      # smoothing lambda
      lambda_mapped[i] = (1-alpha) * lambda_new + alpha * lambda_mapped[i]
    }
    
    # unmap lambda
    lambdas_k = unmap_lambda(lambda_mapped, lambda_map, lambda0)
    
    # minimum of zero for penalty strengths
    lambdas_k[which(lambdas_k < 0)] = 0
    
    # assigning new lambda to global list
    Lambdas[[k+1]] = relist(lambdas_k, Lambdas[[1]])
    
    # updating lambda vector locally for next iteration
    lambda = lambdas_k
    
    # old length of convergence check indices
    oldlength <- length(convInd)
    
    if(k > 2){ # after 2 iterations, check whether any lambda > 1e5 and exclude from check
      convInd = which(lambda <= 1e6)
    }
    
    if(silent < 2){
      cat("outer", k, "-", paste0(psname, ":"), round(lambda, 3), "\n")
      
      # print only if something changes
      if(length(convInd) != oldlength & length(seq_along(lambda)[-convInd]) > 0){
        cat(psname, seq_along(lambda)[-convInd], "excluded from convergence check (> 1e6)", "\n")
      }
    }
    
    # convergence check
    # if(all(abs(lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol)){
    if(max(abs(
      (lambda - unlist(Lambdas[[k]]))[convInd] / unlist(Lambdas[[k]])[convInd]
    )) < tol){
      
      if(silent < 2){
        cat("Converged\n")
      }
      break
    }
    
    if(k == maxiter){
      cat("No convergence\n")
      warning("No convergence\n")
    } 
  }
  
  # final model fit
  lambda = lambda * smoothing # scaling lambda by smoothing factor
  
  if(silent < 2){
    if(any(smoothing != 1)){
      cat("Smoothing factor:", smoothing, "\n")
    }
    cat("Final model fit with", paste0(psname, ":"), round(lambda, 3), "\n")
  }
  
  # fitting the model conditional on lambda: current local lambda will be pulled by f
  opt = stats::optim(newpar, obj$fn, newgrad, 
                     method = "BFGS", hessian = TRUE, # return hessian in the end
                     control = control)
  
  # setting new optimum par for next iteration
  newpar = opt$par 
  
  # reporting to extract penalties
  mod = obj$report() 
  
  # evaluating current Hessian
  # J = obj$he()
  J = opt$hessian
  
  # computing inverse Hessian
  J_inv = MASS::ginv(J) 
  
  # saving entire model object
  if(saveall){
    allmods[[k+1]] = mod
  }
  
  #############################################
  
  # assign RTMB obj to return object
  mod$obj <- obj
  
  # if all intermediate models should be returned, assign
  if(saveall) {
    mod$allmods = allmods
  }
  
  # assign final lambda to return object
  mod[[psname]] = lambda
  
  # assigning all lambdas to return object
  mod[[paste0("all_", psname)]] = Lambdas
  
  # calculating unpenalized log-likelihood at final parameter values
  lambda = rep(0, length(lambda))
  dat[[psname]] = lambda
  
  # format parameter to list
  # skeleton = utils::as.relistable(par)
  # parlist = utils::relist(opt$par, skeleton)
  parlist = obj$env$parList(opt$par)
  mod[[argname_par]] = parlist # and assing to return object
  
  mod[[paste0("relist_", argname_par)]] = obj$env$parList
  
  # assign estimated parameter as vector
  mod[[paste0(argname_par, "_vec")]] = opt$par
  
  # assign log-likelihood at optimum to return object
  mod$llk = -pnll(parlist)
  
  ## calculating effective degrees of freedom for final model
  mod$edf = list()
  for(i in 1:n_re){
    edoF_i = numeric(nrow(re_inds[[i]]))
    for(j in 1:nrow(re_inds[[i]])){
      idx = re_inds[[i]][j,]
      edoF_i[j] = edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]]))
    }
    mod$edf[[i]] = edoF_i
  }
  
  # number of fixed parameters
  # mod$n_fixpar = length(unlist(par[!(names(par) %in% random)]))
  mod$n_fixpar = length(opt$par)
  
  # assing conditinoal Hessian
  mod$Hessian_conditional = J
  
  # removing penalty list from model object
  mod = mod[names(mod) != "Pen"] 
  
  
  if(joint_unc){
    ### constructing joint object
    parlist$loglambda = log(mod[[psname]])
    
    # finding the number of similar random effects for each random effect
    # indvec = rep(1:n_re, times = re_lengths)
    
    # computing log determinants
    logdetS = numeric(length(S))
    for(i in 1:length(S)){
      logdetS[i] = gdeterminant(S[[i]])
    }
    
    ## defining joint negative log-likelihood
    jnll = function(par) {
      
      environment(pnll) = environment()
      
      # overloading assignment operators, currently necessary
      "[<-" <- ADoverload("[<-") 
      "c" <- ADoverload("c")
      "diag<-" <- ADoverload("diag<-")
      
      dat[[psname]] = exp(par$loglambda)
      
      l_p = -pnll(par[names(par) != "loglambda"])
      
      ## computing additive constants (missing from only penalized likelihood)
      const = 0
      for(i in 1:n_re){
        for(j in 1:nrow(re_inds[[i]])){
          k = length(re_inds[[i]][j,])
          
          if(i == 1){
            loglam = par$loglambda[j]
          } else{
            loglam = par$loglambda[re_lengths[i-1] + j]
          }
          
          const = const - k * log(2*pi) + k * loglam + logdetS[i]
        }
      }
      
      l_joint = l_p + 0.5 * const
      -l_joint
    }
    
    if(is.null(map)){
      map = list(loglambda = lambda_map)
    } else{
      map$loglambda = lambda_map
    }
    
    # creating joint AD object
    obj_joint = MakeADFun(jnll, parlist,
                          random = names(par)[names(par) != "loglambda"], # REML, everything random except lambda
                          map = map)
    
    # assigning object to return object
    mod$obj_joint = obj_joint
  }
  
  class(mod) = "qreml model"
  return(mod)
}


