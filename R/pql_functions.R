# helper function for penalty and pql
reshape_lambda <- function(num_elements, lambda) {
  start <- 1
  result <- lapply(num_elements, function(len) {
    # Extract sub-vector from lambda based on the number of elements
    sub_vector <- lambda[start:(start + len - 1)]
    start <<- start + len
    return(sub_vector)
  })
  return(result)
}


#' Computes penalty based on quadratic form
#'
#' This function computes penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i}
#' and is intended to be used inside the penalized negative log-likelihood function when fitting models with splines or simple random effects with \code{RTMB} via penalized quasi-likelihood (PQL) with the \code{pql()} function.
#'
#' @param re_coef Coefficient vector, matrix or list of coefficient vectors/ matrices.\cr\cr
#' Each list entry corresponds to a different smooth/ random effect with its own associated penalty matrix in \code{S}.
#' When several smooths/ random effects of the same kind are present, it is convenient to pass them as a matrix, where each row corresponds to one smooth/ random effect.\cr\cr
#' Caution: The formatting of \code{re_coef} needs to match the structure of the parameter list in your penalized negative log-likelihood function, 
#' i.e. you cannot have two random effect vectors of different names (different list elements in the parameter list), combine them into a matrix inside your likelihood and pass the matrix to \code{penalty}.
#' If these are seperate random effects, each with its own name, they need to be passed as a list to \code{penalty}. Moreover, the ordering of \code{re_coef} needs to match the character vector \code{random} specified in \code{pql}.
#' @param S Penalty matrix or list of penalty matrices matching the structure of \code{re_coef} and also the dimension of the individuals smooths/ random effects.
#' @param lambda Penalty strength parameter. Vector that has a length corresponding to the total number of random effects/ spline coefficients in \code{re_coef}.
#'
#' @return Returns the penalty value and reports to \code{pql()} for a seamless experience.
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
penalty = function(re_coef, S, lambda) {
  # "[<-" <- ADoverload("[<-") # currently necessary
  
  # Convert re_coef to a list of matrices (even if originally a vector)
  if (!is.list(re_coef)) {
    re_coef = list(re_coef)
  }
  
  re_coef = lapply(re_coef, function(x) {
    if (is.vector(x)) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  # Get number of distinct random effects (of the same structure)
  n_re = length(re_coef)
  
  # Get the number of similar random effects for each distinct random effect
  re_lengths = sapply(re_coef, nrow)  # All elements are matrices now
  
  # Precompute start and end indices for lambda
  end = cumsum(re_lengths)
  start = c(1, end[-length(end)] + 1)
  
  # Ensure S is a list of length n_re, replicating it if necessary
  if (!is.list(S)) {
    S = list(S)
  }
  if (length(S) == 1) {
    S = rep(S, n_re)
  }
  
  RTMB::REPORT(S) # Report penalty matrix list
  
  # Initialize penalty variables
  Pen = vector("list", n_re)
  pen = 0
  
  # Loop over distinct random effects - each now a matrix
  for (i in 1:n_re) {
    current_re = re_coef[[i]]  # current_re is always a matrix now
    
    # Vectorized calculation of penalty for each random effect
    quadform = rowSums(current_re %*% S[[i]] * current_re)
    Pen[[i]] = quadform
    
    # Apply lambda directly using precomputed indices
    pen = pen + sum(lambda[start[i]:end[i]] * quadform)
  }
  
  RTMB::REPORT(Pen) # Report the penalty list for pql update
  
  0.5 * pen
}

#' Penalized quasi-likelihood (PQL) algorithm for models with simple random effects
#'
#' This algorithm can be used very flexible to fit statistical models that involves \strong{penalized splines} or simple \strong{i.i.d. random effects} with \code{RTMB} that have penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i}
#' PQL is typically much faster than the full Laplace approximation method, but may be slightly less accurate regarding the estimation of the penalty strength parameters.
#' The user has to specify the penalized negative log-likelihood function \code{pnll} structured as dictated by \code{RTMB} and use the \code{penalty} function contained in \code{LaMa} to compute the penalty inside the likelihood.
#'
#' @param pnll Penalized negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{penalty} function from \code{LaMa} to compute the penalty.
#' @param par Named list of initial parameters. The random effects can be vectors or matrices, the latter summarising several random effects of the same structure, each one being a row in the matrix.
#' @param dat Initial data list, that contains the data used in the likelihood function, hyperparameters, and the initial penalty strength. If initial penalty strength vector is not called \code{lambda}, you need to specify its name in \code{dat}. 
#' Its length needs to match the to the total number of random effects.
#' @param random Vector of names of the random effects in \code{par} that are penalized. Caution: The ordering of \code{random} needs to match the order of the random effects passed to \code{penalty()} inside the likelihood function.
#' @param penalty Name given to the penalty parameter in \code{dat}. Defaults to \code{"lambda"}.
#' @param alpha Optional hyperparamater for exponential smoothing of the penalty strengths. For larger values smoother convergence is to be expected but the algorithm may need more iterations.
#' @param maxiter Maximum number of iterations.
#' @param tol Convergence tolerance for the penalty strength parameters.
#' @param inner_tol Convergence tolerance for the inner optimization.
#' @param silent Integer silencing level: 0 corresponds to full printing of inner and outer iteratinos, 1 to printing of outer iterations only, and 2 to no printing.
#' @param saveall Logical, if TRUE, then all model objects from each iteration are saved in the final model object. Defaults to FALSE.
#'
#' @return Returns a model list influenced by the users report statements in \code{pnll}
#' @export
#'
#' @import RTMB
#'
#' @examples
#'data = elephant[1:1000,] # subset
#'# initial parameter list
#'par = list(logmu = log(c(0.3, 1)), 
#'           logsigma = log(c(0.2, 0.7)),
#'           logkappa = log(c(0.2, 0.7)),
#'           beta0 = c(-2,2),
#'           betaspline = matrix(rep(0, 18), nrow = 2))
#'# data object with initial penalty strength lambda
#'dat = list(step = data$step, angle = data$angle, tod = data$tod, N = 2, lambda = rep(10,2))
#'# model matrices
#'modmat = make_matrices(~ s(tod, bs = "cp"), data = data.frame(tod = 1:24), 
#'    knots = list(tod = c(0,24))) # wrapping points
#'dat$Z = modmat$Z # spline design matrix
#'dat$S = modmat$S # penalty matrix
#'# penalized negative log-likelihood function
#'pnll = function(par) {
#' getAll(par, dat) # makes everything contained available without $
#' Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = TRUE)
#' delta = stationary_p(Gamma, t = 1, ad = TRUE)
#' mu = exp(logmu)
#' sigma = exp(logsigma)
#' kappa = exp(logkappa)
#' # calculating all state-dependent densities
#' allprobs = matrix(1, nrow = length(step), ncol = N)
#' ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
#' for(j in 1:N){
#'   allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])*dvm(angle[ind],0,kappa[j])
#' } 
#' -forward_g(delta, Gamma[,,tod], allprobs, ad = TRUE) +
#'   penalty(betaspline, S, lambda) # this does all the penalization work
#'}
#' # model fitting
#' mod = pql(pnll, par, dat, random = "betaspline")
pql = function(pnll, # penalized negative log-likelihood function
               par, # initial parameter list
               dat, # initial dat object, currently needs to be called dat!
               random, # names of parameters in par that are random effects/ penalized
               penalty = "lambda", # name given to the penalty parameter in dat
               alpha = 0, # exponential smoothing parameter
               maxiter = 100, # maximum number of iterations
               tol = 1e-5, # tolerance for convergence
               inner_tol = 1e-10, # tolerance for inner optimization
               silent = 1, # print level
               saveall = FALSE) # save all intermediate models?
{
  
  # setting the argument name for par because later updated par is returned
  argname_par = as.character(substitute(par))
  argname_dat = as.character(substitute(dat))
  
  # setting the environment for mllk to be the local environment such that it pull the right lambda
  environment(pnll) = environment() 
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re = length(random) 
  
  # list to save all model objects
  allmods = list() 
  
  # initial lambda locally
  lambda = dat[[penalty]]
  
  # experimentally, changing the name of the data object in pnll to dat
  if(argname_dat != "dat"){
    body(pnll) <- parse(text=gsub(argname_dat, "dat", deparse(body(pnll))))
  }
  
  # creating the objective function as wrapper around pnll to pull lambda from local
  f = function(par){
    environment(pnll) = environment()
    
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    getLambda = function(x) lambda
    
    dat[[penalty]] = DataEval(getLambda, rep(advector(1), 0))
    
    pnll(par)
  }
  
  # creating the RTMB objective function
  if(silent %in% 0:1) cat("Creating AD function\n")
  
  obj = MakeADFun(func = f, parameters = par, silent = TRUE) # silent and replacing with own prints
  newpar = obj$par # saving initial paramter value as vector to initialize optimization in loop
  
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
  
  # finding the indices of the random effects to later index Hessian
  re_inds = list() 
  for(i in 1:n_re){
    re_dim = dim(as.matrix(par[[random[i]]]))
    re_inds[[i]] = matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])
    if(dim(re_inds[[i]])[2] == 1) re_inds[[i]] = t(re_inds[[i]]) # if only one column, then transpose
  }
  
  # get number of similar random effects for each distinct random effect (of same structure)
  re_lengths = sapply(re_inds, function(x) if (is.vector(x)) 1 else nrow(x))
  
  # initialize list of penalty strength parameters
  Lambdas = list()
  Lambdas[[1]] = reshape_lambda(re_lengths, lambda) # reshaping to match structure of random effects
  
  if(silent < 2){
    cat("Initializing with", paste0(penalty, ":"), round(lambda, 3), "\n")
  }
  
  # computing rank deficiency for each penalty matrix to use in correction term
  m = numeric(length(S)) 
  for(i in 1:length(m)) {
    m[i] = nrow(S[[i]]) - Matrix::rankMatrix(S[[i]])
  } 
  
  ### updating algorithm
  # loop over outer iterations until convergence or maxiter
  for(k in 1:maxiter){
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    opt = stats::optim(newpar, obj$fn, newgrad, 
                       method = "BFGS", 
                       control = list(reltol = inner_tol, maxit = 1000))
    
    # setting new optimum par for next iteration
    newpar = opt$par 
    
    # reporting to extract penalties
    mod = obj$report() 
    
    # evaluating current Hessian
    J = obj$he()
    
    # computing Fisher information matrix
    J_inv = MASS::ginv(J) 
    
    # saving entire model object
    allmods[[k]] = mod 
    
    ## updating all lambdas
    lambdas_k = list() # temporary lambda list
    
    # looping over distinct random effects (matrices)
    for(i in 1:n_re){
      # initializing lambda vector for i-th random effect
      lambdas_k[[i]] = numeric(nrow(re_inds[[i]]))
      
      # looping over similar random effects
      for(j in 1:nrow(re_inds[[i]])){
        idx = re_inds[[i]][j,] # indices of this random effect
        
        # effective degrees of freedom for this random effect: J^-1_p J
        edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
        
        # calculating new lambda based on updating rule
        lambda_new = as.numeric((edoF - m[i]) / mod$Pen[[i]][j]) # m is correction if S_i does not have full rank
        
        # potentially smoothing new lambda
        lambdas_k[[i]][j] = (1-alpha) * lambda_new + alpha * Lambdas[[k]][[i]][j]
      }
      
      # minimum of zero for penalty strengths
      lambdas_k[[i]][which(lambdas_k[[i]] < 0)] = 0
    }
    
    # assigning new lambda to global list
    Lambdas[[k+1]] = lambdas_k
    
    # updating lambda vector locally for next iteration
    lambda = unlist(lambdas_k) 
    
    if(silent < 2){
      cat("outer", k, "-", paste0(penalty, ":"), round(lambda, 3), "\n")
    }
    
    # convergence check
    if(max(abs(lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol){
      if(silent < 2){
        cat("Converged\n")
      }
      break
    }
    
    if(k == maxiter) warning("No convergence\n")
  }
  
  # assign RTMB obj to return object
  mod$obj <- obj
  
  # if all intermediate models should be returned, assign
  if(saveall) {
    mod$allmods = allmods
  }
  
  # assign final lambda to return object
  mod[[penalty]] = lambda
  
  # assigning all lambdas to return object
  mod[[paste0("all_", penalty)]] = Lambdas
  
  # calculating unpenalized log-likelihood at final parameter values
  lambda = rep(0, length(lambda))
  dat[[penalty]] = lambda
  
  # format parameter to list
  skeleton = utils::as.relistable(par)
  parlist = utils::relist(opt$par, skeleton)
  mod[[argname_par]] = parlist # and assing to return object
  
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
  mod$n_fixpar = length(unlist(par[!(names(par) %in% random)]))
  
  # assing conditinoal Hessian
  mod$Hessian_conditional = J
  
  # removing penalty list from model object
  mod = mod[names(mod) != "Pen"] 
  
  #############################
  ### constructing joint object
  parlist$loglambda = log(mod[[penalty]])
  
  # finding the number of similar random effects for each random effect
  # indvec = rep(1:n_re, times = re_lengths)
  
  # computing log determinants
  logdetS = numeric(length(S))
  for(i in 1:length(S)){
    logdetS[i] = determinant(S[[i]])$modulus
  }
  
  ## defining joint negative log-likelihood
  jnll = function(par) {
    
    environment(pnll) = environment()
    
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    dat[[penalty]] = exp(par$loglambda)
    
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
  
  # creating joint AD object
  obj_joint = MakeADFun(jnll, parlist,
                        random = names(par)[names(par) != "loglambda"]) # REML, everything random except lambda
  
  # assigning object to return object
  mod$obj_joint = obj_joint
  
  return(mod)
}
