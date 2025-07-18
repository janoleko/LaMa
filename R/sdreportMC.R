#' Monte Carlo version of \code{sdreport}
#' 
#' After optimisation of an AD model, \code{sdreportMC} can be used to calculate samples of confidence intervals of all model parameters and \code{REPORT()}ed quantities
#' including nonlinear functions of random effects and parameters.
#' 
#' @details
#' This function simply samples from the approximate multivariate normal distribution of the maximum likelihood estimate (MLE) of the parameters
#' \deqn{ \hat{\theta} \sim N(\theta, H^{-1}), }
#' where \eqn{H} is the Hessian matrix of the negative log-likelihood function at the MLE. 
#' It then returns either the sampled parameters or \code{REPORT()}ed transformations of them. If \code{CI = TRUE}, it does not return the samples but directly computes confidence intervals.
#' 
#' If you are interested in several quantities, calling \code{sdreportMC} once with a vector \code{what} will generally be faster than calling it several times with single elements of \code{what}.
#'
#' @param obj object returned by \code{MakeADFun()} after optimisation or model of class \code{qremlModel} as returned by \code{\link{qreml}}.
#' @param what vector of strings with names of parameters and/or \code{REPORT()}ed quantities to be reported
#' @param nSamples number of samples to draw from the multivariate normal distribution of the MLE
#' @param Hessian optional Hessian matrix. If not provided, it will be computed from the object
#' @param CI logical. If \code{TRUE}, only confidence intervals instead of samples will be returned
#' @param probs vector of probabilities for the confidence intervals (ignored if no CIs are computed)
#'
#' @return named list corresponding to the elements of \code{what}. Each element has the structure of the corresponding quantity with an additional dimension added for the samples.
#' For example, if a quantity is a vector, the list contains a matrix. If a quantity is a matrix, the list contains an array. If quantity is an array, the list contains an array with one extra dimension.
#' @export
#' 
#' @import RTMB
#' @importFrom mgcv rmvn
#' @importFrom MASS ginv
#'
#' @examples
#' # fitting an HMM to the trex data and running sdreportMC
#' ## negative log-likelihood function
#' nll = function(par) {
#'   getAll(par, dat) # makes everything contained available without $
#'   Gamma = tpm(eta) # computes transition probability matrix from unconstrained eta
#'   delta = stationary(Gamma) # computes stationary distribution
#'   # exponentiating because all parameters strictly positive
#'   mu = exp(logmu)
#'   sigma = exp(logsigma)
#'   kappa = exp(logkappa)
#'   # reporting statements for sdreportMC
#'   REPORT(mu)
#'   REPORT(sigma)
#'   REPORT(kappa)
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
#'   for(j in 1:N){
#'     allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])*dvm(angle[ind],0,kappa[j])
#'   }
#'   -forward(delta, Gamma, allprobs) # simple forward algorithm
#' }
#'
#' ## initial parameter list
#' par = list(
#'  logmu = log(c(0.3, 1)),       # initial means for step length (log-transformed)
#'   logsigma = log(c(0.2, 0.7)), # initial sds for step length (log-transformed)
#'   logkappa = log(c(0.2, 0.7)), # initial concentration for turning angle (log-transformed)
#'   eta = rep(-2, 2)             # initial t.p.m. parameters (on logit scale)
#' )   
#' ## data and hyperparameters
#' dat = list(
#'   step = trex$step[1:500],   # hourly step lengths
#'   angle = trex$angle[1:500], # hourly turning angles
#'   N = 2
#' )
#'
#' ## creating AD function
#' obj = MakeADFun(nll, par, silent = TRUE) # creating the objective function
#'
#' ## optimising
#' opt = nlminb(obj$par, obj$fn, obj$gr) # optimization
#'
#' ## running sdreportMC
#' # `mu` has report statement, `delta` is automatically reported by `forward()`
#' sdrMC = sdreportMC(obj, 
#'                    what = c("mu", "delta"), 
#'                    nSamples = 50)
#' dim(sdrMC$delta)
#' # now a matrix with 50 samples (rows)
sdreportMC = function(obj, 
                      what, 
                      nSamples = 1000,
                      Hessian = NULL, 
                      CI = FALSE, 
                      probs = c(0.025, 0.975)){
  n <- nSamples
  
  # check whether obj is either a qremlModel or an RTMB object
  if(inherits(obj, "qremlModel")){
    obj <- obj$obj # qremlModel has RTMB obj as entry named obj
  }
  if(is.null(obj$env)){
    stop("'obj' seems not to be generated by 'RTMB::MakeADFun()'")
  }
  
  if(is.null(Hessian)){
    # check if Hessian can be computed
    # if not, obj is a marginal llk -> no Hessian implemented
    H = tryCatch(obj$he(), error = function(e) "no Hessian")
    
    # if marginal llk, compute joint Hessian
    if(!is.matrix(H)){
      cat("Computing joint Hessian...\n")
      
      sdr = sdreport(obj, getJointPrecision = TRUE) # takes time
      H = sdr$jointPrecision
    } else{ # else just run sdreport to get point estimate from obj
      sdr = sdreport(obj, ignore.parm.uncertainty = TRUE)
    }
  } else{
    H = Hessian
    sdr = sdreport(obj, ignore.parm.uncertainty = TRUE)
  }
  
  # extract relisting function
  relist_par = obj$env$parList
  # extract parameter list from sdreport
  parlist = as.list(sdr, "Estimate")
  # save parameter list skeleton to reshape later
  # skeleton = utils::as.relistable(parlist)
  # save parameter names
  parnames = names(parlist)
  # turn estimated parameter list into vector
  # parvec <- obj$env$par
  parvec <- sdr$par.fixed
  # parvec = unlist(parlist)
  # compute Fisher information
  I = ginv(as.matrix(H))
  
  # sample from multivariate normal distribution
  #simpars = mvtnorm::rmvnorm(n, parvec, I)
  simpars <- rmvn(n, parvec, I)
  
  ## check which elements of what are parameters and which are reported
  whatpar = what[what %in% names(parlist)]
  
  report0 = obj$report()
  whatreport = what[what %in% names(report0)]
  
  ## first deal with parameters
  if(length(whatpar) > 0){
    if(length(whatpar) == 1){
      Samplepars1 = lapply(1:n, function(t){
        # samplelist = utils::relist(simpars[t,], skeleton)
        samplelist = relist_par(simpars[t,])
        samplelist[[whatpar]]})
      
      Samplepars = stats::setNames(vector("list", length(whatpar)), whatpar)
      Samplepars[[1]] = Samplepars1
    } else{
      Samplepars1 = lapply(1:n, function(t){
        # samplelist = utils::relist(simpars[t,], skeleton)
        samplelist = relist_par(simpars[t,])
        samplelist[whatpar]})
      
      # create list skeleton, first layer: what, second layer: samples
      Samplepars = stats::setNames(vector("list", length(whatpar)), whatpar)
      
      # reshaping
      for (name in whatpar) {
        Samplepars[[name]] = lapply(Samplepars1, function(sublist) sublist[[name]])
      }
    }
    
    # reshaping elements in Samplepars
    for(i in 1:length(Samplepars)){
      thisSamplepars = Samplepars[[i]]
      # check if vector etc
      if(is.vector(thisSamplepars[[1]])){
        if(length(thisSamplepars[[1]]) == 1){
          # if only scalars, create vector
          thisSamplepars = unlist(thisSamplepars)
        } else{
          # if vectors, create matrix
          thisSamplepars = do.call(rbind, thisSamplepars)
        }
      } else if(is.array(thisSamplepars[[1]])){
        # if arrays, create array
        thisSamplepars = array(unlist(thisSamplepars), dim = c(dim(thisSamplepars[[1]]), length(thisSamplepars)))
      }
      Samplepars[[i]] = thisSamplepars
    }
  } else{
    Samplepars = NULL
  }
  
  ## then deal with reported quantities
  if(length(whatreport) > 0){
    cat("Sampling reported quantities...\n")
    if(length(whatreport) == 1){
      Samplereports1 = lapply(1:n, function(t) obj$report(simpars[t,])[[whatreport]])
      Samplereports = stats::setNames(vector("list", length(whatreport)), whatreport)
      Samplereports[[1]] = Samplereports1
    } else {
      Samplereports1 = lapply(1:n, function(t){
        reportlist = obj$report(simpars[t,])
        reportlist[whatreport]})
      
      # create list skeleton, first layer: what, second layer: samples
      Samplereports = stats::setNames(vector("list", length(whatreport)), whatreport)
      
      # reshaping
      for (name in whatreport) {
        Samplereports[[name]] = lapply(Samplereports1, function(sublist) sublist[[name]])
      }
    }
    
    # reshaping elements in Samplereports
    for(i in 1:length(Samplereports)){
      thisSamplereports = Samplereports[[i]]
      # check if vector etc
      if(is.vector(thisSamplereports[[1]])){
        if(length(thisSamplereports[[1]]) == 1){
          # if only scalars, create vector
          thisSamplereports = unlist(thisSamplereports)
        } else{
          # if vectors, create matrix
          thisSamplereports = do.call(rbind, thisSamplereports)
        }
      } else if(is.array(thisSamplereports[[1]])){
        # if arrays, create array
        thisSamplereports = array(unlist(thisSamplereports), dim = c(dim(thisSamplereports[[1]]), length(thisSamplereports)))
      }
      Samplereports[[i]] = thisSamplereports
    }
  } else{
    Samplereports = NULL
  }
  
  # if 
  if(CI == FALSE){
    if(is.null(Samplepars)){
      return(Samplereports)
    } else if(is.null(Samplereports)){
      return(Samplepars)
    } else{
      return(list(par = Samplepars, report = Samplereports))
    }
    
  } else{
    if(!is.null(Samplepars)){
      CIpar = list()
      for(i in 1:length(Samplepars)){
        thisSamplepars = Samplepars[[i]]
        if(is.vector(thisSamplepars)){
          qthisSamplepars = stats::quantile(thisSamplepars, probs)
        } else if(is.matrix(thisSamplepars)){
          qthisSamplepars = apply(thisSamplepars, 2, stats::quantile, probs)
        } else if(is.array(thisSamplepars) & !is.matrix(thisSamplepars)){
          qthisSamplepars = apply(thisSamplepars, 1:(length(dim(thisSamplepars))-1), stats::quantile, probs)
        }
        CIpar[[i]] = qthisSamplepars
      }
      names(CIpar) = names(Samplepars)
    } else{
      CIpar = NULL
    }
    
    if(!is.null(Samplereports)){
      CIreport = list()
      for(i in 1:length(Samplereports)){
        thisSamplereports = Samplereports[[i]]
        if(is.vector(thisSamplereports)){
          qthisSamplereports = stats::quantile(thisSamplereports, probs)
        } else if(is.matrix(thisSamplereports)){
          qthisSamplereports = apply(thisSamplereports, 2, stats::quantile, probs)
        } else if(is.array(thisSamplereports) & !is.matrix(thisSamplereports)){
          qthisSamplereports = apply(thisSamplereports, 1:(length(dim(thisSamplereports))-1), stats::quantile, probs)
        }
        CIreport[[i]] = qthisSamplereports
      }
      names(CIreport) = names(Samplereports)
    } else{CIreport = 0}
    
    if(is.null(CIpar)){
      return(CIreport)
    } else if(is.null(CIreport)){
      return(CIpar)
    } else{
      return(list(CIpar = CIpar, CIreport = CIreport))
    }
  }
}
