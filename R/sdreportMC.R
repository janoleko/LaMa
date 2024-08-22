#' Monte Carlo version of sdreport
#' 
#' After optimization of an AD model, sdreportMC can be used to calculate samples of confidence intervals of all model parameters and \code{REPORT()}ed quantities
#' including non linear functions of random effects and parameters.
#'
#' @param obj Object returned by \code{MakeADFun()} after optimization
#' @param what Vector of strings with names of parameters and \code{REPORT()}ed quantities to be reported
#' @param Hessian Optional Hessian matrix. If not provided, it will be computed from the object
#' @param CI Logical. If \code{TRUE}, only confidence intervals instead of samples will be returned
#' @param n Number of samples to draw from the multivariate normal distribution of the MLE
#' @param probs Vector of probabilities for the confidence intervals (ignored if no CIs are computed)
#'
#' @return Named list corresponding to the elements of \code{what}. Each element has the structure of the corresponding quantity with an additional dimension added for the samples.
#' For example, if a quantity is a vector, the list contains a matrix. If a quantity is a matrix, the list contains an array.
#' @export
#' 
#' @import RTMB
#'
#' @examples
#' # currently no examples
sdreportMC = function(obj, what, Hessian = NULL, CI = FALSE, n = 1000, probs = c(0.025, 0.975)){
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
  
  # extract parameter list from sdreport
  parlist = as.list(sdr, "Estimate")
  # save parameter list skeleton to reshape later
  skeleton = as.relistable(parlist)
  # save parameter names
  parnames = names(parlist)
  # turn estimated parameter list into vector
  parvec = unlist(parlist)
  # compute Fisher information
  I = MASS::ginv(as.matrix(H))
  
  # sample from multivariate normal distribution
  simpars = mvtnorm::rmvnorm(n, parvec, I)
  
  ## check which elements of what are parameters and which are reported
  whatpar = what[what %in% names(parlist)]
  
  report0 = obj$report()
  whatreport = what[what %in% names(report0)]
  
  ## first deal with parameters
  if(length(whatpar) > 0){
    if(length(whatpar) == 1){
      Samplepars1 = lapply(1:n, function(t){
        samplelist = relist(simpars[t,], skeleton)
        samplelist[[whatpar]]})
      
      Samplepars = setNames(vector("list", length(whatpar)), whatpar)
      Samplepars[[1]] = Samplepars1
    } else{
      Samplepars1 = lapply(1:n, function(t){
        samplelist = relist(simpars[t,], skeleton)
        samplelist[whatpar]})
      
      # create list skeleton, first layer: what, second layer: samples
      Samplepars = setNames(vector("list", length(whatpar)), whatpar)
      
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
      Samplereports = setNames(vector("list", length(whatreport)), whatreport)
      Samplereports[[1]] = Samplereports1
    } else {
      Samplereports1 = lapply(1:n, function(t){
        reportlist = obj$report(simpars[t,])
        reportlist[whatreport]})
      
      # create list skeleton, first layer: what, second layer: samples
      Samplereports = setNames(vector("list", length(whatreport)), whatreport)
      
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
          qthisSamplepars = quantile(thisSamplepars, probs)
        } else if(is.matrix(thisSamplepars)){
          qthisSamplepars = apply(thisSamplepars, 2, quantile, probs)
        } else if(is.array(thisSamplepars) & !is.matrix(thisSamplepars)){
          qthisSamplepars = apply(thisSamplepars, 1:(length(dim(thisSamplepars))-1), quantile, probs)
        }
        CIpar[[i]] = qthisSamplepars
      }
      names(CIpar) = names(Samplepars)
    } else{CIpar = NULL}
    if(!is.null(Samplereports)){
      CIreport = list()
      for(i in 1:length(Samplereports)){
        thisSamplereports = Samplereports[[i]]
        if(is.vector(thisSamplereports)){
          qthisSamplereports = quantile(thisSamplereports, probs)
        } else if(is.matrix(thisSamplereports)){
          qthisSamplereports = apply(thisSamplereports, 2, quantile, probs)
        } else if(is.array(thisSamplepars) & !is.matrix(thisSamplepars)){
          qthisSamplereports = apply(thisSamplereports, 1:(length(dim(thisSamplereports))-1), quantile, probs)
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
