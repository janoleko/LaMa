#' Calculate pseudo-residuals
#' 
#' @description
#' For HMMs, pseudo-residuals are used to assess the goodness-of-fit of the model. 
#' These are based on the cumulative distribution function (CDF)
#' \deqn{F_{X_t}(x_t) = F(x_t \mid x_1, \dots, x_{t-1}, x_{t+1}, \dots, x_T)}
#' and can be used to quantify whether an observation is extreme relative to its model-implied distribution.
#' 
#' This function calculates such residuals via probability integral transform, based on the local state probabilities obtained by \code{\link{stateprobs}} or \code{\link{stateprobs_g}} and the respective parametric family.
#'
#' @details
#' When used for discrete pseudo-residuals, this function is just a wrapper for \code{\link{pseudo_res_discrete}}.
#'
#' @seealso \code{\link{plot.LaMaResiduals}} for plotting pseudo-residuals.
#' 
#' @param obs vector of continuous-valued observations (of length n)
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson) or CDF function to evaluate directly.
#' If a discrete CDF is passed, the \code{discrete} argument needs to be set to \code{TRUE} because this cannot determined automatically.
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(n,N), where N is the number of states) as computed by \code{\link{stateprobs}}, \code{\link{stateprobs_g}} or \code{\link{stateprobs_p}}
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}}, \code{\link{forward_g}} or \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function and avoid calculating local state proabilities manually.
#' In this case, a call should look like \code{pseudo_res(obs, "norm", par, mod = mod)}.
#' @param normal logical, if \code{TRUE}, returns Gaussian pseudo residuals
#'
#' These will be approximately standard normally distributed if the model is correct.
#' @param discrete logical, if \code{TRUE}, computes discrete pseudo residuals (which slightly differ in their definition)
#'
#' By default, will be determined using \code{dist} argument, but only works for standard discrete distributions.
#' When used with a special discrete distribution, set to \code{TRUE} manually. See \code{\link{pseudo_res_discrete}} for details.
#' @param randomise for discrete pseudo residuals only. Logical, if \code{TRUE}, return randomised pseudo residuals. Recommended for discrete observations.
#' @param seed for discrete pseudo residuals only. Integer, seed for random number generation
#'
#' @return vector of pseudo residuals
#' @export
#'
#' @examples
#' ## continuous-valued observations
#' obs = rnorm(100)
#' stateprobs = matrix(0.5, nrow = 100, ncol = 2)
#' par = list(mean = c(1,2), sd = c(1,1))
#' pres = pseudo_res(obs, "norm", par, stateprobs)
#'
#' ## discrete-valued observations
#' obs = rpois(100, lambda = 1)
#' par = list(lambda = c(1,2))
#' pres = pseudo_res(obs, "pois", par, stateprobs)
#' 
#' ## custom CDF function
#' obs = rnbinom(100, size = 1, prob = 0.5)
#' par = list(size = c(0.5, 2), prob = c(0.4, 0.6))
#' pres = pseudo_res(obs, pnbinom, par, stateprobs, 
#'                   discrete = TRUE)
#' # if discrete CDF function is passed, 'discrete' needs to be set to TRUE
#' 
#' ## no CDF available, only density (artificial example)
#' obs = rnorm(100)
#' par = list(mean = c(1,2), sd = c(1,1))
#' cdf = function(x, mean, sd) integrate(dnorm, -Inf, x, mean = mean, sd = sd)$value
#' pres = pseudo_res(obs, cdf, par, stateprobs)
#' 
#' ## full example with model object
#' step = trex$step[1:200]
#' 
#' nll = function(par){
#'   getAll(par)
#'   Gamma = tpm(logitGamma)
#'   delta = stationary(Gamma)
#'   mu = exp(logMu); REPORT(mu)
#'   sigma = exp(logSigma); REPORT(sigma)
#'   allprobs = matrix(1, length(step), 2)
#'   ind = which(!is.na(step))
#'   for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' par = list(logitGamma = c(-2,-2), 
#'            logMu = log(c(0.3, 2.5)), 
#'            logSigma = log(c(0.3, 0.5)))
#'            
#' obj = MakeADFun(nll, par)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = obj$report()
#' 
#' pres = pseudo_res(step, "gamma2", list(mean = mod$mu, sd = mod$sigma),
#'                   mod = mod)
#'
#' plot(pres)
pseudo_res = function(obs, 
                      dist, 
                      par,
                      stateprobs = NULL,
                      mod = NULL,
                      normal = TRUE,
                      discrete = NULL, 
                      randomise = TRUE,
                      seed = NULL) {
  
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$type)){
      stop("'mod' contains no type.")
    }
    if(!(mod$type) %in% c("homogeneous", "inhomogeneous", "periodic")){
      stop("'mod' contains invalid type.")
    }
    
    if(is.null(mod$delta)){
      stop("'mod' contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("'mod' contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("'mod' contains no state-dependent probabilities.")
    }
    
    if(mod$type == "periodic"){
      if(is.null(mod$tod)){
        stop("'mod' contains no cyclic indexing variable.")
      }
    }
    
    # calculate state probabilities based on model object
    if(mod$type == "homogeneous"){
      stateprobs = stateprobs(mod = mod)
    }
    if(mod$type == "inhomogeneous"){
      stateprobs = stateprobs_g(mod = mod)
    }
    if(mod$type == "periodic"){
      stateprobs = stateprobs_p(mod = mod)
    }
  }
  
  # if discrete is not specified, try to determine
  if(is.null(discrete)){
    if(is.character(dist)){
      discrete = dist %in% c("pois", "binom", "geom", "nbinom")
    } else if(is.function(dist)){
      discrete = FALSE
      message("Assuming 'dist' evaluates a continuous CDF. If discrete, please set 'discrete = TRUE'.")
    }
  }
  
  if(discrete){
    cat("Calculating discrete pseudo-residuals\n")
    
    residuals <- pseudo_res_discrete(obs, 
                                     dist,
                                     par, 
                                     stateprobs,
                                     normal,
                                     randomise,
                                     seed)
  } else{
    # Number of observations and number of states
    nObs <- length(obs)              # Length of observations
    N <- ncol(stateprobs)            # Number of states (columns in stateprobs)
    ind <- which(!is.na(obs))
    
    # Check that the number of rows in `stateprobs` matches the length of `obs`
    if (nrow(stateprobs) != nObs) {
      stop("The number of rows in 'stateprobs' must match the length of 'obs'.")
    }
    
    # Construct the CDF function name dynamically, e.g., "pnorm" for "norm"
    if(is.character(dist)){
      cdf_name <- paste0("p", dist)
      cdf_func <- get(cdf_name, mode = "function")
    } else if(is.function(dist)){
      cdf_func <- Vectorize(dist)
    } else{
      stop("'dist' must be a character string or a function.")
    }
    
    
    # Initialize a matrix to store CDF values for each observation and state
    cdf_values <- matrix(0, nrow = nObs, ncol = N)
    
    # Loop over each state to compute the CDF values for each observation
    for (state in 1:N) {
      
      # Extract the parameters for the current state
      # can be either vector
      current_par <- lapply(par, function(param){
        if(is.matrix(param)){
          if(ncol(param) != N | nrow(param) != nObs){
            stop("Parameter matrix dimensions must match number of observations and states")
          } else{
            return(param[, state])
          }
        } else if(is.vector(param)){
          if(length(param) == 1){
            return(param)
          } else{
            if(length(param) != N){
              stop("Parameter vector must have the same length as the number of states")
            } else{
              return(param[state])
            }
          }
        }
      })
      
      # evaluate the CDF function at each observation with these parameters
      cdf_values[ind, state] <- do.call(cdf_func, args = c(list(obs[ind]), current_par))
    }
    
    # Compute pseudo-residuals by weighting CDF values with state probabilities
    residuals <- rowSums(cdf_values * stateprobs)
    
    if(normal){
      residuals <- qnorm(residuals)
    }
  }
  
  # handle infinite value
  residuals[is.infinite(residuals)] <- NA
  
  class(residuals) <- "LaMaResiduals"
  
  return(residuals)
}


#' Calculate pseudo-residuals for discrete-valued observations
#' 
#' @description
#' For HMMs, pseudo-residuals are used to assess the goodness-of-fit of the model. 
#' These are based on the cumulative distribution function (CDF)
#' \deqn{F_{X_t}(x_t) = F(x_t \mid x_1, \dots, x_{t-1}, x_{t+1}, \dots, x_T)}
#' and can be used to quantify whether an observation is extreme relative to its model-implied distribution.
#' 
#' This function calculates such residuals for \strong{discrete-valued} observations, based on the local state probabilities obtained by \code{\link{stateprobs}} or \code{\link{stateprobs_g}} and the respective parametric family.
#'
#' @details
#' For discrete observations, calculating pseudo residuals is slightly more involved, as the CDF is a step function.
#' Therefore, one can calculate the lower and upper CDF values for each observation. 
#' By default, this function does exactly that and then randomly samples the interval in between to give approximately Gaussian psuedo-residuals.
#' If \code{randomise} is set to \code{FALSE}, the lower, upper and mean pseudo-residuasl are returned.
#'
#' @param obs vector of discrete-valued observations (of length n)
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson) or CDF function to evaluate directly.
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(n,N), where N is the number of states)
#' @param normal logical, if \code{TRUE}, returns Gaussian pseudo residuals
#'
#' These will be approximately standard normally distributed if the model is correct.
#' @param randomise logical, if \code{TRUE}, return randomised pseudo residuals. Recommended for discrete observations.
#' @param seed integer, seed for random number generation
#'
#' @return vector of pseudo residuals
#' @export
#'
#' @examples
#' obs = rpois(100, lambda = 1)
#' stateprobs = matrix(0.5, nrow = 100, ncol = 2)
#' par = list(lambda = c(1,2))
#' pres = pseudo_res_discrete(obs, "pois", par, stateprobs)
pseudo_res_discrete <- function(obs, 
                                dist,
                                par,
                                stateprobs,
                                normal = TRUE,
                                randomise = TRUE,
                                seed = NULL) {
  
  # Number of observations and number of states
  nObs <- length(obs)              # Length of observations
  N <- ncol(stateprobs)            # Number of states (columns in stateprobs)
  ind <- which(!is.na(obs))
  
  # Check that the number of rows in `stateprobs` matches the length of `obs`
  if (nrow(stateprobs) != nObs) {
    stop("The number of rows in 'stateprobs' must match the length of 'obs'.")
  }
  
  # Check that each parameter in `par` has the correct length
  if (!all(sapply(par, length) == N)) {
    stop("Each entry in 'par' must have a length equal to the number of columns in 'stateprobs'.")
  }
  
  if(is.character(dist)){
    # Construct the CDF function name dynamically, e.g., "ppois" for "pois"
    cdf_name <- paste0("p", dist)
    cdf_func <- get(cdf_name, mode = "function")
  } else if(is.function(dist)){
    cdf_func <- dist # if function, use this as CDF
  } else{
    stop("'dist' must be a character string or a function.")
  }
  
  # Initialize a matrix to store CDF values for each observation and state
  cdf_values_lower <- cdf_values_upper <- matrix(0, nrow = nObs, ncol = N)
  cdf_values_random <- matrix(0, nrow = nObs, ncol = N)
  
  # Set the random seed if specified #
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Loop over each state to compute the CDF values for each observation
  for (state in 1:N) {
    
    # Extract the parameters for the current state
    # can be either vector
    current_par <- lapply(par, function(param){
      if(is.matrix(param)){
        if(ncol(param) != N | nrow(param) != nObs){
          stop("Parameter matrix dimensions must match number of observations and states")
        } else{
          return(param[, state])
        }
      } else if(is.vector(param)){
        if(length(param) == 1){
          return(param)
        } else{
          if(length(param) != N){
            stop("Parameter vector must have the same length as the number of states")
          } else{
            return(param[state])
          }
        }
      }
    })
    
    
    # Use `mapply` to evaluate the CDF function at each observation with these parameters
    # safe evaluation of CDF at `obs - 1` in case of errors
    cdf_values_lower[ind, state] <- tryCatch(
      do.call(cdf_func, args = c(list(obs[ind] - 1), current_par)),
      error = function(e) do.call(cdf_func, args = c(list(pmax(obs[ind] - 1, 0)), current_par))
    )
    
    cdf_values_upper[ind, state] <- do.call(cdf_func, args = c(list(obs[ind]), current_par))
    
    if (randomise) {
      cdf_values_random[ind, state] <- sapply(1:length(ind), function(i) {
        stats::runif(1, cdf_values_lower[ind[i], state], cdf_values_upper[ind[i], state])
      })
    }
  }
  
  # Either return randomised pseudo-residuals or lower, upper, and mean
  if (randomise) {
    cat("Randomised between lower and upper\n")
    
    # Compute pseudo-residuals by weighting random CDF values with state probabilities
    residuals <- rowSums(cdf_values_random * stateprobs)
    
    if (normal) {
      return(qnorm(residuals))
    } else {
      return(residuals)
    }
  } else {
    # Compute pseudo-residuals by weighting lower and upper CDF values with state probabilities
    residuals_lower <- rowSums(cdf_values_lower * stateprobs)
    residuals_upper <- rowSums(cdf_values_upper * stateprobs)
    
    if (normal) {
      return(list(
        lower = qnorm(residuals_lower),
        upper = qnorm(residuals_upper),
        mean = qnorm((residuals_lower + residuals_upper) / 2)
      ))
    } else {
      return(list(
        lower = residuals_lower,
        upper = residuals_upper,
        mean = (residuals_lower + residuals_upper) / 2
      ))
    }
  }
}


#' Plot pseudo-residuals
#' 
#' @description
#' Plot pseudo-residuals computed by \code{\link{pseudo_res}}.
#' 
#' @param x pseudo-residuals as returned by \code{\link{pseudo_res}}
#' @param hist logical, if \code{TRUE}, adds a histogram of the pseudo-residuals
#' @param col character, color for the QQ-line (and density curve if \code{histogram = TRUE})
#' @param lwd numeric, line width for the QQ-line (and density curve if \code{histogram = TRUE})
#' @param main optional character vector of main titles for the plots of length 2 (or 3 if \code{histogram = TRUE})
#' @param ... currently ignored. For method consistency
#' 
#' @returns NULL, plots the pseudo-residuals in a 2- or 3-panel layout
#' @export
#'
#' @importFrom graphics par lines hist abline
#' @importFrom stats acf na.pass qqnorm
#' @importFrom RTMB dnorm
#'
#' @examples
#' ## pseudo-residuals for the trex data
#' step = trex$step[1:200]
#' 
#' nll = function(par){
#'   getAll(par)
#'   Gamma = tpm(logitGamma)
#'   delta = stationary(Gamma)
#'   mu = exp(logMu); REPORT(mu)
#'   sigma = exp(logSigma); REPORT(sigma)
#'   allprobs = matrix(1, length(step), 2)
#'   ind = which(!is.na(step))
#'   for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' par = list(logitGamma = c(-2,-2), 
#'            logMu = log(c(0.3, 2.5)), 
#'            logSigma = log(c(0.3, 0.5)))
#'            
#' obj = MakeADFun(nll, par)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = obj$report()
#' 
#' pres = pseudo_res(step, "gamma2", list(mean = mod$mu, sd = mod$sigma),
#'                   mod = mod)
#'                   
#' plot(pres)
plot.LaMaResiduals <- function(
    x, 
    hist = TRUE,
    col = "darkblue", 
    lwd = 1.5,
    main = NULL,
    ...
    ) {
  
  # Extract mean if residuals is a list
  res <- if (is.list(x)) x$mean else x
  res_clean <- na.omit(res)
  
  columns <- if (hist) 3 else 2
  
  if (is.null(main)) main <- rep("", columns)
  if (length(main) < columns) {
    main <- c(main, rep("", columns - length(main)))
  }
  
  old_par <- par(mfrow = c(1, columns))
  on.exit(par(old_par))
  
  # QQ Plot
  qqnorm(res_clean, main = main[1], 
         xlab = "theoretical quantiles", ylab = "sample quantiles",
         bty = "n", pch = 16, col = "#00000070")
  # qqline(res_clean, col = col, lwd = lwd)
  abline(a = 0, b = 1, col = col, lwd = lwd)
  
  # Histogram with normal curve
  if (hist) {
    r <- range(res_clean)
    dr <- diff(r)
    xgrid <- seq(r[1] - dr / 4, r[2] + dr / 4, length.out = 200)
    dens <- dnorm(xgrid)
    ylim <- c(0, max(dens) * 1.1)
    hist(res_clean, main = main[2], border = "white", ylim = ylim,
         prob = TRUE, xlab = "pseudo-residuals", ylab = "density")
    # curve(dnorm(x), add = TRUE, col = col, lwd = lwd)
    lines(xgrid, dens, col = col, lwd = lwd)
  }
  
  # ACF
  acf(res, na.action = na.pass, main = main[columns],
      xlab = "lag", bty = "n")
}


