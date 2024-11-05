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
#' @param obs vector of continuous-valued observations (of length n)
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(n,N), where N is the number of states) as computed by \code{\link{stateprobs}}, \code{\link{stateprobs_g}} or \code{\link{stateprobs_p}}
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson)
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
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
#' pres = pseudo_res(obs, stateprobs, "norm", par)
#'
#' ## discrete-valued observations
#' obs = rpois(100, lambda = 1)
#' stateprobs = matrix(0.5, nrow = 100, ncol = 2)
#' par = list(lambda = c(1,2))
#' pres = pseudo_res(obs, stateprobs, "pois", par)
pseudo_res = function(obs, 
                      stateprobs, 
                      dist, 
                      par,
                      normal = TRUE,
                      discrete = NULL, 
                      randomise = TRUE,
                      seed = NULL) {
  
  # if discrete is not specified, try to determine
  if(is.null(discrete)){
    discrete = dist %in% c("pois", "binom", "geom", "nbinom")
  }
  
  if(discrete){
    cat("Discrete pseudo-residuals are calculated\n")
    
    return(pseudo_res_discrete(obs, 
                               stateprobs, 
                               dist, 
                               par, 
                               normal,
                               randomise,
                               seed))
  } else{
    # Number of observations and number of states
    nObs <- length(obs)              # Length of observations
    N <- ncol(stateprobs)            # Number of states (columns in stateprobs)
    
    # Check that the number of rows in `stateprobs` matches the length of `obs`
    if (nrow(stateprobs) != nObs) {
      stop("The number of rows in `stateprobs` must match the length of `obs`.")
    }
    
    # Construct the CDF function name dynamically, e.g., "pnorm" for "norm"
    cdf_name <- paste0("p", dist)
    cdf_func <- get(cdf_name, mode = "function")
    
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
      cdf_values[, state] <- do.call(cdf_func, args = c(list(obs), current_par))
    }
    
    # Compute pseudo-residuals by weighting CDF values with state probabilities
    residuals <- rowSums(cdf_values * stateprobs)
    
    if(normal){
      return(qnorm(residuals))
    } else{
      return(residuals)
    }
  }
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
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(n,N), where N is the number of states)
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson)
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
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
#' pres = pseudo_res_discrete(obs, stateprobs, "pois", par)
pseudo_res_discrete <- function(obs, 
                                stateprobs, 
                                dist, 
                                par,
                                normal = TRUE,
                                randomise = TRUE,
                                seed = NULL) {
  # Number of observations and number of states
  nObs <- length(obs)              # Length of observations
  N <- ncol(stateprobs)            # Number of states (columns in stateprobs)
  
  # Check that the number of rows in `stateprobs` matches the length of `obs`
  if (nrow(stateprobs) != nObs) {
    stop("The number of rows in `stateprobs` must match the length of `obs`.")
  }
  
  # Check that each parameter in `par` has the correct length
  if (!all(sapply(par, length) == N)) {
    stop("Each entry in `par` must have a length equal to the number of columns in `stateprobs`.")
  }
  
  # Construct the CDF function name dynamically, e.g., "pnorm" for "norm"
  cdf_name <- paste0("p", dist)
  cdf_func <- get(cdf_name, mode = "function")
  
  # Initialize a matrix to store CDF values for each observation and state
  cdf_values_lower <- cdf_values_upper <- matrix(0, nrow = nObs, ncol = N)
  cdf_values_random <- matrix(0, nrow = nObs, ncol = N)
  
  # Set the random seed if specified
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
    cdf_values_lower[, state] <- tryCatch(
      do.call(cdf_func, args = c(list(obs-1), current_par)),
      error = function(e) do.call(cdf_func, args = c(list(pmax(obs-1, 0)), current_par))
    )
    
    cdf_values_upper[, state] <- do.call(cdf_func, args = c(list(obs), current_par))
    
    if (randomise) {
      cdf_values_random[, state] <- sapply(1:nObs, function(i) {
        stats::runif(1, cdf_values_lower[i, state], cdf_values_upper[i, state])
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
