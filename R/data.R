#' Elephant Movement Data
#'
#' Synthetic data set of hourly step lengths and turning angles of an elephant.
#'
#' @format A data frame with 10.000 rows and 4 variables:
#' \describe{
#'   \item{tod}{time of day variable ranging from 1 to 24}
#'   \item{step}{hourly step lengths in kilometres}
#'   \item{angle}{hourly turning angles in radians}
#'   \item{state}{hidden state variable}
#' }
#' @source Generated for example purposes.
"elephant"

#' Shark Acceleration Data
#'
#' Synthetic data set of minutely overall dynamic body acceleration (ODBA) of a shark.
#'
#' @format A data frame with 5.000 rows and 3 variables:
#' \describe{
#'   \item{ODBA}{overall dynamci body acceleration}
#'   \item{logODBA}{logarithm of overall dynamic body acceleration}
#'   \item{state}{hidden state variable}
#' }
#' @source Generated for example purposes.
"shark"