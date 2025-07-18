% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_matrix_functions.R
\name{process_hid_formulas}
\alias{process_hid_formulas}
\title{Process and standardise formulas for the state process of hidden Markov models}
\usage{
process_hid_formulas(formulas, nStates, ref = NULL)
}
\arguments{
\item{formulas}{formulas for the transition process of a hidden Markov model, either as a single formula, a list of formulas, or a matrix.}

\item{nStates}{number of states of the Markov chain}

\item{ref}{optional vector of reference categories for each state, defaults to \code{1:nStates}. 
If provided, must be of length \code{nStates} and contain valid state indices. 
If a formula matrix is provided, this cannot be specified because reference categries are specified by one \code{"."} entry in each row.}
}
\value{
named list of formulas of length \code{nStates * (nStates - 1)}, where each formula corresponds to a transition from state \eqn{i} to state \eqn{j}, excluding transitions from \code{i} to \code{ref[i]}.
}
\description{
Process and standardise formulas for the state process of hidden Markov models
}
\examples{
# single formula for all non-reference category elements
formulas = process_hid_formulas(~ s(x), nStates = 3)
# now a list of length 6 with names tr.ij, not including reference categories

# different reference categories
formulas = process_hid_formulas(~ s(x), nStates = 3, ref = c(1,1,1))

# different formulas for different entries (and only for 2 of 6)
formulas = list(tr.12 ~ s(x), tr.23 ~ s(y))
formulas = process_hid_formulas(formulas, nStates = 3, ref = c(1,1,1))
# also a list of length 6, remaining entries filled with tr.ij ~ 1

# matrix input with reference categories
formulas = matrix(c(".", "~ s(x)", "~ s(y)",
                    "~ g", ".", "~ I(x^2)",
                    "~ y", "~ 1", "."), 
                    nrow = 3, byrow = TRUE)
# dots define reference categories
formulas = process_hid_formulas(formulas, nStates = 3)
}
