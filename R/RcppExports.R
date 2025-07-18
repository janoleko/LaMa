# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

forward_cpp_h <- function(allprobs, delta, Gamma) {
    .Call(`_LaMa_forward_cpp_h`, allprobs, delta, Gamma)
}

forward_cpp_g <- function(allprobs, delta, Gamma) {
    .Call(`_LaMa_forward_cpp_g`, allprobs, delta, Gamma)
}

forward_cpp_g_tracks <- function(allprobs, Delta, Gamma, trackInd) {
    .Call(`_LaMa_forward_cpp_g_tracks`, allprobs, Delta, Gamma, trackInd)
}

forward_cpp_h_tracks <- function(allprobs, Delta, Gamma, trackInd) {
    .Call(`_LaMa_forward_cpp_h_tracks`, allprobs, Delta, Gamma, trackInd)
}

forward_cpp_p <- function(allprobs, delta, Gamma, tod) {
    .Call(`_LaMa_forward_cpp_p`, allprobs, delta, Gamma, tod)
}

forward_cpp_s <- function(allprobs, delta, Gamma, agsizes) {
    .Call(`_LaMa_forward_cpp_s`, allprobs, delta, Gamma, agsizes)
}

forward_cpp_sp <- function(allprobs, delta, Gamma, agsizes, tod) {
    .Call(`_LaMa_forward_cpp_sp`, allprobs, delta, Gamma, agsizes, tod)
}

logalpha_cpp <- function(allprobs, delta, Gamma) {
    .Call(`_LaMa_logalpha_cpp`, allprobs, delta, Gamma)
}

logbeta_cpp <- function(allprobs, Gamma) {
    .Call(`_LaMa_logbeta_cpp`, allprobs, Gamma)
}

viterbi_g_cpp <- function(allprobs, delta, Gamma) {
    .Call(`_LaMa_viterbi_g_cpp`, allprobs, delta, Gamma)
}

rep_times <- function(x, times) {
    .Call(`_LaMa_rep_times`, x, times)
}

tpm_g_cpp <- function(Z, beta, N, byrow) {
    .Call(`_LaMa_tpm_g_cpp`, Z, beta, N, byrow)
}

tpm_g2_cpp <- function(Eta, N, byrow, ref) {
    .Call(`_LaMa_tpm_g2_cpp`, Eta, N, byrow, ref)
}

semigroup_cpp <- function(Q, times) {
    .Call(`_LaMa_semigroup_cpp`, Q, times)
}

tpm_thinned_t_cpp <- function(Gamma, t) {
    .Call(`_LaMa_tpm_thinned_t_cpp`, Gamma, t)
}

