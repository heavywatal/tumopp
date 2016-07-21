#' Calculate MATH score from variant allele freqneucies
#' @param freqs numeric vector
#' @return numeric
#' @rdname frequency
#' @export
math_score = function(freqs) {
    med = stats::median(freqs)
    mad = stats::median(abs(freqs - med))
    mad / med
}
