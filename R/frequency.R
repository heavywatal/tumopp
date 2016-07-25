#' Calculate MATH (mutant-allele tumor heterogeneity) score from variant allele freqneucies
#' @inheritParams stats::mad
#' @return numeric
#' @rdname frequency
#' @export
math_score = function(x, constant=1.4826, na.rm=FALSE) {
    med = stats::median(x)
    mad = stats::mad(x, center=med, constant=constant, na.rm=FALSE)
    mad / med
}
