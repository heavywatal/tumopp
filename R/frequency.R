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

#' Make true spectrum from which allele frequencies are generated
#' @param descendants column of raw population including ancestors
#' @param lowerbound to discard peripheral branches
#' @return numeric vector
#' @rdname frequency
#' @export
descendants_fs = function(descendants, lowerbound=0) {
    fs = descendants / max(descendants)
    fs[(lowerbound < fs) & (fs < 1)]
}
