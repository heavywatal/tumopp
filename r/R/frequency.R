#' Calculate MATH (mutant-allele tumor heterogeneity) score from variant allele frequencies
#' @inheritParams stats::mad
#' @return numeric
#' @rdname frequency
#' @export
math_score = function(x, constant=1.4826, na.rm=FALSE) {
  med = stats::median(x, na.rm = na.rm)
  mad = stats::mad(x, center = med, constant = constant, na.rm = na.rm)
  mad / med
}

#' Count number of extant descendants
#' @inheritParams igraph::ego_size
#' @return integer
#' @rdname frequency
#' @export
count_descendants = function(graph, nodes) {
  x = igraph::ego_size(graph, order = 1073741824L, nodes = nodes, mode = "out")
  as.integer(x + 1L) %/% 2L
}

#' Translate descendant numbers to frequencies
#' @param descendants column of raw population including ancestors
#' @param lowerbound to discard peripheral branches
#' @return numeric vector applicable to math_score
#' @rdname frequency
#' @export
descfreqs = function(descendants, lowerbound=0) {
  frac = descendants / max(descendants, na.rm = TRUE)
  dplyr::if_else(frac < lowerbound, NA_real_, frac)
}
