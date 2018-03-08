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

#' Calculate expected allele frequencies in the extant cells
#' @inheritParams igraph::ego_size
#' @return integer
#' @rdname frequency
#' @export
allele_freqs = function(graph, nodes) {
  stopifnot("1" %in% nodes)
  x = igraph::ego_size(graph, order = 1073741824L, nodes = nodes, mode = "out")
  x = as.integer(x + 1L) %/% 2L
  x / max(x)
}

#' Extract cells whose expected allele frequencies are above threshold
#' @param population tibble with descendants and id column
#' @param threshold lowerbound of detectable allele frequency
#' @return integer IDs
#' @rdname frequency
#' @export
detectable_mutants_all = function(population, threshold) {
  population$id[population$allelefreq > threshold]
}

#' @rdname frequency
#' @export
detectable_mutants = function(graph, nodes, threshold) {
  n = length(nodes)
  counts = paths_to_origin(graph, nodes) %>%
    purrr::flatten_int() %>%
    table()
  counts[(counts / n) > threshold] %>% names() %>% as.integer()
}
