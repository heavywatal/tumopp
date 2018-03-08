#' Make ms-like neutral variation matrix
#' @param population tibble of raw population
#' @param nsam number of cells to sample
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @return int matrix [nsam x segsites]
#' @rdname ms
#' @export
sprinkle_mutations = function(population, nsam, mu=NULL, segsites=NULL) {
  graph = make_igraph(population)
  samples = sample(filter_extant(population)$id, nsam) %>% as.character()
  lineages = igraph::neighborhood(graph, order = 1073741824L, nodes = samples, mode = "in") %>%
    purrr::map(~ .x$name)
  common_ancs = purrr::reduce(lineages, intersect)
  nodes = purrr::flatten_chr(lineages) %>% unique() %>% setdiff(common_ancs)
  if (is.null(segsites)) {
    if (is.null(mu)) stop("specify either mu or segsites")
    segsites = stats::rpois(1L, length(nodes) * mu)
  } else if (!is.null(mu)) warning("mu is ignored if segsites is given")
  mutant_ids = sample(nodes, segsites, replace = TRUE) %>% sort()
  lineages %>%
    purrr::map(~ mutant_ids %in% .) %>%
    purrr::flatten_int() %>%
    matrix(nrow = nsam, byrow = TRUE)
}

#' Summary statistics of neutral mutations
#' @param msout tibble from `sprinkle_mutations()`
#' @param lowerbound to discard peripheral branches
#' @return tibble
#' @rdname ms
#' @export
sample_stats = function(msout, lowerbound=0) {
  vaf = colMeans(msout)
  tibble::tibble(math = math_score(vaf[vaf > lowerbound]))
}
