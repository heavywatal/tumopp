#' Sample bulk of cells near the specified coord
#' @param population tibble of extant cells
#' @param center named (x, y, z) vector, list, or tibble
#' @param size number of cells
#' @return tibble
#' @rdname sample
#' @export
sample_bulk = function(population, center, size=100L) {
  if (is.list(center)) {
    return(purrr::pmap_dfr(center, function(x, y, z, ...) {
      sample_bulk(population, center = c(x = x, y = y, z = z), size = size)
    }))
  }
  population %>%
    dplyr::mutate(dist = dist_euclidean(., center)) %>%
    dplyr::arrange(.data$dist) %>%
    utils::head(size)
}

#' Make ms-like neutral variation matrix from genealogy column of extant cells
#' @param genealogy list of int vectors
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @return int matrix [nsam x segsites]
#' @rdname sample
#' @export
sprinkle_mutations = function(genealogy, mu=NULL, segsites=NULL) {
  common_ancs = purrr::reduce(genealogy, intersect)
  nodes = purrr::flatten_int(genealogy) %>% unique() %>% setdiff(common_ancs)
  if (is.null(segsites)) {
    if (is.null(mu)) stop("specify either mu or segsites")
    segsites = stats::rpois(1L, length(nodes) * mu)
  } else if (!is.null(mu)) warning("mu is ignored if segsites is given")
  mutant_ids = sample(nodes, segsites, replace = TRUE) %>% sort()
  genealogy %>%
    purrr::map(~ mutant_ids %in% .) %>%
    purrr::flatten_int() %>%
    matrix(nrow = length(genealogy), byrow = TRUE)
}

#' Summary statistics of neutral mutations
#' @param msout tibble from `sprinkle_mutations()`
#' @param lowerbound to discard peripheral branches
#' @return tibble
#' @rdname summarize
#' @export
sample_stats = function(msout, lowerbound=0) {
  vaf = colMeans(msout)
  tibble::tibble(math = math_score(vaf[vaf > lowerbound]))
}
