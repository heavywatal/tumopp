#' Sample bulk of cells near the specified coord
#' @param population tibble of extant cells
#' @param center named vector or single-row tbl
#' @param size number of cells
#' @return tibble
#' @rdname sample
#' @export
sample_bulk = function(population, center, size=100L) {
    population %>%
    dplyr::mutate(dist= dist_euclidean(., center)) %>%
    dplyr::arrange(.data$dist) %>%
    head(size)
}

#' Make ms-like matrix from genealogy column of population
#' @param genealogy list of int vectors
#' @param nsam number of samples
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @return int matrix [nsam x segsites]
#' @rdname sample
#' @export
make_samples = function(genealogy, nsam=length(genealogy), mu=NULL, segsites=NULL) {
    samples = sample(genealogy, nsam)
    common_ancs = purrr::reduce(samples, intersect)
    nodes = purrr::flatten_int(samples) %>% unique() %>% setdiff(common_ancs)
    if (is.null(segsites)) {
        if (is.null(mu)) stop('specify either mu or segsites')
        segsites = stats::rpois(1L, length(nodes) * mu)
    } else if (!is.null(mu)) warning('mu is ignored if segsites is given')
    mutant_ids = sample(nodes, segsites, replace=TRUE)
    samples %>%
        purrr::map(~ mutant_ids %in% .) %>%
        purrr::flatten_int() %>%
        matrix(nrow=nsam, byrow=TRUE)
}

#' Summary statistics of neutral mutations
#' @param msout tibble from `make_samples()`
#' @param lowerbound to discard peripheral branches
#' @return tibble
#' @rdname summarize
#' @export
sample_stats = function(msout, lowerbound=0) {
    vaf = colMeans(msout)
    tibble::tibble(math= math_score(vaf[vaf > lowerbound]))
}
