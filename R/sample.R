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
    nodes = purrr::flatten_int(samples) %>>% unique() %>>% setdiff(common_ancs)
    if (is.null(segsites)) {
        segsites = stats::rpois(1L, length(nodes) * mu)
    } else if (!is.null(mu)) warning('mu is ignored if segsites is given')
    mutant_ids = sample(nodes, segsites, replace=TRUE)
    samples %>>%
        purrr::map(~ mutant_ids %in% .) %>>%
        purrr::flatten_int() %>>%
        matrix(nrow=nsam, byrow=TRUE)
}

#' Summary statistics of neutral mutations
#' @param msout tibble from `make_samples()`
#' @return tibble
#' @rdname summarize
#' @export
sample_stats = function(msout) {
    vaf = colMeans(msout)
    tibble::tibble(math= math_score(vaf))
}
