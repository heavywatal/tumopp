#' Make ms-like matrix from survivors tibble
#' @param survivors tibble
#' @param nsam number of samples
#' @param mu mutation rate per cell division
#' @return logical tibble
#' @rdname sample
#' @export
make_samples = function(survivors, nsam=nrow(survivors), mu=1e-1) {
    max_id = max(survivors$id)
    num_mut_events = stats::rpois(1L, (max_id - 1L) * mu)
    mutant_ids = sample.int(max_id, num_mut_events)
    survivors %>>%
        dplyr::sample_n(nsam) %>>%
        purrr::by_row(~{mutant_ids %in% .$genealogy[[1]]},
            .to='site', .labels=FALSE, .collate='cols') %>>%
        stats::setNames(paste0('s', mutant_ids)) %>>%
        dplyr::select_if(any)
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
