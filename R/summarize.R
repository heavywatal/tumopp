#' Summary statistics of genealogy
#' @param extant tibble
#' @return tibble
#' @rdname summarize
#' @export
genetic_stats = function(extant) {
    dplyr::summarize_(extant,
       mean_age= ~mean(age), sd_age= ~sd(age),
       max_age= ~max(age), min_age= ~min(age),
       evenness= ~wtl::evenness(clade)
    )
}

#' Sample two specimens and measure mean branch lengths
#' @param population tibble
#' @param ncell number of cells per specimen
#' @param npair number of pairs and output rows
#' @return tibble
#' @rdname summarize
#' @export
within_between_samples = function(population, ncell=100L, npair=1L) {
    if (npair > 1L) {
        purrr::map_df(seq_len(npair), ~within_between_samples(population, ncell))
    } else {
        .extant = filter_extant(population)
        .o1 = dplyr::sample_n(.extant, 1)
        .o2 = dplyr::sample_n(.extant, 1)
        .nodes1 = sample_bulk(.extant, .o1, ncell)$id
        .nodes2 = sample_bulk(.extant, .o2, ncell)$id
        .within1 = mean_branch_length(population, .nodes1)
        .within2 = mean_branch_length(population, .nodes2)
        tibble::tibble(
            euclidean= dist_euclidean(.o1, .o2),
            within= mean(c(.within1, .within2)),
            between= mean_branch_length(population, .nodes1, .nodes2))
    }
}

#' Summary statistics of morphology
#' @return tibble
#' @rdname summarize
#' @export
morphological_stats = function(extant) {
    extant %>>%
        dplyr::filter_(~surface) %>>%
        dplyr::summarise_(phi_mean= ~mean(phi), phi_sd= ~stats::sd(phi),
                          r_mean= ~mean(r), r_sd= ~stats::sd(r)) %>>%
        dplyr::mutate(surface=sum(extant$surface) / nrow(extant))
}
