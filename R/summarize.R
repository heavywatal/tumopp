#' Summary statistics of genealogy
#' @param extant tibble
#' @return tibble
#' @rdname summarize
#' @export
genetic_stats = function(extant) {
  dplyr::summarize(
    extant,
    mean_age = mean(.data$age), sd_age = stats::sd(.data$age),
    max_age = max(.data$age), min_age = min(.data$age),
    evenness = wtl::evenness(.data$clade)
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
    purrr::map_dfr(seq_len(npair), ~within_between_samples(population, ncell))
  } else {
    .extant = filter_extant(population)
    .o1 = dplyr::sample_n(.extant, 1L)
    .o2 = dplyr::sample_n(.extant, 1L)
    .nodes1 = sample_bulk(.extant, .o1, ncell)$id
    .nodes2 = sample_bulk(.extant, .o2, ncell)$id
    .within1 = mean_branch_length(population, .nodes1)
    .within2 = mean_branch_length(population, .nodes2)
    tibble::tibble(
      euclidean = dist_euclidean(.o1, .o2),
      within = mean(c(.within1, .within2)),
      between = mean_branch_length(population, .nodes1, .nodes2)
    )
  }
}

#' Summary statistics of morphology
#' @return tibble
#' @rdname summarize
#' @export
morphological_stats = function(extant) {
  extant %>%
    dplyr::filter(.data$surface)
  dplyr::summarise(
    phi_mean = mean(.data$phi), phi_sd = stats::sd(.data$phi),
    r_mean = mean(.data$r), r_sd = stats::sd(.data$r)
  )
  dplyr::mutate(surface = sum(extant$surface) / nrow(extant))
}
