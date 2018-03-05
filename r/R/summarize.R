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
