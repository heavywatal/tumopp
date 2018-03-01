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

#' Sample specimens and measure mean branch lengths
#' @param population tibble
#' @param ncell number of cells per specimen
#' @param nsam number of samples
#' @return tibble
#' @rdname summarize
#' @export
within_between_samples = function(population, ncell=100L, nsam=2L) {
  stopifnot(nsam > 1L)
  .extant = filter_extant(population)
  .centers = .extant %>%
    dplyr::select(.data$x, .data$y, .data$z) %>%
    dplyr::sample_n(nsam) %>%
    dplyr::mutate(
      sample_id = seq_along(.data$x),
      samples = purrr::pmap(., function(x, y, z, ...) {
        sample_bulk(.extant, center = c(x = x, y = y, z = z), size = ncell)$id
      }),
      within = purrr::map_dbl(.data$samples, ~mean_branch_length(population, .x))
    )
  tidyr::crossing(lhs = .centers$sample_id, rhs = .centers$sample_id) %>%
    dplyr::filter(.data$lhs < .data$rhs) %>%
    dplyr::mutate(res = purrr::pmap(., ~{
      row_x = .centers[.x, ]
      row_y = .centers[.y, ]
      tibble::tibble(
        euclidean = dist_euclidean(row_x, row_y),
        within = mean(row_x$within, row_y$within),
        between = mean_branch_length(population, row_x$samples[[1]], row_y$samples[[1]])
      )
    })) %>%
    tidyr::unnest() %>%
    dplyr::select(-.data$lhs, -.data$rhs)
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
