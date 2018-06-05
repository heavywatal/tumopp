#' Sample cells from a population
#'
#' @description
#' `sample_bulk` samples bulk of cells near the specified coord
#' @param extant tibble of extant cells
#' @param center named (x, y, z) vector, list, or tibble
#' @param ncell number of cells per specimen
#' @return tibble
#' @rdname sample
#' @export
sample_bulk = function(extant, center, ncell=100L) {
  extant %>%
    dplyr::mutate(r_sample = dist_euclidean(., center)) %>%
    dplyr::top_n(-ncell, .data$r_sample) %>%
    dplyr::arrange(.data$r_sample) %>%
    utils::head(ncell)
}

#' `sample_uniform_regions` sample uniformly distributed regions
#' @param nsam number of regions to sample
#' @rdname sample
#' @export
sample_uniform_regions = function(extant, nsam=2L, ncell=10L) {
  extant = dplyr::select(extant, .data$id, .data$x, .data$y, .data$z)
  clusters = extant %>%
    dplyr::select(-.data$id) %>%
    stats::kmeans(centers = nsam, iter.max = 32L)
  extant %>%
    dplyr::mutate(cluster = clusters$cluster) %>%
    tidyr::nest(-.data$cluster) %>%
    dplyr::arrange(.data$cluster) %>%
    dplyr::bind_cols(tibble::as_tibble(clusters$centers)) %>%
    purrr::pmap_dfr(function(x, y, z, data, ...) {
      sampled = sample_bulk(data, center = c(x = x, y = y, z = z), ncell = ncell)
      utils::head(sampled, 1L) %>%
        dplyr::select(.data$id, .data$x, .data$y, .data$z) %>%
        dplyr::mutate(samples = list(sampled$id))
    })
}

#' `sample_random_regions` samples multiple regions at random
#' @rdname sample
#' @export
sample_random_regions = function(extant, nsam=2L, ncell=10L) {
  stopifnot(nsam > 1L)
  extant %>%
    dplyr::select(.data$id, .data$x, .data$y, .data$z) %>%
    dplyr::sample_n(nsam) %>%
    dplyr::mutate(samples = purrr::pmap(., function(x, y, z, ...) {
      sample_bulk(extant, center = c(x = x, y = y, z = z), ncell = ncell)$id
    }))
}
