#' Sample bulk of cells near the specified coord
#' @param extant tibble of extant cells
#' @param center named (x, y, z) vector, list, or tibble
#' @param size number of cells
#' @return tibble
#' @rdname sample
#' @export
sample_bulk = function(extant, center, size=100L) {
  extant %>%
    dplyr::mutate(r_sample = dist_euclidean(., center)) %>%
    dplyr::top_n(-size, .data$r_sample) %>%
    dplyr::arrange(.data$r_sample) %>%
    utils::head(size)
}

#' Sample multiple regions from extant cells
#' @param nsam number of regions to sample
#' @param ncell number of cells per specimen
#' @return nested tibble
#' @rdname sample
#' @export
sample_random_regions = function(extant, nsam=2L, ncell=10L) {
  stopifnot(nsam > 1L)
  extant %>%
    dplyr::select(.data$id, .data$x, .data$y, .data$z) %>%
    dplyr::sample_n(nsam) %>%
    dplyr::mutate(samples = purrr::pmap(., function(x, y, z, ...) {
      sample_bulk(extant, center = c(x = x, y = y, z = z), size = ncell)$id
    }))
}

# Make union of integers from list
union_int = function(x) {
  unique(purrr::flatten_int(x))
}

# Shortcut with different defaults
combn_int_list = function(x, m, FUN=union_int, simplify=FALSE, ...) {
  utils::combn(x, m, FUN = FUN, simplify = simplify, ...)
}

#' Make combinations of cell IDs for various number of samples
#' @param ids list of ID vectors
#' @return nested tibble
#' @rdname sample
#' @export
combn_sample_ids = function(ids, nsam=seq_along(ids)) {
  tibble::tibble(
    nsam = nsam,
    nodes = purrr::map(nsam, ~combn_int_list(ids, .x))
  ) %>%
    tidyr::unnest()
}

#' Calculate expected allele capture rate on various combinations of samples
#' @param combinations nested tibble from combn_sample_ids()
#' @inheritParams detectable_mutants_all
#' @return tibble
#' @rdname sample
#' @export
summarize_capture_rate = function(combinations, population, threshold = 0.01) {
  ids = detectable_mutants_all(population, threshold)
  len_ids = length(ids)
  dplyr::transmute(
    combinations,
    threshold,
    .data$nsam,
    capture_rate = purrr::map_dbl(.data$nodes, ~sum(.x %in% ids)) / len_ids
  )
}

#' Sample specimens and measure mean branch lengths
#' @param graph igraph
#' @param regions output of sample_random_regions()
#' @return tibble
#' @rdname sample
#' @export
within_between_samples = function(graph, regions) {
  regions %>%
    dplyr::mutate(
      within = purrr::map_dbl(.data$samples, ~mean_branch_length(graph, as.character(.x)))
    ) %>%
    purrr::transpose() %>%
    purrr::cross2(.x = ., .y = ., .filter = ~.x$id >= .y$id) %>%
    purrr::map_dfr(~{
      row_i = .x[[1L]]
      row_j = .x[[2L]]
      tibble::tibble(
        region_i = row_i$id,
        region_j = row_j$id,
        within_i = row_i$within,
        within_j = row_j$within,
        euclidean = dist_euclidean(row_i, row_j),
        between = mean_branch_length(.graph, as.character(row_i$samples), as.character(row_j$samples))
      )
    }) %>%
    dplyr::mutate(
      within = 0.5 * (.data$within_i + .data$within_j),
      fst = wtl::fst_HBK(.data$within, .data$between)
    )
}
