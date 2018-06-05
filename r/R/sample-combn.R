#' Sample cells from a population
#'
#' @description
#' `combn_sample_ids` makes combinations of cell IDs for various number of samples
#' @param ids list of ID vectors
#' @param nsam number of regions to sample
#' @return tibble
#' @rdname sample-combn
#' @export
combn_sample_ids = function(ids, nsam=seq_along(ids)) {
  tibble::tibble(
    nsam = nsam,
    nodes = purrr::map(nsam, ~combn_int_list(ids, .x))
  ) %>%
    tidyr::unnest()
}

# Make union of integers from list
union_int = function(x) {
  unique(purrr::flatten_int(x))
}

# Shortcut with different defaults
combn_int_list = function(x, m, FUN=union_int, simplify=FALSE, ...) {
  utils::combn(x, m, FUN = FUN, simplify = simplify, ...)
}

#' `summarize_capture_rate` calculates expected allele capture rate on various combinations of samples
#' @param combinations nested tibble from `combn_sample_ids()`
#' @inheritParams detectable_mutants_all
#' @rdname sample-combn
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

#' `within_between_samples` summarizes branch lengths
#' @param graph igraph
#' @param regions output of `sample_uniform_regions()`
#' @rdname sample-combn
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
