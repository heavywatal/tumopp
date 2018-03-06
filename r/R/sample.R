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

#' Sample specimens and measure mean branch lengths
#' @param regions output of sample_random_regions()
#' @param graph igraph
#' @return tibble
#' @rdname sample
#' @export
within_between_samples = function(regions, graph) {
  regions %>%
    dplyr::mutate(
      within = purrr::map_dbl(.data$samples, ~mean_branch_length.igraph(graph, as.character(.x)))
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
        between = mean_branch_length.igraph(.graph, as.character(row_i$samples[[1L]]), as.character(row_j$samples[[1L]]))
      )
    }) %>%
    dplyr::mutate(within = 0.5 * (.data$within_i + .data$within_j))
}

#' Traceback ancestors and calculate union of the IDs
#' @param nodes igraph vertices
#' @return union of IDs
#' @rdname sample
#' @export
make_ancestor_union = function(graph, nodes) {
  igraph::all_simple_paths(graph, from = "1", to = nodes, mode = "out") %>%
    purrr::map(~as.integer(.x$name)) %>%
    purrr::flatten_int() %>%
    unique()
}

#' Make ms-like neutral variation matrix
#' @param population tibble of raw population
#' @param nsam number of cells to sample
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @return int matrix [nsam x segsites]
#' @rdname sample
#' @export
sprinkle_mutations = function(population, nsam, mu=NULL, segsites=NULL) {
  graph = make_igraph(population)
  samples = sample(filter_extant(population)$id, nsam) %>% as.character()
  lineages = igraph::neighborhood(graph, order = 1073741824L, nodes = samples, mode = "in") %>%
    purrr::map(~ .x$name)
  common_ancs = purrr::reduce(lineages, intersect)
  nodes = purrr::flatten_chr(lineages) %>% unique() %>% setdiff(common_ancs)
  if (is.null(segsites)) {
    if (is.null(mu)) stop("specify either mu or segsites")
    segsites = stats::rpois(1L, length(nodes) * mu)
  } else if (!is.null(mu)) warning("mu is ignored if segsites is given")
  mutant_ids = sample(nodes, segsites, replace = TRUE) %>% sort()
  lineages %>%
    purrr::map(~ mutant_ids %in% .) %>%
    purrr::flatten_int() %>%
    matrix(nrow = nsam, byrow = TRUE)
}

#' Summary statistics of neutral mutations
#' @param msout tibble from `sprinkle_mutations()`
#' @param lowerbound to discard peripheral branches
#' @return tibble
#' @rdname summarize
#' @export
sample_stats = function(msout, lowerbound=0) {
  vaf = colMeans(msout)
  tibble::tibble(math = math_score(vaf[vaf > lowerbound]))
}
