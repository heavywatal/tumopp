#' Make parent-child list from raw population tbl or igraph
#' @param obj population tbl or igraph
#' @return tibble
#' @rdname graph
#' @export
make_edgelist = function(obj) {
  if (igraph::is_igraph(obj)) {
    igraph::as_data_frame(obj, "edges") %>%
      tibble::as_tibble()
  } else {
    dplyr::filter(obj, .data$age > 0L) %>%
      dplyr::transmute(
        from = purrr::map2_int(.data$genealogy, .data$age, `[[`),
        to = .data$id
      ) %>%
      dplyr::arrange(.data$to) %>%
      dplyr::mutate_all(as.character)
  }
}

#' Make igraph from raw population tbl
#' @param population tbl
#' @return tibble
#' @rdname graph
#' @export
make_igraph = function(population) {
  igraph::graph_from_data_frame(make_edgelist(population))
}

#' Mean branch length within/between sub-graphs
#' @param from,to node id
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length = function(population, from=integer(0L), to=from) {
  population %>%
    filter_connected(union(from, to)) %>%
    make_igraph() %>%
    mean_branch_length.igraph(as.character(from), as.character(to))
}

#' Mean branch length within/between sub-graphs
#' @param graph igraph
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length.igraph = function(graph, from=igraph::V(graph), to=from) {
  .d = igraph::distances(graph, from, to, mode = "all", weights = NA)
  .n = length(from) * length(to) - length(intersect(from, to))
  sum(.d) / .n
}

# Slow; do not use
mean_branch_length.R = function(genealogy) {
  .len = length(genealogy)
  .sum_lengths = sum(lengths(genealogy, FALSE)) * (.len - 1L)
  .sum_intersects = tidyr::crossing(from = seq_len(.len), to = seq_len(.len)) %>%
    dplyr::filter(.data$from < .data$to) %>%
    purrr::pmap_int(function(from, to) {
      length(intersect(genealogy[[from]], genealogy[[to]]))
    }) %>%
    sum()
  (.sum_lengths - 2 * .sum_intersects) / choose(.len, 2L)
}

#' Set coordinates of nodes and edges for plotting
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(population) {
  graph = make_igraph(population)
  nodes = igraph::layout_as_tree(graph, flip.y = FALSE) %>%
    tibble::as_tibble() %>%
    stats::setNames(c("pos", "age")) %>%
    dplyr::mutate(id = igraph::V(graph)$name)
  .extant = population %>%
    dplyr::transmute(id = as.character(.data$id), extant = .data$death == 0, .data$clade)
  .children = dplyr::select(nodes, posend = .data$pos, ageend = .data$age, .data$id) %>%
    dplyr::left_join(.extant, by = "id")
  make_edgelist(graph) %>%
    dplyr::left_join(nodes, by = c(from = "id")) %>%
    dplyr::left_join(.children, by = c(to = "id"))
}
