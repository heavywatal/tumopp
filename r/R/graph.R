#' Make igraph from raw population tbl
#' @param population tbl
#' @return tibble
#' @rdname graph
#' @export
make_igraph = function(population) {
  population %>%
    dplyr::transmute(
      from = .data$ancestor,
      to = .data$id
    ) %>%
    dplyr::filter(.data$from > 0L) %>%
    dplyr::arrange(.data$to) %>%
    dplyr::mutate_all(as.character) %>%
    igraph::graph_from_data_frame()
}

distances_from_origin = function(graph, nodes=character(0L)) {
  igraph::distances(graph, "1", nodes, mode = "out", weights = NA, algorithm = "unweighted") %>%
    as.integer()
}

paths_to_origin = function(graph, nodes=character(0L)) {
  igraph::shortest_paths(graph, from = "1", to = nodes, mode = "out", weights = NA, output = "vpath")$vpath %>%
    purrr::map(~as.integer(.x$name))
}

#' Mean branch length within/between sub-graphs
#' @param graph igraph
#' @param from,to igraph vertices
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length = function(graph, from=igraph::V(graph), to=from) {
  .d = igraph::distances(graph, from, to, mode = "all", weights = NA, algorithm = "unweighted")
  .n = length(from) * length(to) - sum(from %in% to)
  sum(.d) / .n
}

#' Set coordinates of nodes and edges for plotting
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(population) {
  extra_cols = population %>%
    dplyr::transmute(name = as.character(.data$id), extant = .data$death == 0, .data$clade)
  make_igraph(population) %>%
    wtl::igraph_layout(igraph::as_tree(flip.y = FALSE)) %>%
    dplyr::rename(pos = "x", age = "y", posend = "xend", ageend = "yend") %>%
    dplyr::left_join(extra_cols, by = c(to = "name"))
}
