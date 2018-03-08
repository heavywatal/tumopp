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
  graph = make_igraph(population)
  nodes = igraph::layout_as_tree(graph, flip.y = FALSE) %>%
    tibble::as_tibble() %>%
    stats::setNames(c("pos", "age")) %>%
    dplyr::mutate(id = igraph::V(graph)$name)
  .extant = population %>%
    dplyr::transmute(id = as.character(.data$id), extant = .data$death == 0, .data$clade)
  .children = dplyr::select(nodes, posend = .data$pos, ageend = .data$age, .data$id) %>%
    dplyr::left_join(.extant, by = "id")
  igraph::as_data_frame(graph, "edges") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(nodes, by = c(from = "id")) %>%
    dplyr::left_join(.children, by = c(to = "id"))
}
