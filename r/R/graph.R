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

#' Mean branch length within/between sub-graphs
#' @param from,to node id
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length = function(population, from=integer(0L), to=from) {
  population %>%
    make_igraph() %>%
    mean_branch_length.igraph(as.character(from), as.character(to))
}

#' Mean branch length within/between sub-graphs
#' @param graph igraph
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length.igraph = function(graph, from=igraph::V(graph), to=from) {
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
