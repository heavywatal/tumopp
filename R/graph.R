#' Make parent-child list from raw population tbl or igraph
#' @param obj population tbl or igraph
#' @return tibble
#' @rdname graph
#' @export
make_edgelist = function(obj) {
    if (igraph::is_igraph(obj)) {
        igraph::as_data_frame(obj, 'edges') %>%
        tibble::as_tibble()
    } else {
        dplyr::filter_(obj, ~age > 0L) %>%
        dplyr::transmute_(
          from= ~purrr::map2_int(genealogy, age, `[[`),
          to= ~id
        ) %>%
        dplyr::arrange_(~to) %>%
        dplyr::mutate_all(as.character)
    }
}

#' Make igraph from raw population tbl
#' @param population tbl
#' @param ids integer vector
#' @return tibble
#' @rdname graph
#' @export
make_igraph = function(population, ids=integer(0L)) {
    if (length(ids) > 0L) {
        population = filter_connected(population, ids)
    }
    igraph::graph_from_data_frame(make_edgelist(population))
}

#' Mean branch length within/between sub-graphs
#' @param from,to node id
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length = function(population, from=integer(0L), to=from) {
    graph = make_igraph(population, union(from, to))
    mean_branch_length.igraph(graph, as.character(from), as.character(to))
}

#' Mean branch length within/between sub-graphs
#' @param graph igraph
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length.igraph = function(graph, from=igraph::V(graph), to=from) {
    .d = igraph::distances(graph, from, to, mode='all', weights=NA)
    .n = length(from) * length(to) - length(intersect(from, to))
    sum(.d) / .n
}

# Slow; do not use
mean_branch_length.R = function(genealogy) {
    .len = length(genealogy)
    .sum_lengths = sum(lengths(genealogy, FALSE)) * (.len - 1)
    .sum_intersects = tidyr::crossing(from=seq_len(.len), to=seq_len(.len)) %>%
        dplyr::filter_(~from < to) %>%
        purrr::pmap_int(function(from, to) {
            length(intersect(genealogy[[from]], genealogy[[to]]))
        }) %>% sum()
    (.sum_lengths - 2 * .sum_intersects) / choose(.len, 2)
}

#' Set coordinates of nodes and edges for plotting
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(population) {
    graph = make_igraph(population)
    nodes = igraph::layout_as_tree(graph, flip.y=FALSE) %>%
        tibble::as_tibble() %>%
        stats::setNames(c('pos', 'age')) %>%
        dplyr::mutate_(id= ~igraph::V(graph)$name)
    .extant = population %>%
        dplyr::transmute_(id= ~as.character(id), extant= ~death == 0)
    .children = dplyr::select_(nodes, posend= ~pos, ageend= ~age, ~id) %>%
        dplyr::left_join(.extant, by='id')
    make_edgelist(graph) %>%
        dplyr::left_join(nodes, by=c(from='id')) %>%
        dplyr::left_join(.children, by=c(to='id'))
}
