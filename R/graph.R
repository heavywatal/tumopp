#' Make parent-child list from raw population tbl or igraph
#' @param obj population tbl or igraph
#' @return tibble
#' @rdname graph
#' @export
make_edgelist = function(obj) {
    if (igraph::is_igraph(obj)) {
        igraph::as_data_frame(obj, 'edges') %>>%
        tibble::as_tibble()
    } else {
        dplyr::filter_(obj, ~age > 0L) %>>%
        dplyr::transmute_(
          from= ~purrr::map2_int(genealogy, age, `[[`),
          to= ~id
        ) %>>%
        dplyr::arrange_(~to) %>>%
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
    .sum_intersects = tidyr::crossing(from=seq_len(.len), to=seq_len(.len)) %>>%
        dplyr::filter_(~from < to) %>>%
        purrr::pmap_int(function(from, to) {
            length(intersect(genealogy[[from]], genealogy[[to]]))
        }) %>>% sum()
    (.sum_lengths - 2 * .sum_intersects) / choose(.len, 2)
}

#' Set coordinates of nodes and edges for plotting
#' @param graph igraph
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(graph) {
    lo = igraph::layout_as_tree(graph, flip.y=FALSE)
    nodes = tibble::as_tibble(lo) %>>%
        stats::setNames(c('pos', 'age')) %>>%
        dplyr::mutate_(
          id= ~igraph::V(graph)$name,
          extant= ~(igraph::degree(graph, mode='out') == 0L))
    make_edgelist(graph) %>>%
    dplyr::left_join(dplyr::select_(nodes, ~-extant), by=c(from='id')) %>>%
    dplyr::left_join(dplyr::select_(nodes, posend= ~pos, ageend= ~age, ~extant, to= ~id), by='to')
}

#' Plot genealogy
#' @param .data tbl from layout_genealogy()
#' @param label character or expression
#' @param xmax numeric
#' @param hist logical: whether add age_histogram() or not.
#' @return gg
#' @rdname graph
#' @export
ggplot_genealogy = function(.data, label='', xmax=20, hist=FALSE) {
    age_lim=c(0, max(.data$ageend, xmax))
    .tree = ggplot2::ggplot(.data)+
    ggplot2::geom_segment(ggplot2::aes_(~age, ~pos, xend=~ageend, yend=~posend), alpha=0.3, size=0.3)+
    ggplot2::geom_point(data=dplyr::filter_(.data, ~extant),
        ggplot2::aes_(x=~ageend, y=~posend),
        size=0.8, colour='dodgerblue', alpha=0.2)+
    ggplot2::coord_cartesian(xlim=age_lim)+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title.y=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        panel.grid.major.y=ggplot2::element_blank())
    if (hist) {
        .tree = .tree + ggplot2::theme(
            axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank())
        .hist = age_histogram(.data, age_lim)
        .top = grid::textGrob(label, x=grid::unit(0.1, 'npc'), just=c('left', 'top'))
        gridExtra::arrangeGrob(.tree, .hist, nrow=2, heights=c(3, 1), top=.top)
    } else {.tree}
}

#' Plot age histogram
#' @param age_lim numeric vector
#' @return gg
#' @rdname graph
age_histogram = function(.data, age_lim) {
    dplyr::filter_(.data, ~extant) %>>%
    ggplot2::ggplot(ggplot2::aes_(~ageend))+
    ggplot2::geom_histogram(binwidth=1, center=0)+
    ggplot2::coord_cartesian(xlim=age_lim)+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        panel.grid.major.y=ggplot2::element_blank())
}
