#' Make parent-child list from raw population tibble
#' @param population tibble
#' @return tibble
#' @rdname graph
#' @export
make_edgelist = function(population) {
    dplyr::filter_(population, ~age > 0L) %>>%
    dplyr::transmute_(
      from= ~purrr::map2_int(genealogy, age, `[[`),
      to= ~id
    ) %>>%
    dplyr::mutate_all(as.integer) %>>%
    dplyr::arrange_(~to)
}

#' Mean branch length within/between sub-graphs
#' @param graph igraph
#' @param from,to node names
#' @return numeric
#' @rdname graph
#' @export
mean_branch_length = function(graph, from=igraph::V(graph), to=from) {
    .d = igraph::distances(graph, from, to, mode='all', weights=NA)
    sum(.d) / sum(.d > 0)
}

#' Calculate Fst by Hudson, Slatkin, and Maddison (1992)
#' @param within,between mean branch length or diversity
#' @return numeric
#' @rdname graph
#' @export
fst_HSM = function(within, between) {
    1.0 - within / between
}

#' Calculate Kst by Hudson, Boos, and Kaplan (1992)
#' @param n number of subpopulations
#' @return numeric
#' @rdname graph
#' @export
fst_HBK = function(within, between, n=2) {
    (between - within) / (between + within / (n - 1))
}

#' Set coordinates of nodes and edges for plotting
#' @param edgelist tibble from make_edgelist()
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(edgelist) {
    g = igraph::graph_from_data_frame(edgelist)
    lo = igraph::layout_as_tree(g, flip.y=FALSE)
    nodes = tibble::as_tibble(lo) %>>%
        stats::setNames(c('pos', 'age')) %>>%
        dplyr::mutate_(
          id= ~igraph::V(g)$name %>>% as.integer(),
          extant= ~(igraph::degree(g, mode='out') == 0L))
    edgelist %>>%
    dplyr::left_join(dplyr::select_(nodes, ~-extant), by=c(from='id')) %>>%
    dplyr::left_join(dplyr::select_(nodes, posend= ~pos, ageend= ~age, ~extant, to= ~id), by='to')
}

#' Plot genealogy
#' @param .data tibble from layout_genealogy()
#' @param label character or expression
#' @param xmax numeric
#' @param hist logical: whether add age_histogram() or not.
#' @return gg
#' @rdname graph
#' @export
ggplot_genealogy = function(.data, label='', xmax=20, hist=TRUE) {
    age_lim=c(0, max(.data$ageend, xmax))
    .tree = ggplot2::ggplot(.data)+
    ggplot2::geom_segment(ggplot2::aes_(~age, ~pos, xend=~ageend, yend=~posend), alpha=0.3, size=0.3)+
    ggplot2::geom_point(data=dplyr::filter_(.data, ~extant),
        ggplot2::aes_(x=~ageend, y=~posend),
        size=0.8, colour='dodgerblue', alpha=0.2)+
    ggplot2::coord_cartesian(xlim=age_lim)+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title=ggplot2::element_blank(),
        axis.text=ggplot2::element_blank(),
        axis.ticks=ggplot2::element_blank(),
        panel.grid.major.y=ggplot2::element_blank())
    if (hist) {
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
