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

#' Set coordinates of nodes and edges for plotting
#' @param edgelist tibble from make_edgelist()
#' @return tibble
#' @rdname graph
#' @export
layout_genealogy = function(edgelist) {
    g = igraph::graph_from_data_frame(edgelist)
    lo = igraph::layout.reingold.tilford(g, circular=FALSE)
    nodes = tibble::as_tibble(lo) %>>%
        stats::setNames(c('x', 'y')) %>>%
        dplyr::mutate_(
          id= ~igraph::V(g)$name %>>% as.integer(),
          extant= ~(igraph::degree(g, mode='out') == 0L),
          y= ~max(y) - y)
    edgelist %>>%
    dplyr::left_join(dplyr::select_(nodes, ~-extant), by=c(from='id')) %>>%
    dplyr::left_join(dplyr::select_(nodes, xend= ~x, yend= ~y, ~extant, to= ~id), by='to')
}

#' Plot genealogy
#' @param .data tibble from layout_genealogy()
#' @param age_lim numeric vector
#' @param hist logical: whether add age_histogram() or not.
#' @return gg
#' @rdname graph
#' @export
ggplot_genealogy = function(.data, age_lim=c(0, max(.data$yend)), hist=TRUE) {
    .tree = ggplot2::ggplot(.data)+
    ggplot2::geom_segment(ggplot2::aes_(~x, ~y, xend=~xend, yend=~yend), alpha=0.3, size=0.3)+
    ggplot2::geom_point(data=dplyr::filter_(.data, ~extant),
        ggplot2::aes_(x=~xend, y=~yend),
        size=1, colour='dodgerblue', alpha=0.2)+
    ggplot2::scale_y_continuous(trans='reverse')+
    ggplot2::coord_cartesian(ylim=age_lim)+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title=ggplot2::element_blank(),
        axis.text=ggplot2::element_blank(),
        axis.ticks=ggplot2::element_blank(),
        panel.grid.major.x=ggplot2::element_blank())
    if (hist) {
        .hist = age_histogram(.data, age_lim)
        gridExtra::arrangeGrob(.tree, .hist, nrow=1, widths=c(3, 1))
    } else {.tree}
}

#' Plot age histogram
#' @return gg
#' @rdname graph
age_histogram = function(.data, age_lim) {
    dplyr::filter_(.data, ~extant) %>>%
    ggplot2::ggplot(ggplot2::aes_(~yend))+
    ggplot2::geom_histogram(binwidth=1, center=0)+
    ggplot2::scale_x_continuous(trans='reverse')+
    ggplot2::coord_flip(xlim=age_lim)+
    ggplot2::labs(x='Number of cell divisions')+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank())
}
