#' ggplot for 2D lattice
#' @param .data a data.frame
#' @param colour a column name for colour
#' @param limit for value range
#' @rdname plot
#' @export
gglattice2D = function(.data, colour='ancestors', limit=maxabs(.data)) {
    ggplot2::ggplot(.data, aes(x, y))+
    ggplot2::geom_point(aes_string(colour=colour), alpha=0.66, size=80/limit)+
    ggplot2::scale_colour_hue(na.value='white', drop=FALSE)+
    ggplot2::expand_limits(x=limit * c(-1, 1), y=limit * c(-1, 1))+
    ggplot2::theme_grey()+
    ggplot2::theme(panel.background=element_rect(fill='grey90'))+
    ggplot2::theme(panel.grid=element_blank())+
    ggplot2::theme(axis.title=element_blank())
}
