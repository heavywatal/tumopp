#' ggplot for 2D lattice
#' @param .data a data.frame
#' @param colour a column name for colour
#' @param limit for value range
#' @rdname plot
#' @export
gglattice2D = function(.data, colour='ancestors', limit=maxabs(.data)) {
    ggplot2::ggplot(.data, ggplot2::aes(x, y))+
    ggplot2::geom_point(ggplot2::aes_string(colour=colour), alpha=0.66, size=80/limit)+
    ggplot2::scale_colour_hue(na.value='white', drop=FALSE)+
    ggplot2::expand_limits(x=limit * c(-1, 1), y=limit * c(-1, 1))+
    ggplot2::theme_grey()+
    ggplot2::theme(panel.background=ggplot2::element_rect(fill='grey90'))+
    ggplot2::theme(panel.grid=ggplot2::element_blank())+
    ggplot2::theme(axis.title=ggplot2::element_blank())
}

#' Write animation GIF for serial sections of 3D tumor
#' @param filename a string
#' @param .data a data.frame
#' @param width an integer
#' @param height an integer
#' @rdname plot
#' @export
save_serial_section = function(filename='serial_section.gif', .data, width=720, height=640) {
    .lim = maxabs(.data)
    section_plots = dplyr::group_by(.data, z) %>>%
        dplyr::do(plt={
            gglattice2D(., limit=.lim)+
            ggplot2::geom_hline(yintercept=.$z[1])
        })
    animation::saveGIF({for (p in section_plots$plt) {print(p)}},
        filename, outdir=getwd(), interval=0.15,
        ani.width=width, ani.height=height, loop=TRUE, autobrowse=FALSE)
}
