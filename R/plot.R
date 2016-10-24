#' ggplot for frequency spectrum
#' @param freqs numeric vector
#' @rdname plot
#' @export
ggfreqspec = function(freqs) {
    tibble::tibble(x=freqs) %>>%
    ggplot2::ggplot(ggplot2::aes_(~x, ~..density..))+
    ggplot2::geom_histogram(bins=25)+
    ggplot2::coord_cartesian(xlim=c(0, 1))+
    ggplot2::labs(x='frequency of alleles (or living descendants)')+
    wtl::theme_wtl()
}

#' ggplot for 2D lattice
#' @param .data data.frame
#' @param colour column name to colorcode
#' @param size relative size of points
#' @param alpha opacity [0, 1]
#' @param limit for value range
#' @rdname plot
#' @export
gglattice2d = function(.data, colour='clade', alpha=0.66, size=1, limit=max_abs_xyz(.data)) {
    ggplot2::ggplot(.data, ggplot2::aes_(~x, ~y))+
    ggplot2::geom_point(ggplot2::aes_string(colour=colour), alpha=alpha, size=size*80/limit)+
    ggplot2::coord_equal(xlim=limit * c(-1, 1), ylim=limit * c(-1, 1))+
    wtl::theme_wtl()+
    ggplot2::theme(axis.title=ggplot2::element_blank())
}

#' Write animation GIF for serial sections of 3D tumor
#' @param filename string
#' @param ... passed to gglattice2d
#' @param width integer
#' @param height integer
#' @rdname plot
#' @export
save_serial_section = function(.data, filename='serial_section.gif', ..., width=720, height=640) {
    if (!requireNamespace('animation', quietly=TRUE)) {
        stop('ERROR: animation is not installed')
    }
    .lim = max_abs_xyz(.data)
    section_plots = dplyr::group_by_(.data, ~z) %>>%
        dplyr::do(plt={
            gglattice2d(., ..., limit=.lim)+
            ggplot2::geom_hline(yintercept=.$z[1])+
            wtl::theme_wtl()
        })
    animation::saveGIF({for (p in section_plots$plt) {print(p)}},
        filename, outdir=getwd(), interval=0.15,
        ani.width=width, ani.height=height, loop=TRUE, autobrowse=FALSE)
}
