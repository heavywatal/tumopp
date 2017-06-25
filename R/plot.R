#' ggplot for frequency spectrum
#' @param freqs numeric vector
#' @rdname plot
#' @export
histogram_freqspec = function(freqs) {
    tibble::tibble(x=freqs) %>%
    ggplot2::ggplot(ggplot2::aes_(~x, ~..density..))+
    ggplot2::geom_histogram(bins=25)+
    ggplot2::coord_cartesian(xlim=c(0, 1))+
    ggplot2::labs(x='frequency of alleles (or living descendants)')+
    wtl::theme_wtl()
}

#' ggplot for 2D lattice
#' @param .tbl tbl with extant cells
#' @param colour column name to colorcode
#' @param size relative size of points
#' @param alpha opacity [0, 1]
#' @param limit for value range
#' @rdname plot
#' @export
plot_lattice2d = function(.tbl, colour='clade', alpha=0.66, size=1, limit=max_abs_xyz(.tbl)) {
    ggplot2::ggplot(.tbl, ggplot2::aes_(~x, ~y))+
    ggplot2::geom_point(ggplot2::aes_string(colour=colour), alpha=alpha, size=size*80/limit)+
    ggplot2::coord_equal(xlim=limit * c(-1, 1), ylim=limit * c(-1, 1))+
    wtl::theme_wtl()+
    ggplot2::theme(axis.title=ggplot2::element_blank())
}

#########1#########2#########3#########4#########5#########6#########7#########

#' Plot genealogy
#' @param .tbl tbl from layout_genealogy()
#' @param xmax numeric
#' @param colour character
#' @return gg
#' @rdname plot-igraph
#' @export
plot_genealogy = function(.tbl, xmax=max(.tbl$ageend), colour='dodgerblue') {
    ggplot2::ggplot(.tbl)+
    ggplot2::geom_segment(ggplot2::aes_(~age, ~pos, xend=~ageend, yend=~posend), alpha=0.3, size=0.3)+
    ggplot2::geom_point(data=dplyr::filter_(.tbl, ~extant),
        ggplot2::aes_(x=~ageend, y=~posend), size=0.8, colour=colour, alpha=0.2)+
    ggplot2::coord_cartesian(xlim=c(0, xmax), expand=FALSE)+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title.y=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        panel.grid.major.y=ggplot2::element_blank())
}

#' Plot age histogram
#' @param alpha opacity [0, 1]
#' @param ... passed to aes_()
#' @return gg
#' @rdname plot-igraph
#' @export
plot_bar_age = function(.tbl, xmax=max(.tbl$ageend), alpha=1.0, ...) {
    dplyr::filter_(.tbl, ~extant) %>%
    ggplot2::ggplot(ggplot2::aes_(~ageend, ...))+
    ggplot2::geom_bar(alpha=alpha)+
    ggplot2::coord_cartesian(xlim=c(0, xmax))+
    wtl::theme_wtl()+
    ggplot2::theme(
        axis.title=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        axis.ticks.y=ggplot2::element_blank(),
        panel.grid.major.y=ggplot2::element_blank())
}

#########1#########2#########3#########4#########5#########6#########7#########

#' Write animation GIF for serial sections of 3D tumor
#' @param .tbl tbl with extant cells
#' @param filename string
#' @param ... passed to gglattice2d
#' @param width integer
#' @param height integer
#' @rdname plot-animation
#' @export
save_serial_section = function(.tbl, filename='serial_section.gif', ..., width=720, height=640) {
    if (!requireNamespace('animation', quietly=TRUE)) {
        stop('ERROR: animation is not installed')
    }
    .lim = max_abs_xyz(.tbl)
    section_plots = dplyr::group_by_(.tbl, ~z) %>%
        dplyr::do(plt={
            plot_lattice2d(., ..., limit=.lim)+
            ggplot2::geom_hline(yintercept=.$z[1])+
            wtl::theme_wtl()
        })
    animation::saveGIF({for (p in section_plots$plt) {print(p)}},
        filename, outdir=getwd(), interval=0.15,
        ani.width=width, ani.height=height, loop=TRUE, autobrowse=FALSE)
}
