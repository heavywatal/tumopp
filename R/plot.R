#' ggplot for 2D lattice
#' @param .data a data.frame
#' @param .color a column name to colorcode
#' @param .palette name for RColorBrewer::brewer.pal()
#' @param limit for value range
#' @rdname plot
#' @export
gglattice2d = function(.data, .color='ancestor', .palette='Spectral', limit=maxabs(.data)) {
    ggplot2::ggplot(.data, ggplot2::aes(x, y))+
    ggplot2::geom_point(ggplot2::aes_string(colour=.color), alpha=0.66, size=80/limit)+
    ggplot2::scale_colour_brewer(palette=.palette)+
    ggplot2::expand_limits(x=limit * c(-1, 1), y=limit * c(-1, 1))+
    ggplot2::theme_grey()+
    ggplot2::theme(panel.background=ggplot2::element_rect(fill='grey90'))+
    ggplot2::theme(panel.grid=ggplot2::element_blank())+
    ggplot2::theme(axis.title=ggplot2::element_blank())
}

#' Write animation GIF for serial sections of 3D tumor
#' @param filename a string
#' @param width an integer
#' @param height an integer
#' @rdname plot
#' @export
save_serial_section = function(filename='serial_section.gif', .data, width=720, height=640) {
    .lim = maxabs(.data)
    section_plots = dplyr::group_by(.data, z) %>>%
        dplyr::do(plt={
            gglattice2d(., limit=.lim)+
            ggplot2::geom_hline(yintercept=.$z[1])
        })
    animation::saveGIF({for (p in section_plots$plt) {print(p)}},
        filename, outdir=getwd(), interval=0.15,
        ani.width=width, ani.height=height, loop=TRUE, autobrowse=FALSE)
}

#' plot tumor in 3d with rgl
#' @rdname plot
#' @export
plot_tumor3d = function(.data, .color='ancestor', .palette='Spectral') {
    .data[.color] = as.factor(.data[[.color]])
    num_colors = length(levels(.data[[.color]]))
    .palette = RColorBrewer::brewer.pal(num_colors, .palette)
    thres = sphere_radius(nrow(.data)) * 0.6
    .data %>>%
        dplyr::filter(sqrt(x^2 + y^2 + z^2) > thres) %>>%
        with(rgl::spheres3d(x, y, z, color=.palette[.[[.color]]],
                            radius=1, alpha=1))
    rgl::box3d()
    rgl::view3d(15, 15, 15, 0.9)
}
