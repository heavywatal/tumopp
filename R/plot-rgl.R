#' plot tumor in 3d with rgl
#' @param .data data.frame
#' @param colour column name to colorcode
#' @param .palette name for RColorBrewer::brewer.pal()
#' @rdname plot-rgl
#' @export
plot_tumor3d = function(.data, colour='clade', .palette='Spectral') {
    if (!requireNamespace('rgl', quietly=TRUE)) {
        stop('ERROR: rgl is not installed')
    }
    colcol = .data[[colour]]
    if (!is.factor(colcol)) {
        colcol = as.factor(colcol)
    }
    num_colors = length(levels(colcol))
    .colors = wtl::brewer_palette(.palette, num_colors)
    with(.data, {
        rgl::spheres3d(x, y, z, color=.colors[colcol], radius=1, alpha=1)})
    rgl::box3d()
    rgl::view3d(15, 15, 15, 0.9)
}

#' Shortcut of plot_tumor3d() %>>% snapshot3d()
#' @param filename string
#' @param ... passed to plot_tumor3d()
#' @return filename
#' @rdname plot-rgl
#' @export
snapshot_surface = function(.data, filename=tempfile('rgl_', fileext='.png'), ...) {
    on.exit(rgl::rgl.close())
    rgl::open3d(useNULL=FALSE)
    dplyr::filter_(.data, ~surface) %>>% plot_tumor3d(...)
    rgl::snapshot3d(filename)
    filename
}
