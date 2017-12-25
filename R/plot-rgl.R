#' plot tumor in 3d with rgl
#' @param .tbl data.frame
#' @param colour column name to colorcode
#' @param .palette name for RColorBrewer::brewer.pal()
#' @param .reverse logical for order of color vector
#' @param .min minimum limit of axes
#' @rdname plot-rgl
#' @export
plot_tumor3d = function(.tbl, colour="clade", .palette="Spectral", .reverse=FALSE, .min=NULL) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("ERROR: rgl is not installed")
  }
  colcol = .tbl[[colour]]
  if (!is.factor(colcol)) {
    colcol = as.factor(colcol)
  }
  num_colors = length(levels(colcol))
  .palette = wtl::brewer_palette(.palette, num_colors)
  if (.reverse) .palette = rev(.palette)
  .col = .palette[colcol]
  .lim = if (is.null(.min)) {
    NULL
  } else {
    max(max_abs_xyz(.tbl), .min) %>% {
      c(-., .)
    }
  }
  with(.tbl, {
    rgl::plot3d(
      x, y, z,
      xlab = "", ylab = "", zlab = "", axes = FALSE,
      type = "s", col = .col, alpha = 1, radius = 1,
      aspect = TRUE, xlim = .lim, ylim = .lim, zlim = .lim
    )
  })
  rgl::box3d()
  rgl::view3d(15, 15, 15, 0.9)
}

#' Shortcut of plot_tumor3d() %>% snapshot3d()
#' @param filename string
#' @param ... passed to plot_tumor3d()
#' @return filename
#' @rdname plot-rgl
#' @export
snapshot_surface = function(.tbl, filename=tempfile("rgl_", fileext = ".png"), ...) {
  on.exit(rgl::rgl.close())
  rgl::open3d(useNULL = FALSE)
  dplyr::filter(.tbl, .data$surface) %>% plot_tumor3d(...)
  rgl::snapshot3d(filename)
  filename
}
