#' ggplot for frequency spectrum
#' @param freqs numeric vector
#' @rdname plot
#' @export
histogram_freqspec = function(freqs) {
  tibble::tibble(x = freqs) %>%
    ggplot2::ggplot(ggplot2::aes_(~x, ~..density..)) +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "frequency of alleles (or living descendants)")
}

#' ggplot for 2D lattice
#' @param .tbl tbl with extant cells
#' @param colour column name to colorcode
#' @param size relative size of points
#' @param alpha opacity `[0, 1]`
#' @param limit for value range
#' @rdname plot
#' @export
plot_lattice2d = function(.tbl, colour="clade", alpha=1, size=1, limit=max_abs_xyz(.tbl)) {
  size = size * 96 / (limit - 0.5)
  ggplot2::ggplot(.tbl, ggplot2::aes_(~x, ~y)) +
    ggplot2::geom_point(ggplot2::aes_string(colour = colour), alpha = alpha, size = size) +
    ggplot2::coord_equal(xlim = limit * c(-1, 1), ylim = limit * c(-1, 1))
}

# #######1#########2#########3#########4#########5#########6#########7#########

#' Plot genealogy
#' @param .tbl tbl from layout_genealogy()
#' @param xmax numeric
#' @param colour character
#' @return gg
#' @rdname plot-igraph
#' @export
plot_genealogy = function(.tbl, xmax=max(.tbl$ageend), colour="dodgerblue") {
  ggplot2::ggplot(.tbl) +
    ggplot2::geom_segment(ggplot2::aes_(~age, ~pos, xend = ~ageend, yend = ~posend), alpha = 0.3, size = 0.3) +
    ggplot2::geom_point(
      data = dplyr::filter(.tbl, .data$extant),
      ggplot2::aes_(x = ~ageend, y = ~posend), size = 0.8, colour = colour, alpha = 0.2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), expand = FALSE)
}

#' Plot age histogram
#' @param alpha opacity `[0, 1]`
#' @param ... passed to aes_()
#' @return gg
#' @rdname plot-igraph
#' @export
plot_bar_age = function(.tbl, xmax=max(.tbl$ageend), alpha=1.0, ...) {
  dplyr::filter(.tbl, .data$extant) %>%
    ggplot2::ggplot(ggplot2::aes_(~ageend, ...)) +
    ggplot2::geom_bar(alpha = alpha) +
    ggplot2::coord_cartesian(xlim = c(0, xmax))
}

# #######1#########2#########3#########4#########5#########6#########7#########

#' Plot capture_rate ~ nsam of biopsy
#' @param data tbl from summarize_capture_rate()
#' @rdname plot-biopsy
#' @export
plot_capture_rate = function(data) {
  ggplot2::ggplot(data, ggplot2::aes_(~nsam, ~capture_rate)) +
    ggplot2::stat_summary(fun.y = mean, geom = "bar", alpha = 0.6) +
    ggplot2::geom_jitter(size = 2, alpha = 0.3, width = 0.25, height = 0, colour = "dodgerblue") +
    ggplot2::stat_summary(fun.data = wtl::mean_sd, geom = "errorbar", width = 0.2) +
    ggplot2::coord_cartesian(ylim = c(0, 1))
}

# #######1#########2#########3#########4#########5#########6#########7#########

#' Plot serial sections of 3D tumor
#' @param .tbl tbl with extant cells
#' @inheritParams ggplot2::ggsave
#' @param ... passed to plot_lattice2d
#' @rdname plot-section
#' @export
save_serial_section = function(.tbl, filename="png/section_%03d.png", scale=6, dpi=72, ...) {
  .lim = max_abs_xyz(.tbl)
  tidyr::nest(.tbl, -.data$z) %>%
    dplyr::arrange(.data$z) %>%
    dplyr::mutate(i = seq_len(nrow(.))) %>%
    purrr::pwalk(function(z, data, i) {
      .outfile = sprintf(filename, i)
      # TODO: fix color for each lineage
      .p = plot_lattice2d(data, ..., limit = .lim) +
        ggplot2::geom_hline(yintercept = z[1L], colour = "#999999", size = 1.5) +
        ggplot2::labs(title = sprintf("z =%4.1f", z)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        )
      ggplot2::ggsave(.outfile, .p, width = 1, height = 1, scale = scale, dpi = dpi)
    })
}
