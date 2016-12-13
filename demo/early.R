#!/usr/bin/env Rscript
library(pipeR)
library(tidyverse)
library(animation)
library(wtl)
library(tumorr)
#########1#########2#########3#########4#########5#########6#########7#########

(.args = wtl::command_args())
indir = .args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

snapshots = tumorr::read_snapshots()$snapshots[[1]]
norigins = snapshots %>>% dplyr::filter(time == 0) %>>% nrow()
nclades = max(4L, norigins)
founders = snapshots %>>% group_by(time) %>>% dplyr::filter(n() == nclades) %>>% (id)
roots = snapshots %>>% dplyr::filter(id < max(founders)) %>>% (id) %>>% setdiff(founders)

.data = snapshots %>>% dplyr::mutate_(
    clade= ~purrr::map_int(genealogy, ~{setdiff(.x, roots)[1L]}),
    clade= ~factor(clade, levels=founders)) %>>% (?.)

.lim = max_abs_xyz(snapshots)

plot_snapshot = function(.data) {
    .t = .data$time[1]
    .N = nrow(.data)
    plot_lattice2d(.data, 'clade', limit=.lim)+
    labs(title=sprintf('t = %.5f, N =%4d', .t, .N))+
    theme(legend.position='none')
}

if (FALSE) {
    .data %>>%
    dplyr::filter(time==sample(unique(time), 1)) %>>%
    plot_snapshot()
}

.out = .data %>>%
    dplyr::group_by(time) %>>%
    dplyr::do(plt={plot_snapshot(.)})

animation::saveGIF({for (.p in .out$plt[-1]) {print(.p)}},
    'earlysteps.gif', outdir=getwd(),
    interval=0.15, ani.width=400, ani.height=400, loop=TRUE, autobrowse=FALSE)

if (FALSE) {
    grob = gridExtra::arrangeGrob(grobs=time_plots$plt, nrow=1)
    ggsave('earlysteps.png', grob, width=54, height=1, scale=3, limitsize=FALSE)
}
