#!/usr/bin/env Rscript
command_args = function() {
    .argv = commandArgs(trailingOnly=FALSE)
    l = list()
    l$file = sub('^--file=', '', grep('^--file=', .argv, value=TRUE))
    l$srcdir = dirname(normalizePath(l$file))
    l$args = grep('^[^-]', .argv[-1], value=TRUE)
    return(l)
}
(args = command_args())

library(pipeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(animation)

library(tumorr)
#load_all('~/git/tumorr')

indir = args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

conf = tumorr::read_conf() %>>% (?.)
snapshots = tumorr::read_snapshots(conf)
nzero = snapshots %>>% dplyr::filter(time == 0) %>>% nrow()
anc_colours = max(4, nzero)
anc_ids = exclusive_ancestors(raw, anc_colours)

snapshots = tumorr::filter_ancestors(snapshots, anc_ids)

.lim = maxabs(snapshots)

plot_snapshot = function(.data) {
    .t = .data$time[1]
    .N = nrow(.data)
    gglattice2D(.data, 'ancestors', .lim)+
    labs(title=sprintf('t = %.5f, N =%4d', .t, .N))+
    theme(legend.position='none')
}

if (FALSE) {
    snapshots %>>%
    dplyr::filter(time==sample(unique(time), 1)) %>>%
    plot_snapshot()
}

time_plots = snapshots %>>%
    group_by(time) %>>%
    do(plt={plot_snapshot(., .lim)})

animation::saveGIF({for (.p in time_plots$plt) {print(.p)}},
    'earlysteps.gif',
    , outdir=getwd(), interval=0.15, ani.width=400, ani.height=400, loop=TRUE, autobrowse=FALSE)

if (FALSE) {
    grob = gridExtra::arrangeGrob(grobs=time_plots$plt, nrow=1)
    ggsave('earlysteps.png', grob, width=54, height=1, scale=3, limitsize=FALSE)
}
