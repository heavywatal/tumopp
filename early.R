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

indir = args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

raw = read_tsv('snapshots.tsv.gz') %>>% (?.)
nzero = raw %>>% dplyr::filter(time == 0) %>>% nrow()

snapshots = raw %>>%
    mutate(ancestors= ifelse(is.na(ancestors), id, ancestors),
           ancestors= tumorr::first_ancestors(ancestors, nzero),
           ancestors= as.factor(ancestors)) %>>% (?.)

theme2D =
    theme(panel.background=element_rect(fill='grey90'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())


gglattice2D = function(.data, .colour=~ancestors, .size=1.6, hex=FALSE) {
    .p = ggplot(.data, aes(x, y))
    .p = if (hex) {.p+
        geom_point(aes_(colour=.colour), alpha=0.66, size=.size)+
        scale_colour_hue(na.value='white')
    } else {.p+
        geom_raster(aes_(fill=.colour), alpha=0.66)+
        scale_fill_hue(na.value='white')
    }
    .p + theme2D
}

.lim = rep(maxabs(snapshots)) * c(-1, 1)

plot_snapshot = function(.data, .lim) {
    .t = .data$time[1]
    .N = nrow(.data)
    gglattice2D(.data, .colour= ~ancestors)+
    coord_cartesian(x=.lim, y=.lim)+
    labs(title=sprintf('t = %.5f, N =%4d', .t, .N))+
    theme(legend.position='none')
}

if (FALSE) {
    snapshots %>>%
    dplyr::filter(time==0) %>>%
    plot_snapshot(.lim)
}

time_plots = snapshots %>>%
    group_by(time) %>>%
    do(plt={plot_snapshot(., .lim)})

animation::saveGIF({for (.p in time_plots) {print(.p)}},
    'earlysteps.gif', loop=1, interval=0.1, outdir=getwd(), ani.width=400, ani.height=400)

if (FALSE) {
    grob = gridExtra::arrangeGrob(grobs=time_plots$plt, nrow=1)
    ggsave('earlysteps.png', grob, width=54, height=1, scale=3, limitsize=FALSE)
}
