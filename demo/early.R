#!/usr/bin/env Rscript
library(tidyverse)
library(animation)
library(wtl)
library(tumorr)
#########1#########2#########3#########4#########5#########6#########7#########

# tumopp -D2 -Cmoore -O4 -R128 -N128 -w 0 0

(.args = wtl::command_args())
indir = .args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

snapshots = tumorr::read_snapshots()
norigins = snapshots %>% dplyr::filter(time == 0) %>% nrow()
nclades = max(4L, norigins)
founders = head(snapshots, nclades)$id %>% print()
roots = head(snapshots, nclades)$genealogy %>%
    flatten_int() %>% unique() %>% setdiff(founders) %>% print()

.tbl = snapshots %>% dplyr::mutate(
    clade= purrr::map_int(.data$genealogy, ~{setdiff(.x, roots)[1L]}),
    clade= factor(.data$clade, levels=founders)) %>% print()

.lim = max_abs_xyz(snapshots)

plot_snapshot = function(data, time) {
    .N = nrow(data)
    plot_lattice2d(data, 'clade', limit=.lim)+
    labs(title=sprintf('t = %.5f, N =%4d', time, .N))+
    theme(legend.position='none')
}

if (FALSE) {
    .tbl %>%
    tidyr::nest(-time) %>%
    dplyr::sample_n(1) %>%
    {plot_snapshot(.$data[[1]], .$time)}
}

.out = .tbl %>%
    tidyr::nest(-time) %>%
    dplyr::mutate(plt=purrr::pmap(., plot_snapshot)) %>%
    print()

animation::saveGIF({for (.p in .out$plt) {print(.p)}},
    'earlysteps.gif', outdir=getwd(),
    interval=0.15, ani.width=400, ani.height=400, loop=TRUE, autobrowse=FALSE)

if (FALSE) {
    grob = gridExtra::arrangeGrob(grobs=time_plots$plt, nrow=1)
    ggsave('earlysteps.png', grob, width=54, height=1, scale=3, limitsize=FALSE)
}
