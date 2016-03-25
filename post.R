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

library(tumorr)

indir = args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

conf = wtl::read.conf('program_options.conf') %>>% (?.)
population = read_tsv('population.tsv.gz') %>>% (?.)
is_hex = conf[['coord']] == 'hex'

(n0 = sum(is.na(population$ancestors)))

final_cells = population %>>%
    dplyr::filter(death == 0) %>>%
    mutate(ancestors= ifelse(is.na(ancestors), id, ancestors),
           ancestors= tumorr::first_ancestors(ancestors, n0),
           ancestors= as.factor(ancestors)) %>>% (?.)


#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] == 2) {

if (is_hex) {
    population = trans_coord_hex(population, hex=is_hex)
}  # fi hex

.p = final_cells %>>% gglattice2D(~ancestors) + theme2D
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

#########1#########2#########3#########4#########5#########6#########7#########
} else {  # 3D
#########1#########2#########3#########4#########5#########6#########7#########

library(animation)
library(rgl)

if (is_hex) {
    population = trans_coord_fcc(final_cells)
}

.lim = rep(maxabs(final_cells)) * c(-1, 1)

plot_section = function(.data) {.data %>>%
    gglattice2D(~ancestors, 5, hex=is_hex)+
    geom_hline(yintercept=.data$z[1])+
    coord_cartesian(x=.lim, y=.lim)+
    theme2D+
    theme(legend.position='right')
}
plot_section(final_cells %>>% dplyr::filter(z == 0))

section_plots = final_cells %>>%
    group_by(z) %>>%
    do(plt={plot_section(.)})

animation::saveGIF({
   for (p in section_plots) {
       print(p)
   }},
   'serial_section.gif', loop=TRUE, interval=0.15, outdir=getwd())

#########1#########2#########3#########4#########5#########6#########7#########

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(-25, 15, 40, zoom=0.9)
clear3d()
axes3d()
thres = sphere_radius(nrow(final_cells)) * 0.6
final_cells %>>%
    dplyr::filter(sqrt(x^2 + y^2 + z^2) > thres) %>>%
    with(spheres3d(x, y, z,
                   radius=1, col=ancestors, alpha=0.8))
title3d('', '', 'x', 'y', 'z')
writeWebGL('.', 'rgl.html', snapshot=FALSE, width=600, height=600)


#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 2D/3D
#########1#########2#########3#########4#########5#########6#########7#########
