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
#load_all('~/git/tumorr')

indir = args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

conf = tumorr::read_conf() %>>% (?.)
snapshots = tumorr::read_snapshots(conf)
nzero = snapshots %>>% dplyr::filter(time == 0) %>>% nrow()
anc_colours = max(4, nzero)
anc_ids = exclusive_ancestors(snapshots, anc_colours)

population = tumorr::read_population(conf) %>>% (?.)
#nzero = sum(is.na(population$ancestors))

final_cells = population %>>%
    dplyr::filter(death == 0) %>>%
    tumorr::filter_ancestors(anc_ids) %>>% (?.)


#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] == 2) {

.p = final_cells %>>%
    gglattice2D('ancestors')
.p
ggsave('ancestors.png', .p, width=8, height=7)

#########1#########2#########3#########4#########5#########6#########7#########
} else {  # 3D
#########1#########2#########3#########4#########5#########6#########7#########

library(rgl)
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
