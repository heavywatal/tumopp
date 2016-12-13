#!/usr/bin/env Rscript
library(pipeR)
library(tidyverse)
library(wtl)
library(tumorr)
#########1#########2#########3#########4#########5#########6#########7#########

tumopp('-h')
result = tumopp('-D3 -Chex -k1e9 -Lstep -a0.2 -m1.0' %>>% wtl::split_chr())

result$population[[1]] %>>%
    filter_extant() %>>%
    dplyr::filter(z == 0) %>>%
    plot_lattice2d()

result$population[[1]] %>>%
    filter_extant() %>>%
    dplyr::filter(z == 0) %>>%
    dplyr::filter(surface) %>>%
    plot_lattice2d()

#########1#########2#########3#########4#########5#########6#########7#########
if (result[['dimensions']] > 2L) {

library(rgl)

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
result$population[[1]] %>>%
    dplyr::filter(surface) %>>%
    plot_tumor3d()
title3d('', '', 'x', 'y', 'z')

.outfile = snapshot_surface(result$population[[1]])
system(sprintf('open %s', .outfile))

writeWebGL('.', 'rgl.html', snapshot=FALSE, width=600, height=600)

#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 3D
#########1#########2#########3#########4#########5#########6#########7#########
