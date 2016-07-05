#!/usr/bin/env Rscript
library(pipeR)
library(dplyr)
library(ggplot2)

library(wtl)
library(tumorr)
#load_all('~/git/tumorr')
(.args = wtl::command_args())
#########1#########2#########3#########4#########5#########6#########7#########

indir = .args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

conf = read_conf()
population = read_population(conf) %>>% colorcode_survivors()

population %>>%
    dplyr::filter(z == 0) %>>%
    gglattice2d('ancestor') %>>%
    (ggsave('section_z0.png', ., width=8, height=7))

population %>>%
    dplyr::filter(z == 0) %>>%
    dplyr::filter(surface) %>>%
    gglattice2d('ancestor') %>>%
    (ggsave('section_z0_surface.png', ., width=8, height=7))

#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] > 2) {

library(rgl)
if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
population %>>%
    dplyr::filter(surface) %>>%
    plot_tumor3d()
title3d('', '', 'x', 'y', 'z')
writeWebGL('.', 'rgl.html', snapshot=FALSE, width=600, height=600)
rgl::snapshot3d('colorcode3d.png')

#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 3D
#########1#########2#########3#########4#########5#########6#########7#########
