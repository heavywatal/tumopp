#!/usr/bin/env Rscript
library(pipeR)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

(.argv = commandArgs(trailingOnly=FALSE))
(..file.. = sub('--file=', '', grep('--file=', .argv, value=TRUE)))
.project = dirname(normalizePath(..file..))

.argv = commandArgs(trailingOnly=TRUE)
indir = .argv[1]
indir = ifelse(is.na(indir), '.', indir)
setwd(indir)

conf = wtl::read.conf('program_options.conf') %>>% (?.)
history = read_tsv('mutation_history.tsv.gz') %>>% (?.)
population = read_tsv('population.tsv.gz') %>>% (?.)
evolution = read_tsv('evolution_history.tsv.gz',
                col_types=list(sites=col_character())) %>>% (?.)

unnested = population %>>%
    mutate(sites=strsplit(sites, ':')) %>>%
    unnest(sites) %>>%
    dplyr::select(-size) %>>%
    full_join(history %>>% add_rownames('sites'), by='sites') %>>%
    dplyr::select(-sites) %>>%
    arrange(time, x, y) %>>% (?.)

unnested_evolution = evolution %>>%
    mutate(sites=strsplit(sites, ':')) %>>%
    unnest(sites) %>>% (?.)

#source(file.path(.project, 'sample.R'))

#########1#########2#########3#########4#########5#########6#########7#########

trans_coord_hex = function(mtrx) {mtrx %>>%
    mutate(y= y + x * 0.5) %>>%
    mutate(x= x * sqrt(3.0 / 4.0))
}

#  hexagonal close packed
trans_coord_hcc = function(mtrx) {
    trans_coord_hex(mtrx) %>>%
    mutate(x= x + ifelse(z %% 2 == 1, sqrt(3) / 3, 0)) %>>%
    mutate(z= z * sqrt(2.0 / 3.0))
}

#  face centered cubic (cubic close packed)
trans_coord_fcc = function(mtrx) {
    trans_coord_hex(mtrx) %>>%
    mutate(x= x + z / sqrt(3.0)) %>>%
    mutate(z= z * sqrt(2.0 / 3.0))
}

gglattice2D = function(.data, .colour=~size, .size=1.6) {
    .p = ggplot(.data, aes(x, y))
    .p = if (conf[['coord']] == 'hex') {.p+
        geom_point(aes_(colour=.colour), alpha=0.66, size=.size)+
        scale_colour_hue(na.value='white')
    } else {.p+
        geom_raster(aes_(fill=.colour), alpha=0.66)+
        scale_fill_hue(na.value='white')
    }
    .p
}

theme2D =
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())

#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] == 2) {

#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['coord']] == 'hex') {

population = trans_coord_hex(population)
evolution = trans_coord_hex(evolution)
unnested = trans_coord_hex(unnested)
unnested_evolution = trans_coord_hex(unnested_evolution)

}  # fi hex
#########1#########2#########3#########4#########5#########6#########7#########

.p = population %>>%
    mutate(ancestor=as.factor(ancestor)) %>>%
    gglattice2D(~ancestor) + theme2D
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

.p = unnested %>>%
    mutate(size=ifelse(60 <= size & size <= 100, size, NA)) %>>%
    dplyr::select(-starts_with('origin_'), -effect, -mutant) %>>%
    dplyr::filter(!duplicated(.)) %>>%
    group_by(x, y) %>>%
    dplyr::filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
    mutate(size=as.factor(size)) %>>%
    gglattice2D + theme2D
#.p
ggsave('early_mid_mutations.png', .p, width=7, height=7)


.heat_colours = c('#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')
plot_mutation_history_2d = function(.data) {
    .p = ggplot(.data, aes(x, y))
    .p = if (conf[['coord']] == 'hex') {.p+
        geom_point(aes(colour=size), alpha=0.5)+
        scale_colour_gradientn(colours=.heat_colours, na.value='white')
    } else {.p +
        geom_raster(aes(fill=size), alpha=0.5)+
        scale_fill_gradientn(colours=.heat_colours, na.value='white')
    }
    .p + theme2D
}
.p = plot_mutation_history_2d(unnested)
#.p
ggsave('gradient.png', .p, width=7, height=7)


maxabs = with(evolution %>>% dplyr::filter(size <= 256), max(abs(c(x, y)))) %>>% (?.)

plot_early_evolution_2d = function(.time) {
    .data = unnested_evolution %>>%
        dplyr::filter(time==.time) %>>%
        mutate(size=as.factor(sites))
    .data %>>%
        bind_rows(data_frame(x=maxabs+seq_len(4), y=maxabs+seq_len(4),
                          fitness=-1, size=as.factor(seq_len(4)))) %>>%
    gglattice2D(~ancestor, 10)+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    labs(title=paste0('t = ', .time, ', N = ', nrow(.data)))+
    theme2D+
    theme(legend.position='none')
}
#print(plot_early_evolution_2d(10))

if (FALSE) {
animation::saveGIF({
    for (.t in unique(evolution$time)) {
        print(plot_early_evolution_2d(.t))
    }},
    'evolution.gif', loop=1, interval=0.1, outdir=getwd())
}

#########1#########2#########3#########4#########5#########6#########7#########
} else {  # 3D
#########1#########2#########3#########4#########5#########6#########7#########

library(animation)
library(rgl)

if (conf[['coord']] == 'hex') {
    population = population %>>% trans_coord_fcc
}
maxabs = with(population, max(abs(c(x, y, z)))) %>>% (?.)

plot_early_mutations_3d_section = function(.z=0, .data) {.data %>>%
    dplyr::filter(z==.z) %>>%
    mutate(ancestor= as.factor(ancestor)) %>>%
    gglattice2D(~ancestor, 5)+
    geom_hline(yintercept=.z)+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    theme2D+
    theme(legend.position='right')
}
#plot_early_mutations_3d_section(0, population)

#animation::saveGIF({
#    for (i in unique(population %>>% arrange(z) %>>% (z))) {
#        print(plot_early_mutations_3d_section(i, population))
#    }},
#    'serial_section.gif', loop=TRUE, interval=0.15, outdir=getwd())

#########1#########2#########3#########4#########5#########6#########7#########

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(-25, 15, 40, zoom=0.9)
clear3d()
axes3d()
with(population %>>% dplyr::filter(sqrt(x^2 + y^2 + z^2)>12), spheres3d(x, y, z,
                radius=1, col=ancestor, alpha=0.8))
title3d('', '', 'x', 'y', 'z')
writeWebGL()


if (FALSE) {  # Thinking about hex 3D neighbors

## transformation

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
clear3d()
axes3d()
.r = 2
seq(-.r, .r) %>>%
    (expand.grid(x=., y=., z=.)) %>>%
#    trans_coord_hcc %>>%
    trans_coord_fcc %>>%
    mutate(r= sqrt(x*x+y*y+z*z)) %>>%
#    dplyr::filter(r <= 3) %>>% (?.) %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


## cells, volume, radius

.r = 60
.cells = seq(-.r, .r) %>>%
    (expand.grid(grid='regular', x=., y=., z=.)) %>>%
    tbl_df %>>%
    bind_rows(trans_coord_fcc(.) %>>% mutate(grid='hex')) %>>%
    mutate(r= sqrt(x*x+y*y+z*z)) %>>%
    dplyr::filter(r < 0.75 * .r) %>>%
    group_by(grid) %>>%
    mutate(n = order(order(r))) %>>%
    (?.)

.curve = data_frame(grid='regular', r=seq_len(100) / 2) %>>%
    mutate(n=pi * 4 * r^3 / 3) %>>%
    bind_rows(mutate(., grid='hex', n= sqrt(2) * n)) %>>%
    (?.)

.p = .cells %>>%
    ggplot(aes(n, r, colour=grid))+
    geom_point()+
    geom_path(data=.curve)
.p
ggsave('radius-nodes.png', .p, width=7, height=7)

.cells %>>% group_by(grid) %>>% tally
summary(.cells)

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
clear3d()
axes3d()
.cells %>>%
    dplyr::filter(r > 47) %>>% (?.) %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


## minimum

.hex_xy = read_csv('x,y,z
0,0,0
1,0,0
0,1,0
1,0,-1
')

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
clear3d()
axes3d()
.hex_xy %>>%
    trans_coord_fcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')

## neighbors

.hex_xy = read_csv('x,y,z
0,0,0
0,1,0
0,-1,0
-1,0,0
-1,1,0
1,0,0
1,-1,0')

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(20, 10, 60)
clear3d()
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x > 0 | (x==0 & y==0)) %>>% mutate(z=-1)) %>>%
    bind_rows(.hex_xy %>>% dplyr::filter(x < 0 | (x==0 & y==0)) %>>% mutate(z=1)) %>>% (?.) %>>%
    trans_coord_fcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(40, 20, 60)
clear3d()
axes3d()
.hex_xy %>>%
    bind_rows(.hex_xy %>>% mutate(z=-1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=1) %>>% dplyr::filter(x < 0 | (x==0 & y==0))) %>>% (?.) %>>%
    bind_rows(.hex_xy %>>% mutate(z=5)) %>>%
    bind_rows(.hex_xy %>>% mutate(z=4) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=6) %>>% dplyr::filter(x > 0 | (x==0 & y==0))) %>>%
    trans_coord_hcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


}  # fi 3D neighbors

#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 2D/3D
#########1#########2#########3#########4#########5#########6#########7#########
