#!/usr/bin/Rscript
library(pipeR)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(animation)

(.argv = commandArgs(trailingOnly=FALSE))
(..file.. = sub('--file=', '', grep('--file=', .argv, value=TRUE)))

.argv = commandArgs(trailingOnly=TRUE)
indir = .argv[1]
indir = ifelse(is.na(indir), '.', indir)

conf = wtl::read.conf(file.path(indir, 'program_options.conf')) %>>% (?.)
history = read_tsv(file.path(indir, 'mutation_history.tsv.gz')) %>>% (?.)
population = read_tsv(file.path(indir, 'population.tsv.gz')) %>>% (?.)
evolution = read_tsv(file.path(indir, 'evolution_history.tsv.gz'),
                col_types=list(sites=col_character())) %>>% (?.)

unnested = population %>>%
    mutate(sites=strsplit(sites, '\\|')) %>>%
    unnest(sites) %>>%
    bind_rows(population %>>% filter(sites=='')) %>>%
    dplyr::select(-size) %>>%
    full_join(history %>>% name_rows, by=c(sites='.rownames')) %>>%
    dplyr::select(-sites) %>>%
    arrange(time, x, y) %>>% (?.)

unnested_evolution = evolution %>>%
    mutate(sites=strsplit(sites, '\\|')) %>>%
    unnest(sites) %>>%
    filter(extract_numeric(sites) <= 4) %>>% (?.)

source(file.path(dirname(..file..), 'sample.R'))


#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] == 2) {

#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['coord']] == 'hex') {

unnested = unnested %>>%
    mutate(y= y + x * 0.5, x= x * sqrt(3) * 0.5)

unnested_evolution = unnested_evolution %>>%
    mutate(y= y + x * 0.5, x= x * sqrt(3) * 0.5)

}  # fi hex
#########1#########2#########3#########4#########5#########6#########7#########


gglattice2D = function(.data, .size=1.6) {
    .p = ggplot(.data, aes(x, y))
    .p = if (conf[['coord']] == 'hex') {.p+
        geom_point(aes(colour=size), alpha=0.66, size=.size)+
        scale_colour_hue(na.value='white')
    } else {.p+
        geom_raster(aes(fill=size), alpha=0.66)+
        scale_fill_hue(na.value='white')
    }
    .p
}

theme2D =
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())

.p = unnested %>>%
    mutate(size=ifelse(size <= 4, size, NA)) %>>%
    dplyr::select(-starts_with('origin_'), -effect) %>>%
    filter(!duplicated(.)) %>>%
    group_by(x, y) %>>%
    filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
    mutate(size=as.factor(size)) %>>%
    gglattice2D + theme2D
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

.p = unnested %>>%
    mutate(size=ifelse(60 <= size & size <= 100, size, NA)) %>>%
    dplyr::select(-starts_with('origin_'), -effect) %>>%
    filter(!duplicated(.)) %>>%
    group_by(x, y) %>>%
    filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
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


maxabs = with(evolution %>>% filter(size <= 256), max(abs(c(x, y)))) %>>% (?.)

plot_early_evolution_2d = function(.time) {
    .data = unnested_evolution %>>%
        filter(time==.time) %>>%
        mutate(size=as.factor(sites))
    .data %>>%
        bind_rows(data_frame(x=maxabs+seq_len(4), y=maxabs+seq_len(4),
                          fitness=-1, size=as.factor(seq_len(4)))) %>>%
    gglattice2D(10)+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    labs(title=paste0('t = ', .time, ', N = ', nrow(.data)))+
    theme2D+
    theme(legend.position='none')
}
#print(plot_early_evolution_2d(10))

animation::saveGIF({
    for (.t in unique(evolution$time)) {
        print(plot_early_evolution_2d(.t))
    }},
    'evolution.gif', loop=1, interval=0.1, outdir=getwd())

#########1#########2#########3#########4#########5#########6#########7#########
} else {  # 3D
#########1#########2#########3#########4#########5#########6#########7#########

#  hexagonal close packed
trans_coord_hcc = function(mtrx) {mtrx %>>%
    mutate(y= y + x * 0.5) %>>%
    mutate(x= x * sqrt(3) * 0.5) %>>%
    mutate(x= x + ifelse(z %% 2 == 1, sqrt(3) / 3, 0)) %>>%
    mutate(z= z * 0.816497)  # TODO
}

#  face centered cubic (cubic close packed)
trans_coord_fcc = function(mtrx) {mtrx %>>%
    mutate(y= y + x * 0.5) %>>%
    mutate(x= x * sqrt(3) * 0.5) %>>%
    mutate(x= x + z * sqrt(3) / 3) %>>%
    mutate(z= z * 0.816497)  # TODO
}

plot_early_mutations_3d_section = function(.z=0, .data) {.data %>>%
    filter(z==.z) %>>%
    bind_rows(data_frame(x=maxabs+seq_len(8), y=maxabs+seq_len(8), z=.z,
                          fitness=-1, size=as.factor(seq_len(8)))) %>>%
    #(?.) %>>% (?str(.)) %>>%
    gglattice2D(5)+
    geom_hline(yintercept=.z)+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    theme2D+
    theme(legend.position='right')
}
#plot_early_mutations_3d_section(0, early3d)

early3d = unnested %>>%
    filter(size <= 13) %>>%
    dplyr::select(-effect, -starts_with('origin_')) %>>%
    filter(!duplicated(.)) %>>%
    trans_coord_fcc %>>%
    mutate(size=as.factor(size)) %>>% (?.)

maxabs = with(population, max(abs(c(x, y, z)))) %>>% (?.)

animation::saveGIF({
    for (i in unique(early3d %>>% arrange(z) %>>% (z))) {
        print(plot_early_mutations_3d_section(i, early3d))
    }},
    'serial_section.gif', loop=TRUE, interval=0.15, outdir=getwd())


#########1#########2#########3#########4#########5#########6#########7#########

library(rgl)

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(-25, 15, 40, zoom=0.9)
clear3d()
axes3d()
with(early3d %>>% filter(sqrt(x^2 + y^2 + z^2)>12), spheres3d(x, y, z,
                radius=1, col=marker, alpha=0.8))
title3d('', '', 'x', 'y', 'z')
writeWebGL()


if (FALSE) {  # Thinking about hex 3D neighbors
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
    bind_rows(.hex_xy %>>% filter(x > 0 | (x==0 & y==0)) %>>% mutate(z=-1)) %>>%
    bind_rows(.hex_xy %>>% filter(x < 0 | (x==0 & y==0)) %>>% mutate(z=1)) %>>% (?.) %>>%
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
    bind_rows(.hex_xy %>>% mutate(z=-1) %>>% filter(x < 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=1) %>>% filter(x < 0 | (x==0 & y==0))) %>>% (?.) %>>%
    bind_rows(.hex_xy %>>% mutate(z=5)) %>>%
    bind_rows(.hex_xy %>>% mutate(z=4) %>>% filter(x > 0 | (x==0 & y==0))) %>>%
    bind_rows(.hex_xy %>>% mutate(z=6) %>>% filter(x > 0 | (x==0 & y==0))) %>>%
    trans_coord_hcc %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')


if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
rgl::clear3d()
rgl::view3d(70, 5, 60)
clear3d()
axes3d()
.r = 2; seq(-.r, .r) %>>%
    (expand.grid(x=., y=., z=.)) %>>%
#    trans_coord_hcc %>>%
    trans_coord_fcc %>>%
    mutate(r= sqrt(x*x+y*y+z*z)) %>>%
#    filter(r <= 3) %>>% (?.) %>>%
with(spheres3d(x, y, z, color='#009999', radius=0.51, alpha=0.6))
title3d('', '', 'x', 'y', 'z')

#  even
max(z, d + z/2)

#  odd


}  # fi 3D neighbors

#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 2D/3D
#########1#########2#########3#########4#########5#########6#########7#########
