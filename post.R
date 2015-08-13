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

unnested = population %>>%
    mutate(sites=strsplit(sites, '\\|')) %>>%
    unnest(sites) %>>%
    bind_rows(population %>>% filter(sites=='')) %>>%
    full_join(history %>>% name_rows, by=c(sites='.rownames')) %>>%
    dplyr::select(-sites) %>>%
    arrange(x, y) %>>% (?.)

source(file.path(dirname(..file..), 'sample.R'))

#########1#########2#########3#########4#########5#########6#########7#########
if (conf[['dimensions']] == 2) {

.heat_colours = c('#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')

plot_early_mutations_2d = function(.data, num_tracing=4) {.data %>>%
    mutate(size=ifelse(size <= num_tracing, size, NA), effect=NULL) %>>%
    dplyr::select(-starts_with('origin_')) %>>%
    filter(!duplicated(.)) %>>%
    group_by(x, y) %>>%
    filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
    mutate(size=as.factor(size)) %>>%
    ggplot(aes(x, y, fill=size))+
    geom_raster(alpha=0.66)+
    scale_fill_hue(na.value='white')+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())
}

plot_mutation_history_2d = function(.data) {.data %>>%
    ggplot(aes(x, y, fill=size))+
    geom_raster(alpha=0.5)+
    scale_fill_gradientn(colours=.heat_colours, na.value='white')+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())
}

.threshold = history$size %>>% nth(4)
.p = plot_early_mutations_2d(unnested, .threshold)
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

.p = plot_mutation_history_2d(unnested)
#.p
ggsave('gradient.png', .p, width=7, height=7)

#########1#########2#########3#########4#########5#########6#########7#########
} else {
#########1#########2#########3#########4#########5#########6#########7#########

plot_early_mutations_3d = function(.z=0, .data) {.data %>>%
    filter(z==.z) %>>%
    bind_rows(data_frame(x=maxabs+seq_len(8), y=maxabs+seq_len(8), z=.z,
                          fitness=-1, marker=as.factor(seq_len(8)))) %>>%
    #(?.) %>>% (?str(.)) %>>%
    ggplot(aes(x, y, fill=marker))+
    geom_raster(alpha=0.66, stat='identity', position='identity', interpolate=FALSE)+
    geom_hline(yintercept=.z)+
    scale_fill_hue(na.value='white')+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())+
    theme(legend.position='right')
}
#plot_early_mutations_3d(20, early3d)

early3d = unnested %>>%
    filter(size <= 8) %>>%
    dplyr::select(-effect, -starts_with('origin_')) %>>%
    filter(!duplicated(.)) %>>%
    mutate(marker=as.factor(size), size=NULL) %>>% (?.)

maxabs = with(population, max(abs(c(x, y, z)))) %>>% (?.)

animation::saveGIF({
    for (i in seq(-maxabs, maxabs)) {
        print(plot_early_mutations_3d(i, early3d))
    }},
    'animation.gif', loop=TRUE, interval=0.15, outdir=getwd())


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
writeWebGL()

#clear3d()
#axes3d()
#with(early3d, spheres3d(early3d,# %>>% filter(sqrt(x^2 + y^2 + z^2)>14),
#                radius=1, col=marker, alpha=ifelse(marker==3, 0.8, 0.06)))

#########1#########2#########3#########4#########5#########6#########7#########
}  # fi 2D/3D
#########1#########2#########3#########4#########5#########6#########7#########
