#!/usr/bin/Rscript
library(pipeR)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(animation)

.argv = commandArgs(trailingOnly=TRUE)
indir = .argv[1]
indir = ifelse(is.na(indir), '.', indir)

conf = wtl::read.conf(file.path(indir, 'program_options.conf')) %>>% (?.)
history = read_tsv(file.path(indir, 'mutation_history.tsv.gz')) %>>% (?.)
#population = read_tsv(file.path(indir, 'population.tsv.gz'), col_types='iicd') %>>% (?.)
population = read_tsv(file.path(indir, 'population.tsv.gz')) %>>% (?.)

unnested = population %>>%
    mutate(sites=strsplit(sites, '\\|')) %>>%
    unnest(sites) %>>%
    bind_rows(population %>>% filter(sites=='')) %>>%
    full_join(history %>>% name_rows, by=c(sites='.rownames')) %>>%
    dplyr::select(-sites) %>>%
    arrange(x, y) %>>% (?.)

#########1#########2#########3#########4#########5#########6#########7#########

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

#########1#########2#########3#########4#########5#########6#########7#########

plot_early_mutations_3d = function(.z=20, .data) {.data %>>%
    filter(z==.z) %>>%
    bind_rows(expand.grid(x=maxabs+2, y=maxabs+2, z=.z,
                          fitness=-1, size=as.factor(seq_len(8)))) %>>%
    ggplot(aes(x, y, fill=size))+
    geom_raster(alpha=0.66)+
    geom_hline(yintercept=.z)+
    scale_fill_hue(na.value='white')+
    coord_cartesian(x=c(-maxabs, maxabs), y=c(-maxabs, maxabs))+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())+
    theme(legend.position='right')
}
#plot_early_mutations_3d(15)

#########1#########2#########3#########4#########5#########6#########7#########

if (conf[['dimensions']] == 2) {

.threshold = history$size %>>% nth(4)
.p = plot_early_mutations_2d(unnested, .threshold)
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

.p = plot_mutation_history_2d(unnested)
#.p
ggsave('gradient.png', .p, width=7, height=7)

} else {

early3d = unnested %>>%
    mutate(size=ifelse(size <= 8, size, NA), effect=NULL) %>>%
    dplyr::select(-starts_with('origin_')) %>>%
    filter(!duplicated(.)) %>>%
    group_by(x, y, z) %>>%
    filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
    mutate(size=as.factor(size)) %>>% (?.)

maxabs = with(population, max(abs(c(x, y, z)))) %>>% (?.)

animation::saveGIF({
    for (i in seq(-maxabs, maxabs)) {
        print(plot_early_mutations_3d(i, early3d))
    }},
    'animation.gif', loop=TRUE, interval=0.15, outdir=getwd())

}

