#!/usr/bin/Rscript
library(pipeR)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

.argv = commandArgs(trailingOnly=TRUE)
indir = .argv[1]
indir = ifelse(is.na(indir), '.', indir)

history = read_tsv(file.path(indir, 'mutation_history.tsv.gz')) %>>% (?.)

population = read_tsv(file.path(indir, 'population.tsv.gz'), col_types='iicd') %>>% (?.)

.data = population %>>%
    mutate(sites=strsplit(sites, '\\|')) %>>%
    unnest(sites) %>>%
    bind_rows(population %>>% filter(sites=='')) %>>%
    full_join(history %>>% name_rows, by=c(sites='.rownames')) %>>%
    dplyr::select(-sites) %>>%
    arrange(x, y) %>>% (?.)

.threshold = history$size %>>% nth(4)

.p = .data %>>%
    mutate(size=ifelse(size <= .threshold, size, NA), effect=NULL) %>>%
    dplyr::select(-starts_with('origin_')) %>>%
    filter(!duplicated(.)) %>>%
    group_by(x, y) %>>%
    filter(!(n() > 1 & is.na(size))) %>>% ungroup %>>%
    mutate(size=as.factor(size)) %>>%
    ggplot(aes(x, y, fill=size))+
    geom_tile(alpha=0.66)+
    scale_fill_hue(na.value='white')+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())
#.p
ggsave('early_mutations.png', .p, width=7, height=7)

.heat_colours = c('#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')

.p = .data %>>%
    ggplot(aes(x, y, fill=size))+
    geom_tile(alpha=0.5)+
    scale_fill_gradientn(colours=.heat_colours, na.value='white')+
    theme(panel.background=element_rect(fill='grey80'))+
    theme(panel.grid=element_blank())+
    theme(axis.title=element_blank())
#.p
ggsave('gradient.png', .p, width=7, height=7)
