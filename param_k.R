#!/usr/bin/env Rscript
library(pipeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

library(wtl)
#########1#########2#########3#########4#########5#########6#########7#########

(.argv = commandArgs(trailingOnly=TRUE))
(..file.. = sub('--file=', '', grep('--file=', .argv, value=TRUE)))
.project = dirname(normalizePath(..file..))

indirs = .argv
(indirs = if (length(indirs) == 0) {list.files(pattern='tumopp_k*')} else {indirs})
stopifnot(length(indirs) > 0)
#########1#########2#########3#########4#########5#########6#########7#########

read_conf = function(infiles) {
    data_frame(.path=infiles) %>%
        group_by(.path) %>% do({
            wtl::read.conf(.$.path) %>% as.data.frame()
        }) %>% ungroup()
}

conf = read_conf(file.path(indirs, 'program_options.conf')) %>>%
    dplyr::rename(indir=.path) %>>%
    mutate(indir=dirname(indir)) %>>% (?.)

demography = conf %>>%
    group_by(indir, shape) %>>% do({
        read_tsv(file.path(.$indir, 'population.tsv.gz'))
    }) %>>%
    dplyr::select(birth, death) %>>%
    gather(event, time, birth, death) %>>%
    dplyr::filter(!(time == 0 & event == 'death')) %>>%  # alive
    mutate(event= factor(event, levels=c('death', 'birth'))) %>>%
    arrange(time, event) %>>%
    mutate(dn = ifelse(event == 'birth', 1, -1),
           size = cumsum(dn)) %>>%
    group_by(time) %>>%
    summarise(size=last(size)) %>>% (?.)

.breaks = c(1, 10, 1000, 100000)
.labels = c(1, 10, 1000, expression(infinity))

font_size = 20
theme_set(theme_bw(font_size, 'sans'))
theme_update(panel.grid=element_blank())

.p = demography %>>%
    group_by(shape, time) %>>%
    summarise(size=max(size)) %>>%
    ungroup() %>>%
    dplyr::filter(shape %in% .breaks) %>>%
    ggplot(aes(time, size, group=shape, colour=shape))+
    scale_colour_gradient(low='dodgerblue', high='tomato', trans='log10',
        breaks=.breaks, labels=.labels,
        guide=guide_legend(reverse=TRUE, title=expression(italic(k)),
            override.aes=list(size=3)))+
    geom_line(size=1, alpha=0.7)+
    labs(x='Time', y='Number of tumor cells')+
    theme(legend.position=c(0.02, 0.98), legend.justification=c(0, 1),
          legend.key=element_blank())
#.p

ggsave('param_k.pdf', .p, width=5, height=5)
