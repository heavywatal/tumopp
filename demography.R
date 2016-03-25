#!/usr/bin/env Rscript
(.argv = commandArgs(trailingOnly=FALSE))
(..file.. = sub('^--file=', '', grep('^--file=', .argv, value=TRUE)))
(.argv = grep('^[^-]', .argv[-1], value=TRUE))
.project = dirname(normalizePath(..file..))

library(pipeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

indir = .argv[1]
if (!is.na(indir)) {
    setwd(indir)
}

population = read_tsv('population.tsv.gz') %>>% (?.)
demography = population %>>%
    dplyr::select(birth, death) %>>%
    gather(event, time, birth, death) %>>%
    dplyr::filter(!(time == 0 & event == 'death')) %>>%  # alive
    mutate(event= factor(event, levels=c('death', 'birth'))) %>>%
    arrange(time, event) %>>%
    mutate(dn = ifelse(event == 'birth', 1, -1),
           size = cumsum(dn)) %>>%
    group_by(time) %>>%
    summarise(size=last(size)) %>>% (?.)

n0 = min(demography$size)
nmax = max(demography$size)
tmax = log2(nmax / n0)

exponential_growth = function(time, n0=8, r=1) {
    n0 * exp(r * time)
}

surface_growth = function(time, n0=8, r=log(2)) {
    (n0^(1/3) + r * (pi * 4 / 3)^(1/3) * time) ^ 3
}

.p = demography %>>%
    #sample_frac(0.3) %>>%
    group_by(time) %>>%
    summarise(size=max(size)) %>>%
    ggplot(aes(time, size))+
    geom_line()+
    stat_function(fun=exponential_growth, args=c(n0, log(2)), colour='#FF0000')+
    stat_function(fun=exponential_growth, args=c(n0, log(2 - 0.2)), colour='#990000')+
    stat_function(fun=exponential_growth, args=c(n0, 1.0), colour='#00FF00')+
    stat_function(fun=surface_growth, args=c(n0, log(2)), colour='blue')+
    theme_bw()+
    coord_cartesian(ylim=range(demography$size))
#.p

ggsave('demography.png', .p, width=7, height=7)
