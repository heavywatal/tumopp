#!/usr/bin/env Rscript
library(pipeR)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)

.argv_all = commandArgs(trailingOnly=FALSE)
..file.. = sub('--file=', '', grep('--file=', .argv_all, value=TRUE))
.project = dirname(normalizePath(..file..))

(.argv = commandArgs(trailingOnly=TRUE))
indir = .argv[1]
indir = ifelse(is.na(indir), '.', indir)
setwd(indir)

population = read_tsv('population.tsv.gz') %>>% (?.)
demography = population %>>%
    dplyr::select(birth, death) %>>%
    gather(event, time, birth, death) %>>%
    dplyr::filter(event=='birth' | time > 0) %>>%
    arrange(time) %>>%
    mutate(dn = ifelse(event == 'birth', 1, -1),
           size = cumsum(dn)) %>>% (?.)

n0 = sum(demography$time == 0)
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

ggsave('demography.png', .p, width=7, height=7)
