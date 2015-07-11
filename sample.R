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

population %>>%
    mutate(half=ifelse(x<1, 'left', 'right')) %>>%
    group_by(half) %>>%
    sample_n(5)
