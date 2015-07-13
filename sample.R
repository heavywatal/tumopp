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

.sample = population %>>%
    mutate(half=ifelse(x<1, 'left', 'right')) %>>%
    group_by(half) %>>%
    sample_n(5) %>>%
    mutate(serial=seq_len(n())) %>>% (?.)


# summary statistics
## the number of distinct CNAs (different strings)

.unnested = .sample %>>%
    mutate(site=strsplit(sites, '\\|'), sites=NULL) %>>%
    unnest(site) %>>%
    mutate(site=as.integer(site)) %>>%
    (?.)

.sites = unique(.unnested$site)
length(.sites)

## Shannon index of binary patterns

.binary = .unnested %>>%
    group_by(half, serial) %>>%
    do(data.frame(site=.sites, exists=.sites %in% .$site)) %>>%
    ungroup %>>%
    arrange(site) %>>%
    (?.)

.binary %>>%
    group_by(site) %>>%
    summarise(p= sum(exists) / n()) %>>%
    summarise(entropy=sum(-p * log(p) - (1-p) * log(1 - p)))

## the total number of alterations

nrow(.unnested)

## the number of variegated alterations
## the number of side-variegated alterations

.binary %>>% as.data.frame

.binary %>>%
    group_by(site, half) %>>%
    summarise(any=any(exists), all=all(exists))

#  lhs rhs
#  --- ---: ERROR
#  --- any: regional
#  --- all: side-specific
#  any any: variegated
#  any all: side-variegated
#  all all: public

