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

library(wtl)
library(tumorr)
#load_all('../tumorr')
#########1#########2#########3#########4#########5#########6#########7#########

indirs = if (length(.argv) > 0) {.argv} else {list.files(pattern='tumopp_.+')}
stopifnot(length(indirs) > 0)
print(indirs)
#########1#########2#########3#########4#########5#########6#########7#########

conf = tumorr::read_conf(indirs)
altered = tumorr::altered_params(conf)

demography = tumorr::read_population(conf, altered) %>%
    tumorr::extract_demography() %>>% (?.)

font_size = 20
theme_set(theme_bw(font_size, 'sans'))
theme_update(panel.grid=element_blank())

n0 = min(demography$size)
nmax = max(demography$size)
tmax = log2(nmax / n0)

#########1#########2#########3#########4#########5#########6#########7#########
if (altered == 'shape') {

    .breaks = c(1, 10, 1000, 100000)
    .labels = c(1, 10, 1000, expression(infinity))

    .p = demography %>>%
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

#########1#########2#########3#########4#########5#########6#########7#########
} else if (altered == 'packing') {
    .p = demography %>>%
        ggplot(aes(time, size, group=packing, colour=packing))+
        geom_line(size=1, alpha=0.7)+
        labs(x='Time', y='Number of tumor cells')+
        theme(legend.position=c(0.02, 0.98), legend.justification=c(0, 1),
              legend.key=element_blank())

#########1#########2#########3#########4#########5#########6#########7#########
} else {

    exponential_growth = function(time, n0=8, r=1) {
        n0 * exp(r * time)
    }

    surface_growth = function(time, n0=8, r=log(2)) {
        (n0^(1/3) + r * (pi * 4 / 3)^(1/3) * time) ^ 3
    }

    .p = demography %>>%
        #sample_frac(0.3) %>>%
        group_by_(.dots=c(altered, 'time')) %>>%
        summarise(size=max(size)) %>>% ungroup() %>>%
        ggplot(aes(time, size))+
        geom_path(aes_string(group=altered, colour=altered))+
        stat_function(fun=exponential_growth, args=c(n0, log(2)), linetype='dashed', colour='#009900')+
        stat_function(fun=exponential_growth, args=c(n0, log(2 - 0.3)), linetype='dashed', colour='#990000')+
        stat_function(fun=exponential_growth, args=c(n0, 1.0), linetype='dashed', colour='#000099')+
        stat_function(fun=surface_growth, args=c(n0, log(2)), linetype='dashed', colour='#000000')+
        theme_bw()+
        coord_cartesian(ylim=range(demography$size))
}
#########1#########2#########3#########4#########5#########6#########7#########
#.p
ggsave(sprintf('demography_%s.pdf', altered), .p, width=5, height=5)
