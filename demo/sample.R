#!/usr/bin/env Rscript
library(tidyverse)

#########1#########2#########3#########4#########5#########6#########7#########

sample_glands = function(popultion) {population %>%
    mutate(half=ifelse(x<1, 'left', 'right')) %>%
    group_by(half) %>%
    sample_n(5) %>%
    mutate(serial=seq_len(n()))
}

unnest_ = function(samples) {samples %>%
    mutate(site=strsplit(sites, ':'), sites=NULL) %>%
    unnest(site) %>%
    mutate(site=as.integer(site))
}

binarize = function(unnested) {unnested %>%
    (~ .sites = unique(.$site)) %>%
    group_by(half, serial) %>%
    do(data.frame(site=.sites, exists=.sites %in% .$site)) %>%
    ungroup() %>%
    arrange(site)
}

# summary statistics
##  the number of distinct CNAs (different strings)
##  the total number of alterations
##  Shannon index of binary patterns
##  the number of variegated alterations
##  the number of side-variegated alterations

summarize_alt = function(binary) {binary %>%
    group_by(site) %>%
    summarise(alterations=sum(exists), p=alterations / n()) %>%
    summarise(distinct=n(), alterations=sum(alterations), entropy=sum(-p * log(p) - (1-p) * log(1 - p)))
}

#  lhs rhs
#  --- ---: ERROR
#  --- any: regional
#  --- all: side-specific
#  any any: variegated
#  any all: side-variegated
#  all all: public
summarize_pattern = function(binary) {binary %>%
    group_by(site, half) %>%
    summarise(any=any(exists), all=all(exists)) %>%
    group_by(site) %>%
    summarise(public=all(all),
        side_variegated=(!public) & any(all) & all(any),
        variegated=(!any(all)) & all(any),
        side_specific=(!public) & (!any(any)) & any(all),
        regional=(!any(all)) & (!variegated) & any(any)) %>%
    summarise_at(vars(-site), sum)
}

summarize_glands = function(samples) {
    .binary = samples %>%
        unnest_() %>%
        binarize()
    bind_cols(summarize_alt(.binary), summarize_pattern(.binary))
}

population %>%
    sample_glands() %>%
    summarize_glands()

#########1#########2#########3#########4#########5#########6#########7#########

.repeated = rdply(100, {
    population %>%
        sample_glands() %>%
        summarize_glands()
}, .id='sampling') %>% as_tibble() %>% print()
write_tsv(.repeated, "samples.tsv.gz")

.p = .repeated %>%
    gather(variable, value, -sampling) %>%
    ggplot(aes(value))+
    geom_histogram()+
    facet_wrap(~variable)+
    theme_bw()+
    theme(panel.grid.minor=element_blank())
#.p
ggsave('repeated_sampling.png', .p, width=7, height=7)
