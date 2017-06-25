library(tidyverse)
library(wtl)
library(tumorr)

refresh('tumorr')

(.result = tumopp(str_split('-D3 -Chex -k24 -Lstep', ' ')[[1]]))
(.population = .result$population[[1]])
(.extant = .population %>% filter_extant())

.n = 100
.o1 = dplyr::sample_n(.extant %>% dplyr::filter(z == 0), 1) %>% print()
.o2 = dplyr::sample_n(.extant %>% dplyr::filter(z == 0), 1) %>% print()
.specimen1 = sample_bulk(.extant, .o1, .n)
.specimen2 = sample_bulk(.extant, .o2, .n)

# T_S
(.ts1 = mean_branch_length(.population, .specimen1$id))
(.ts2 = mean_branch_length(.population, .specimen2$id))
(.ts = mean(c(.ts1, .ts2)))

# T_B
(.tb = mean_branch_length(.population, .specimen1$id, .specimen2$id))

# T_D
(.td = .tb - .ts)

# T_T
(.tt = (.tb + .ts) / 2)

# Fst
fst_HBK(.ts, .tb)
fst_HSM(.ts, .tb)

#########1#########2#########3#########4#########5#########6#########7#########

loadNamespace('cowplot')

.plot_genealogy = function(population, .nodes1, .nodes2) {
    .layout = layout_genealogy(population) %>%
        dplyr::mutate(specimen=
            ifelse(to %in% .nodes1, 1L,
            ifelse(to %in% .nodes2, 2L, NA)) %>% as.factor())
    .specimens = dplyr::filter(.layout, !is.na(specimen))
    plot_genealogy(.layout)+
        ggplot2::geom_point(data= .specimens,
            ggplot2::aes_(x=~ageend, y=~posend, colour=~specimen),
            size=2, alpha=0.8)+
        ggplot2::theme(legend.position='top',
            panel.border=ggplot2::element_blank(), panel.grid=element_blank(),
            axis.title=element_blank(), axis.text=element_blank(),
            axis.ticks=element_blank())
}

.plot_specimens = function(population, o1, o2, size=100L) {
    .extant = filter_extant(population)
    .specimen1 = sample_bulk(.extant, o1, size)
    .specimen2 = sample_bulk(.extant, o2, size)
    .plat = gglattice2d(.extant %>% dplyr::filter(-1 <= z, z < 1), 'z', , 0.5)+
        geom_point(data=.specimen1 %>% dplyr::filter(-1 <= z, z < 1), colour='tomato', alpha=0.7)+
        geom_point(data=.specimen2 %>% dplyr::filter(-1 <= z, z < 1), colour='turquoise', alpha=0.7)
    .nodes1 = as.character(.specimen1$id)
    .nodes2 = as.character(.specimen2$id)
    .pgen = .plot_genealogy(population, .nodes1, .nodes2)
    list(lattice=.plat, genealogy=.pgen)
}
.plts = .plot_specimens(.population, .o1, .o2)
.gg = cowplot::plot_grid(plotlist=.plts)
gprint(.gg)
ggsave('sample_phylogeny_Pifany.png', .gg, width=2, height=1, scale=6)

.n = 40
purrr::map_df(seq_len(.n), ~within_between_samples(.population)) %>%
    dplyr::mutate(fst=fst_HBK(within, between)) %>%
    ggplot(aes(euclidean, fst))+
    geom_point()+
    theme_wtl()


#########1#########2#########3#########4#########5#########6#########7#########
## Draw random tree

.draw_tree = function() {
    .result = tumopp(str_split('-D3 -Chex -Lstep -N60', ' ')[[1]])
    .tbl = .result$population[[1]] %>% layout_genealogy()
    .tree = ggplot2::ggplot(.tbl)+
        ggplot2::geom_segment(ggplot2::aes_(~age, ~pos, xend=~ageend, yend=~posend), alpha=1, size=1)+
        ggplot2::geom_point(data=dplyr::filter(.tbl, .data$extant),  ggplot2::aes_(x=~ageend, y=~posend),
            size=2, colour='dodgerblue', alpha=1)+
        ggplot2::geom_point(data=dplyr::filter(.tbl, !.data$extant), ggplot2::aes_(x=~ageend, y=~posend),
            size=1.2, colour='black', alpha=1)+
        wtl::theme_wtl()+
        ggplot2::theme(
            axis.title=ggplot2::element_blank(),
            axis.text=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank(),
            panel.grid.major=ggplot2::element_blank())
    .outfile = sprintf('tree_%s.png', as.character(.result$seed))
    message(.outfile)
    ggsave(.outfile, .tree, width=1, height=1, scale=4)
}
.draw_tree()
