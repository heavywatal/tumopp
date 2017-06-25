library(tidyverse)
library(wtl)
library(tumorr)

.o = tibble(x=0, y=0, z=0)
.n = 20
.x = seq(-.n, .n)
.extant = tidyr::crossing(x=.x, y=.x, z=.x, coord='rect') %>%
    dplyr::bind_rows(trans_coord_hex(.) %>% dplyr::mutate(coord='hex')) %>% print()

.tbl = .extant %>%
    dplyr::mutate(dist=dist_euclidean(.extant, .o)) %>%
    dplyr::filter(dist <= .n * sqrt(0.5)) %>%
    dplyr::group_by(coord) %>%
    dplyr::arrange(dist) %>%
    dplyr::mutate(n=seq_along(dist)) %>% print()

ggplot(.tbl, aes(dist, n))+
stat_function(fun=volume_function('rect', 3), colour='dodgerblue', size=3, alpha=0.6)+
stat_function(fun=volume_function('hex', 3), colour='tomato', size=3, alpha=0.6)+
geom_path(aes(linetype=coord))+
wtl::theme_wtl()+
theme(legend.position='top')

ggplot(.tbl, aes(n, dist))+
stat_function(fun=radius_function('rect', 3), colour='dodgerblue', size=3, alpha=0.6)+
stat_function(fun=radius_function('hex', 3), colour='tomato', size=3, alpha=0.6)+
geom_path(aes(linetype=coord))+
wtl::theme_wtl()+
theme(legend.position='top')
