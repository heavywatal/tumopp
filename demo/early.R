#!/usr/bin/env Rscript
library(tidyverse)
library(wtl)
library(tumorr)
#########1#########2#########3#########4#########5#########6#########7#########

# tumopp -D2 -k100 -Cmoore -Lconst -O4 -R256 -N256 -w 0 0 -o Cmoore_Lconst
# tumopp -D2 -k100 -Cmoore -Lstep -O4 -R256 -N256 -w 0 0 -o Cmoore_Lstep
# tumopp -D2 -k100 -Cmoore -Llinear -O4 -R256 -N256 -w 0 0 -o Cmoore_Llinear
# tumopp -D2 -k100 -Chex -Lconst -O4 -R256 -N256 -w 0 0 -o Chex_Lconst
# tumopp -D2 -k100 -Chex -Lstep -O4 -R256 -N256 -w 0 0 -o Chex_Lstep
# tumopp -D2 -k100 -Chex -Llinear -O4 -R256 -N256 -w 0 0 -o Chex_Llinear

(.args = wtl::command_args())
indir = .args$args[1]
if (!is.na(indir)) {
    setwd(indir)
}

snapshots = tumorr::read_snapshots()
norigins = snapshots %>% dplyr::filter(time == 0) %>% nrow()
nclades = max(4L, norigins)
founders = head(snapshots, nclades)$id %>% print()
roots = head(snapshots, nclades)$genealogy %>%
    flatten_int() %>% unique() %>% setdiff(founders) %>% print()

.tbl = snapshots %>% dplyr::mutate(
    clade= purrr::map_int(.data$genealogy, ~{setdiff(.x, roots)[1L]}),
    clade= factor(.data$clade, levels=founders)) %>% print()

.lim = tumorr::max_abs_xyz(snapshots)

plot_snapshot = function(data, time) {
    .N = nrow(data)
    tumorr::plot_lattice2d(data, 'clade', alpha=1.0, limit=.lim)+
    scale_color_brewer(palette='Spectral')+
    # labs(title=sprintf('t = %.5f, N =%4d', time, .N))+
    # theme_void()+
    theme(axis.text=element_blank(), axis.ticks=element_blank())+
    theme(legend.position='none')
}

.out = .tbl %>%
    tidyr::nest(-time) %>%
    dplyr::mutate(plt=purrr::pmap(., plot_snapshot)) %>%
    print()

dir.create('png', mode='0755')
purrr::iwalk(.out$plt, ~{
    .outfile = sprintf('png/snapshot_%d.png', .y)
    message(.outfile)
    ggsave(.outfile, .x, width=1, height=1, scale=6, dpi=72)
  })

# magick -loop 1 -delay 8 png/*.png snapshot-noloop.gif
# magick snapshot-noloop.gif -layers Optimize snapshot-noloop-opt.gif

if (FALSE) {
    # tiled-png for CSS sprite
    grob = gridExtra::arrangeGrob(grobs=time_plots$plt, nrow=1)
    ggsave('earlysteps.png', grob, width=54, height=1, scale=3, limitsize=FALSE)
}
