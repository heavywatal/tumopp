library(tidyverse)
library(wtl)
library(tumopp)

refresh("tumopp/r")

(.result = tumopp(str_split("-N40000 -D3 -Chex -k24 -Lstep", " ")[[1]]))
(.population = .result$population[[1]])
(.extant = .population %>% filter_extant())
(.graph = make_igraph(.population))

# #######1#########2#########3#########4#########5#########6#########7#########

.regions = sample_random_regions(.extant, nsam=8L, ncell=120L) %>% print()

.threshold = 0.05

.combn_biopsy = .regions$samples %>%
  purrr::map(~detectable_mutants(.graph, .x, .threshold)) %>%
  combn_sample_ids() %>%
  print()

summarize_capture_rate(.combn_biopsy, .population, .threshold) %>%
  print() %>%
  plot_capture_rate()

c(0.01, 0.03, 0.05) %>%
  purrr::map_dfr(~summarize_capture_rate(.combn_biopsy, .population, .x)) %>%
  tidyr::unnest() %>%
  print() %>%
  plot_capture_rate() +
  facet_wrap(~threshold)

detectable_mutants_all(.population, .threshold)

.plot_allelefreq_biopsy = function(population, sample_ancestors, threshold) {
  .df = population %>%
    dplyr::filter(.data$allelefreq >= threshold) %>%
    dplyr::transmute(
      .data$id,
      .data$allelefreq,
      captured = .data$id %in% sample_ancestors
    )
  ggplot(.df, aes(allelefreq, group=captured))+
    geom_histogram(aes(fill=captured), binwidth=0.02)+
    theme_bw()
}
.plot_allelefreq_biopsy(.population, .combn_biopsy$nodes[[3]], .threshold)

# #######1#########2#########3#########4

Rprof()
.distances = within_between_samples(.graph, .regions) %>% print()
Rprof(NULL)
summaryRprof()

.xmax = max(.distances$euclidean)
.ymax = max(max(.distances$fst), 0.6)
.distances %>%
  ggplot(aes(euclidean, fst)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ x + 0) +
  coord_cartesian(xlim = c(0, .xmax), ylim = c(0, .ymax)) +
  theme_wtl()
