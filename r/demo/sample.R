library(tidyverse)
library(wtl)
library(tumopp)

refresh("tumopp/r")

.result = tumopp::read_results('tumopp_20180501_182513') %>% print()

(.result = tumopp(str_split("-N40000 -D2 -Chex -k24 -Lstep", " ")[[1]]))
(.population = .result$population[[1]])
(.extant = .population %>% filter_extant())
(.regions = sample_uniform_regions(.extant, 8L, 100L))

.sampled = .regions %>%
  dplyr::transmute(cluster = id, id = samples) %>%
  tidyr::unnest() %>%
  dplyr::left_join(.extant, by = "id") %>%
  print()

.extant %>%
  plot_lattice2d() +
  geom_point(data = .sampled, aes(x, y)) +
  scale_colour_brewer(palette = "Spectral", guide = FALSE)

(.graph = make_igraph(.population))

# #######1#########2#########3#########4#########5#########6#########7#########

.threshold = 0.01

.detectable = .regions$samples %>% purrr::map(~detectable_mutants(.graph, .x, .threshold))
.combn_biopsy = .detectable %>% combn_sample_ids() %>% print()
.tbl = summarize_capture_rate(.combn_biopsy, .population, .threshold) %>% print()
.tbl %>% plot_capture_rate()

.do = function(threshold, samples = .regions$samples, graph=.graph, population=.population) {
  samples %>%
    purrr::map(~detectable_mutants(graph, .x, threshold)) %>%
    combn_sample_ids() %>%
    summarize_capture_rate(population, threshold)
}
.tbl_capture = c(0.01, 0.03, 0.05) %>% wtl::mcmap_dfr(.do) %>% print()

.tbl_capture %>%
  plot_capture_rate() +
  facet_wrap(~threshold) +
  theme_bw()

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
