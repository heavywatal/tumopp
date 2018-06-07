library(tidyverse)
library(wtl)
library(tumopp)

refresh("tumopp/r")

add_cluster = function(extant, regions) {
  regions %>%
    dplyr::transmute(cluster = .data$id, id = .data$samples) %>%
    tidyr::unnest() %>%
    {dplyr::left_join(extant, ., by = "id")}
}

(.result = tumopp(str_split("-N40000 -D2 -Chex -k24 -Lconst", " ")[[1]]))
(.population = .result$population[[1]])

(.extant = .population %>% filter_extant())
(.regions = sample_uniform_regions(.extant, 8L, 100L))
(.sampled = add_cluster(.extant, .regions))
(.graph = make_igraph(.population))

.sampled %>%
  plot_lattice2d(size = 0.3) +
  geom_point(data = function(x) {dplyr::filter(x, !is.na(cluster))}, aes(x, y), size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Spectral", guide = FALSE) +
  theme(axis.title = element_blank())

# #######1#########2#########3#########4#########5#########6#########7#########

.threshold = 0.01
.detectable = .regions$samples %>% purrr::map(~detectable_mutants(.graph, .x, .threshold))
.combn_biopsy = .detectable %>% combn_sample_ids() %>% print()
.tbl = summarize_capture_rate(.combn_biopsy, .population, .threshold) %>% print()
.tbl %>% plot_capture_rate()

.combn_capture_rate = function(population, regions, graph, thr_range = c(0.01, 0.03, 0.05), ...) {
  wtl::mcmap_dfr(thr_range, function(thr) {
    regions$samples %>%
      purrr::map(~detectable_mutants(graph, .x, thr)) %>%
      combn_sample_ids() %>%
      summarize_capture_rate(population, thr)
  })
}
.tbl_capture = .combn_capture_rate(.population, .regions, .graph) %>% print()
.tbl_capture %>% plot_capture_rate() + facet_wrap(~threshold)

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

# #######1#########2#########3#########4#########5#########6#########7#########

.tr_P = c(random='Push-1', roulette='Push-2', mindrag='Push-3', minstraight='Push-4')
.tr_L = c(const='Constant-rate', step='Step-function', linear='Linear-function')

setwd("~/working/tumopp/spatial")
result_dirs = fs::dir_ls()
raw_results = tumopp::read_results(result_dirs)

results = raw_results %>%
  select_if(~ n_distinct(.x) > 1L) %>%
  dplyr::select(-outdir, -seed, -drivers) %>%
  print()

.population = results$population[[1]]

df_sampled = results %>%
  dplyr::mutate(
    local = factor(.tr_L[local], levels=.tr_L),
    path = factor(.tr_P[path], levels=.tr_P),
    extant = purrr::map(population, filter_extant),
    regions = purrr::map(extant, sample_uniform_regions, nsam=8L, ncell=100L),
    extant = purrr::map2(extant, regions, add_cluster),
    graph = purrr::map(population, make_igraph)
  ) %>% print()

df_extant = df_sampled %>%
  dplyr::select(local, path, shape, extant) %>%
  dplyr::mutate(extant = purrr::map(extant, ~{
      dplyr::mutate(.x, clade = as.factor(as.integer(clade)))
  })) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  print()

.p_sampled = df_extant %>%
  sample_n(40000) %>%
  plot_lattice2d(size = 0.3) +
  geom_point(data = function(x) {dplyr::filter(x, !is.na(cluster))}, aes(x, y), size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Spectral", guide = FALSE) +
  facet_grid(local + path ~ shape + replicate) +
  wtl::erase(axis.title, axis.text, axis.ticks) +
  theme(panel.spacing = unit(0, 'mm'))
.p_sampled
ggsave('samples.png', .p_sampled, width = 12, height = 12)

df_capture = df %>%
  dplyr::mutate(
    capture_tbl = purrr::pmap(., .combn_capture_rate),
    capture_plt = purrr::map(capture_tbl, ~plot_capture_rate(.x) + facet_wrap(~threshold))
  ) %>% print()

.p_capture = df_capture %>%
  dplyr::transmute(
    local = factor(.tr_L[local], levels=.tr_L),
    path = factor(.tr_P[path], levels=.tr_P),
    shape, capture_tbl) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  plot_capture_rate() +
  facet_grid(local + path ~ shape + replicate + threshold) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(panel.spacing = unit(0, 'mm'))
.p_capture
ggsave("capture_rate.png", .p_capture, width=15, height=12)
