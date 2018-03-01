#' Modify population table
#' @param population tibble
#' @param coord string
#' @param dimensions integer
#' @param ... ignored
#' @return tibble
#' @rdname population
modify_population = function(population, coord, dimensions, ..., num_clades=4L) {
  extant = filter_extant(population)
  strelem = get_se(coord, dimensions)
  col_surface = detect_surface(extant, strelem) %>%
    dplyr::select(.data$id, .data$surface)
  if (coord == "hex") {
    population = trans_coord_hex(population)
  }
  max_phi = c(hex = 12L, moore = 27L, neumann = 6L)[coord]
  population %>%
    set_graph_property() %>%
    dplyr::mutate(r = dist_euclidean(.), phi = .data$phi / max_phi) %>%
    dplyr::left_join(col_surface, by = "id")
}

#' Filter extant cells
#' @return tibble
#' @rdname population
#' @export
filter_extant = function(population) {
  dplyr::filter(population, .data$death == 0)
}

#' Add age and clade column
#' @return tibble
#' @rdname population
set_graph_property = function(population) {
  .graph = make_igraph(population)
  .nodes = as.character(population$id)
  .out = dplyr::mutate(
    population,
    age = igraph::distances(.graph, .nodes, "1", mode = "in", algorithm = "unweighted")[, 1L],
    age = as.integer(.data$age),
    descendants = igraph::neighborhood.size(.graph, order = 1073741824L, nodes = .nodes, mode = "out"),
    descendants = as.integer(.data$descendants)
  )
  founders = list_clade_founders(.out, 4L)
  clade_data = founders %>%
    as.character() %>%
    rlang::set_names() %>%
    purrr::map_dfr(~{
      tibble::tibble(id = igraph::subcomponent(.graph, .x, mode = "out")$name)
    }, .id = "clade") %>%
    dplyr::mutate(
      id = as.integer(.data$id),
      clade = factor(.data$clade, levels = as.character(founders))
    )
  dplyr::left_join(.out, clade_data, by = "id")
}

#' @param num_clades integer
#' @return ids of clade founders
#' @rdname population
list_clade_founders = function(population, num_clades) {
  origin = sum(population$age == 0L)
  stopifnot(num_clades >= origin)
  num_divisions = num_clades - origin
  roots = utils::head(population$id, num_divisions)
  seq_len(num_divisions + num_clades) %>% setdiff(roots)
}

#' Extract demography from raw population data
#' @return tibble
#' @rdname population
#' @export
extract_demography = function(population) {
  population %>%
    dplyr::select(.data$birth, .data$death) %>%
    tidyr::gather("event", "time", c("birth", "death")) %>%
    dplyr::filter(!(.data$time == 0 & .data$event == "death")) %>% # alive
    dplyr::mutate(event = factor(.data$event, levels = c("death", "birth"))) %>%
    dplyr::arrange(.data$time, .data$event) %>%
    dplyr::mutate(dn = ifelse(.data$event == "birth", 1, -1)) %>%
    dplyr::group_by(.data$time) %>%
    dplyr::summarise(dn = sum(.data$dn)) %>%
    dplyr::mutate(size = cumsum(.data$dn))
}
