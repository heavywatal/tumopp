#' Modify population table
#' @param ... a row of nested tibble
#' @return tibble
#' @rdname population
modify_population = function(..., num_clades=4L) {
    result = list(...)
    population = result$population %>% set_id()
    extant = filter_extant(population)
    strelem = get_se(result$coord, result$dimensions)
    col_surface = detect_surface(extant, strelem) %>%
        dplyr::select(.data$id, .data$surface)
    if (result$coord == 'hex') {
        population = trans_coord_hex(population)
    }
    result$max_phi = c(hex=12, moore=27, neumann=6)[result$coord]
    result$population = population %>%
        dplyr::mutate(r= dist_euclidean(.), phi= .data$phi / result$max_phi) %>%
        set_clades(num_clades) %>%
        dplyr::left_join(count_descendants(extant), by='id') %>%
        dplyr::left_join(col_surface, by='id') %>%
        list()
    result
}

#' Filter extant cells
#' @param population tibble including ancestors
#' @return tibble
#' @rdname population
#' @export
filter_extant = function(population) {
    dplyr::filter(population, .data$death == 0)
}

#' Filter connected cells
#' @param ids integer vector of extant/sample cells
#' @return tibble
#' @rdname population
#' @export
filter_connected = function(population, ids) {
    ids = dplyr::filter(population, .data$id %in% ids)$genealogy %>%
        purrr::flatten_int() %>%
        unique()
    dplyr::filter(population, .data$id %in% ids)
}

#' Extract age and id from genealogy column
#' @return tibble
#' @rdname population
set_id = function(population) {
    dplyr::mutate(population,
        genealogy= stringr::str_split(.data$genealogy, ':') %>% purrr::map(as.integer),
        age= lengths(.data$genealogy) - 1L,
        id= purrr::map2_int(.data$genealogy, .data$age + 1L, `[`))
}

#' Add a column of ancestor ids by which branches are classified
#' @param num_clades integer
#' @return tibble
#' @rdname population
set_clades = function(population, num_clades) {
    origin = sum(population$age == 0L)
    stopifnot(num_clades >= origin)
    num_divisions = num_clades - origin
    roots = utils::head(population$id, num_divisions)
    founders = seq_len(num_divisions + num_clades) %>% setdiff(roots)
    dplyr::mutate(population,
        clade= purrr::map_int(.data$genealogy, ~{setdiff(.x, roots)[1L]}),
        clade= factor(.data$clade, levels=founders))
}

#' Add a column of living descendants number
#' @return tibble with $id and $discendants
#' @rdname population
count_descendants = function(population) {
    if (nrow(population) == 0L) {
        return(tibble::tibble(id=integer(0L), descendants=integer(0L)))
    }
    population$genealogy %>%
        purrr::flatten_int() %>%
        table() %>%
        tibble::as_tibble() %>%
        stats::setNames(c('id', 'descendants')) %>%
        dplyr::mutate(id= as.integer(.data$id))
}

#' Extract demography from raw population data
#' @return tibble
#' @rdname population
#' @export
extract_demography = function(population) {population %>%
    dplyr::select(.data$birth, .data$death) %>%
    tidyr::gather_('event', 'time', c('birth', 'death')) %>%
    dplyr::filter(!(.data$time == 0 & .data$event == 'death')) %>%  # alive
    dplyr::mutate(event= factor(.data$event, levels=c('death', 'birth'))) %>%
    dplyr::arrange(.data$time, .data$event) %>%
    dplyr::mutate(dn= ifelse(.data$event == 'birth', 1, -1)) %>%
    dplyr::group_by(.data$time) %>%
    dplyr::summarise(dn= sum(.data$dn)) %>%
    dplyr::mutate(size= cumsum(.data$dn))
}
