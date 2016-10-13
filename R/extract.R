#' extract param names from conf
#' @param conf a data.frame
#' @return a string vector
#' @rdname extract
#' @export
altered_params = function(conf) {
    dplyr::select_(conf, ~-path, ~-out_dir, ~-seed) %>>%
    dplyr::summarise_each(dplyr::funs(length(unique(.)))) %>>%
    unlist() %>>% (.[. > 1L]) %>>% names()
}

#' Extract demography from raw population data
#' @return tibble
#' @rdname extract
#' @export
extract_demography = function(population) {population %>>%
    dplyr::select_(~birth, ~death) %>>%
    tidyr::gather_('event', 'time', c('birth', 'death')) %>>%
    dplyr::filter_(~!(time == 0 & event == 'death')) %>>%  # alive
    dplyr::mutate_(event= ~factor(event, levels=c('death', 'birth'))) %>>%
    dplyr::arrange_(~time, ~event) %>>%
    dplyr::mutate_(dn= ~ifelse(event == 'birth', 1, -1)) %>>%
    dplyr::group_by_(~time) %>>%
    dplyr::summarise_(dn= ~sum(dn)) %>>%
    dplyr::mutate_(size= ~cumsum(dn))
}

#' Filter extant cells
#' @param population tibble including ancestors
#' @return tibble
#' @rdname extract
#' @export
filter_extant = function(population) {
    dplyr::filter_(population, ~death == 0)
}

#' Filter connected cells
#' @param ids integer vector of extant/sample cells
#' @return tibble
#' @rdname extract
#' @export
filter_connected = function(population, ids) {
    ids = dplyr::filter_(population, ~ id %in% ids)$genealogy %>>%
        purrr::flatten_int() %>>%
        unique()
    dplyr::filter_(population, ~ id %in% ids)
}

#' Extract age and id from genealogy column
#' @return tibble
#' @rdname extract
set_id = function(population) {
    dplyr::mutate_(population,
        genealogy= ~stringr::str_split(genealogy, ':') %>>% purrr::map(as.integer),
        age= ~lengths(genealogy) - 1L,
        id= ~purrr::map2_int(genealogy, age + 1L, `[`))
}

#' Modify population table
#' @param result a row of nested tibble
#' @return tibble
#' @rdname extract
modify_population = function(result, num_clades=4L) {
    population = result$population[[1L]] %>>% set_id()
    extant = filter_extant(population)
    strelem = get_se(result$coord, result$dimensions)
    col_surface = detect_surface(extant, strelem) %>>%
        dplyr::select_(~id, ~surface)
    if (result$coord == 'hex') {
        population = trans_coord_hex(population)
    }
    result$population[[1L]] = population %>>%
        set_clades(num_clades) %>>%
        dplyr::left_join(count_descendants(extant), by='id') %>>%
        dplyr::left_join(col_surface, by='id')
    result
}

#' Add a column of ancestor ids by which branches are classified
#' @param num_clades integer
#' @return tibble
#' @rdname extract
set_clades = function(population, num_clades) {
    origin = sum(population$age == 0L)
    stopifnot(num_clades >= origin)
    num_divisions = num_clades - origin
    roots = head(population$id, num_divisions)
    founders = seq_len(num_divisions + num_clades) %>>% setdiff(roots)
    dplyr::mutate_(population,
        clade= ~purrr::map_int(genealogy, ~{setdiff(.x, roots)[1L]}),
        clade= ~factor(clade, levels=founders))
}

#' Add a column of living descendants number
#' @return tibble with $id and $discendants
#' @rdname extract
count_descendants = function(population) {
    if (nrow(population) == 0L) {
        return(tibble::tibble(id=integer(0L), descendants=integer(0L)))
    }
    population$genealogy %>>%
        purrr::flatten_int() %>>%
        table() %>>%
        tibble::as_tibble() %>>%
        stats::setNames(c('id', 'descendants')) %>>%
        dplyr::mutate_(id= ~as.integer(id))
}
