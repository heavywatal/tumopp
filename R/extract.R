#' extract param names from conf
#' @param conf a data.frame
#' @return a string vector
#' @rdname extract
#' @export
altered_params = function(conf) {
    dplyr::select_(conf, ~-path, ~-out_dir, ~-seed) %>>%
    dplyr::summarise_each(dplyr::funs(length(unique(.)))) %>>%
    unlist() %>>% (.[. > 1]) %>>% names()
}

#' extract demography from population data
#' @param grouped_df a grouped_df
#' @return a grouped data.frame
#' @rdname extract
#' @export
extract_demography = function(grouped_df) {grouped_df %>>%
    dplyr::select_(~birth, ~death) %>>%
    tidyr::gather_(~event, ~time, ~birth, ~death) %>>%
    dplyr::filter_(~!(time == 0 & event == 'death')) %>>%  # alive
    dplyr::mutate_(event= ~factor(event, levels=c('death', 'birth'))) %>>%
    dplyr::arrange_(~time, ~event) %>>%
    dplyr::mutate_(dn= ~ifelse(event == 'birth', 1, -1),
                   size= ~cumsum(dn)) %>>%
    dplyr::group_by_(~time, add=TRUE) %>>%
    dplyr::summarise_(size= ~dplyr::last(size))
}

#' Filter extant cells
#' @return tibble
#' @rdname extract
#' @export
filter_extant = function(raw_population) {
    dplyr::filter_(raw_population, ~death == 0)
}

#' Modify population table
#' @param result a row of nested tibble
#' @return tibble
#' @rdname extract
modify_population = function(result, num_clades=4L) {
    population = result$population[[1]] %>>%
        dplyr::mutate_(
          genealogy= ~stringr::str_split(genealogy, ':') %>>% purrr::map(as.integer),
          age= ~lengths(genealogy) - 1L,
          id= ~purrr::map2_int(genealogy, age + 1L, `[`)
        ) %>>%
        set_clades(num_clades) %>>%
        count_descendants()
    extant = population %>>%
        filter_extant() %>>%
        detect_surface(get_se(result$coord, result$dimensions)) %>>%
        dplyr::select_(~id, ~surface)
    population = dplyr::left_join(population, extant, by='id')
    if (result$coord == 'hex') {
        population = trans_coord_hex(population)
    }
    result$population[[1]] = population
    result
}

#' Add a column of ancestor ids by which branches are classified
#' @param raw_population tibble including ancestors
#' @param num_clades integer
#' @return tibble
#' @rdname extract
set_clades = function(raw_population, num_clades) {
    origin = sum(raw_population$age == 0L)
    stopifnot(num_clades >= origin)
    num_divisions = num_clades - origin
    roots = head(raw_population$id, num_divisions)
    founders = seq_len(num_divisions + num_clades) %>>% setdiff(roots)
    dplyr::mutate_(raw_population,
      clade= ~purrr::map_int(genealogy, ~{setdiff(.x, roots)[1L]}),
      clade= ~factor(clade, levels=founders)
    )
}

#' Add a column of living descendants number
#' @return tibble with $id and $discendants
#' @rdname extract
count_descendants = function(raw_population) {
    .counts = raw_population %>>%
        filter_extant() %>>%
        (purrr::flatten_int(.$genealogy)) %>>%
        table() %>>%
        tibble::as_tibble() %>>%
        stats::setNames(c('id', 'descendants')) %>>%
        dplyr::mutate_(id= ~as.integer(id))
    dplyr::left_join(raw_population, .counts, by='id')
}
