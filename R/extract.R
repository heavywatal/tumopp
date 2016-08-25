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
    dplyr::mutate_(event=~ factor(event, levels=c('death', 'birth'))) %>>%
    dplyr::arrange_(~time, ~event) %>>%
    dplyr::mutate_(dn =~ ifelse(event == 'birth', 1, -1),
           size =~ cumsum(dn)) %>>%
    dplyr::group_by_(~time, add=TRUE) %>>%
    dplyr::summarise_(size=~ dplyr::last(size))
}

#' cell ids that existed at the specified tumor size
#' @param raw_population a data.frame including ancestors
#' @param size an integer
#' @param origin an integer
#' @return an integer vector
#' @rdname extract
#' @export
exclusive_ancestors = function(raw_population, size, origin=1) {
    stopifnot(size >= origin)
    mothers = head(raw_population, size - origin)$id
    seq_len(length(mothers) + size) %>>%
    setdiff(mothers)
}

#' cell ids that existed at the specified tumor size (from snapshots)
#' @param snapshots a data.frame
#' @return an integer vector
#' @rdname extract
#' @export
exclusive_ancestors_ss = function(snapshots, size) {
    dplyr::group_by_(~snapshots, ~time) %>>%
    dplyr::filter_(~ n() == size) %>>%
    (stringr::str_extract(.$genealogy, '\\d+$'))
}

#' split ancestor string and extract one matched
#' @param genealogy list of integer vectors
#' @param ids integers
#' @return factor vector of extracted ancestors
#' @rdname extract
#' @export
extract_ancestor = function(genealogy, ids) {
    purrr::map_int(genealogy, ~ max(.[. %in% ids])) %>>%
    factor(levels=ids)
}

#' Extract survivors and add some columns
#' @param result a row of a nested tibble of conf and list(pop)
#' @param n an integer number of exclusive ancestors
#' @return a data.frame
#' @rdname extract
#' @export
extract_survivors = function(result, n=4) {
    population = result$population[[1]]
    num_founders = sum(lengths(population$genealogy) == 1L)
    anc_ids = exclusive_ancestors(population, n, num_founders)
    population = population %>>%
        dplyr::filter_(~ death == 0) %>>%
        dplyr::mutate_(ancestor=~ extract_ancestor(genealogy, anc_ids)) %>>%
        detect_surface(get_se(result$coord, result$dimensions))
    if (result$coord == 'hex') {
        population = trans_coord_hex(population)
    }
    population
}

#' Add a column of living descendants number
#' @return tibble with $id and $discendants
#' @rdname extract
#' @export
count_descendants = function(raw_population) {
    .counts = raw_population %>>%
        dplyr::filter_(~death == 0) %>>%
        (purrr::flatten_int(.$genealogy)) %>>%
        table() %>>%
        tibble::as_tibble() %>>%
        stats::setNames(c('id', 'descendants')) %>>%
        dplyr::mutate_(id=~ as.integer(id))
    dplyr::left_join(raw_population, .counts, by='id')
}
