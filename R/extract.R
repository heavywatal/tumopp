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
#' @param population a raw data.frame
#' @param size an integer
#' @param origin an integer
#' @return an integer vector
#' @rdname extract
#' @export
exclusive_ancestors = function(population, size, origin=1) {
    stopifnot(size >= origin)
    mothers = head(population, size - origin)$genealogy %>>%
        stringr::str_extract('\\d+$')
    seq_len(length(mothers) + size) %>>% setdiff(mothers)
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
#' @param genealogy a string vector of colon-separated ancestors
#' @param ids a string vector of exclusive ancestor ids
#' @return a factor vector of extracted ancestors
#' @rdname extract
#' @export
#' @examples
#' extract_ancestor(c('1:2:9', '1:2:4:6'), c('2', '5', '6'))
extract_ancestor = function(genealogy, ids) {
    stringr::str_split(genealogy, ':', length(ids) + 1) %>>%
    purrr::map_chr(~ max(.[. %in% ids])) %>>%
    factor(levels=ids)
}

#' E ancestor column
#' @param raw_population a data.frame including ancestors
#' @param n an integer number of exclusive ancestors
#' @return a data.frame
#' @rdname extract
#' @export
colorcode_survivors = function(raw_population, n=4) {
    num_founders = length(grep(':', raw_population$genealogy, invert=TRUE))
    anc_ids = exclusive_ancestors(raw_population, n, num_founders)
    raw_population %>>%
        dplyr::filter_(~ death == 0) %>>%
        dplyr::mutate_(ancestor=~ extract_ancestor(genealogy, anc_ids))
}
