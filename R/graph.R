#' cell ids that existed at the specified tumor size
#' @param population a raw data.frame
#' @return a data.frame
#' @rdname graph
#' @export
make_edgelist = function(population) {
    dplyr::transmute(population,
      from=stringr::str_extract(genealogy, '\\d+(?=:\\d+$)'),
      to=stringr::str_extract(genealogy, '\\d+$')
    ) %>>%
    dplyr::filter(!is.na(from)) %>>%
    dplyr::mutate_all(as.integer) %>>%
    dplyr::arrange(to)
}
