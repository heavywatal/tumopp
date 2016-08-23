#' cell ids that existed at the specified tumor size
#' @param population a raw data.frame
#' @return a data.frame
#' @rdname graph
#' @export
make_edgelist = function(population) {
    dplyr::filter_(population, ~divs > 0L) %>>%
    dplyr::transmute_(
      from= ~map2_int(genealogy, divs, `[[`),
      to= ~id
    ) %>>%
    dplyr::mutate_all(as.integer) %>>%
    dplyr::arrange_(~to)
}
