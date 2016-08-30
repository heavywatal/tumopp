#' Read config files
#' @param indirs a string vector
#' @return tibble
#' @rdname read
#' @export
read_confs = function(indirs='.') {
    stats::setNames(,indirs) %>>%
    file.path('program_options.conf') %>>%
    purrr::map_df(wtl::read_boost_ini, .id='path')
}

#' Read populations
#' @return list of tibbles
#' @rdname read
#' @export
read_populations = function(indirs='.') {
    file.path(indirs, 'population.tsv.gz') %>>%
    purrr::map(readr::read_tsv, col_types=readr::cols(
        beta= readr::col_double(),
        delta= readr::col_double(),
        rho= readr::col_double()))
}

#' Read confs and populations as a nested tibble
#' @return nested tibble
#' @rdname read
#' @export
read_results = function(indirs='.') {
    read_confs(indirs) %>>%
    dplyr::mutate(population= read_populations(indirs)) %>>%
    modify_population()
}

#' read snapshots
#' @param conf tibble
#' @return a grouped data.frame
#' @rdname read
#' @export
read_snapshots = function(conf) {
    dplyr::group_by_(conf, ~path) %>>%
    dplyr::do({
        x = readr::read_tsv(file.path(.$path, 'snapshots.tsv.gz'))
        if (.$coord == 'hex') {
            x = trans_coord_hex(x)
        }
        x
    })
}
