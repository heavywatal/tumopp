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
#' @return a grouped data.frame
#' @rdname read
#' @export
read_snapshots = function(indirs='.') {
    read_confs(indirs) %>>%
    dplyr::mutate(path=indirs) %>>%
    purrr::by_row(~{
        .d = readr::read_tsv(file.path(.x$path, 'snapshots.tsv.gz'))
        if (.x$coord == 'hex') .d = trans_coord_hex(.d)
        .d %>>% set_id()
    }, .to='snapshots')
}
