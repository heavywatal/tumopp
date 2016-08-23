#' read config files
#' @param indirs a string vector
#' @return tibble
#' @rdname read
#' @export
read_conf = function(indirs='.') {
    stats::setNames(,indirs) %>>%
    purrr::map_df(~wtl::read_boost_ini(file.path(., 'program_options.conf')),
                  .id='path')
}

#' read a population
#' @param conf tibble
#' @return tibble
#' @rdname read
#' @export
read_population = function(conf) {
    stopifnot(nrow(conf) == 1L)
    .cols = readr::cols(
        beta= readr::col_double(),
        delta= readr::col_double(),
        rho= readr::col_double())
    file.path(conf$path, 'population.tsv.gz') %>>%
        readr::read_tsv(col_types=.cols) %>>%
        modify_population()
}

#' Modify population table
#' @param raw_population tibble including ancestors
#' @return tibble
#' @rdname read
#' @export
modify_population = function(raw_population) {
    population = raw_population %>>%
        dplyr::mutate_(
          genealogy= ~stringr::str_split(genealogy, ':') %>>% purrr::map(as.integer),
          divs= ~lengths(genealogy) - 1L,
          id= ~purrr::map2_int(genealogy, divs + 1L, `[`)
        ) %>>%
        count_descendants()
    population
}

#' read snapshots
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
