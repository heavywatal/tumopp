#' read config files
#' @param indirs a string vector
#' @return a data.frame
#' @rdname read
#' @export
read_conf = function(indirs='.') {
    stats::setNames(,indirs) %>>%
    purrr::map_df(~wtl::read_boost_ini(file.path(., 'program_options.conf')),
                  .id='path')
}


#' read a population
#' @param conf a data.frame
#' @return a data.frame
#' @rdname read
#' @export
read_population = function(conf) {
    stopifnot(nrow(conf) == 1L)
    .types = list(
        beta= readr::col_double(),
        delta= readr::col_double(),
        rho= readr::col_double()
    )
    se = get_se(conf$coord, conf$dimensions)
    .path = file.path(conf$path, 'population.tsv.gz')
    x = readr::read_tsv(.path, col_types=.types)
    x = detect_surface(x, se)
    if (conf$coord == 'hex') {
        x = trans_coord_hex(x)
    }
    x
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

#' get commandline arguments
#' @return a list
#' @rdname read
#' @export
command_args = function() {
    .argv = commandArgs(trailingOnly=FALSE)
    l = list()
    l$file = sub('^--file=', '', grep('^--file=', .argv, value=TRUE))
    l$srcdir = dirname(normalizePath(l$file))
    l$args = grep('^[^-]', .argv[-1], value=TRUE)
    return(l)
}
