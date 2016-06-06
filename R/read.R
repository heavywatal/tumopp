#' read a config file
#' @param filename a string
#' @rdname read
#' @return a data.frame
.read_conf = function(filename) {
    readr::read_delim(filename, '=', col_names=c('key', 'val'), comment='#') %>>%
    dplyr::summarise_each(dplyr::funs(paste0(., collapse='\t'))) %>>%
    {paste(.$key, .$val, sep='\n')} %>>%
    readr::read_tsv()
}

#' read config files
#' @param indirs a string vector
#' @return a data.frame
#' @rdname read
#' @export
read_conf = function(indirs='.') {
    stats::setNames(,indirs) %>>%
    purrr::map_df(~.read_conf(file.path(., 'program_options.conf')),
                  .id='path')
}


#' read a population
#' @return a data.frame
#' @rdname read
.read_population = function(conf) {
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

#' read populations
#' @param conf a data.frame
#' @param params a string vector; columns to be preserved
#' @return a nested data.frame
#' @rdname read
#' @export
read_population = function(conf, params=NULL) {
    x = purrr::by_row(conf, .read_population, .to='population')
    if (is.null(params)) {
        x
    } else {
        x[c('path', params, 'population')]
    }
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
