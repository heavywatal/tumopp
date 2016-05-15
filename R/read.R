#' read a config file
#' @param filename a string
#' @return a data.frame
.read_conf = function(filename) {
    readr::read_delim(filename, '=', col_names=c('key', 'val'), comment='#') %>>%
    dplyr::summarise_each(dplyr::funs(paste0(., collapse='\t'))) %>>%
    {paste(.$key, .$val, sep='\n')} %>>%
    readr::read_tsv()
}

#' read config files
#' @param indir a string vector
#' @return a data.frame
#' @rdname read
#' @export
read_conf = function(indir='.') {
    dplyr::data_frame(path=indir) %>>%
    dplyr::group_by_(~path) %>>% dplyr::do({
        .read_conf(file.path(.$path, 'program_options.conf'))
    }) %>>%
    dplyr::ungroup()
}

#' read population
#' @param conf a data.frame
#' @param params a string vector; columns to be preserved
#' @return a grouped data.frame
#' @rdname read
#' @export
read_population = function(conf, params=NULL) {
    dplyr::group_by_(conf, .dots=c('path', params)) %>>%
    dplyr::do({
        x = readr::read_tsv(file.path(.$path, 'population.tsv.gz')) %>>%
            dplyr::mutate_(surface=~ genealogy %in% extract_surface(.))
        if (.$coord == 'hex') {
            x = trans_coord_hex(x)
        }
        x
    })
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
