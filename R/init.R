#' tumorr: R interface to tumopp
#' @docType package
#' @name tumorr
NULL


#' read config files
#' @param indir a string vector
#' @return a data.frame
#' @rdname read
#' @export
read_conf = function(indir) {
    dplyr::data_frame(path=indir) %>%
    dplyr::group_by(path) %>% dplyr::do({
        wtl::read.conf(file.path(.$path, 'program_options.conf'))
    }) %>%
    dplyr::ungroup()
}

#' read population
#' @param conf a data.frame
#' @param params a string vector; columns to be preserved
#' @return a grouped data.frame
#' @rdname read
#' @export
read_population = function(conf, params) {
    dplyr::group_by_(conf, .dots=c('path', params)) %>>%
    dplyr::do({
        read_tsv(file.path(.$path, 'population.tsv.gz'))
    })
}

#' extract demography from population data
#' @param grouped_df a grouped_df
#' @return a grouped data.frame
#' @rdname extract
#' @export
extract_demography = function(grouped_df) {grouped_df %>%
    dplyr::select(birth, death) %>%
    tidyr::gather(event, time, birth, death) %>%
    dplyr::filter(!(time == 0 & event == 'death')) %>%  # alive
    dplyr::mutate(event= factor(event, levels=c('death', 'birth'))) %>%
    dplyr::arrange(time, event) %>%
    dplyr::mutate(dn = ifelse(event == 'birth', 1, -1),
           size = cumsum(dn)) %>%
    dplyr::group_by(time, add=TRUE) %>%
    dplyr::summarise(size=last(size))
}

#' extract param names from conf
#' @param conf a data.frame
#' @return a string vector
#' @rdname extract
#' @export
altered_params = function(conf) {
    dplyr::select(conf, -path, -out_dir, -seed) %>%
    dplyr::summarise_each(funs(length(unique(.)))) %>%
    unlist() %>>% (.[. > 1]) %>% names()
}

#' split ancestor string and extract an integer
#' @param string ancestors column
#' @param n consider first n ancestors
#' @return an integer vector
#' @rdname extract
#' @export
#' @examples
#' first_ancestors(c('1:2:9', '2:4:6'), 4)
first_ancestors = function(string, n) {
    sapply(stringr::str_split(string, ':', n + 1), function(x) {
        x = suppressWarnings(as.integer(x))
        max(x[x <= n], na.rm=TRUE)
    })
}
