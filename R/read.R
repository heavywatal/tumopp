#' read a config file
#' @param filename a string
#' @return a data.frame
#' @rdname read
.read_conf = function(filename) {
    dict = list()
    for (line in scan(filename, what=character(), sep="\n", quiet=TRUE)) {
        if (substr(line, 1, 1) == "#") {next}
        pair = unlist(strsplit(line, " = "))
        dict[pair[1]] = if (grepl('[^[:digit:].-]', pair[2])) {
            gsub("^'|'$", '', gsub('^"|"$', '', pair[2]))
        } else {as.numeric(pair[2])}
    }
    data.frame(dict, stringsAsFactors=FALSE)
}

#' read config files
#' @param indir a string vector
#' @return a data.frame
#' @rdname read
#' @export
read_conf = function(indir='.') {
    dplyr::data_frame(path=indir) %>%
    dplyr::group_by(path) %>% dplyr::do({
        .read_conf(file.path(.$path, 'program_options.conf'))
    }) %>%
    dplyr::ungroup()
}

#' read snapshots
#' @param conf a data.frame
#' @return a grouped data.frame
#' @rdname read
#' @export
read_snapshots = function(conf) {
    dplyr::group_by_(conf, .dots=c('path')) %>%
    dplyr::do({
        x = readr::read_tsv(file.path(.$path, 'snapshots.tsv.gz'))
        if (.$coord == 'hex') {
            x = trans_coord_hex(x)
        }
        x
    })
}

#' read population
#' @param conf a data.frame
#' @param params a string vector; columns to be preserved
#' @return a grouped data.frame
#' @rdname read
#' @export
read_population = function(conf, params=NULL) {
    dplyr::group_by_(conf, .dots=c('path', params)) %>%
    dplyr::do({
        x = readr::read_tsv(file.path(.$path, 'population.tsv.gz'))
        if (.$coord == 'hex') {
            x = trans_coord_hex(x)
        }
        x
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

#' cell ids at the size
#' @param snapshots a data.frame
#' @param size an integer
#' @return an integer vector
#' @rdname extract
#' @export
exclusive_ancestors = function(snapshots, size) {
    dplyr::group_by(snapshots, time) %>>%
    dplyr::filter(n()==size) %>>% (id)
}

#' split ancestor string and extract an integer
#' @param string ancestors column
#' @param anc_ids an integer vector
#' @return an integer vector
#' @rdname extract
#' @export
#' @examples
#' first_ancestors(c('1:2:9', '1:2:4:6'), c(2, 5, 6))
first_ancestors = function(string, anc_ids) {
    suppressWarnings(
    sapply(stringr::str_split(string, ':', max(anc_ids)), function(x) {
        x = as.integer(x)
        max(x[x %in% anc_ids])
    })) %>%
    pmax(min(anc_ids))
}

#' filter first ancesters for coloring
#' @param mtrx a data.frame
#' @param anc_ids an integer vector
#' @return a modified data.frame
#' @rdname extract
#' @export
filter_ancestors = function(mtrx, anc_ids) {
    dplyr::mutate(mtrx,
        ancestors= ifelse(is.na(ancestors), id, ancestors),
        ancestors= first_ancestors(ancestors, anc_ids),
        ancestors= as.factor(ancestors))
}
