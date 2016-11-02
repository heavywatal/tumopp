#' Run C++ simulation
#' @param args command line arguments as a string vector or list of strings
#' @param nsam number of samples to measure genetic and physical distance
#' @return nested tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0L), nsam=0L) {
    if (is.list(args)) {
        purrr::map_df(args, tumopp, .id='args')
    } else {
        message(paste(args, collapse=' '))
        result = cpp_tumopp(args, nsam)
        .out = wtl::read_boost_ini(result[1L]) %>>%
            dplyr::mutate(population=list(readr::read_tsv(result[2L]))) %>>%
            modify_population()
        if (nsam > 0L) {
            .dist = readr::read_tsv(result[3L])
            .out = .out %>>% dplyr::mutate(distances=list(.dist))
        }
        .out
    }
}

#' Make argment list for tumopp()
#' @param alt named list of altered arguments
#' @param const unnamed vector of constant arguments
#' @param nreps number of repeats
#' @return list of character vectors
#' @rdname tumopp
#' @export
make_args = function(alt, const=NULL, nreps=1L) {
    altered = purrr::invoke(expand.grid, alt, stringsAsFactors=FALSE)
    prefix = format(Sys.time(), '%Y%m%d_%H%M_')
    paste0(prefix, seq_len(nreps)) %>>%
        purrr::map_df(~dplyr::mutate(altered, o=.x)) %>>%
        purrr::by_row(~{
            alt_par = paste0('-', names(.), .)
            c(const, alt_par)
        }, .labels=FALSE) %>>%
        `[[`('.out') %>>%
        stats::setNames(purrr::map_chr(., paste, collapse=' '))
}
