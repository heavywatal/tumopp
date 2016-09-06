#' Run C++ simulation
#' @param args command line arguments as a string vector or list of strings
#' @return nested tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0)) {
    if (is.list(args)) {
        purrr::map_df(args, ~{
            message(paste(.x, collapse=' '))
            tumopp(.x)
        }, .id='args')
    } else {
        results = cpp_tumopp(args)
        wtl::read_boost_ini(results[1]) %>>%
        dplyr::mutate(population=list(readr::read_tsv(results[2]))) %>>%
        modify_population()
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
