#' Run C++ simulation
#' @param args command line arguments as a string vector or list of strings
#' @param npair number of samples to measure genetic and physical distance
#' @return nested tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0L), npair=0L) {
    if (is.list(args)) {
        purrr::map_dfr(args, tumopp, .id='args')
    } else {
        message(paste(args, collapse=' '))
        nsam_nrep = c('0', '0')
        result = cpp_tumopp(c(nsam_nrep, args), npair)
        .out = wtl::read_boost_ini(result[1L]) %>%
            dplyr::mutate(population=list(readr::read_tsv(result[2L]))) %>%
            purrr::pmap_dfr(modify_population)
        .out$drivers = list(readr::read_tsv(result[3L]))
        if (npair > 0L) {
            .dist = readr::read_tsv(result[4L])
            .out = .out %>% dplyr::mutate(distances=list(.dist))
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
    paste0(prefix, seq_len(nreps)) %>%
    purrr::map_dfr(~dplyr::mutate(altered, o=.x)) %>%
    purrr::pmap(function(...) {
        .params = c(...)
        .names = names(.params)
        .template = ifelse(nchar(.names) > 1L, '--%s=%s', '-%s%s')
        c(const, sprintf(.template, .names, .params))
    }) %>%
    stats::setNames(purrr::map_chr(., paste, collapse=' '))
}
