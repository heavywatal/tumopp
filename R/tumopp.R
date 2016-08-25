#' Run C++ simulation
#' @param args command line arguments as a string vector
#' @return tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0)) {
    results = cpp_tumopp(args)
    wtl::read_boost_ini(results[1]) %>>%
    dplyr::mutate(population=list(readr::read_tsv(results[2]))) %>>%
    modify_population()
}
