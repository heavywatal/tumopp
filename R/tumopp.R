#' Run C++ simulation
#' @param args command line arguments as a string
#' @return tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0)) {
    results = cpp_tumopp(args)
    output = wtl::read_boost_ini(results[1])
    pop = readr::read_tsv(results[2])
    output = dplyr::mutate(output, population=list(pop))
    output
}
