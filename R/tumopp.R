#' Run C++ simulation
#' @param args command line arguments as a string vector
#' @return tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0)) {
    results = cpp_tumopp(args)
    conf = wtl::read_boost_ini(results[1])
    pop = readr::read_tsv(results[2]) %>>% modify_population()
    dplyr::mutate(conf, population=list(pop))
}
