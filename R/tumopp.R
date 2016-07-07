#' Run C++ simulation
#' @param args command line arguments as a string
#' @return tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0)) {
    cpp_tumopp(args) %>>% readr::read_tsv()
}
