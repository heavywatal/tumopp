#' tumorr: R interface to tumopp
#' @docType package
#' @name tumorr
NULL



#' split ancestor string and extract an integer
#' @param string ancestors column
#' @param n consider first n ancestors
#' @return an integer vector
#' @export
#' @examples
#' first_ancestors(c('1:2:9', '2:4:6'), 4)
first_ancestors = function(string, n) {
    sapply(stringr::str_split(string, ':', n + 1), function(x) {
        x = suppressWarnings(as.integer(x))
        max(x[x <= n], na.rm=TRUE)
    })
}
