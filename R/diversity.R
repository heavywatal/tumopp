#' H / H_max
#' @param species a factor vector
#' @return a number
#' @rdname diversity
#' @export
evenness = function(species) {
    freqs = table(species)
    wtl::shannon_index(freqs) / log(length(freqs))
}
