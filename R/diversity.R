#' H / H_max
#' @param species a factor vector
#' @return a number
#' @rdname diversity
#' @export
evenness = function(species) {
    shannon_index(species) / log(length(levels(species)))
}

#' Shannon index
#' @return a number
#' @rdname diversity
#' @export
shannon_index = function(species) {
    freqs = table(species)
    freqs = freqs[freqs > 0] / sum(freqs)
    -sum(freqs * log(freqs))
}

#' Simpson index
#' @return a number
#' @rdname diversity
#' @export
simpson_index = function(species) {
    freqs = table(species)
    freqs = freqs[freqs > 0]
    sum(freqs ^ 2) / (sum(freqs) ^ 2)
}
