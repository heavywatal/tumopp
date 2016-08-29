#' Summary statistics of genealogy
#' @param survivors tibble
#' @return tibble
#' @rdname summarize
#' @export
genetic_stats = function(survivors) {
    dplyr::summarize_(survivors,
       mean_age= ~mean(age), sd_age= ~sd(age),
       max_age= ~max(age), min_age= ~min(age),
       evenness= ~wtl::evenness(clade)
    )
}

#' Summary statistics of morphology
#' @param coord string
#' @return tibble
#' @rdname summarize
#' @export
morphological_stats = function(survivors, coord=c('moore', 'hex', 'neumann')) {
    .phi_max = c(hex=12, moore=27, neumann=6)[match.arg(coord)]
    survivors %>>%
        dplyr::filter_(~surface) %>>%
        dplyr::mutate_(r= ~sqrt(x^2 + y^2 + z^2), phi= ~phi / .phi_max) %>>%
        dplyr::summarise_(phi_mean= ~mean(phi), phi_sd= ~stats::sd(phi),
                          r_mean= ~mean(r), r_sd= ~stats::sd(r)) %>>%
        dplyr::mutate(surface=sum(survivors$surface) / nrow(survivors))
}
