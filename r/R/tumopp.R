#' R interface of tumopp, tumor growth simulator in C++
#' @docType package
#' @name tumopp
#' @useDynLib tumopp
#' @importFrom magrittr %>%
#' @importFrom rlang .data
NULL
# to suppress NOTE
utils::globalVariables(c(".", "n"))

#' Run C++ simulation
#' @param args command line arguments as a string vector or list of strings
#' @param npair number of samples to measure genetic and physical distance
#' @param nsam number of samples for ms-like output
#' @return nested tibble
#' @rdname tumopp
#' @export
tumopp = function(args=character(0L), npair=0L, nsam=0L) {
  if (is.list(args)) {
    purrr::map_dfr(args, tumopp, .id = "args")
  } else {
    if (length(args) == 1L) {
      args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
    }
    message(paste(args, collapse = " "))
    nrep = as.integer(nsam > 0L)
    result = cpp_tumopp(c(nsam, nrep, args), npair = npair, nsam = nsam)
    if (length(result) == 0L) return(invisible(NULL))
    .out = wtl::read_boost_ini(result["config"]) %>%
      dplyr::mutate(population = list(readr::read_tsv(result["specimens"]))) %>%
      dplyr::mutate(population = purrr::pmap(., modify_population)) %>%
      dplyr::mutate(drivers = list(readr::read_tsv(result["drivers"])))
    if (npair > 0L) {
      .dist = readr::read_tsv(result["distances"])
      .out = .out %>% dplyr::mutate(distances = list(.dist))
    }
    if (nsam > 0L) {
      .out$ms = wtl::parse_ms(strsplit(result["ms"], "\n")[[1]])
    }
    .out
  }
}

#' Make argument list for tumopp()
#' @param alt named list of altered arguments
#' @param const unnamed vector of constant arguments
#' @param nreps number of repeats
#' @return list of character vectors
#' @rdname tumopp
#' @export
make_args = function(alt, const=NULL, nreps=1L) {
  altered = purrr::invoke(expand.grid, alt, stringsAsFactors = FALSE)
  prefix = format(Sys.time(), "%Y%m%d_%H%M_")
  paste0(prefix, seq_len(nreps)) %>%
    purrr::map_dfr(~dplyr::mutate(altered, o = .x)) %>%
    purrr::pmap(function(...) {
      .params = c(...)
      .names = names(.params)
      .template = ifelse(nchar(.names) > 1L, "--%s=%s", "-%s%s")
      c(const, sprintf(.template, .names, .params))
    }) %>%
    stats::setNames(purrr::map_chr(., paste, collapse = " "))
}
