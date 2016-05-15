#' calculate radius from volume
#' @param volume a numeric vector
#' @return radius a numeric vector
#' @rdname sphere
#' @export
#' @examples
#' sphere_radius(42)
sphere_radius = function(volume) {
    (volume * 3 / (4 * pi)) ^ (1/3)
}

#' calculate volume from radius
#' @param radius a numeric vector
#' @return volume a numeric vector
#' @rdname sphere
#' @export
#' @examples
#' sphere_volume(42)
sphere_volume = function(radius) {
    4 * pi * (radius^3) / 3
}
