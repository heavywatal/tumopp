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

#' get max value for plot limits
#' @param mtrx a data.frame with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @rdname coord
#' @export
#' @examples
#' maxabs(data.frame(x=2, y=-3, z=4))
maxabs = function(mtrx) {
    max(abs(c(mtrx$x, mtrx$y, mtrx$z)))
}

#' 2D transformation into hexagonal grid
#' @return transformed matrix
#' @rdname coord
trans_coord_hex_xy = function(mtrx) {
    dplyr::mutate_(mtrx, y=~ y + x * 0.5) %>>%
    dplyr::mutate_(x=~ x * sqrt(3.0 / 4.0))
}

#' 3D transformation (hexagonal close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_hcc = function(mtrx) {
    trans_coord_hex_xy(mtrx) %>>%
    dplyr::mutate_(x=~ x + ifelse(z %% 2 == 1, sqrt(3) / 3, 0)) %>>%
    dplyr::mutate_(z=~ z * sqrt(2.0 / 3.0))
}

#' 3D transformation (face centered cubic, cubic close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_fcc = function(mtrx) {
    trans_coord_hex_xy(mtrx) %>>%
    dplyr::mutate_(x=~ x + z / sqrt(3.0)) %>>%
    dplyr::mutate_(z=~ z * sqrt(2.0 / 3.0))
}

#' coordinate transformation
#' @return transformed matrix
#' @rdname coord
#' @export
trans_coord_hex = function(mtrx) {
    if ('z' %in% names(mtrx)) {
        trans_coord_fcc(mtrx)
    } else {
        trans_coord_hex_xy(mtrx)
    }
}
