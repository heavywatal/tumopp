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
#' @param mtrx data.framw with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @export
#' @examples
#' maxabs(population)
maxabs = function(mtrx) {
    max(abs(c(mtrx$x, mtrx$y, mtrx$z)))
}

#' 2D transformation into hexagonal grid
#' @param mtrx 2D matrix
#' @return transformed matrix
#' @rdname hex
#' @export
trans_coord_hex = function(mtrx) {mtrx %>>%
    mutate(y= y + x * 0.5) %>>%
    mutate(x= x * sqrt(3.0 / 4.0))
}

#' 3D transformation (hexagonal close packed)
#' @param mtrx 3D matrix
#' @return transformed matrix
#' @rdname hex
#' @export
#
trans_coord_hcc = function(mtrx) {
    trans_coord_hex(mtrx) %>>%
    mutate(x= x + ifelse(z %% 2 == 1, sqrt(3) / 3, 0)) %>>%
    mutate(z= z * sqrt(2.0 / 3.0))
}

#' 3D transformation (face centered cubic, cubic close packed)
#' @param mtrx 3D matrix
#' @return transformed matrix
#' @rdname hex
#' @export
trans_coord_fcc = function(mtrx) {
    trans_coord_hex(mtrx) %>>%
    mutate(x= x + z / sqrt(3.0)) %>>%
    mutate(z= z * sqrt(2.0 / 3.0))
}
