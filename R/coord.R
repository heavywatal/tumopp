#' get max value for plot limits
#' @param mtrx a data.frame with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @rdname coord
#' @export
#' @examples
#' max_abs_xyz(data.frame(x=2, y=-3, z=4))
max_abs_xyz = function(mtrx) {
    max(abs(mtrx[c('x', 'y', 'z')]))
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

#' Rotate
#' @param theta radian angle
#' @param axis a string
#' @return modified data.frame
#' @rdname coord
#' @export
rotate = function(mtrx, theta, axis=c('z', 'x', 'y')) {
    axis = match.arg(axis)
    .x = mtrx[['x']]
    .y = mtrx[['y']]
    .z = mtrx[['z']]
    .sin = sin(theta)
    .cos = cos(theta)
    if (axis == 'z') {
        dplyr::mutate(mtrx,
            x= .x * .cos - .y * .sin,
            y= .x * .sin + .y * .cos)
    } else if (axis == 'x') {
        dplyr::mutate(mtrx,
            y= .y * .cos - .z * .sin,
            z= .y * .sin + .z * .cos)
    } else {  # y
        dplyr::mutate(mtrx,
            x= .x * .cos + .z * .sin,
            z= - .x * .sin + .z * .cos)
    }
}
