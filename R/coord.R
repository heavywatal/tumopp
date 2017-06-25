#' get max value for plot limits
#' @param .tbl data.frame with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @rdname coord
#' @export
#' @examples
#' max_abs_xyz(data.frame(x=2, y=-3, z=4))
max_abs_xyz = function(.tbl) {
    max(abs(.tbl[c('x', 'y', 'z')]))
}

#' 2D transformation into hexagonal grid
#' @return transformed matrix
#' @rdname coord
trans_coord_hex_xy = function(.tbl) {
    dplyr::mutate(.tbl, y= .data$y + .data$x * 0.5) %>%
    dplyr::mutate(x= .data$x * sqrt(3.0 / 4.0))
}

#' 3D transformation (hexagonal close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_hcc = function(.tbl) {
    trans_coord_hex_xy(.tbl) %>%
    dplyr::mutate(x= .data$x + ifelse(.data$z %% 2 == 1, sqrt(3) / 3, 0)) %>%
    dplyr::mutate(z= .data$z * sqrt(2.0 / 3.0))
}

#' 3D transformation (face centered cubic, cubic close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_fcc = function(.tbl) {
    trans_coord_hex_xy(.tbl) %>%
    dplyr::mutate(x= .data$x + .data$z / sqrt(3.0)) %>%
    dplyr::mutate(z= .data$z * sqrt(2.0 / 3.0))
}

#' 2D/3D transformation into hexagonal grid
#' @return transformed matrix
#' @rdname coord
#' @export
trans_coord_hex = function(.tbl) {
    trans_coord_fcc(.tbl)
}

#' Translate cells to minimize the range of symmetric axes
#' @return modified data.frame
#' @rdname coord
#' @export
centering = function(.tbl) {
    .offset = dplyr::filter(.tbl, .data$surface) %>%
        dplyr::summarise_at(c('x', 'y', 'z'), dplyr::funs((min(.) + max(.)) * 0.5))
    dplyr::mutate(.tbl,
        x= .data$x - .offset$x,
        y= .data$y - .offset$y,
        z= .data$z - .offset$z)
}

#' Rotate
#' @param theta radian angle
#' @param axis a string
#' @return modified data.frame
#' @rdname coord
#' @export
rotate = function(.tbl, theta, axis=c('z', 'x', 'y')) {
    axis = match.arg(axis)
    .x = .tbl[['x']]
    .y = .tbl[['y']]
    .z = .tbl[['z']]
    .sin = sin(theta)
    .cos = cos(theta)
    if (axis == 'z') {
        dplyr::mutate(.tbl,
            x= .x * .cos - .y * .sin,
            y= .x * .sin + .y * .cos)
    } else if (axis == 'x') {
        dplyr::mutate(.tbl,
            y= .y * .cos - .z * .sin,
            z= .y * .sin + .z * .cos)
    } else {  # y
        dplyr::mutate(.tbl,
            x= .x * .cos + .z * .sin,
            z= - .x * .sin + .z * .cos)
    }
}

#' Calculate distance from a specified cell
#' @param point named vector or tibble
#' @return numeric vector
#' @rdname coord
#' @export
dist_euclidean = function(.tbl, point=c(x=0, y=0, z=0)) {
    sqrt((.tbl[['x']] - point[['x']])^2 + (.tbl[['y']] - point[['y']])^2 + (.tbl[['z']] - point[['z']])^2)
}
