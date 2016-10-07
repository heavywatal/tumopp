#' get max value for plot limits
#' @param .data data.frame with (x, y, z) columns
#' @return max(abs(x, y, z))
#' @rdname coord
#' @export
#' @examples
#' max_abs_xyz(data.frame(x=2, y=-3, z=4))
max_abs_xyz = function(.data) {
    max(abs(.data[c('x', 'y', 'z')]))
}

#' 2D transformation into hexagonal grid
#' @return transformed matrix
#' @rdname coord
trans_coord_hex_xy = function(.data) {
    dplyr::mutate_(.data, y=~ y + x * 0.5) %>>%
    dplyr::mutate_(x=~ x * sqrt(3.0 / 4.0))
}

#' 3D transformation (hexagonal close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_hcc = function(.data) {
    trans_coord_hex_xy(.data) %>>%
    dplyr::mutate_(x=~ x + ifelse(z %% 2 == 1, sqrt(3) / 3, 0)) %>>%
    dplyr::mutate_(z=~ z * sqrt(2.0 / 3.0))
}

#' 3D transformation (face centered cubic, cubic close packed)
#' @return transformed matrix
#' @rdname coord
trans_coord_fcc = function(.data) {
    trans_coord_hex_xy(.data) %>>%
    dplyr::mutate_(x=~ x + z / sqrt(3.0)) %>>%
    dplyr::mutate_(z=~ z * sqrt(2.0 / 3.0))
}

#' 2D/3D transformation into hexagonal grid
#' @return transformed matrix
#' @rdname coord
#' @export
trans_coord_hex = function(.data) {
    trans_coord_fcc(.data)
}

#' Rotate
#' @param theta radian angle
#' @param axis a string
#' @return modified data.frame
#' @rdname coord
#' @export
rotate = function(.data, theta, axis=c('z', 'x', 'y')) {
    axis = match.arg(axis)
    .x = .data[['x']]
    .y = .data[['y']]
    .z = .data[['z']]
    .sin = sin(theta)
    .cos = cos(theta)
    if (axis == 'z') {
        dplyr::mutate(.data,
            x= .x * .cos - .y * .sin,
            y= .x * .sin + .y * .cos)
    } else if (axis == 'x') {
        dplyr::mutate(.data,
            y= .y * .cos - .z * .sin,
            z= .y * .sin + .z * .cos)
    } else {  # y
        dplyr::mutate(.data,
            x= .x * .cos + .z * .sin,
            z= - .x * .sin + .z * .cos)
    }
}

#' Calculate distance from a specified cell
#' @param center single row in .data
#' @return numeric vector
#' @rdname coord
#' @export
dist_euclidean = function(.data, center) {
    with(.data, sqrt((x - center$x)^2 + (y - center$y)^2 + (z - center$z)^2))
}

#' Return volume function of specified coord and dims
#' @param coord hex or not
#' @param dimensions integer
#' @return function with radius argument
#' @rdname coord
#' @export
volume_function = function(coord='hex', dimensions=3L) {
    switch(dimensions, NULL,
        switch(coord, hex=volume_hex2d, volume_rect2d),
        switch(coord, hex=volume_hex3d, volume_rect3d))
}

#' Return radius function of specified coord and dims
#' @return function with volume argument
#' @rdname coord
#' @export
radius_function = function(coord='hex', dimensions=3L) {
    switch(dimensions, NULL,
        switch(coord, hex=radius_hex2d, radius_rect2d),
        switch(coord, hex=radius_hex3d, radius_rect3d))
}

volume_rect2d = function(radius) {
    pi * (radius^2)
}
volume_rect3d = function(radius) {
    4 * pi * (radius^3) / 3
}
radius_rect2d = function(area) {
    sqrt(area / pi)
}
radius_rect3d = function(volume) {
    ((volume / pi) * 0.75) ^ (1.0 / 3.0)
}
volume_hex2d = function(radius) volume_rect2d(radius) * 2 / sqrt(3)
volume_hex3d = function(radius) volume_rect3d(radius) * sqrt(2)
radius_hex2d = function(volume) radius_rect2d(volume) * sqrt(0.75)
radius_hex3d = function(volume) radius_rect3d(volume) * sqrt(0.5)
