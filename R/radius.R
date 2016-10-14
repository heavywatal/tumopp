#' Return volume function of specified coord and dims
#' @param coord hex or not
#' @param dimensions integer
#' @return function with radius argument
#' @rdname radius
#' @export
volume_function = function(coord='hex', dimensions=3L) {
    switch(dimensions, NULL,
        switch(coord, hex=volume_hex2d, volume_rect2d),
        switch(coord, hex=volume_hex3d, volume_rect3d))
}

#' Return radius function of specified coord and dims
#' @return function with volume argument
#' @rdname radius
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
    4.0 * pi * (radius^3) / 3.0
}
radius_rect2d = function(volume) {
    sqrt(volume / pi)
}
radius_rect3d = function(volume) {
    ((volume / pi) * 0.75) ^ (1.0 / 3.0)
}
volume_hex2d = function(radius) volume_rect2d(radius) * 2.0 / sqrt(3.0)
volume_hex3d = function(radius) volume_rect3d(radius) * sqrt(2.0)
radius_hex2d = function(volume) radius_rect2d(volume * sqrt(0.75))
radius_hex3d = function(volume) radius_rect3d(volume * sqrt(0.5))
