#' Convert data.frame to cimg
#' @param mtrx a data.frame with (x, y, z) columns
#' @return cimg
#' @rdname morphology
#' @export
df2img = function(mtrx) {
    vars = c('x', 'y', 'z')
    mtrx = mtrx[vars]
    .grid = dplyr::summarise_each_(mtrx, dplyr::funs(min, max), vars=vars) %>>%
        {expand.grid(x= seq(.$x_min, .$x_max),
                     y= seq(.$y_min, .$y_max),
                     z= seq(.$z_min, .$z_max))} %>>% tibble::as_tibble()
    joined = dplyr::left_join(.grid, mtrx %>>% dplyr::mutate(v=1), by=vars)
    joined = tidyr::replace_na(joined, list(v=0))
    arr = reshape2::acast(joined, x ~ y ~ z, `[`, 1, value.var='v', fill=0)
    dim(arr) = c(dim(arr), 1)
    imager::as.cimg(arr)
}

#' Convert cimg to data.frame
#' @param img a cimg object
#' @return a data.frame with (x, y, z) columns
#' @rdname morphology
#' @export
img2df = function(img) {
    as.array(img) %>>%
    {dim(.) = head(dim(.), 3); .} %>>%
    reshape2::melt(c('x', 'y', 'z')) %>>%
    tibble::as_tibble()
}

#' Structuring element
#' @param coord a string
#' @param dimensions an integer
#' @return cimg structuring element
#' @rdname morphology
#' @export
get_se = function(coord=c('moore', 'neumann', 'hex'), dimensions=3) {
     coord = match.arg(coord)
     df = tibble::tibble(x=c(-1, 0, 1)) %>>%
         dplyr::mutate_(y=~x, z=~x) %>>%
         tidyr::expand_(dots=c('x', 'y', 'z'))
     if (coord == 'neumann') {
         df = dplyr::filter_(df, ~abs(x) + abs(y) + abs(z) < 2)
     } else if (coord == 'hex') {
         idx = trans_coord_hex(df) %>>% {.$x^2 + .$y^2 + .$z^2 < 1.1}
         df = dplyr::filter(df, idx)
     }
     if (dimensions < 3) {
         df = dplyr::filter_(df, ~z == 0)
     }
     df2img(df)
}

#' Filter img by surface
#' @param se a cimg object
#' @return cimg
#' @rdname morphology
#' @export
filter_surface = function(img, se) {
    img - imager::erode(img, se)
}

#' extract cells on suface
#' @return a logical vector
#' @rdname morphology
#' @export
detect_surface = function(mtrx, se) {
    axes = c('x', 'y', 'z')
    mins = dplyr::summarise_each_(mtrx, dplyr::funs(min), vars=axes)
    img = df2img(mtrx) %>>% filter_surface(se)
    product = img2df(img) %>>% dplyr::transmute_(
        x=~ x + mins$x - 1,
        y=~ y + mins$y - 1,
        z=~ z + mins$z - 1,
        surface=~ value > 0)
    dplyr::left_join(mtrx, product, by=axes)
}
