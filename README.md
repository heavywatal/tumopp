# tumorr

`tumorr` is R interface of [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Installation

1.  Install [tumopp](https://github.com/heavywatal/tumopp) to `~/local/` or `/usr/local/`

2.  Install [devtools](https://github.com/hadley/devtools) in R:
    `install.packages('devtools')`

3.  Create `~/.R/Makevars`:

    ```
    CXX1XSTD = -std=c++14
    LDFLAGS = -L${HOME}/local/lib -L/usr/local/lib
    ```

    [`CXX_STD = CXX14` in `src/Makevars` will be supported in R 3.4.0]
    (http://gallery.rcpp.org/articles/rcpp-and-c++11-c++14-c++17/)

4. `devtools::install_github('heavywatal/tumorr')`
