# tumorr

`tumorr` is R interface of [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Installation

1.  Install [tumopp](https://github.com/heavywatal/tumopp) to `~/local/` or `/usr/local/`

2.  Install [devtools](https://github.com/hadley/devtools) in R:
    `install.packages('devtools')`

3.  Create `~/.R/Makevars` to specify the location of tumopp and boost libraries:
    ```
    CPPFLAGS = -I${HOME}/local/include -L/usr/local/include
    LDFLAGS = -L${HOME}/local/lib -L/usr/local/lib
    ```

4. `devtools::install_github('heavywatal/tumorr')`
