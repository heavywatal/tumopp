# tumorr

`tumorr` is R interface of [tumopp](https://github.com/heavywatal/tumopp),
a tumor growth simulator in C++.

## Installation

```r
devtools::install_github('heavywatal/tumorr')
```

Before that, you may need to add the following line to your `~/.R/Makevars`:

```
CXX1XSTD = -std=c++14
```
