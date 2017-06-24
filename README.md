# tumopp
Tumor growth simulator in C++.

## Dependencies

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler
- [CMake](https://cmake.org/)
- [Boost C++ Libraries](http://www.boost.org/)
- [sfmt-class](//github.com/heavywatal/sfmt-class)
- [cxxwtl](//github.com/heavywatal/cxxwtl)

## Installation

```sh
git clone https://github.com/heavywatal/tumopp.git
mkdir build-tumopp
cd build-tumopp/
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/local ../tumopp
make
make install
```

## Usage

Run tumopp in R via [`tumorr` package](//github.com/heavywatal/tumorr)

Alternatively, it can be executed as a command-line program:
```sh
tumopp -h
tumopp -N20000 -D3 -Chex -k100 -d0.1 -m0.5 -w -o OUTPUT
```

## Reference

Watal M. Iwasaki and Hideki Innan (2017)
"Simulation Framework for Generating Intratumor Heterogeneity Patterns in a Cancer Cell Population"
[*bioRxive preprint*](https://doi.org/10.1101/109801)
