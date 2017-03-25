# tumopp
Tumor growth simulator in C++.

## Dependencies

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler
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

## Reference

Watal M. Iwasaki and Hideki Innan (2017)
"Simulation Framework for Generating Intratumor Heterogeneity Patterns in a Cancer Cell Population"
[*bioRxive preprint*](https://doi.org/10.1101/109801)
