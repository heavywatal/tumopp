# tumopp

Tumor growth simulator in C++.

## Dependencies

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/)
- [Boost C++ Libraries](http://www.boost.org/) (>= 1.64.0)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)

## Installation

Use `BOOST_ROOT` environment variable so that CMake can find your boost library.

```sh
git clone https://github.com/heavywatal/tumopp.git
mkdir build-tumopp
cd build-tumopp/
export BOOST_ROOT=${HOME}/local
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/local ../tumopp
make
make install
```

## Usage

Run tumopp in R via [`tumorr` package](https://github.com/heavywatal/tumorr)

Alternatively, it can be executed as a command-line program:
```sh
tumopp -h
tumopp -N20000 -D3 -Chex -k100 -d0.1 -m0.5 -w -o OUTPUT
```

- [Online documentation generated with doxygen](https://heavywatal.github.io/tumopp/)
- @ref params


## Reference

Watal M. Iwasaki and Hideki Innan (2017)
"Simulation Framework for Generating Intratumor Heterogeneity Patterns in a Cancer Cell Population"
[*PLOS ONE* 12(9): e0184229](https://doi.org/10.1371/journal.pone.0184229)

[Project page on GitHub](https://github.com/heavywatal/tumopp)
