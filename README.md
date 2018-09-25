# tumopp

[![Build Status](https://travis-ci.org/heavywatal/tumopp.svg?branch=master)](https://travis-ci.org/heavywatal/tumopp)

Tumor growth simulator in C++.

## Dependencies

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/)
- [clippson](https://github.com/heavywatal/clippson)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)

## Installation

The easiest way is to use [Homebrew](https://brew.sh/)/[Linuxbrew](http://linuxbrew.sh/).
The following command installs tumopp and all the dependencies:
```sh
brew install heavywatal/tap/tumopp
```

Alternatively, you can get the source code from GitHub manually:
```sh
git clone https://github.com/heavywatal/tumopp.git
cd tumopp/
mkdir build
cd build/
DESTINATION=${HOME}/local
cmake -DCMAKE_INSTALL_PREFIX=$DESTINATION ..
make -j2
make install
```

If needed, set `CMAKE_PREFIX_PATH` variable so that CMake can find prerequisite libraries,
e.g., `cmake -DCMAKE_PREFIX_PATH=$(brew --prefix) ..`


## Usage

Run tumopp via [R package](https://github.com/heavywatal/rtumopp/)

Alternatively, it can be executed as a command-line program:
```sh
tumopp -h
tumopp -N20000 -D3 -Chex -k100 -d0.1 -m0.5 -o OUTPUT_DIR
```

- [Online documentation generated with doxygen](https://heavywatal.github.io/tumopp/)
- @ref params


## Reference

Watal M. Iwasaki and Hideki Innan (2017)
"Simulation Framework for Generating Intratumor Heterogeneity Patterns in a Cancer Cell Population"
[*PLOS ONE* 12(9): e0184229](https://doi.org/10.1371/journal.pone.0184229)

[Project page on GitHub](https://github.com/heavywatal/tumopp)
