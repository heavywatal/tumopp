# tumopp

[![Build Status](https://travis-ci.org/heavywatal/tumopp.svg?branch=master)](https://travis-ci.org/heavywatal/tumopp)

*tumopp* is a fast and flexible tumor growth simulator.
The core library is written in C++.
It can be installed and run via [R package](https://github.com/heavywatal/rtumopp/).


## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.4.0)

The following libraries are automatically installed or optional:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)


## Installation

See [tumopp R package](https://github.com/heavywatal/rtumopp/).


## Alternative installation for command-line execution

The easiest way is to use [Homebrew](https://brew.sh/)/[Linuxbrew](http://linuxbrew.sh/).
The following command installs tumopp with its dependencies:
```sh
brew install heavywatal/tap/tumopp
```

You can manually install the latest version from source code to an arbitrary `DESTINATION`:
```sh
git clone https://github.com/heavywatal/tumopp.git
cd tumopp/
mkdir build
cd build/
DESTINATION=${HOME}/local
cmake -DCMAKE_INSTALL_PREFIX=$DESTINATION ..
make -j2
make install
PATH=${DESTINATION}/bin:$PATH
```

Example:
```sh
tumopp -h
tumopp -N20000 -D3 -Chex -k100 -d0.1 -m0.5 -o OUTPUT_DIR
```


## API Document

- [Online documentation generated with doxygen](https://heavywatal.github.io/tumopp/)
- @ref params


## Reference

Watal M. Iwasaki and Hideki Innan (2017)
"Simulation Framework for Generating Intratumor Heterogeneity Patterns in a Cancer Cell Population"
[*PLOS ONE* 12(9): e0184229](https://doi.org/10.1371/journal.pone.0184229)

[Project page on GitHub](https://github.com/heavywatal/tumopp)
