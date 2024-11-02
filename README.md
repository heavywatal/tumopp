# tumopp

[![Build status](https://github.com/heavywatal/tumopp/workflows/build/badge.svg)](https://github.com/heavywatal/tumopp/actions)

*tumopp* is a fast and flexible tumor growth simulator.
The core library is written in C++.
It can be installed and run via [R package](https://github.com/heavywatal/rtumopp/).


## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++17 compiler (clang++ >= Apple LLVM 12, g++ >= 8)
- [CMake](https://cmake.org/) (>= 3.15.0)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [pcglite](https://github.com/heavywatal/pcglite)
- [zlib](https://zlib.net)


## Installation

See [tumopp R package](https://heavywatal.github.io/rtumopp/).


## Alternative installation for command-line execution

The easiest way is to use [Homebrew](https://brew.sh/).
The following command installs tumopp with its dependencies:
```sh
brew install heavywatal/tap/tumopp
```

You can manually install the latest version from source code to an arbitrary `DESTINATION`:
```sh
git clone https://github.com/heavywatal/tumopp.git
cd tumopp/
DESTINATION=${HOME}/local
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$DESTINATION
cmake --build build -j 2
cmake --install build -j 2
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
