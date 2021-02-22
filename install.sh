#!/bin/bash

set -e

# Install bamtools
cd bamtools
mkdir -p build
cd build
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PWD/../ ..
make
make install

# Install LRez
cd ../../
make
