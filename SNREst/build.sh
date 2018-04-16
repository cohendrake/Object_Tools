#! /bin/bash

build_type="Debug"

[ $# -ge 1 ] && build_type=$1

cmake -DCMAKE_BUILD_TYPE=${build_type} .

make
