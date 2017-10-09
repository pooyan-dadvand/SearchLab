#!/bin/zsh

# Kratos OmniCompile
clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt 2> /dev/null
rm *.cmake 2> /dev/null
rm -rf CMakeFiles\ 2> /dev/null

export CC=clang
export CXX=clang++

cmake ..                                                                                                            \
-DCMAKE_C_COMPILER=${CC}                                                                                            \
-DCMAKE_CXX_COMPILER=${CXX}                                                                                         \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -fopenmp=libomp -Wsign-compare -std=c++11 -Wno-overloaded-virtual -mavx2 -g"  \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -fopenmp=libomp -mavx2 -g"                                                        \
-DKRATOS_PATH="/home/roigcarlo/Kratos/"                                                                             \
-DCMAKE_BUILD_TYPE=Release
