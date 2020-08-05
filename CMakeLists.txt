cmake_minimum_required (VERSION 3.2)
project (TripleSphere CXX)

# Debug Build
set(CMAKE_BUILD_TYPE DEBUG)

include_directories ("include")

# Project TestMC: validate with Monte Carlo result
add_executable (TestMC  tests/TestMC.cpp  include/TripleSphere.h)