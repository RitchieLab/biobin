#!/bin/bash

########################## preload cmake, sqlite, and boost
#module load cmake/3.17.3
#module load sqlite/3.34.1
#module load boost/1.65.1 

cmake -S . -B build -D CMAKE_BUILD_TYPE=Debug