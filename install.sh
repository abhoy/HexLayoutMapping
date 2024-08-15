#!/bin/bash

if [[ "$1" == "help" || $# -ne 1 ]]; then
  echo "Requires Command Line Arguments"
  echo "help: This Message"
  echo "cmake: Create build directory, compile and install"
  echo "make: Recompile to reflect the update"
  echo "remake: Compile and install"
fi

if [ "$1" == "cmake" ]; then 
  #Checking for existance of build directory
  if [ ! -d "./build" ]; then 
    mkdir ./build
  fi
  cd build && cmake .. && cmake --build . && cd ..
fi
 
if [ "$1" == "make" ]; then 
    cd build && make && cd ..
fi

if [ "$1" == "remake" ]; then 
    cd build && cmake .. && cmake --build . && cd ..
fi


