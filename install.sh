#!/bin/bash

echo
echo "Building..."
echo

mkdir build
cd build
cmake ..
make

# Only run make install if the --no-install argument wasn't provided
if [[ "$1" != "--no-install" ]]; then
  echo
  echo "Installing..."
  echo
  make install
fi
