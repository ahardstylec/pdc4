#!/bin/sh
cd build
rm -R *
cmake ..
make
cp ../extrinsics.yml ./extrinsics.yml
cp ../intrinsics.yml ./intrinsics.yml
