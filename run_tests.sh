#!/bin/bash

# rm -rf build
mkdir build
cd build
cmake ..
make -j1

MPIRUN="mpirun -np"
PARTIES=3
dirlist=`ls test_*.out`

for test in ${dirlist}
do
	echo "Running test:" $test
    # make $test;
    $MPIRUN $PARTIES ./$test
done

cd ..