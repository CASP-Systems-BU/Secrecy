#!/bin/bash

MPIRUN="mpirun -np"
PARTIES=3
tests=( test-primitives test-share test-rsz test-and-all test-equality test-inequality test-rca test-join test-select test-relational test-conversion test-sort test-distinct test-group test-group-odd-even test-in test-mask test-arithmetic-to-binary test-ageq test-q1 test-q2 test-q3 test-group-join test-group-min-max test-tpch-q4 test-tpch-q4-baseline test-tpch-q6 test-tpch-q6-baseline test-tpch-q13 test-tpch-q13-baseline test-qpwd test-qpwd-baseline test-qcredit test-qcredit-baseline)

make clean
for test in "${tests[@]}"
do
	echo "Running test:" $test
    make $test; $MPIRUN $PARTIES ./$test
done
