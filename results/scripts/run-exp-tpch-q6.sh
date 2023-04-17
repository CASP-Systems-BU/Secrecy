#!/bin/bash

MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-tpch-q6"

make clean
make exp-tpch-q6  > /dev/null 2>&1
for ROWS in 1000 10000 100000
  do
    echo "Running experiment exp-tpch-q6 with ROWS=$ROWS"
    $MPIRUN_LOCAL $TEST $ROWS
done