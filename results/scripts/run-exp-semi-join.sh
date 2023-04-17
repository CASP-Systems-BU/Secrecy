#!/bin/bash

### ### ###                ### ### ###

### ### ### INITIALIZATION ### ### ###

### ### ###                ### ### ###

### paths configuration ###
MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-semi-join"

for ROWS1 in 1000 10000 100000 1000000 10000000 100000000
  do
    for ROWS2 in 16 64 256 1024
    do
      echo "Running experiment with ROWS_LEFT=$ROWS1 and ROWS_RIGHT=$ROWS2"
      $MPIRUN_LOCAL $TEST $ROWS1 $ROWS2
    done
done