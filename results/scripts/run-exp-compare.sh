#!/bin/bash

### ### ###                ### ### ###

### ### ### INITIALIZATION ### ### ###

### ### ###                ### ### ###

### paths configuration ###
MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-cmp-swap"

for ROWS in 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608
  do
    for COLS in 1 2 4 8
    do
      echo "Running experiment with ROWS=$ROWS and COLUMNS=$COLS"
      $MPIRUN_LOCAL $TEST $ROWS $COLS
    done
done