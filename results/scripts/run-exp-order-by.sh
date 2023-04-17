#!/bin/bash

### ### ###                ### ### ###

### ### ### INITIALIZATION ### ### ###

### ### ###                ### ### ###

### paths configuration ###
MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-order-by"

for ROWS in 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608
do
  for COLS in 1 2 4 8
  do
    for ATT in 1 2 3
    do
      if [ "$ATT" -le "$COLS" ]; then
        echo "Running experiment with ROWS=$ROWS, COLUMNS=$COLS and ATT=$ATT"
        $MPIRUN_LOCAL $TEST $ROWS $COLS $ATT
      fi
    done
  done
done