#!/bin/bash

### ### ###                ### ### ###

### ### ### INITIALIZATION ### ### ###

### ### ###                ### ### ###

### paths configuration ###
MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-join"

for INPUT_SIZE in 100 1000 10000 20000 30000
do
  echo "Running experiment with INPUT_SIZE=$INPUT_SIZE"
  $MPIRUN_LOCAL $TEST $INPUT_SIZE
done