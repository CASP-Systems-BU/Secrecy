#!/bin/bash

### ### ###                ### ### ###

### ### ### INITIALIZATION ### ### ###

### ### ###                ### ### ###

### paths configuration ###
MPIRUN_LOCAL="mpirun -n 3"
MPIRUN="mpirun -n 3 -f host_file"
TEST="./exp-equality"

### LOCAL RUN (3 MPI Processes on the same node)
for INPUT_SIZE in 1000 10000 100000 1000000 10000000 #100000000
do
  echo "Running experiment with INPUT_SIZE=$INPUT_SIZE"
  $MPIRUN_LOCAL $TEST $INPUT_SIZE
done