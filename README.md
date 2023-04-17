# SECRECY: Secure collaborative analytics in untrusted clouds
This repository includes the implementation of the Secrecy relational Multi-Party Computation framework described in the [USENIX NSDI'23 paper](https://www.usenix.org/conference/nsdi23/presentation/liagouris) by John Liagouris, Vasiliki Kalavri, Muhammad Faisal, and Mayank Varia.

You can cite the paper using the BibTeX below:

```
@inproceedings {285183,
author = {John Liagouris and Vasiliki Kalavri and Muhammad Faisal and Mayank Varia},
title = {{SECRECY}: Secure collaborative analytics in untrusted clouds},
booktitle = {20th USENIX Symposium on Networked Systems Design and Implementation (NSDI 23)},
year = {2023},
isbn = {978-1-939133-33-5},
address = {Boston, MA},
pages = {1031--1056},
url = {https://www.usenix.org/conference/nsdi23/presentation/liagouris},
publisher = {USENIX Association},
month = apr,
}
```

**NOTICE**: This is an academic proof-of-concept prototype and has not received careful code review. This implementation is NOT ready for production use.

## Building and running Secrecy

### Repository organization
This repository is organized as follows:
- The `deployment` folder contains automation scripts to deploy secrecy in AWS in either same region or cross region setups.
- The `src` folder contains the core functionality of Secrecy, including the implementation of MPC primitives, relational oblivious operators, and party communication.
- The `examples` folder contains the implementation of example queries with the Secrecy API.
- The `src/test` folder contains various unit and end-to-end tests.
- The `src/experiments` folder contains the implementation of various microbenchmarks and performance experiments.
- Plotting scripts and other helper utilies are located in the `results/scripts` folder.


### Dependencies
To build Secrecy, you will need to install:
- CMake
- [Libsodium](https://libsodium.gitbook.io/doc/installation)
- an MPI implementation, such as [OpenMPI](https://www.open-mpi.org/software/ompi/v4.0/) or [MPICH](https://www.mpich.org/downloads/).
- Create a new directory `include/external-lib` and clone `https://github.com/mfaisal97/sql-parser` inside it.

<!-- Build and run the tests
------------
Change to the `tests` directory.

1. Build and run all tests: 
   - Run `make tests`. 

2. Build and run an individual test: 
   - Run `make test-xyz` to build a test, where `xyz` is the test name. For instance, run `make test-equality` to build the binary equality test. 
   - Execute the test with `mpirun -np 3 test-xyz`.

Run an example
---------
Change to the `examples` directory.

1. Build all examples: 
   - Run `make all`. 

2. Build and run an individual example, e.g. the comorbidity query: 
   - Build the example with `make comorbidity`.
   - Run the example with `mpirun -np 3 comorbidity <NUM_ROWS_1> <NUM_ROWS_2>`.

Run the experiments
---------
Change to the `experiments` directory.

1. Build all experiments: 
   - Run `make all`. 

2. Build and run an individual experiment, e.g. the equality microbenchmark: 
   - Build the experiment with `make exp-equality`.
   - Run it with `mpirun -np 3 exp-equality <INPUT_SIZE>`. -->

<!-- Specifying dependencies on Linux
-------------
To build and run Secrecy on linux, edit the provided `Makefile` as follows:
- Use the variables `CFLAGS= -03 -Wall` and `DEP= -lsodium -lm`
- Specify the dependency in the end of the target, for example:

    `exp-equality:   exp_equality.c $(PRIMITIVES) $(MPI) $(CFLAGS) -o exp-equality exp_equality.c $(PRIMITIVES) $(DEP)` -->
    
### Building and running using the `CMakeLists.txt` file:

- Make sure you have pkg-config installed. 

   On Linux:
   ``` sudo apt update
   sudo apt install pkg-config 
   ```    
   On OSX:
   ```
   brew install pkg-config
   ```
- To compile and run all test cases:
   ```
   ./run_tests.sh
   ```
- To run a specific test file or experiment: 
   - First use cmake to create the build dir and make file:
      ```
      mkdir build
      cd build
      cmake ..
      ```
   - Build and run a test file:
      ```
      cd build
      make test_planner_q1.out
      mpirun -np 3 ./test_planner_q1.out [X] [X]  // Will compile and run Q1 with X rows per input
      ```
      
## License
Secrecy is distributed under the terms of the Apache License (Version 2.0).
