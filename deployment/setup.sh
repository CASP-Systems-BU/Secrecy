#!/bin/bash

sudo apt update
sudo apt install git cmake pkg-config build-essential manpages-dev gfortran wget -y

mkdir -p include/external-lib 
cd include/external-lib && git clone https://github.com/mfaisal97/sql-parser && cd ../..

sudo mkdir -p /src
sudo chmod -R 777 /src

# Setup libsoduim
cd /src
wget https://download.libsodium.org/libsodium/releases/libsodium-1.0.18-stable.tar.gz
tar -xzf libsodium-1.0.18-stable.tar.gz
cd /src/libsodium-stable
./configure
make -j16
make
make check
sudo make install
sudo ldconfig



# # Setup MPICH
# cd /src
# wget https://www.mpich.org/static/downloads/4.1/mpich-4.1.tar.gz
# tar -xzf mpich-4.1.tar.gz
# cd /src/mpich-4.1
# ./configure
# make -j4
# make
# sudo make install
# sudo ldconfig


# Setup OPEN-MPI
cd /src
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz
tar -xzf openmpi-4.1.5.tar.gz
cd /src/openmpi-4.1.5
./configure
make -j16
sudo make install
sudo ldconfig