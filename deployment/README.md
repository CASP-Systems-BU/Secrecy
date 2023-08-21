# Deployment Machine
We use a deployment machine to create the computing cluster and deploy Secrecy. The deployment machine can be used to run queries and provide input paramters for experiments. 

## Setup
1. Create a new Ubuntu machine.
2. Generate a new ssh pair key with the default location and with name 'id_rsa'.
3. Clone this repository in the home directory of the new machine:
` git clone https://github.com/CASP-Systems-BU/secrecy `
4. Install Ansible: 
` sudo apt update `
` sudo apt install ansible python3-pip awscli`
` pip3 install boto3 `
5. Update the `secrets.sh` script to reflect your AWS secrets.

## Adding an experiment 
Check the experiment example under `experiments` directory

## Running an experiment
We have two modes to run the experiments: same-region (LAN setup) and cross-region (WAN setup). You can run an existing experiment using the script `deploy.sh`:
1. Make sure the AWS secrets are loaded as indicated in setup.
2. Make sure you are in the `secrecy/deployment` directory.
3. Run the `deploy.sh` script which has the following parameters:
```
# Assign the provided arguments to variables or use default values
exp_name=${1}                       # Described by yaml file.
run_name=${2}                       # Used as a prefix for the spawn machines.
deploy_type=${3:-lan}               # Either `lan` or `wan`.
machine_type=${4:-r5.xlarge}        # What type of AWS instance to use.
mpi_type=${5:-openmpi}              # Either `openmpi` or `mpich`.
machines_num=${6:-3}                # Number of machines in the MPI cluster.
```
For example, run from `secrecy/deployment`, the following command.
```
./deploy.sh experiment-example friday-run lan r5.xlarge openmpi
```
