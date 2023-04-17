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
6. Run ` aws configure ` to add your secrets (or Run `./secrets.sh` to the AWS secerets). You need to do this step for each new terminal.

## Running an experiment
We have two modes to run the experiments: same-region (LAN setup) and cross-region (WAN setup). You can run an existing experiment using either the script `run.yaml` in either directories as follows:
1. Make sure the AWS secrets are loaded as indicated in setup.
2. modify the `experiment-examples.yaml` file to show the experiment to execute and its paramters.
3. Run `ansible-playbook run.yaml` from the corresponding setup directory.
