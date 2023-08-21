#!/bin/bash

# Check the number of arguments
if [ $# -lt 2 ] || [ $# -gt 6 ]; then
    echo "Usage: $0 exp_name run_name [deploy_type] [machine_type] [mpi_type] [machines_num]"
    exit 1
fi

# Assign the provided arguments to variables or use default values
exp_name=${1}                       # Described by yaml file.
run_name=${2}                       # Used as a prefix for the spawn machines.
deploy_type=${3:-lan}               # Either `lan` or `wan`.
machine_type=${4:-r5.xlarge}        # What type of AWS instance to use.
mpi_type=${5:-openmpi}              # Either `openmpi` or `mpich`.
machines_num=${6:-3}                # Number of machines in the MPI cluster.


lan_openmpi_images=("ami-07f39960accd8cac3" "ami-07f39960accd8cac3" "ami-07f39960accd8cac3" "ami-07f39960accd8cac3")
lan_mpich_images=("ami-084020ff2b8f2012e" "ami-084020ff2b8f2012e" "ami-084020ff2b8f2012e" "ami-084020ff2b8f2012e")
lan_regions=("us-east-2" "us-east-2" "us-east-2" "us-east-2")


wan_openmpi_images=("ami-07f39960accd8cac3" "ami-0d1a64b4421553fbe" "ami-0692c69cad93f375f" "ami-0470cf319addb98af")
wan_mpich_images=("ami-084020ff2b8f2012e" "ami-0d301312b55985630" "ami-0a19cec8272cba63b" "ami-026deefba960f1d34")
wan_regions=("us-east-2" "us-east-1" "us-west-1" "us-west-2")


# name - region - instance_type - image_id

machine_names=()
machine_regions=()
machine_types=()
machine_images=()

for i in $(seq 1 $machines_num)
do
    machine_names+=("${run_name}_${i}")

    # update the machine_regions, machine_types, machine_images.
    if [ $deploy_type == "lan" ]; then
        machine_regions+=("${lan_regions[($i-1)%4]}")
        machine_types+=("${machine_type}")
        if [ $mpi_type == "openmpi" ]; then
            machine_images+=("${lan_openmpi_images[($i-1)%4]}")
        elif [ $mpi_type == "mpich" ]; then
            machine_images+=("${lan_mpich_images[($i-1)%4]}")
        else
            echo "mpi_type must be either openmpi or mpich"
            exit 1
        fi
    elif [ $deploy_type == "wan" ]; then
        machine_regions+=("${wan_regions[($i-1)%4]}")
        machine_types+=("${machine_type}")
        if [ $mpi_type == "openmpi" ]; then
            machine_images+=("${wan_openmpi_images[($i-1)%4]}")
        elif [ $mpi_type == "mpich" ]; then
            machine_images+=("${wan_mpich_images[($i-1)%4]}")
        else
            echo "mpi_type must be either openmpi or mpich"
            exit 1
        fi
    else
        echo "deploy_type must be either lan or wan"
        exit 1
    fi
done


# echo ${machine_names[@]}
# echo ${machine_regions[@]}
# echo ${machine_types[@]}
# echo ${machine_images[@]}

# TODO: create a json string that has list of object that has the following keys: name, region, instance_type, image_id
#       for example: [{"name": "test_1", "region": "us-east-2", "instance_type": "r5.xlarge", "image_id": "ami-07f39960accd8cac3"}, {"name": "test_2", "region": "us-east-2", "instance_type": "r5.xlarge", "image_id": "ami-07f39960accd8cac3"}]
#       use the variables machine_names, machine_regions, machine_types, machine_images
#       using a for loop and string concatenation
#       use the variable machine_json to store the json string
machine_json="["
for i in $(seq 0 $((${#machine_names[@]}-1)))
do
    machine_json+="{\"name\": \"${machine_names[$i]}\", \"region\": \"${machine_regions[$i]}\", \"instance_type\": \"${machine_types[$i]}\", \"image_id\": \"${machine_images[$i]}\"}"
    if [ $i -lt $((${#machine_names[@]}-1)) ]; then
        machine_json+=","
    fi
done
machine_json+="]"

input_json="{\"exp_name\":../experiments/${exp_name}.yaml, \"machine_json\":${machine_json}}"

source ./secrets.sh
ansible-playbook ./common/run.yaml -e "${input_json}"
