#!/bin/bash

export AWS_ACCESS_KEY="YOUR_ACCESS_KEY"
export AWS_SECRET_KEY="YOUR_SECRET_KEY"
export ANSIBLE_HOST_KEY_CHECKING=False

eval `ssh-agent`
ssh-add PATH_TO_YOUR_SSH_KEY