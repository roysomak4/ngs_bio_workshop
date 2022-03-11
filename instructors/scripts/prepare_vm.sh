#!/bin/bash

VM_IP=$1

# create non root user and configure SSH access
echo "Creating user and configuring SSH access for ${VM_IP}..."
scp create_user.sh root@${VM_IP}:/root/
ssh -t root@${VM_IP} "/root/create_user.sh"
echo "User created and SSH access configured."

# configure VM to install NGS apps and mount NFS share
echo "Installing NGS apps and mounting NFS share for genomic assets on ${VM_IP}..."
scp configure_vm.sh bioseq@${VM_IP}:/home/bioseq/
ssh -t bioseq@${VM_IP} "/home/bioseq/configure_vm.sh"
echo "VM configuration complete."