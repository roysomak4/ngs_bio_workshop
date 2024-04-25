# Instructions to setup the Virtual Machines for the course

## 1. Create the virtual machines on Digitalocean platform

Please refer to the instructions at https://gitlab.com/roysomak4/digitalocean_manage_vms

After cloning the repository, run the following command to create the VMs
```
python3 create_vm.py hosts/<list of vm names>.txt --image ubuntu-22-04-x64 --worker-size s-4vcpu-8gb-intel --region nyc3 --tag-names ngs --ssh-keys <ssh-key-id-from-digitalocean,ssh-key-id-from-digitalocean>
``` 
How to get the SSH key IDs
```
doctl compute ssh-key ls
```

## 2. Add the IP addresses of the newly created VMs to the known_hosts list

The script above generates a output file with list of IP addresses of the newly created VMs (`vm_ips.txt`) in the `output` folder of the repo.

Use the `ssh-keyscan` utility to add the ips to the known_hosts list. Update the path to the location of the vm_ips.txt file on your system

```
ssh-keyscan -f vm_ips.txt >> ~/.ssh/known_hosts
```

## 3. Run the VM prep script

This bash script will help initialize the VM by creating non-root user and install the necessary apps for the NGS hands on workshop. It will also mount the remote NFS server with sequencing assets, such as the human genome reference sequence, dbSNP and Cosmic VCF files, and other database files.  

Change to the script directory inside the repo

```
cd instructors/scripts
```

Create a bash array of the ip addresses. Update the path to the location of the vm_ips.txt file on your system
```
vm_ips=(`cat vm_ips.txt`)
```

Execute the VM prep script using GNU parallel to iterate through the list of VMs in the newly created bash array
```
parallel -j 2 --line-buffer "./prepare_vm.sh {}" ::: "${vm_ips[@]}"
```
The `--line-buffer` argument to GNU parallel outputs the stdout in almost realtime on the console. It is good for development and debug purposes. However, when setting up multiple VMs simultaneously, it is better to not use that argument. All the stdout from the scripts will be printed at the end of the process. 
