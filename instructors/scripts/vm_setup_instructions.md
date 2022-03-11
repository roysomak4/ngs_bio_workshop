# Instructions to setup the Virtual Machines for the course

## 1. Create the virtual machines on Digitalocean platform

Please refer to the instructions at https://gitlab.com/roysomak4/digitalocean_manage_vms

## 2. Add the IP addresses of the newly created VMs to the known_hosts list

The script above generates a output file with list of IP addresses of the newly created VMs. Use that list to complete this step.  
Example of ip list

```
# nano vm_ips.txt

192.168.10.1
192.168.10.2
192.168.10.3
192.168.10.4
192.168.10.5
```

Use the `ssh-keyscan` utility to add the ips to the known_hosts list

```
ssh-keyscan -f vm_ips.txt >> ~/.ssh/known_hosts
```

## 3. Run the VM prep script

This bash script will help initialize the VM by creating non-root user and install the necessary apps for the NGS hands on workshop. It will also mount the remote NFS server with sequencing assets, such as the human genome reference sequence, dbSNP and Cosmic VCF files, and other database files.

```
ssh-keyscan -f ~/projects/digitalocean_manage_vms/output/vm_ips.txt >>~/.ssh/known_host


```
