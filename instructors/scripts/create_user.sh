#! /bin/bash

# update repos and install dependencies
apt update
# apt upgrade -y
apt install -y unzip nfs-common tree
# create user
useradd -m -U bioseq -s /bin/bash
echo -e "bioseq123\nbioseq123" | passwd bioseq
usermod -aG sudo bioseq
sed -i "s/^PasswordAuthentication no/PasswordAuthentication yes/" /etc/ssh/sshd_config
sed -i "s/^KbdInteractiveAuthentication no/KbdInteractiveAuthentication yes/" /etc/ssh/sshd_config
systemctl restart sshd
cp -R /root/.ssh /home/bioseq/
chown -R bioseq:bioseq /home/bioseq/.ssh
echo -e "%sudo ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/passwd_override
