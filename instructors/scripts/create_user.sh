#! /bin/bash

useradd -m -U bioseq -s /bin/bash
echo -e "bioseq123\nbioseq123" | passwd bioseq
usermod -aG sudo bioseq
sed -i "s/^PasswordAuthentication no/PasswordAuthentication yes/" /etc/ssh/sshd_config
systemctl restart sshd
cp -R /root/.ssh /home/bioseq/
chown -R bioseq:bioseq /home/bioseq/.ssh
echo -e "%sudo ALL=(ALL:ALL) NOPASSWD:ALL" > /etc/sudoers.d/passwd_override
