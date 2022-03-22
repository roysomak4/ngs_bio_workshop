#!/bin/bash
# update repos and install dependencies
sudo apt update
sudo apt upgrade -y
sudo apt install -y unzip nfs-common tree

# mount nfs assets share on vm
mkdir assets
sudo bash -c 'echo -e "159.203.87.19:/exports/assets\t/home/bioseq/assets\tnfs\trw,soft,intr,noacl,noatime,timeo=900,retrans=3,proto=tcp,vers=4\t0\t0" >> /etc/fstab'
sudo mount -a -v
sudo chown -R bioseq assets

# install ngs apps
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
chmod +x Mambaforge-Linux-x86_64.sh 
./Mambaforge-Linux-x86_64.sh -b 
rm -r Mambaforge-Linux-x86_64.sh
mambaforge/bin/mamba install -y -c bioconda bwa samtools fastp fastqc sambamba varscan==2.4.4 abra2

# Download annovar (refer to download link in private repo - only for administrator)
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz
rm -r annovar.latest.tar.gz
rm -r annovar/humandb

# # Install google drive download utility
# wget https://raw.githubusercontent.com/roysomak4/gdrive_download_file/master/gdrive_download.sh
# chmod +x gdrive_download.sh 
# sudo mv gdrive_download.sh /usr/local/bin/

# prepare sample data folder
mkdir data
# gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
cp assets/sample_data.zip data/
unzip data/sample_data.zip -d data/
mv data/grip_course_hd701_R1.fq.gz data/hd701_R1.fq.gz
mv data/grip_course_hd701_R2.fq.gz data/hd701_R2.fq.gz

# clean up
rm data/sample_data.zip data/grip_course.bed
rm -r Downloads

# update environment
mambaforge/bin/mamba init -q
sudo reboot