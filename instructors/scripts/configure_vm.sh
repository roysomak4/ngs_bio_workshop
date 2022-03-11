#!/bin/bash
# update repos and install dependencies
sudo apt update
sudo apt upgrade -y
sudo apt install -y unzip nfs-common

# mount nfs assets share on vm
mkdir assets
sudo mount -t nfs -o rw,soft,intr,noacl,noatime,timeo=900,retrans=3,proto=tcp,vers=4 159.203.87.19:/exports/assets assets
sudo chown -R bioseq assets

# install ngs apps
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh
chmod +x Mambaforge-pypy3-Linux-x86_64.sh 
./Mambaforge-pypy3-Linux-x86_64.sh -b 
rm -r Mambaforge-pypy3-Linux-x86_64.sh
mamba install -y -c bioconda bwa samtools fastp fastqc sambamba varscan==2.4.4 abra2

# Download annovar (refer to download link in private repo - only for administrator)
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz
rm -r annovar.latest.tar.gz
rm -r annovar/humandb

# Install google drive download utility
wget https://raw.githubusercontent.com/roysomak4/gdrive_download_file/master/gdrive_download.sh
chmod +x gdrive_download.sh 
sudo mv gdrive_download.sh /usr/local/bin/

# prepare sample data folder
mkdir data
gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
mv Downloads/sample_data.zip data/
unzip data/sample_data.zip -d data/

# clean up
rm data/sample_data.zip data/grip_course.bed
rm -r Downloads