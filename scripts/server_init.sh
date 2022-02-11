#!/bin/bash

sudo apt update
sudo apt upgrade -y
sudo apt install -y unzip pigz
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh
chmod +x Mambaforge-pypy3-Linux-x86_64.sh 
./Mambaforge-pypy3-Linux-x86_64.sh 
mamba install -qy -c bioconda bwa samtools fastp sambamba varscan==2.4.4 abra2
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
wget https://raw.githubusercontent.com/roysomak4/gdrive_download_file/master/gdrive_download.sh
chmod +x gdrive_download.sh 
sudo mv gdrive_download.sh /usr/local/bin/

gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
gdrive_download.sh 1yANVV31SMqwzRNnEFiLSktxj8eQL5CSp annotation_db_and_vcfs.zip
gdrive_download.sh 1mgEkHKU7hPSRkFlPcoRVE1S-QKoHvzwf gatk_b37.tar.gz
gdrive_download.sh 1pq5zpiXsOAEklZAoNF9gB0M4Nlsacpqj gatk_b37.tar.gz.md5
cd Downloads
md5sum gatk_b37.tar.gz
cat gatk_b37.tar.gz.md5 
tar -zxvf gatk_b37.tar.gz 
unzip sample_data.zip
unzip annotation_db_and_vcfs.zip
tar -zxvf grip_vcfs.tar.gz
tar -zxvf grip_course_annovar_db.tar.gz