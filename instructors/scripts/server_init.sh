#!/bin/bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y unzip pigz pv
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
wget https://raw.githubusercontent.com/roysomak4/gdrive_download_file/master/gdrive_download.sh
chmod +x gdrive_download.sh 
sudo mv gdrive_download.sh /usr/local/bin/

# prepare sample data and assets folder
mkdir data assets
gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
mv Downloads/sample_data.zip data/
unzip data/sample_data.zip -d data/
rm data/sample_data.zip
##############
gdrive_download.sh 1yANVV31SMqwzRNnEFiLSktxj8eQL5CSp annotation_db_and_vcfs.zip
mv Downloads/annotation_db_and_vcfs.zip assets/
unzip assets/annotation_db_and_vcfs.zip -d assets/
rm assets/annotation_db_and_vcfs.zip
tar -zxvf assets/grip_vcfs.tar.gz -C assets/
tar -zxvf assets/grip_course_annovar_db.tar.gz -C assets/
mv assets/grip_course_annovar_db assets/human_db
rm assets/grip_vcfs.tar.gz
rm assets/grip_course_annovar_db.tar.gz
rm assets/*.md5
#################
gdrive_download.sh 1mgEkHKU7hPSRkFlPcoRVE1S-QKoHvzwf gatk_b37.tar.gz
gdrive_download.sh 1pq5zpiXsOAEklZAoNF9gB0M4Nlsacpqj gatk_b37.tar.gz.md5
mv Downloads/gatk_b37.tar.gz assets/
tar -zxvf assets/gatk_b37.tar.gz -C assets/
rm assets/gatk_b37.tar.gz

rm -r Downloads
