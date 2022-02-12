#!/bin/bash

gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
gdrive_download.sh 1yANVV31SMqwzRNnEFiLSktxj8eQL5CSp annotation_db_and_vcfs.zip
gdrive_download.sh 1mgEkHKU7hPSRkFlPcoRVE1S-QKoHvzwf gatk_b37.tar.gz
gdrive_download.sh 1pq5zpiXsOAEklZAoNF9gB0M4Nlsacpqj gatk_b37.tar.gz.md5

md5sum -c Downloads/gatk_b37.tar.gz.md5
mkdir assets data

mv Downloads/gatk_b37.tar.gz assets/
tar -zxvf assets/gatk_b37.tar.gz
rm assets/gatk_b37.tar.gz

mv Downloads/annotation_db_and_vcfs.zip assets/
unzip assets/annotation_db_and_vcfs.zip
tar -zxvf assets/grip_vcfs.tar.gz
tar -zxvf assets/grip_course_annovar_db.tar.gz
rm -r assets/annotation_db_and_vcfs.zip assets/grip_vcfs.tar.gz assets/grip_course_annovar_db.tar.gz

mv Downloads/sample_data.zip data/
unzip data/sample_data.zip
rm data/sample_data.zip

