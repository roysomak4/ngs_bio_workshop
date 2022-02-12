# Sequencing Informatics Hands On Lab
## GRIP course 2019

### Instructions

#### Steps to be executed by the course instructor
1. Create Virtual Machines (VMs) using predefined images for hands on practical
2. Once the VMs are created in Digitalocean, distribute the public IP addresses to the students. Each student is assigned one IP address.

#### Steps to be followed by the students
##### Preparing the system for accessing the VM
> **Installing PuTTy (Windows OS)**
1. Download PuTTy on your system from the following link: https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html
2. Click on `putty-64bit-0.72-installer.msi` for a 64 bit windows installer.
3. Once downloaded, double click on the installer file and follow the setup instructions.
4. This will install `putty.exe`, `pagent.exe` and `pscp.exe` on your windows computer in the folder `C:\Program Files\PuTTY`.
5. Download SSH keys for accessing the VMs from: https://drive.google.com/open?id=1mWcMYyyhXGIU15udYzdczzLBevrZGXv2
6. Open `pagent.exe` located in `C:\Program Files\PuTTY`
7. The `pagent` starts in the status bar (lower right).
8. Right click and choose the option `Add key`.
9. Choose the `id_grip_rsa.ppk` file that was downloaded.
10. Start `putty.exe`. Enter the IP address provided by your instructor in in the `hostname (IP address)` textbox.
11. Click `Open`. This should open up the terminal requesting to enter the username for the VM to login in to.

> **Linux or Mac**
1. Download SSH keys for accessing the VMs from: https://drive.google.com/open?id=1mWcMYyyhXGIU15udYzdczzLBevrZGXv2
2. Copy the ssh keys (public and private) into you SSH folder
   ```
   cp <source folder>/id_grip_rsa* ~/.ssh/
   ``` 

##### Begin NGS data analysis
1. The students will login to their respective VMs with username `bioseq`. This user is a non root user with root previlages enabled.  
**Linux / Mac users**
   ```
   ssh bioseq@<ip address>
   ```
**Windows users**  
   >Follow instructions above from steps 6 to 10. When prompted for username, enter `bioseq`
 
2. There should not be any password prompt. If prompted for password, type `grip123`
3. A login screen on the command line should be seen.
4. Download the NGS assets and sample data (reference sequence, bed files, vcf files, databases, etc)
    ```
    gdrive_download.sh 1p7vZebLxPkvlUrHAsyrh44WyS5Q4iV2I sample_data.zip
    gdrive_download.sh 1yANVV31SMqwzRNnEFiLSktxj8eQL5CSp annotation_db_and_vcfs.zip
    gdrive_download.sh 1mgEkHKU7hPSRkFlPcoRVE1S-QKoHvzwf gatk_b37.tar.gz
    gdrive_download.sh 1pq5zpiXsOAEklZAoNF9gB0M4Nlsacpqj gatk_b37.tar.gz.md5
    ```
5. Change into `Downloads` directory
   ```
   cd Downloads
   ``` 
6. Verify that the downloaded file are intact
   ```
   cat gatk_b37.tar.gz.md5
   md5sum gatk_b37.tar.gz
   ```
7. extract the downloaded files
   ```
   tar -zxvf gatk_b37.tar.gz
   unzip sample_data.zip
   unzip annotation_db_and_vcfs.zip
   tar -zxvf grip_vcfs.tar.gz
   tar -zxvf grip_course_annovar_db.tar.gz

   ```
7. Create `sample_data` directory. Move sample data and bed file to `sample_data` and change into that directory
   ```
   mkdir sample_data
   mv *.fq.gz sample_data/
   mv *bed sample_data/
   cd sample_data
   ```
8. Rename directory for annotation database
   ```
   mv grip_course_annovar_db humandb
   ```
8. Perform FASTQC step
   ```
   fastp -i grip_course_hd701_R1.fq.gz -o hd701_R1_qc.fq.gz -I grip_course_hd701_R2.fq.gz -O hd701_R2_qc.fq.gz -z 4 -w 2 -h hd701_fq_qc.html
   ```
9. Align sequence to GRCh37
   ```
   bwa mem -M -v 1 -t 2 -R "@RG\\tID:hd701\\tSM:hd701\\tPL:ILLUMINA_NEXTSEQ\\tPI:150\\tCN:roy_lab" /home/bioseq/Downloads/gatk_b37/human_g1k_v37.fasta hd701_R1_qc.fq.gz hd701_R2_qc.fq.gz >hd701_raw.sam
   ```
10. Encode the SAM file into a BAM file
    ```
    samtools view -Shu hd701_raw.sam >hd701_raw.bam
    ```
11. Sort and index raw BAM file
    ```
    sambamba sort -p -t 2 -o hd701_sorted.bam hd701_raw.bam
    ```
12. Mark PCR duplicates
    ```
    sambamba markdup -t 2 -p hd701_sorted.bam hd701_dedup.bam
    ```
13. Create target regions for indel realignment
    ```
    java -Xmx3g -jar /home/bioseq/abra/abra2-2.20.jar --in hd701_dedup.bam --out hd701_realigned.bam --ref /home/bioseq/Downloads/gatk_b37/human_g1k_v37.fasta --threads 2 --targets grip_course.bed --tmpdir /tmp/ >hd701_abra.log
    ```
14. Index the realigned bam
    ```
    sambamba index hd701_realigned.bam
    ``` 
15. Call variant
    ```
    samtools mpileup -BA -q 20 -Q 30 -d 4000 -l grip_course.bed -f /home/bioseq/Downloads/gatk_b37/human_g1k_v37.fasta hd701_realigned.bam | varscan -Xmx4G mpileup2vcf --min-coverage 8 --min-var-freq 0.05 --p-value 0.05 --min-avg-qual 30 --strand-filter 1 --output-vcf 1 --variants >hd701_raw.vcf
    ```
16. Compress and index VCF file
    ```
    bgzip hd701_raw.vcf
    tabix -p vcf hd701_raw.vcf.gz
    ```
17. Annotate variants
    ```
    perl /home/bioseq/annovar/table_annovar.pl hd701_raw.vcf.gz /home/bioseq/Downloads/humandb/ -out hd701_annotated -buildver hg19 -remove -nastring . -otherinfo -vcfinput -thread 2 -maxgenethread 2 -protocol refGene,cytoBand,cosmic85,clinvar_20150330 -operation g,r,f,f
    ```