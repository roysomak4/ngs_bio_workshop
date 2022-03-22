# NGS Bioinformatics Hands On Lab

1. Examine content of the FASTQ file
   ```
   zcat data/hd701_R1.fq.gz | less -N
   ```
2. Perform FASTQ processing prior to alignment
   ```
   fastp -i data/hd701_R1.fq.gz -o data/hd701_R1_processed.fq.gz -I data/hd701_R2.fq.gz -O data/hd701_R2_processed.fq.gz -z 4 -w 2 -h data/hd701_fq_qc.html
   ```
3. Generate FASTQ quality control data
   ```
   fastqc -o data -f fastq data/*.gz
   ```
4. :checkered_flag: 
   Download FASTQ quality control data from server (Windows)
   ```
   scp 'bioseq@<your-ip-address>:~/data/*.html' ./Downloads/
   ```
   Download FASTQ quality control data from server (Linux / macOS)
   ```
   scp 'bioseq@<your-ip-address>:~/data/*.html' ~/Downloads/
   ```

5. Align sequence to GRCh37
   ```
   bwa mem -M -v 1 -t 2 -R "@RG\\tID:hd701\\tSM:hd701\\tPL:ILLUMINA_NEXTSEQ\\tPI:150\\tCN:roy_lab" assets/gatk_b37/human_g1k_v37.fasta data/hd701_R1_processed.fq.gz data/hd701_R2_processed.fq.gz >data/hd701_raw.sam
   ```
6. Convert the SAM file into a BAM file
   ```
   samtools view -Shu data/hd701_raw.sam >data/hd701_raw.bam
   ```
7. Sort and index raw BAM file
   ```
   sambamba sort -p -t 2 -o data/hd701_sorted.bam data/hd701_raw.bam
   ```
8. Mark PCR duplicates in BAM file
   ```
   sambamba markdup -t 2 -p data/hd701_sorted.bam data/hd701_dedup.bam
   ```
9. Perform indel realignment
   ```
   abra2 --in data/hd701_dedup.bam --out data/hd701_realigned.bam --ref assets/gatk_b37/human_g1k_v37.fasta --threads 2 --targets assets/ngs_workshop.bed --index --tmpdir /tmp/ >data/hd701_abra.log
   ```
10. Inspect aligned sequences
    ```
    samtools view -h data/hd701_realigned.bam | less -N
    ```
11. Call variant
    ```
    samtools mpileup -BA -q 20 -Q 30 -d 4000 -l assets/ngs_workshop.bed -f assets/gatk_b37/human_g1k_v37.fasta data/hd701_realigned.bam | varscan -Xmx4G mpileup2vcf --min-coverage 8 --min-var-freq 0.05 --p-value 0.05 --min-avg-qual 30 --strand-filter 1 --output-vcf 1 --variants >data/hd701_raw.vcf
    ```
12. Compress and index VCF file
    ```
    bgzip data/hd701_raw.vcf
    tabix -p vcf data/hd701_raw.vcf.gz
    ```
13. Inspect variants in a raw VCF file
    ```
    zcat data/hd701_raw.vcf.gz | less -N
    ```
14. Annotate variants
    ```
    annovar/table_annovar.pl data/hd701_raw.vcf.gz assets/humandb/ -out data/hd701_annotated -buildver hg19 -remove -nastring . -otherinfo -vcfinput -thread 2 -maxgenethread 2 -protocol refGene,cytoBand,cosmic85,clinvar_20150330 -operation g,r,f,f
    ```
15. Inspect annotated VCF file
    ```
    less data/hd701_annotated.hg19_multianno.vcf
    ```
16. :checkered_flag: 
    Download BAM files for visualization (Windows)
    ```
    scp bioseq@<your-ip-address>:~/data/hd701_realigned.ba* ./Downloads/
    ```
    Download BAM files for visualization (Linux / macOS)
    ```
    scp 'bioseq@<your-ip-address>:~/data/hd701_realigned.ba*' ~/Downloads/
    ```
