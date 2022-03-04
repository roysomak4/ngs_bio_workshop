# Sequencing Informatics Hands On Lab

1. Visualize the FASTQ file
   ```
   zcat data/grip_course_hd701_R1.fq.gz | less -N
   ```
2. Perform FASTQC step
   ```
   fastp -i data/grip_course_hd701_R1.fq.gz -o data/hd701_R1_qc.fq.gz -I data/grip_course_hd701_R2.fq.gz -O data/hd701_R2_qc.fq.gz -z 4 -w 2 -h data/hd701_fq_qc.html
   ```
3. Align sequence to GRCh37
   ```
   bwa mem -M -v 1 -t 2 -R "@RG\\tID:hd701\\tSM:hd701\\tPL:ILLUMINA_NEXTSEQ\\tPI:150\\tCN:roy_lab" assets/gatk_b37/human_g1k_v37.fasta data/hd701_R1_qc.fq.gz data/hd701_R2_qc.fq.gz >data/hd701_raw.sam
   ```
4. Encode the SAM file into a BAM file
   ```
   samtools view -Shu data/hd701_raw.sam >data/hd701_raw.bam
   ```
5. Sort and index raw BAM file
   ```
   sambamba sort -p -t 2 -o data/hd701_sorted.bam data/hd701_raw.bam
   ```
6. Mark PCR duplicates
   ```
   sambamba markdup -t 2 -p data/hd701_sorted.bam data/hd701_dedup.bam
   ```
7. Create target regions for indel realignment
   ```
   abra2 --in data/hd701_dedup.bam --out data/hd701_realigned.bam --ref assets/gatk_b37/human_g1k_v37.fasta --threads 2 --targets assets/grip_course.bed --index --tmpdir /tmp/ >data/hd701_abra.log
   ```
8. Index the realigned bam
   ```
   sambamba index data/hd701_realigned.bam
   ```
9. Visualize the sequence alignment
   ```
   samtools view -h data/hd701_realigned.bam | less -N
   ```
10. Call variant
    ```
    samtools mpileup -BA -q 20 -Q 30 -d 4000 -l assets/grip_course.bed -f assets/gatk_b37/human_g1k_v37.fasta data/hd701_realigned.bam | varscan -Xmx4G mpileup2vcf --min-coverage 8 --min-var-freq 0.05 --p-value 0.05 --min-avg-qual 30 --strand-filter 1 --output-vcf 1 --variants >data/hd701_raw.vcf
    ```
11. Compress and index VCF file

```
bgzip data/hd701_raw.vcf
tabix -p vcf data/hd701_raw.vcf.gz
```

12. Annotate variants

```
annovar/table_annovar.pl data/hd701_raw.vcf.gz assets/human_db/ -out hd701_annotated -buildver hg19 -remove -nastring . -otherinfo -vcfinput -thread 2 -maxgenethread 2 -protocol refGene,cytoBand,cosmic85,clinvar_20150330 -operation g,r,f,f
```
