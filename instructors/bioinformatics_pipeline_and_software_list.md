## pipeline for workshop
1. fastqc (QC only)
2. fastp (QC and filter)
3. bwa (alignment)
4. samtools (SAM -> BAM, sort and index)
5. sambamba (dedup)
6. abra2 (realignment)
7. Variant calling
   1. Round 1
      1. varscan (variant call) -> this will miss <5% VAF variants and FLT3 in first pass
      2.  fastvep (annotation)
      3.  benchmark VCF and demonstrate missing low VAF variant and FLT3 ITD
   2.  Round 2
       1. varscan run at lower VAF filter (to catch low VAF variants)
       2. fastvep (annotation)
       3. benchmark VCF and demonstrate the impact of tweaking VAF filter to identify low VAF variant
   3.  Round 3
       1.  Run vardict (variant call) to pick FLT3 ITD
       2.  fastvep (annotation)
8.  View the BAM file, VCF file and BED file in IGV
9.  Generate clinical report

## List of bioinformatics software tools
1. fastp
2. fastqc
3. bwa
4. samtools
5. sambamba
6. fastvep
7. varscan==2.4.4
8. abra2
9. vardict
10. mosdepth
11. multiqc

### Human genome version
GRCh38
