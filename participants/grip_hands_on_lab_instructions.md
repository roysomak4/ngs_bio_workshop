# NGS Bioinformatics Hands-On Lab

> A web version with copy buttons and progress tracking is in
> [`docs/index.html`](../docs/index.html), served by GitHub Pages. If you change
> a command here, change it there too.

You will take raw sequencing reads from a tumor reference sample all the way to
annotated variant calls, one command at a time. Then you will score your calls
against the known answer and see what each change to the pipeline gains you.

- **Sample:** a commercial tumor reference standard (OncoSpan) with 18 known
  variants in 9 cancer genes, including one large FLT3 internal tandem duplication (ITD)
- **Reads:** 49,642 Illumina read pairs, about 140 bp long, from 17 target regions
- **Genome:** GRCh38 (hg38)

```
FASTQ ─ fastqc, fastp ─ bwa ─ samtools ─ sambamba ─ ABRA2 ─ mosdepth, MultiQC
                                                       │
           ┌───────────────────────────────────────────┼──────────────────────────────┐
     Round 1: VarScan, VAF ≥ 5%        Round 2: VarScan, VAF ≥ 2%        Round 3: VarDict
           └──── fastVEP annotation and benchmark against the truth set after each round ┘
                                                       │
                                        IGV: look at the reads behind each call
```

## How this lab works

**Run every command in your terminal.** Copy a command from a grey box, paste
it into the terminal, and press Enter. Wait for the prompt (`$`) to come back
before pasting the next one.

- macOS / Linux: paste with <kbd>Cmd</kbd>+<kbd>V</kbd> or <kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>V</kbd>
- Windows PowerShell: paste with a right-click (or <kbd>Ctrl</kbd>+<kbd>V</kbd> in Windows Terminal)
- Paste one box at a time. A box with lines ending in `\` is a single command.

**View results in your browser.** Open `http://<your-ip>/` and log in with the
same username and password as SSH. The **QC reports** page shows HTML reports,
and the **IGV** page shows your BAM and VCF files. Files appear there as soon
as the command that writes them finishes.

**"Look inside" boxes** are optional peeks at what a step changed. They are
the best part of the lab, so try them if you have time.

---

## 0. Log in and get oriented

Log in with the IP address, username and password on your handout:

```
ssh bioseq@<your-ip>
```

Your account is already set up. `~/data` holds the input files, and you will
write everything into `~/work`. Move there now:

```
cd ~/work
```

Shortcuts for every file you will use are already loaded. Check them:

```
echo "sample=$SAMPLE threads=$THREADS"; echo "reference=$REF"; echo "targets=$BED"
```

> **Look inside:** `ls -lh ~/data` lists the input files. `cat ~/workshop.env`
> shows how the shortcuts are defined. `head -3 $BED` shows the first target regions.

---

## Part 1 — From reads to an analysis-ready BAM

### 1. Quality check the raw reads (fastqc)

FASTQ files store each read as four lines: a name, the bases, a `+` separator,
and one quality character per base. fastqc summarizes quality across all reads.

```
fastqc -t $THREADS -o qc $R1 $R2
```

> **Look inside:** `zcat $R1 | head -4` prints the first read. Each quality
> character encodes a Phred score: `printf '%d\n' "'F"` prints 70, and
> 70 − 33 = Q37, roughly a 1-in-5,000 chance that the base is wrong.

**View:** portal → QC reports → `ngs_ws_R1_fastqc.html`

### 2. Trim and filter the reads (fastp)

fastp removes adapter sequence and low-quality read ends, and drops reads that
become too short.

```
fastp -i $R1 -I $R2 -o $SAMPLE.trim_R1.fq.gz -O $SAMPLE.trim_R2.fq.gz -w $THREADS -h qc/$SAMPLE.fastp.html -j qc/$SAMPLE.fastp.json
```

> **Look inside:** `jq .summary qc/$SAMPLE.fastp.json` compares reads and Q30
> before and after filtering.

**View:** portal → QC reports → `ngs_ws.fastp.html`

### 3. Align the reads to the genome (bwa)

bwa finds where each read pair comes from in the 3-billion-base GRCh38 genome.
The `-R` option tags every read with the sample name. This step takes a couple
of minutes, mostly to load the genome index into memory.

```
bwa mem -v 1 -t $THREADS -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" $REF $SAMPLE.trim_R1.fq.gz $SAMPLE.trim_R2.fq.gz > $SAMPLE.sam
```

> **Look inside:** SAM is plain text with one line per read. Columns 1–9 are
> the read name, flag, chromosome, position, mapping quality (MAPQ), CIGAR
> (how the read lines up, e.g. `140M` = 140 matching bases), and mate information:
>
> ```
> grep -v '^@' $SAMPLE.sam | head -3 | cut -f1-9
> ```

### 4. Convert to BAM, sort and index (samtools)

BAM is the compressed binary form of SAM. Sorting by position and indexing let
tools jump straight to any region.

```
samtools view -b -o $SAMPLE.raw.bam $SAMPLE.sam
```
```
samtools sort -@ $THREADS -o $SAMPLE.sorted.bam $SAMPLE.raw.bam
```
```
samtools index $SAMPLE.sorted.bam
```

> **Look inside:** `ls -lh $SAMPLE.sam $SAMPLE.*.bam` compares file sizes.
> `samtools view -c $SAMPLE.sorted.bam chr12:25245350` counts the reads that
> cover KRAS codon 12.

### 5. Mark PCR duplicates (sambamba)

Library PCR makes identical copies of some fragments. Copies are not
independent evidence, so they are flagged (SAM flag 1024) and ignored by the
variant callers.

```
sambamba markdup -t $THREADS $SAMPLE.sorted.bam $SAMPLE.dedup.bam 2> qc/$SAMPLE.markdup.log
```
```
samtools index $SAMPLE.dedup.bam
```

> **Look inside:** `samtools flagstat $SAMPLE.dedup.bam` shows how many reads
> are now marked as duplicates.

### 6. Realign reads around indels (ABRA2)

Aligners place reads one at a time, so reads carrying an insertion or deletion
are often misaligned near its ends. ABRA2 assembles the reads in each target
region and realigns them together.

```
mkdir -p ~/tmp
```
```
abra2 --in $SAMPLE.dedup.bam --out $SAMPLE.realigned.bam --ref $REF --targets $BED --threads $THREADS --tmpdir ~/tmp > qc/$SAMPLE.abra2.log 2>&1
```
```
samtools index $SAMPLE.realigned.bam
```

From here on, `$BAM` refers to `ngs_ws.realigned.bam`.

> **Look inside (in IGV, step 8):** compare `ngs_ws.dedup.bam` and
> `ngs_ws.realigned.bam` at the **FBXW7 frameshift**. After ABRA2, more reads
> show the 1-base deletion instead of mismatches near their ends.

### 7. Check coverage and combine the QC (mosdepth, samtools, MultiQC)

Low-VAF variants need deep coverage. mosdepth measures depth across each
target region.

```
mosdepth -t $THREADS -n -b $BED --thresholds 100,250,500 qc/$SAMPLE $BAM
```
```
samtools flagstat $BAM > qc/$SAMPLE.flagstat.txt
```

MultiQC collects every QC result so far into one report:

```
multiqc -f -o qc/multiqc qc
```

> **Look inside:** the five targets with the lowest mean depth:
>
> ```
> zcat qc/$SAMPLE.regions.bed.gz | sort -k5,5n | head -5
> ```

**View:** portal → QC reports → `multiqc` → `multiqc_report.html`

---

## Part 2 — Calling variants in three rounds

Each round calls variants from the same BAM, annotates them with fastVEP, and
scores them against the truth set: the 18 variants known to be in this sample.

### Round 1 — VarScan with a 5% VAF cutoff

`samtools mpileup` counts the bases at every position in the targets. VarScan
then reports positions where the alternate allele passes its thresholds:
at least 5% of reads (`--min-var-freq 0.05`), at least 8× coverage, and
support on both strands.

```
samtools mpileup -B -q 20 -Q 30 -d 4000 -l $BED -f $REF $BAM | varscan mpileup2vcf --min-coverage 8 --min-var-freq 0.05 --p-value 0.05 --min-avg-qual 30 --strand-filter 1 --output-vcf 1 --variants > r1.vcf
```

Compress and index the VCF, so IGV and bcftools can read it:

```
bgzip -f r1.vcf && tabix -f -p vcf r1.vcf.gz
```

> **Look inside:** each call with its VarScan allele frequency (`FREQ`):
>
> ```
> bcftools query -f '%CHROM:%POS\t%REF>%ALT\t[%FREQ]\n' r1.vcf.gz
> ```
>
> The read counts behind one call, BRAF V600E:
>
> ```
> samtools mpileup -r chr7:140753336-140753336 -f $REF $BAM
> ```

**Annotate** with fastVEP: genes, transcripts, consequences, HGVS names and
ClinVar. `annotate-vcf` prints the full `fastvep` command before running it.

```
annotate-vcf r1.vcf.gz r1
```

> **Look inside:** gene, consequence and protein change for each call:
>
> ```
> awk -F'\t' -v OFS='\t' '/^##/ {next} !h {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next} $c["HGVSp"] != "-" {print $c["Location"], $c["SYMBOL"], $c["Consequence"], $c["HGVSp"]}' r1.fastvep.tsv | sort -u | column -t
> ```
>
> `r1.fastvep.vcf` holds the same information in the VCF's `CSQ` field.

**Benchmark** against the truth set:

```
benchmark-vcf r1.vcf.gz
```

**What to expect:** most variants are found, but not EGFR L858R (3%), KRAS
G12D (4%) or the FLT3 ITD. BRAF V600E at 6% sits right at the cutoff.

**Discuss:** is a 3% EGFR L858R clinically important? What limits how low the
cutoff can go?

### Round 2 — VarScan with a 2% VAF cutoff

Same command, with `--min-var-freq 0.02`:

```
samtools mpileup -B -q 20 -Q 30 -d 4000 -l $BED -f $REF $BAM | varscan mpileup2vcf --min-coverage 8 --min-var-freq 0.02 --p-value 0.05 --min-avg-qual 30 --strand-filter 1 --output-vcf 1 --variants > r2.vcf
```
```
bgzip -f r2.vcf && tabix -f -p vcf r2.vcf.gz
```
```
annotate-vcf r2.vcf.gz r2
```
```
benchmark-vcf r2.vcf.gz
```

> **Look inside:** how many calls did each round make?
>
> ```
> for f in r1.vcf.gz r2.vcf.gz; do echo "$f: $(bcftools view -H $f | wc -l) calls"; done
> ```

**What to expect:** L858R and G12D are now found. Look at the "not in the
truth set" list: a lower cutoff also lets through sequencing and PCR errors.

**Discuss:** how would a clinical lab decide whether a 2% call is real? Why does
the FLT3 ITD still not appear?

### Round 3 — VarDict for the FLT3 ITD

The FLT3 ITD inserts 279 bases, about twice the read length. No read contains
the whole insertion, so bwa cannot write it as an insertion. Instead it
soft-clips the part of each read that doesn't match the genome. mpileup only
counts aligned bases, so VarScan never sees the ITD at any cutoff.

> **Look inside:** how many reads near the ITD are soft-clipped (an `S` in the CIGAR)?
>
> ```
> samtools view $BAM chr13:28033800-28034200 | awk '$6 ~ /[0-9]+S/' | wc -l
> ```
>
> Compare with all reads there: `samtools view -c $BAM chr13:28033800-28034200`

VarDict assembles soft-clipped reads into the longer sequence they came from,
so it can call the ITD. Its output goes through a strand-bias test
(`teststrandbias.R`) and is converted to VCF (`var2vcf_valid.pl`):

```
vardict-java -G $REF -N $SAMPLE -b $BAM -f 0.02 -c 1 -S 2 -E 3 -g 4 -th $THREADS $BED \
  | teststrandbias.R \
  | var2vcf_valid.pl -N $SAMPLE -E -f 0.02 > r3.vcf
```
```
bgzip -f r3.vcf && tabix -f -p vcf r3.vcf.gz
```
```
annotate-vcf r3.vcf.gz r3
```
```
benchmark-vcf r3.vcf.gz
```

> **Look inside:** the ITD call (the long ALT allele is cut short here):
>
> ```
> bcftools view -H -r chr13:28033500-28034300 r3.vcf.gz | cut -f1-5 | cut -c1-120
> ```

**What to expect:** the FLT3 ITD is found. Its observed VAF is usually lower
than expected, because reads that only partly overlap the insertion are hard
to count.

### Compare the three rounds

One report with all three rounds side by side:

```
benchmark-vcf r1.vcf.gz r2.vcf.gz r3.vcf.gz --html qc/benchmark.html
```

**View:** portal → QC reports → `benchmark.html`. **NEW** marks the round
that first found each variant. Click any position to open it in IGV.

---

## 8. Look at the evidence in IGV

Open portal → **IGV**. Click files on the left to add them as tracks: start with
`ngs_ws.realigned.bam`, `r1.vcf.gz`, `r2.vcf.gz`, `r3.vcf.gz` and
`data/ngs_ws_roi_sorted.bed`. Use the **Jump to** buttons to visit each variant.

1. **KRAS G12D (4%)**: find the few reads carrying the T. Would you trust this
   call? Compare with KRAS G13D (13%) three bases away.
2. **EGFR L858R (3%)**: which VCF tracks show a call here?
3. **FBXW7 frameshift**: add `ngs_ws.dedup.bam` and compare it with the
   realigned BAM (step 6).
4. **FLT3 ITD**: tick **Show soft-clipped bases**, then add the BAM again.
   Coloured soft-clipped ends pile up at the same place: the reads trying to
   cover the inserted sequence.
5. **A call from the "not in the truth set" list** in round 2 (click it in the
   benchmark report): what does the evidence look like compared with a real
   variant?

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `command not found` | Log out and back in, or run `source /etc/profile.d/ngs.sh; source ~/workshop.env` |
| `No such file or directory` | Run `cd ~/work`, then `ls` to check the previous step's output exists |
| A pasted command does nothing or shows `>` | The paste was incomplete. Press <kbd>Ctrl</kbd>+<kbd>C</kbd> and paste the box again |
| The connection dropped | Log in again and `cd ~/work`. Finished steps keep their files, so re-run the step that was interrupted |
| A file is missing in IGV | Click **Refresh list**. BAMs need a `.bai` and VCFs a `.tbi` next to them (`samtools index`, `tabix -p vcf`) |
| Browser asks for a password | Same username and password as SSH |
