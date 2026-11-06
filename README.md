# NGS Hands-on Bioinformatics Workshop

Materials for a hands-on NGS bioinformatics workshop. Participants take a small
targeted-sequencing dataset (49,642 read pairs across 17 regions in 9 cancer
genes) from raw FASTQ to annotated variant calls, then score each calling
strategy against a known truth set.

The pipeline covers QC (fastqc, fastp), alignment (bwa mem), BAM processing
(samtools, sambamba, ABRA2), coverage QC (mosdepth, MultiQC), three rounds of
variant calling (VarScan at 5% VAF, VarScan at 2% VAF, then VarDict for a
279 bp FLT3 ITD), annotation with fastVEP, and review in IGV.

Each participant gets their own DigitalOcean VM, built from a single pre-built
snapshot. Every tool, the hg38 reference and the dataset are already installed.
Participants run commands over SSH and view reports and IGV in a web portal on
their VM.

## For participants

Everything you need is in the participant guide:

**[NGS Hands-On Lab guide](https://roysomak4.github.io/ngs_bio_workshop/)**
(source: [`docs/index.html`](docs/index.html))

The guide walks through logging in to your VM, each pipeline step with
copy-ready commands and "look inside" checks, the three variant-calling rounds,
and the IGV review. It ends with a troubleshooting section. Your instructor will
give you your VM's IP address and the password.

All you need is a laptop with a terminal (Terminal on macOS or Linux,
PowerShell on Windows 10/11) and a web browser. You don't need to install
anything.

## For instructors

Use these docs to run the workshop yourself or adapt its content.

| Doc | What it covers |
|---|---|
| [`instructors/scripts/README.md`](instructors/scripts/README.md) | Building the VM image, provisioning participant VMs with [dovm](https://github.com/roysomak4/do-vm-manager), checking them before class, and teardown |
| [`instructors/testing_the_guide.html`](instructors/testing_the_guide.html) | Previewing the participant guide locally, testing it end to end on a VM, and publishing it with GitHub Pages |
| [`instructors/workshop_deployment_review.html`](instructors/workshop_deployment_review.html) | Design review: pipeline rationale, expected results per round, VM and portal design, image contents |
| [`instructors/bioinformatics_pipeline_and_software_list.md`](instructors/bioinformatics_pipeline_and_software_list.md) | The workshop pipeline and software list |

The `.html` docs show as source on GitHub. Clone the repo and open them in a
browser to read them.

### Repo layout

```
docs/index.html                 participant guide (published with GitHub Pages)
instructors/scripts/
├── dovm/                       VM spec templates: prep VM, course VMs, single-user test VM
├── upload_and_build.sh         builds the image on the prep VM
├── snapshot.sh                 snapshots the prep VM
├── check_vms.sh                checks every participant VM before class
├── build/                      build script, conda environment and lock file,
│                               files baked into the image (portal, helpers, profile)
└── data/ngs_data/              FASTQ, target BED, truth VCF
```

### Quick start

From `instructors/scripts`, with `dovm` and `doctl` logged in to DigitalOcean:

```sh
# First time: copy the templates (don't edit them) and add your SSH key ID(s)
cp dovm/prep_template.yaml dovm/prep.yaml
cp dovm/course_template.yaml dovm/course.yaml

# Once per course: build and snapshot the image
dovm create dovm/prep.yaml
./upload_and_build.sh <prep-ip> pw.txt /path/to/hg38[.tar.gz] [/path/to/Homo_sapiens.GRCh38.<release>.gff3.gz]
./snapshot.sh <prep-ip> [snapshot-name]   # prints the snapshot ID
dovm destroy dovm/prep.yaml

# Class day: put the snapshot ID in dovm/course.yaml under image:, then
dovm validate dovm/course.yaml            # checks the snapshot ID, shows cost
dovm create dovm/course.yaml --login-user bioseq --login-password-file pw.txt
./check_vms.sh pw.txt

# After class
dovm destroy dovm/course.yaml
```

For prerequisites and the details of each step, see
[`instructors/scripts/README.md`](instructors/scripts/README.md).

## License

[MIT](LICENSE)
