# Workshop VM setup

One DigitalOcean droplet per participant, all created from a single pre-built
snapshot. Participants run every pipeline step by pasting commands into an SSH
terminal. Each VM also serves a small password-protected web portal at
`http://<ip>/` for viewing the HTML QC reports and opening BAM, VCF and BED
files in igv.js, so nobody needs `scp` or a desktop IGV install.

The old scripts (NFS share, per-VM conda install, ANNOVAR) are in
`archived_scripts/` for reference only.

## What you need on your laptop

- [dovm](https://github.com/roysomak4/do-vm-manager) and `doctl`, both logged
  in to the DigitalOcean account
- `rsync`, `ssh`, `curl`
- Your indexed hg38 reference, as a directory or `.tar.gz`: one uncompressed
  FASTA with `chr` contig names, its `.fai`, and the bwa index (`.amb .ann .bwt .pac .sa`)
- A password file for the participants, e.g. `pw.txt` (git-ignored). It is
  both the SSH password and the web portal password.

Run every command below from this directory (`instructors/scripts`). dovm
writes droplet IPs to `output/`, which the other scripts read.

## 1. Build the image (once per course year, about an hour)

```sh
dovm create dovm/prep.yaml
ip=$(awk '$1 == "prep-vm" { print $2 }' output/vms_ip_host.txt)
./upload_and_build.sh "$ip" pw.txt ~/projects/assets/genomes/hg38/bwa-mem/bwa-mem.tar.gz
```

The reference can be a directory or a `.tar.gz` with one top-level folder. A
tarball is uploaded as is, checked against `<tarball>.sha256` if that file
exists, and unpacked into `/opt/assets/hg38` on the droplet. On re-runs the
upload is skipped once the reference is unpacked.

`upload_and_build.sh` uploads the reference (resumable), the build bundle and
the workshop data, then runs `build/build_image.sh` on the droplet:

| Step | Result on the image |
|---|---|
| Check the reference | `/opt/assets/hg38`: FASTA, `.fai`, bwa index present; `chr` names; warns if ALT contigs lack an `.alt` file |
| System packages | Caddy, jq, tree; full `apt upgrade` |
| `HOME_MODE 0755` | Participant home directories are readable by the web server |
| Conda environment | `/opt/miniforge/envs/ngs`, including fastVEP, from `build/environment.yml`, or from `build/ngs-workshop.lock.txt` once it exists |
| fastVEP data | Ensembl GFF3 (`ENSEMBL_RELEASE`, default 115) and a ClinVar database in `/opt/assets/fastvep`; set `BUILD_CLINVAR=0` to skip ClinVar |
| Workshop data | `/etc/skel/data`: FASTQ, target BED, cleaned truth VCF (bgzipped and indexed); `/etc/skel/work/qc` |
| Login environment | `/etc/profile.d/ngs.sh` activates the env; `~/workshop.env` sets `SAMPLE`, `REF`, `BED`, `R1`, `R2`, `TRUTH`, ... |
| Commands | `/opt/ngs/bin`: `annotate-vcf`, `benchmark-vcf` (terminal table per round; `--html` writes a side-by-side report of all rounds), `check-vm` |
| Web portal | Caddy on port 80 with the shared password: `/files/` (`~/work`), `/data/` (`~/data`), `/igv.html` |
| Smoke test | `check-vm --full` runs a 10,000-read-pair version of the whole pipeline |

The first build solves the conda environment and copies
`build/ngs-workshop.lock.txt` back to your laptop. **Commit it**: later builds
install exactly those versions. To change a tool version, edit
`environment.yml`, delete the lock file, and rebuild.

If the SSH connection drops, run the same command again.
Finished steps are skipped. The full log is on the droplet at
`/var/log/ngs-image-build.log`.

## 2. Snapshot it

```sh
./snapshot.sh "$ip" ngs-2026-v1
dovm destroy dovm/prep.yaml
```

`snapshot.sh` clears cloud-init state and the build files,
powers the droplet off, snapshots it, and prints the snapshot ID and minimum
disk size. Put the ID in `dovm/course.yaml` under `image:`.

The snapshot's minimum disk size is the prep droplet's disk. Keep both specs on
the same size, and never resize the prep droplet upward.

## 3. Create the participant VMs (class day or the evening before)

```sh
dovm validate dovm/course.yaml                     # checks the snapshot ID, shows cost
dovm create dovm/course.yaml --login-user bioseq --login-password-file pw.txt
./check_vms.sh pw.txt
```

Password logins start about a minute after dovm reports the droplets ready,
because cloud-init creates the `bioseq` user on first boot. `check_vms.sh` runs
`check-vm` on every VM over SSH and tests the portal, including the HTTP range
requests igv.js needs. Hand out `output/vms_ip_host.txt`, one line per
participant.

## 4. During class

- **Participant access:** `ssh bioseq@<ip>` in Terminal (macOS/Linux) or
  PowerShell (Windows 10/11). Browser: `http://<ip>/`, user `bioseq`, same
  password.
- **A VM breaks:** give the participant `spare01` or `spare02` from
  `output/vms_ip_host.txt`.
- **Re-check one VM:** `ssh root@<ip> check-vm` (add `--full` to rerun the mini pipeline).

## 5. After class

```sh
dovm destroy dovm/course.yaml
```

## Files

```
dovm/prep.yaml                  the image-build droplet
dovm/course.yaml                40 students + 2 spares from the snapshot
upload_and_build.sh             laptop: upload reference, bundle and data; run the build
snapshot.sh                     laptop: finalize, power off, snapshot
check_vms.sh                    laptop: check every participant VM
build/build_image.sh            droplet: builds the image (idempotent)
build/environment.yml           conda environment (pinned to major.minor)
build/ngs-workshop.lock.txt     exact versions, written by the first build; commit it
build/files/ngs-profile.sh      -> /etc/profile.d/ngs.sh
build/files/workshop.env        -> ~/workshop.env for each participant
build/files/Caddyfile           -> /etc/caddy/Caddyfile (user and password hash filled in)
build/files/web/                -> /opt/ngs/web: portal landing page and IGV page
build/files/bin/                -> /opt/ngs/bin: annotate-vcf, benchmark-vcf, check-vm
data/ngs_data/                  workshop FASTQ, target BED, truth VCF
```

## Notes

- **Portal security:** plain HTTP with basic auth, by design. The VMs are
  short-lived, hold only synthetic data, and share one password.
- **IGV genome:** `igv.html` uses igv.js's hosted `hg38` genome, which the
  participant's browser downloads from IGV's servers. Only the participant's
  own BAM, VCF and BED files come from the VM.
- **fastVEP and chromosome names:** fastVEP matches `chr1` against Ensembl's
  `1` on its own. Its docs don't say whether the same applies to supplementary
  databases, so the build renames ClinVar's contigs to `chr` form before
  building that database. The build warns if a test annotation comes back
  without ClinVar fields.
- **Target BED:** FLT3 exons 13–15 overlapped after padding and are merged
  into one interval (`exon-13-15`) in `ngs_ws_roi_sorted.bed`, so no target
  region overlaps another. VarDict therefore reports the ITD once.
