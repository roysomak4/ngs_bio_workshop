#!/usr/bin/env bash
# Upload the build bundle to the prep droplet and run build_image.sh there.
#
#   ./upload_and_build.sh <prep-ip> <password-file> [hg38-dir | hg38.tar.gz] [gff3]
#
# password-file  the participants' password (the same file you pass to dovm's
#                --login-password-file); it protects the web portal
# hg38           the indexed reference (FASTA, .fai, bwa index), either as a
#                local directory or as a .tar.gz holding one top-level folder.
#                A tarball is checked against <tarball>.sha256 if present and
#                unpacked on the droplet. Omit it if the reference is already there,
#                or pass "" to skip it while still giving a GFF3.
# gff3           optional local Ensembl GFF3, named Homo_sapiens.GRCh38.<release>.gff3.gz.
#                Uploaded instead of downloading it from the Ensembl FTP on the
#                droplet, which can take 40 minutes. The build uses its release.
#
# Safe to re-run: uploads resume and build_image.sh skips finished steps.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ip=${1:?usage: ./upload_and_build.sh <prep-ip> <password-file> [hg38-dir or .tar.gz] [gff3]}
pwfile=${2:?usage: ./upload_and_build.sh <prep-ip> <password-file> [hg38-dir or .tar.gz] [gff3]}
refdir=${3:-}
gff3=${4:-}
data=$here/data/ngs_data

[[ -s $pwfile ]] || { echo "password file $pwfile is empty or missing" >&2; exit 1; }
build_env=
if [[ -n $gff3 ]]; then
  [[ -s $gff3 ]] || { echo "GFF3 file $gff3 is empty or missing" >&2; exit 1; }
  [[ $(basename "$gff3") =~ ^Homo_sapiens\.GRCh38\.([0-9]+)\.gff3\.gz$ ]] ||
    { echo "GFF3 must be named Homo_sapiens.GRCh38.<release>.gff3.gz, as on the Ensembl FTP" >&2; exit 1; }
  build_env="ENSEMBL_RELEASE=${BASH_REMATCH[1]}"
  gzip -t "$gff3" || { echo "GFF3 file $gff3 is not a valid gzip file" >&2; exit 1; }
fi
mkdir -p "$here/output"
ssh_opts=(-o StrictHostKeyChecking=accept-new -o UserKnownHostsFile="$here/output/known_hosts"
  -o ServerAliveInterval=30)
ssh_cmd="ssh ${ssh_opts[*]}"

if [[ $refdir == *.tar.gz ]]; then
  if ssh "${ssh_opts[@]}" "root@$ip" 'compgen -G "/opt/assets/hg38/*.fa.bwt" >/dev/null'; then
    echo "==> Reference already unpacked on the droplet; skipping upload"
  else
    echo "==> Uploading reference tarball $refdir (resumable)"
    ssh "${ssh_opts[@]}" "root@$ip" mkdir -p /opt/assets/hg38
    rsync -avP --partial -e "$ssh_cmd" "$refdir" "root@$ip:/opt/assets/"
    if [[ -s $refdir.sha256 ]]; then
      rsync -e "$ssh_cmd" "$refdir.sha256" "root@$ip:/opt/assets/"
      ssh "${ssh_opts[@]}" "root@$ip" "cd /opt/assets && sha256sum -c $(basename "$refdir").sha256"
    fi
    # A fresh droplet may still hold the apt lock for unattended-upgrades.
    echo "==> Installing pigz"
    ssh "${ssh_opts[@]}" "root@$ip" \
      "command -v pigz >/dev/null || { apt-get update -q -o DPkg::Lock::Timeout=600 &&
       DEBIAN_FRONTEND=noninteractive apt-get install -yq -o DPkg::Lock::Timeout=600 --no-install-recommends pigz; }"
    echo "==> Unpacking into /opt/assets/hg38"
    ssh "${ssh_opts[@]}" "root@$ip" \
      "tar -I pigz -xf /opt/assets/$(basename "$refdir") -C /opt/assets/hg38 --strip-components=1 --no-same-owner &&
       rm -f /opt/assets/$(basename "$refdir") /opt/assets/$(basename "$refdir").sha256"
  fi
elif [[ -n $refdir ]]; then
  echo "==> Uploading reference from $refdir (resumable)"
  ssh "${ssh_opts[@]}" "root@$ip" mkdir -p /opt/assets/hg38
  rsync -avP --partial -e "$ssh_cmd" "${refdir%/}/" "root@$ip:/opt/assets/hg38/"
fi

if [[ -n $gff3 ]]; then
  dest=/opt/assets/fastvep/$(basename "$gff3")
  if ssh "${ssh_opts[@]}" "root@$ip" "test -s $dest"; then
    echo "==> GFF3 already on the droplet; skipping upload"
  else
    # Upload under .part so an interrupted copy is never mistaken for the
    # finished file, which build_image.sh would then skip downloading.
    echo "==> Uploading GFF3 $gff3 (resumable)"
    ssh "${ssh_opts[@]}" "root@$ip" mkdir -p /opt/assets/fastvep
    rsync -avP --partial -e "$ssh_cmd" "$gff3" "root@$ip:$dest.part"
    ssh "${ssh_opts[@]}" "root@$ip" "gzip -t $dest.part && mv $dest.part $dest"
  fi
fi

echo "==> Uploading build bundle and workshop data"
rsync -av --delete -e "$ssh_cmd" --exclude data/ "$here/build/" "root@$ip:/root/ngs-build/"
rsync -av -e "$ssh_cmd" \
  "$data/ngs_ws_R1.fastq.gz" "$data/ngs_ws_R2.fastq.gz" \
  "$data/ngs_ws_roi_sorted.bed" "$data/ngs_ws_truth_vars.vcf" \
  "root@$ip:/root/ngs-build/data/"
# Not rsync --chmod: macOS ships openrsync, which rejects numeric modes.
ssh "${ssh_opts[@]}" "root@$ip" 'umask 077 && cat >/root/ngs-build/portal_password' <"$pwfile"

echo "==> Running build_image.sh (log on the droplet: /var/log/ngs-image-build.log)"
ssh "${ssh_opts[@]}" "root@$ip" "chmod +x /root/ngs-build/build_image.sh && $build_env /root/ngs-build/build_image.sh"

echo "==> Copying the conda lock file back"
rsync -e "$ssh_cmd" "root@$ip:/root/ngs-build/ngs-workshop.lock.txt" "$here/build/"
echo "Commit build/ngs-workshop.lock.txt so future builds install the same versions."
