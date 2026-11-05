#!/usr/bin/env bash
# Builds the NGS workshop image on a fresh Ubuntu 24.04 prep droplet.
#
# Run as root from /root/ngs-build (upload_and_build.sh copies this directory
# there and runs it). The indexed hg38 reference must already be in
# /opt/assets/hg38. Every step checks what is already in place, so the script
# can be re-run after a failure or a dropped SSH connection.
#
#   ./build_image.sh             build or update the image contents
#   ./build_image.sh --finalize  clean up just before the snapshot (snapshot.sh runs this)
set -euo pipefail

BUILD_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOGIN_USER=${LOGIN_USER:-bioseq}
REF_DIR=/opt/assets/hg38
FASTVEP_DIR=/opt/assets/fastvep
ENSEMBL_RELEASE=${ENSEMBL_RELEASE:-115}
BUILD_CLINVAR=${BUILD_CLINVAR:-1}
CONDA=/opt/miniforge/bin/conda
ENV_BIN=/opt/miniforge/envs/ngs/bin
LOG=/var/log/ngs-image-build.log

step() { printf '\n==> %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

[[ $EUID -eq 0 ]] || die "run as root"

finalize() {
  step "Finalizing for snapshot"
  # cloud-init skips a login user that already exists, so a baked-in account
  # would keep its empty home on every student VM.
  if id "$LOGIN_USER" >/dev/null 2>&1; then
    die "user $LOGIN_USER exists on this droplet; remove it first: userdel -r $LOGIN_USER"
  fi
  rm -rf /tmp/* "$BUILD_DIR"
  "$CONDA" clean -afy >/dev/null
  apt-get clean
  rm -f /root/.bash_history
  # A fresh machine-id and cloud-init state make each droplet from the snapshot
  # run first-boot setup (dovm's login user) as a new instance.
  cloud-init clean --logs --machine-id
  echo "Done. Power off and snapshot now; do not reboot this droplet first."
}

if [[ ${1:-} == --finalize ]]; then
  finalize
  exit 0
fi

exec > >(tee -a "$LOG") 2>&1
echo "Build started $(date -u +%FT%TZ)"

find_reference() {
  step "Checking the uploaded hg38 reference in $REF_DIR"
  shopt -s nullglob
  local fastas=("$REF_DIR"/*.fa "$REF_DIR"/*.fasta "$REF_DIR"/*.fna)
  shopt -u nullglob
  [[ ${#fastas[@]} -eq 1 ]] ||
    die "expected exactly one uncompressed FASTA (.fa/.fasta/.fna) in $REF_DIR, found ${#fastas[@]}"
  REF=${fastas[0]}
  [[ -s $REF.fai ]] || die "missing $REF.fai (samtools faidx)"
  local ext
  for ext in amb ann bwt pac sa; do
    [[ -s $REF.$ext ]] || die "missing bwa index file $REF.$ext"
  done
  grep -qP '^chr1\t' "$REF.fai" ||
    die "contigs are not chr-prefixed; the workshop BED and truth VCF use chr names"
  if grep -qP '_alt\t' "$REF.fai" && [[ ! -s $REF.alt ]]; then
    warn "FASTA has ALT contigs but no $REF.alt; reads in ALT-covered regions will get MAPQ 0"
  fi
  echo "Reference: $REF"
}

install_packages() {
  step "Installing system packages"
  export DEBIAN_FRONTEND=noninteractive
  apt-get update -q
  apt-get upgrade -yq -o Dpkg::Options::=--force-confold
  apt-get install -yq --no-install-recommends caddy jq tree unzip curl ca-certificates nano less
}

set_home_mode() {
  step "Making new home directories readable by the web portal"
  # cloud-init creates the login user with useradd at first boot, which reads
  # HOME_MODE. The default 0750 would stop Caddy from reading ~/work.
  if grep -qE '^#?\s*HOME_MODE' /etc/login.defs; then
    sed -i -E 's/^#?\s*HOME_MODE.*/HOME_MODE\t0755/' /etc/login.defs
  else
    printf 'HOME_MODE\t0755\n' >>/etc/login.defs
  fi
}

install_conda_env() {
  step "Installing Miniforge and the ngs environment"
  if [[ ! -x $CONDA ]]; then
    curl -fsSL -o /tmp/miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash /tmp/miniforge.sh -b -p /opt/miniforge
    rm /tmp/miniforge.sh
  fi
  local lock=$BUILD_DIR/ngs-workshop.lock.txt
  if [[ -d /opt/miniforge/envs/ngs ]]; then
    echo "Environment exists; delete /opt/miniforge/envs/ngs to rebuild it"
  elif [[ -s $lock ]]; then
    echo "Installing from $lock"
    "$CONDA" create -y -n ngs --file "$lock"
  else
    echo "No lock file; solving environment.yml"
    "$CONDA" env create -n ngs -f "$BUILD_DIR/environment.yml"
  fi
  mkdir -p /opt/ngs
  "$CONDA" list -n ngs --explicit --md5 >"$lock"
  cp "$lock" /opt/ngs/ngs-workshop.lock.txt
}

install_fastvep_data() {
  step "Preparing fastVEP GRCh38 data in $FASTVEP_DIR"
  mkdir -p "$FASTVEP_DIR/sa"
  GFF3=$FASTVEP_DIR/Homo_sapiens.GRCh38.$ENSEMBL_RELEASE.gff3.gz
  if [[ ! -s $GFF3 ]]; then
    curl -fL -o "$GFF3.part" \
      "https://ftp.ensembl.org/pub/release-$ENSEMBL_RELEASE/gff3/homo_sapiens/$(basename "$GFF3")"
    mv "$GFF3.part" "$GFF3"
  fi

  if [[ $BUILD_CLINVAR == 1 ]] && ! compgen -G "$FASTVEP_DIR/sa/clinvar*" >/dev/null; then
    local tmp
    tmp=$(mktemp -d)
    curl -fL -o "$tmp/clinvar.vcf.gz" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    # ClinVar uses bare contig names. fastVEP documents chr aliasing for gene
    # models but not for supplementary databases, so match the workshop's chr names.
    # Renamed as text: ClinVar's header has no ##contig lines, so bcftools
    # annotate --rename-chrs has nothing to rename and stops at the first record.
    zcat "$tmp/clinvar.vcf.gz" |
      awk 'BEGIN { FS = OFS = "\t" }
        /^#/ { print; next }
        $1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/ { $1 = "chr" $1 }
        $1 == "MT" { $1 = "chrM" }
        { print }' |
      "$ENV_BIN/bgzip" -c >"$tmp/clinvar.chr.vcf.gz"
    "$ENV_BIN/fastvep" sa-build --source clinvar -i "$tmp/clinvar.chr.vcf.gz" \
      -o "$FASTVEP_DIR/sa/clinvar" --assembly GRCh38
    "$ENV_BIN/bcftools" view -h "$tmp/clinvar.vcf.gz" | grep -m1 '^##fileDate' >"$FASTVEP_DIR/sa/clinvar.release.txt" || true
    rm -rf "$tmp"
  fi
}

install_workshop_data() {
  step "Installing workshop data into /etc/skel"
  local src=$BUILD_DIR/data skel=/etc/skel
  mkdir -p "$skel/data" "$skel/work/qc"
  install -m 644 "$src/ngs_ws_R1.fastq.gz" "$src/ngs_ws_R2.fastq.gz" "$src/ngs_ws_roi_sorted.bed" "$skel/data/"

  # Truth VCF: keep the eight site columns (the sample column holds placeholder
  # values that trigger warnings), fix the swapped INFO types, fill GENE from the BED.
  local tmp
  tmp=$(mktemp -d)
  awk 'BEGIN { FS = OFS = "\t" }
    /^##FORMAT/ { next }
    /^##INFO=<ID=HIAF,/ { print "##INFO=<ID=HIAF,Number=1,Type=Float,Description=\"Expected variant allele fraction\">"; next }
    /^##INFO=<ID=GENE,/ { print "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">"; next }
    /^##/ { sub(/[ \t]+$/, ""); print; next }
    { print $1, $2, $3, $4, $5, $6, $7, $8 }' "$src/ngs_ws_truth_vars.vcf" >"$tmp/truth.vcf"
  awk 'BEGIN { OFS = "\t" } { split($4, a, ";"); print $1, $2, $3, a[1] }' "$src/ngs_ws_roi_sorted.bed" |
    sort -k1,1 -k2,2n | "$ENV_BIN/bgzip" -c >"$tmp/genes.bed.gz"
  "$ENV_BIN/tabix" -p bed "$tmp/genes.bed.gz"
  "$ENV_BIN/bcftools" annotate -a "$tmp/genes.bed.gz" -c CHROM,FROM,TO,INFO/GENE "$tmp/truth.vcf" |
    "$ENV_BIN/bcftools" sort -Oz -o "$skel/data/ngs_ws_truth_vars.vcf.gz"
  "$ENV_BIN/bcftools" index -tf "$skel/data/ngs_ws_truth_vars.vcf.gz"
  rm -rf "$tmp"

  install -m 644 "$BUILD_DIR/files/workshop.env" "$skel/workshop.env"
  grep -q 'workshop.env' "$skel/.bashrc" ||
    printf '\n# NGS workshop variables (SAMPLE, REF, BED, ...)\n[ -f ~/workshop.env ] && . ~/workshop.env\n' >>"$skel/.bashrc"
  chmod -R a+rX "$skel/data" "$skel/work"
}

install_profile_and_tools() {
  step "Installing login profile and workshop commands"
  sed -e "s|@REF@|$REF|" -e "s|@GFF3@|$GFF3|" "$BUILD_DIR/files/ngs-profile.sh" >/etc/profile.d/ngs.sh
  chmod 644 /etc/profile.d/ngs.sh
  mkdir -p /opt/ngs/bin
  install -m 755 "$BUILD_DIR"/files/bin/* /opt/ngs/bin/
}

warm_fastvep_cache() {
  step "Building fastVEP's transcript cache and testing annotation"
  # The first run writes a sidecar cache next to the GFF3. Do it now as root,
  # then confirm an unprivileged user can annotate with the read-only data.
  local tmp out
  tmp=$(mktemp -d)
  chown nobody "$tmp"
  "$ENV_BIN/fastvep" annotate -i /etc/skel/data/ngs_ws_truth_vars.vcf.gz -o "$tmp/root.vcf" \
    --gff3 "$GFF3" --fasta "$REF" --hgvs
  out=$(runuser -u nobody -- bash -c ". /etc/profile.d/ngs.sh && cd $tmp && annotate-vcf /etc/skel/data/ngs_ws_truth_vars.vcf.gz user 2>&1") ||
    die "annotate-vcf failed as an unprivileged user:"$'\n'"$out"
  grep -q 'CSQ=' "$tmp/user.fastvep.vcf" || die "fastVEP output has no CSQ annotations"
  if compgen -G "$FASTVEP_DIR/sa/clinvar*" >/dev/null && ! grep -q 'CLINVAR' "$tmp/user.fastvep.vcf"; then
    warn "ClinVar database built but no ClinVar fields in the test output; check contig naming"
  fi
  rm -rf "$tmp"
}

install_portal() {
  step "Configuring the web portal (Caddy)"
  local pwfile=${PORTAL_PASSWORD_FILE:-$BUILD_DIR/portal_password} hash
  mkdir -p /opt/ngs/web
  install -m 644 "$BUILD_DIR"/files/web/* /opt/ngs/web/
  if [[ -s $pwfile ]]; then
    hash=$(caddy hash-password --plaintext "$(tr -d '\r\n' <"$pwfile")")
    sed -e "s|@USER@|$LOGIN_USER|g" -e "s|@HASH@|$hash|" "$BUILD_DIR/files/Caddyfile" >/etc/caddy/Caddyfile
    shred -u "$pwfile"
  elif ! grep -q "$LOGIN_USER" /etc/caddy/Caddyfile 2>/dev/null; then
    die "no portal password: put the participants' password in $pwfile"
  fi
  caddy validate --adapter caddyfile --config /etc/caddy/Caddyfile
  systemctl enable caddy
  systemctl restart caddy
}

find_reference
install_packages
set_home_mode
install_conda_env
install_fastvep_data
install_workshop_data
install_profile_and_tools
warm_fastvep_cache
install_portal

step "Running the full smoke test"
/opt/ngs/bin/check-vm --full

echo
echo "Build finished $(date -u +%FT%TZ). Next: ./snapshot.sh on your laptop."
