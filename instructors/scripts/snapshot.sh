#!/usr/bin/env bash
# Clean up the prep droplet, power it off and snapshot it.
#
#   ./snapshot.sh <prep-ip> [snapshot-name]
#
# Prints the snapshot ID to put in dovm/course.yaml. Needs doctl (logged in).
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ip=${1:?usage: ./snapshot.sh <prep-ip> [snapshot-name]}
name=${2:-ngs-$(date +%Y)-$(date +%m%d)}
ssh_opts=(-o StrictHostKeyChecking=accept-new -o UserKnownHostsFile="$here/output/known_hosts")

id=$(doctl compute droplet list --format ID,PublicIPv4 --no-header | awk -v ip="$ip" '$2 == ip { print $1 }')
[[ -n $id ]] || { echo "no droplet with public IP $ip" >&2; exit 1; }

echo "==> Finalizing droplet $id ($ip)"
ssh "${ssh_opts[@]}" "root@$ip" /root/ngs-build/build_image.sh --finalize

echo "==> Powering off"
doctl compute droplet-action power-off "$id" --wait

echo "==> Taking snapshot $name (this takes several minutes)"
doctl compute droplet-action snapshot "$id" --snapshot-name "$name" --wait

doctl compute snapshot list --resource droplet --format ID,Name,MinDiskSize,Size,Regions |
  awk -v n="$name" 'NR == 1 || $2 == n'
echo
echo "Set 'image:' in dovm/course.yaml to the ID above, then destroy the prep VM:"
echo "  dovm destroy dovm/prep.yaml"
