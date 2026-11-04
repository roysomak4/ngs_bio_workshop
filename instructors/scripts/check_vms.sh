#!/usr/bin/env bash
# Check every participant VM: run check-vm over SSH and test the web portal with
# the participants' password, including the range requests igv.js needs.
#
#   ./check_vms.sh <password-file> [ip-list]
#
# ip-list defaults to output/vm_ips.txt, written by dovm create. Per-VM logs go
# to output/checks/<ip>.log.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
pwfile=${1:?usage: ./check_vms.sh <password-file> [ip-list]}
iplist=${2:-$here/output/vm_ips.txt}
user=${LOGIN_USER:-bioseq}
pw=$(tr -d '\r\n' <"$pwfile")
logdir=$here/output/checks
mkdir -p "$logdir"
ssh_opts=(-o StrictHostKeyChecking=accept-new -o UserKnownHostsFile="$here/output/known_hosts"
  -o ConnectTimeout=15 -o BatchMode=yes)

hostname_for() {
  awk -v ip="$1" '$2 == ip { print $1 }' "$here/output/vms_ip_host.txt" 2>/dev/null
}

check_one() {
  local ip=$1 vm portal range
  if ssh "${ssh_opts[@]}" "root@$ip" /opt/ngs/bin/check-vm >"$logdir/$ip.log" 2>&1; then vm=PASS; else vm=FAIL; fi
  portal=$(curl -s -m 10 -o /dev/null -w '%{http_code}' -u "$user:$pw" "http://$ip/files/") || true
  range=$(curl -s -m 10 -o /dev/null -w '%{http_code}' -u "$user:$pw" -H 'Range: bytes=0-99' \
    "http://$ip/data/ngs_ws_roi_sorted.bed") || true
  [[ $portal == 200 ]] && portal=PASS || portal="FAIL($portal)"
  [[ $range == 206 ]] && range=PASS || range="FAIL($range)"
  printf '%-12s %-16s %-6s %-10s %s\n' "$(hostname_for "$ip")" "$ip" "$vm" "$portal" "$range"
}

results=$(mktemp)
while read -r ip; do
  [[ -n $ip ]] && check_one "$ip" >>"$results" &
done <"$iplist"
wait

printf '%-12s %-16s %-6s %-10s %s\n' HOST IP VM PORTAL RANGE
sort "$results"
fails=$(grep -c FAIL "$results" || true)
rm -f "$results"
if ((fails)); then
  echo "$fails VM(s) failed. See $logdir/<ip>.log; hand those participants a spare."
  exit 1
fi
echo "All VMs passed."
