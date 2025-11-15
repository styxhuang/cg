#!/usr/bin/env bash
set -euo pipefail
for f in example/ethylene/*.dist.tgt; do
  base="${f%.dist.tgt}"
  awk -f scripts/dist_to_pot.awk "$f" > "${base}.pot"
  echo "generated ${base}.pot"
done
ls -l example/ethylene/*.pot