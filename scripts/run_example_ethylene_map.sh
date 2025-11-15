#!/usr/bin/env bash
set -euo pipefail
TOP=${TOP:-}
TRJ=${TRJ:-}
MAP=${MAP:-example/ethylene/mapping_two_beads.xml}
if [ -z "${TOP}" ] || [ -z "${TRJ}" ]; then
  echo "usage: TOP=topol.tpr TRJ=traj.xtc MAP=${MAP} $0"
  exit 1
fi
if command -v csg_map >/dev/null 2>&1; then
  csg_map --top "${TOP}" --trj "${TRJ}" --cg "${MAP}"
elif command -v conda >/dev/null 2>&1; then
  conda run -n cg bash -lc "csg_map --top '${TOP}' --trj '${TRJ}' --cg '${MAP}'"
else
  echo "csg_map not found"
  exit 1
fi