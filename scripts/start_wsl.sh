#!/usr/bin/env bash
set -euo pipefail
if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then . ~/miniconda3/etc/profile.d/conda.sh; elif [ -f ~/anaconda3/etc/profile.d/conda.sh ]; then . ~/anaconda3/etc/profile.d/conda.sh; elif [ -f /opt/conda/etc/profile.d/conda.sh ]; then . /opt/conda/etc/profile.d/conda.sh; fi
conda activate cg
ROOT=$(cd "$(dirname "$0")"/.. && pwd)
python -m pip install -r "$ROOT/backend/requirements.txt"
cd "$ROOT/frontend" && npm install
cd "$ROOT"
python -m uvicorn backend.app:app --app-dir "$ROOT" --host 0.0.0.0 --port 8000 &
sleep 1
node "$ROOT/frontend/server.js"
