#!/usr/bin/env bash
# BSAP — Docker helper (optional)
# Primary local entry:  streamlit run APP.py
# Headless audit:       python Alpha4.py
# Website only:         python -m http.server 8765  →  /index.html
set -euo pipefail
cd "$(dirname "$0")"

echo "============================================================"
echo "[*] BSAP Docker launcher"
echo "    UI:      http://localhost:8501  (streamlit / APP.py)"
echo "    Website: serve index.html from project root after an audit"
echo "============================================================"

if command -v docker-compose >/dev/null 2>&1; then
  docker-compose up --build
elif command -v docker >/dev/null 2>&1; then
  docker compose up --build
else
  echo "[-] Docker not found. Local run:"
  echo "    streamlit run APP.py"
  echo "    python Alpha4.py"
  exit 1
fi
