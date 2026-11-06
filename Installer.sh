#!/bin/bash
set -e
# --- XILEHXILEH SUITE INSTALLER ---
echo "-------------------------------------------------------"
echo "🧬 Initializing BSAP Setup"
echo "-------------------------------------------------------"

# 1. DETECT SYSTEM COMPILER (Required for forensic C-extensions)
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    if ! command -v gcc >/dev/null 2>&1; then
        echo "[!] WARNING: gcc not found. Attempting to install build-essential..."
        sudo apt-get update && sudo apt-get install -y build-essential
    fi
fi

# 2. DETECT PYTHON EXECUTABLE
if command -v python3 >/dev/null 2>&1; then
    PYTHON_EXE="python3"
elif command -v python >/dev/null 2>&1; then
    PYTHON_EXE="python"
else
    echo "[-] ERROR: Python is not installed. Please install Python 3.8+."
    exit 1
fi

# 3. VERIFY VERSION (3.10+)
if ! $PYTHON_EXE -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
    echo "[-] ERROR: BSAP requires Python 3.10+."
    exit 1
fi
echo "[+] Detected $($PYTHON_EXE --version) using command '$PYTHON_EXE'"

# 4. CREATE VIRTUALENV & INSTALL DEPENDENCIES (single source: Requirements.txt)
if ! $PYTHON_EXE -m venv .venv 2>/dev/null; then
    echo "[!] python3-venv missing. Installing..."
    sudo apt-get update && sudo apt-get install -y python3-venv
    $PYTHON_EXE -m venv .venv
fi
PYTHON_EXE=".venv/bin/python"
$PYTHON_EXE -m pip install --upgrade pip
$PYTHON_EXE -m pip install -r Requirements.txt

# 5. INITIALIZE ENVIRONMENT
if [ -f "Setup_env.py" ]; then
    echo "[*] Running environment scaffolding..."
    $PYTHON_EXE Setup_env.py 
fi

echo "-------------------------------------------------------"
echo "[!] SETUP COMPLETE"
echo "[*] Dashboard:  source .venv/bin/activate && streamlit run APP.py"
echo "[*] Headless:   source .venv/bin/activate && python Alpha4.py"
echo "[*] Website:    python -m http.server 8765  →  http://localhost:8765/index.html"
echo "-------------------------------------------------------"