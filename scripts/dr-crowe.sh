#!/bin/bash
# Dr. Crowe Coder - Bash launcher
# 194 PhD Agent System

# Set environment
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export PYTHONPATH="$SCRIPT_DIR/.."
export CRIOS_HOME="$SCRIPT_DIR/.."
export CRIOS_MODE="dr-crowe-coder"

# Activate virtual environment if it exists
if [ -f "$SCRIPT_DIR/../venv/bin/activate" ]; then
    source "$SCRIPT_DIR/../venv/bin/activate"
fi

# Run Dr. Crowe Coder
python "$SCRIPT_DIR/../src/agents/dr_crowe_coder.py" "$@"