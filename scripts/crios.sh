#!/bin/bash
# CriOS Bash Entry Point
# This script provides a convenient way to run CriOS commands from bash

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Find Python command
PYTHON_CMD=""
for cmd in python3 python py; do
    if command -v $cmd &> /dev/null; then
        PYTHON_CMD=$cmd
        break
    fi
done

if [ -z "$PYTHON_CMD" ]; then
    echo -e "${RED}Error: Python is not installed or not in PATH${NC}"
    echo -e "${YELLOW}Please install Python 3.9+ from https://www.python.org/${NC}"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$($PYTHON_CMD --version 2>&1)
if ! echo "$PYTHON_VERSION" | grep -E "Python 3\.(9|1[0-9])" > /dev/null; then
    echo -e "${YELLOW}Warning: Python 3.9+ is recommended. Found: $PYTHON_VERSION${NC}"
fi

# Check for virtual environment
VENV_PATH="$PROJECT_ROOT/venv"
if [ -d "$VENV_PATH" ]; then
    if [ -f "$VENV_PATH/bin/python" ]; then
        PYTHON_CMD="$VENV_PATH/bin/python"
        echo -e "${GREEN}Using virtual environment at: $VENV_PATH${NC}"
    fi
else
    echo -e "${YELLOW}No virtual environment found. Creating one...${NC}"
    $PYTHON_CMD -m venv "$VENV_PATH"
    if [ $? -eq 0 ]; then
        PYTHON_CMD="$VENV_PATH/bin/python"
        echo -e "${GREEN}Virtual environment created successfully${NC}"
        
        # Install dependencies
        echo -e "${YELLOW}Installing dependencies...${NC}"
        REQUIREMENTS_FILE="$PROJECT_ROOT/requirements.txt"
        if [ -f "$REQUIREMENTS_FILE" ]; then
            $PYTHON_CMD -m pip install --upgrade pip --quiet
            $PYTHON_CMD -m pip install -r "$REQUIREMENTS_FILE" --quiet
            if [ $? -eq 0 ]; then
                echo -e "${GREEN}Dependencies installed successfully${NC}"
            else
                echo -e "${YELLOW}Warning: Some dependencies may not have installed correctly${NC}"
            fi
        fi
    else
        echo -e "${YELLOW}Warning: Could not create virtual environment. Using system Python.${NC}"
    fi
fi

# Set CLI path
CLI_PATH="$PROJECT_ROOT/src/cli/commands.py"

if [ ! -f "$CLI_PATH" ]; then
    echo -e "${RED}Error: CriOS CLI not found at $CLI_PATH${NC}"
    exit 1
fi

# Set PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

# Execute the command
if [ $# -eq 0 ]; then
    # Show help if no arguments provided
    $PYTHON_CMD "$CLI_PATH" --help
else
    # Run with provided arguments
    $PYTHON_CMD "$CLI_PATH" "$@"
fi

# Capture and return exit code
EXIT_CODE=$?
exit $EXIT_CODE