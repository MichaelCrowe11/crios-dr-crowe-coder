# Dr. Crowe Coder - PowerShell launcher
# 194 PhD Agent System

[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
$ErrorActionPreference = "Stop"

# Enable ANSI colors for Windows Terminal/PowerShell 7+
if ($PSVersionTable.PSVersion.Major -ge 7) {
    try { $PSStyle.OutputRendering = 'ANSI' } catch {}
}

# Set environment
$env:PYTHONPATH = "$PSScriptRoot\.."
$env:CRIOS_HOME = "$PSScriptRoot\.."
$env:CRIOS_MODE = "dr-crowe-coder"

# Activate virtual environment if it exists
$venvPath = "$PSScriptRoot\..\venv\Scripts\Activate.ps1"
if (Test-Path $venvPath) {
    & $venvPath
}

# Run Dr. Crowe Coder
python "$PSScriptRoot\..\src\agents\dr_crowe_coder.py" $args