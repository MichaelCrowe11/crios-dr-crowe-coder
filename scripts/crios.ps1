# CriOS PowerShell Entry Point
# Usage: .\scripts\crios.ps1 <crios args...>

param([Parameter(ValueFromRemainingArguments=$true)][string[]]$Arguments)

# Force UTF-8 to avoid encoding crashes on ✓/• etc.
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
# Try to set ANSI rendering if available (PS 7+)
if ($PSVersionTable.PSVersion.Major -ge 7) {
    try { $PSStyle.OutputRendering = 'ANSI' } catch {}
}

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectRoot = Split-Path -Parent $ScriptDir

# Try local editable install first (python -m)
$py = Get-Command python -ErrorAction SilentlyContinue
if (-not $py) { Write-Error "Python not found in PATH"; exit 1 }

Push-Location $ProjectRoot
try {
  python -m src.cli "${Arguments}"
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} finally {
  Pop-Location
}