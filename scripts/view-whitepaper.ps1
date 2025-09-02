# Script to view Dr. Crowe Coder White Paper
param(
    [string]$Format = "browser"
)

$whitePaperPath = "$PSScriptRoot\..\docs\DrCroweCoder_WhitePaper.md"

switch ($Format) {
    "browser" {
        # Try to open in default browser/markdown viewer
        Start-Process $whitePaperPath
    }
    "vscode" {
        # Open in VS Code if available
        code $whitePaperPath
    }
    "notepad" {
        # Open in Notepad++
        notepad $whitePaperPath
    }
    "print" {
        # Display first 100 lines in console
        Get-Content $whitePaperPath -Head 100
    }
    default {
        Write-Host "Opening white paper at: $whitePaperPath"
        Start-Process $whitePaperPath
    }
}

Write-Host ""
Write-Host "White Paper Location: $whitePaperPath"
Write-Host "File Size: $((Get-Item $whitePaperPath).Length / 1KB) KB"
Write-Host ""
Write-Host "View options:"
Write-Host "  .\view-whitepaper.ps1           - Open in default app"
Write-Host "  .\view-whitepaper.ps1 -Format vscode    - Open in VS Code"
Write-Host "  .\view-whitepaper.ps1 -Format notepad   - Open in Notepad"
Write-Host "  .\view-whitepaper.ps1 -Format print     - Show in console"