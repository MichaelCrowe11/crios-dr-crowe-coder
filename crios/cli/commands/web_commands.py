"""
CriOS CLI - Web Commands
Web interface and server management
"""

import logging
from typing import Optional

import typer
from rich.console import Console

console = Console()

# Create sub-app for web commands
app = typer.Typer(name="web", help="🌐 Web interface commands")


@app.command("serve")
def serve_web(
    host: str = typer.Option("0.0.0.0", "--host", help="Host to bind to"),
    port: int = typer.Option(8000, "--port", "-p", help="Port to serve on"),
    reload: bool = typer.Option(False, "--reload", help="Auto-reload on changes"),
    debug: bool = typer.Option(False, "--debug", help="Enable debug mode")
):
    """🚀 Start the CriOS web interface"""
    try:
        import uvicorn
        from ...web.main import app as web_app
        
        console.print(f"🚀 Starting CriOS web interface on {host}:{port}")
        console.print("🌐 Access the interface at http://localhost:8000")
        
        uvicorn.run(
            "crios.web.main:app",
            host=host,
            port=port,
            reload=reload,
            log_level="debug" if debug else "info"
        )
        
    except ImportError:
        console.print("[red]Error:[/red] Web dependencies not installed")
        console.print("Install with: pip install 'crios[web]'")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[red]Error starting web server:[/red] {e}")
        raise typer.Exit(1)


@app.command("ide")  
def launch_ide(
    port: int = typer.Option(3000, "--port", "-p", help="IDE port"),
    open_browser: bool = typer.Option(True, "--open", help="Open browser automatically")
):
    """💻 Launch the immersive IDE"""
    console.print(f"💻 Launching CriOS Immersive IDE on port {port}")
    console.print("🔧 IDE launch functionality coming soon...")