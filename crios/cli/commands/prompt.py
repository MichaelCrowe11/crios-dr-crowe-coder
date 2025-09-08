import json
import typer
from pathlib import Path

app = typer.Typer(help="Emit orchestration prompt templates (system/user/playbook)")

PROMPTS_DIR = Path(__file__).resolve().parents[3] / "prompts"

FILES = {
    "system": "master_system_prompt.md",
    "developer": "developer_playbook.md",
    "user-templates": "user_prompt_templates.md",
    "examples": "mini_examples.md",
    "readme": "README-ORCHESTRATION.md",
    "skeleton": "response_skeleton.md",
    "acceptance": "acceptance_criteria.md",
}

@app.command()
def emit(
    kind: str = typer.Argument(..., help=f"One of: {', '.join(FILES.keys())}"),
    json_out: bool = typer.Option(False, "--json", help="Emit JSON with name+content"),
):
    """Emit the selected prompt file to stdout."""
    if kind not in FILES:
        raise typer.Exit(code=2)
    path = PROMPTS_DIR / FILES[kind]
    if not path.exists():
        typer.echo(f"Missing file: {path}", err=True)
        raise typer.Exit(code=1)
    content = path.read_text(encoding="utf-8")
    if json_out:
        typer.echo(json.dumps({"name": FILES[kind], "content": content}))
    else:
        typer.echo(content)

@app.command()
def list():  # noqa: A003 - CLI verb
    """List available prompt artifact keys."""
    typer.echo("\n".join(FILES.keys()))


@app.callback()
def main(
    emit: str = typer.Option(None, "--emit", metavar="KEY", help=f"Emit prompt by key ({', '.join(FILES.keys())})"),
    list_keys: bool = typer.Option(False, "--list", help="List available prompt keys and exit"),
    json_out: bool = typer.Option(False, "--json", help="When used with --emit, output JSON structure"),
):
    """Root handler enabling `crios prompt --emit system` usage.

    If no options provided, falls through to subcommands (emit, list).
    """
    if list_keys:
        typer.echo("\n".join(FILES.keys()))
        raise typer.Exit()
    if emit:
        if emit not in FILES:
            typer.echo(f"Invalid key: {emit}. Use --list to see options.", err=True)
            raise typer.Exit(code=2)
        path = PROMPTS_DIR / FILES[emit]
        if not path.exists():
            typer.echo(f"Missing file: {path}", err=True)
            raise typer.Exit(code=1)
        content = path.read_text(encoding="utf-8")
        if json_out:
            typer.echo(json.dumps({"name": FILES[emit], "content": content}))
        else:
            typer.echo(content)
        raise typer.Exit()
