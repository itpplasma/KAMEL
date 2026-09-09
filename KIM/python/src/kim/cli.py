"""Command-line interface for KIM."""

import typer

app = typer.Typer(
    help="Run and inspect KIM plasma simulations.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Run and inspect KIM plasma simulations."""
