def test_package_import_and_cli_help() -> None:
    import kim
    from kim.cli import app
    from typer.testing import CliRunner

    assert kim.__version__ == "0.1.0"

    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "Run and inspect KIM plasma simulations." in result.stdout
