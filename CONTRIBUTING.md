# Contributing to KAMEL

Thank you for your interest in contributing to KAMEL! This is an academic research project maintained by the Computational Physics group at TU Graz.

## Reporting Issues

If you encounter bugs or have feature requests, please open a [GitHub Issue](https://github.com/itpplasma/KAMEL/issues). Include:

- A clear description of the problem or suggestion
- Steps to reproduce (for bugs)
- Your platform and compiler versions

## Proposing Changes

1. Fork the repository
2. Create a feature branch from `main`
3. Make your changes
4. Ensure the code builds and tests pass:
   ```bash
   make all
   make test
   ```
5. Open a Pull Request against `main`

## Code Style

This project uses automated formatting tools. Please install the pre-commit hooks before committing:

```bash
pip install pre-commit
pre-commit install
```

The following formatters are enforced:
- **Fortran**: fprettify
- **C/C++**: clang-format
- **Python**: black, isort

## Build and Test

See the [README](README.md) for build instructions and dependencies.

## Questions

For questions about the physics or the code, feel free to open a GitHub Issue or contact us at plasma.itp@tugraz.at.

## Note

KAMEL is primarily a research tool. We welcome contributions but maintainer bandwidth is limited. Please be patient with reviews and responses.
