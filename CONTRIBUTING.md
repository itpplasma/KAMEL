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

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

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
