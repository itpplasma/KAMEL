"""Creation of standalone KIM example cases from packaged resources."""

from __future__ import annotations

import os
import shutil
from importlib import resources
from pathlib import Path
from typing import Iterator

from kim.errors import ConfigurationError


def create_example(destination: Path | str, *, name: str = "periodic") -> Path:
    """Copy a packaged example case into a new directory.

    The destination must not exist, including as a dangling symbolic link.  The
    request returned by this function refers to the copied ``./profiles``
    directory and is intended to be used after changing into the case directory.
    """

    if name != "periodic":
        raise ConfigurationError(f"unknown example: {name!r}")

    target = Path(destination)
    if os.path.lexists(target):
        raise ConfigurationError(f"destination already exists: {target}")

    source_root = resources.files("kim.example_data").joinpath(name)
    created_entries: list[Path] = []
    try:
        target.mkdir(parents=True, exist_ok=False)
        created_entries.append(target)
        for relative_path, source in _resource_files(source_root):
            copied = target / relative_path
            if not copied.parent.exists():
                copied.parent.mkdir(parents=True, exist_ok=False)
                created_entries.append(copied.parent)
            with source.open("rb") as input_file, copied.open("xb") as output_file:
                created_entries.append(copied)
                shutil.copyfileobj(input_file, output_file)
    except BaseException:
        for entry in reversed(created_entries):
            try:
                if entry.is_dir() and not entry.is_symlink():
                    entry.rmdir()
                else:
                    entry.unlink()
            except FileNotFoundError:
                continue
            except OSError:
                # Leave entries that were not created by this call in place.
                continue
        raise

    return target / "request.json"


def _resource_files(
    root: resources.abc.Traversable,
) -> Iterator[tuple[Path, resources.abc.Traversable]]:
    """Yield regular files below a packaged resource directory."""

    def walk(
        resource: resources.abc.Traversable, relative: Path
    ) -> Iterator[tuple[Path, resources.abc.Traversable]]:
        if resource.is_dir():
            for child in sorted(resource.iterdir(), key=lambda item: item.name):
                yield from walk(child, relative / child.name)
        else:
            yield relative, resource

    yield from walk(root, Path())
