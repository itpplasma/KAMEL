"""Validation and immutable staging for KIM radial profile files."""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from kim.config import ProfileConfig, SimulationConfig
from kim.errors import ProfileError
from numpy.typing import NDArray


@dataclass(frozen=True)
class ProfileData:
    """Validated contents of one two-column radial profile."""

    name: str
    units: str
    path: Path
    radius: NDArray[np.float64]
    values: NDArray[np.float64]


@dataclass(frozen=True)
class ProfileValidation:
    """Summary of a successfully validated profile set."""

    point_count: int
    radial_minimum: float
    radial_maximum: float
    resonance_radius: float


@dataclass(frozen=True)
class ProfileCopy:
    """Provenance for one copied profile file."""

    name: str
    source: Path
    destination: Path
    sha256: str


@dataclass(frozen=True)
class ProfileSet:
    """Named source files forming one KIM `r_eff` profile set."""

    density: Path
    electron_temperature: Path
    ion_temperature: Path
    safety_factor: Path
    radial_electric_field: Path
    toroidal_velocity: Path

    @classmethod
    def from_config(cls, config: ProfileConfig) -> ProfileSet:
        """Resolve profile filenames relative to a validated profile directory."""

        directory = config.directory
        return cls(
            density=directory / config.density_file,
            electron_temperature=directory / config.electron_temperature_file,
            ion_temperature=directory / config.ion_temperature_file,
            safety_factor=directory / config.safety_factor_file,
            radial_electric_field=directory / config.radial_electric_field_file,
            toroidal_velocity=directory / config.toroidal_velocity_file,
        )

    @classmethod
    def from_simulation(cls, config: SimulationConfig) -> ProfileSet:
        """Resolve the profile files referenced by a complete simulation request."""

        return cls.from_config(config.profiles)

    def validate_for(self, config: SimulationConfig) -> ProfileValidation:
        """Validate the profile set against a simulation's domain and mode numbers."""

        return self.validate(
            radial_minimum=config.grid.radial_minimum,
            plasma_radius=config.grid.plasma_radius,
            m_mode=config.setup.m_mode,
            n_mode=config.setup.n_mode,
        )

    def validate(
        self,
        *,
        radial_minimum: float,
        plasma_radius: float,
        m_mode: int,
        n_mode: int,
    ) -> ProfileValidation:
        """Validate syntax, units, grids, domain coverage, and q resonance."""

        if radial_minimum >= plasma_radius:
            raise ProfileError("radial_minimum must be smaller than plasma_radius")
        if m_mode == 0 or n_mode == 0:
            raise ProfileError("m_mode and n_mode must be nonzero for q-resonance validation")
        profiles = self._read_existing()
        density = profiles["n"]
        for name in ("Te", "Ti", "q"):
            self._require_matching_grid(profiles[name], density)
        if "Vz" in profiles:
            self._require_matching_grid(profiles["Vz"], density)

        for profile in profiles.values():
            self._require_coverage(profile, radial_minimum, plasma_radius)
        self._validate_positive(profiles["n"])
        self._validate_positive(profiles["Te"])
        self._validate_positive(profiles["Ti"])
        self._validate_density_units(density)
        resonance_radius = self._find_resonance(
            profiles["q"], radial_minimum, plasma_radius, m_mode, n_mode
        )
        return ProfileValidation(
            point_count=density.radius.size,
            radial_minimum=float(density.radius[0]),
            radial_maximum=float(density.radius[-1]),
            resonance_radius=resonance_radius,
        )

    def copy_to(self, destination: Path | str) -> tuple[ProfileCopy, ...]:
        """Copy present inputs into a new directory and record content digests."""

        target = Path(destination)
        self._require_distinct_paths()
        sources = []
        for specification in self._specifications():
            name, source, _units, required = specification
            if not source.is_file():
                if required:
                    raise ProfileError(f"{source}: required profile file is missing")
                continue
            sources.append((name, source))
        if target.exists():
            if not target.is_dir() or any(target.iterdir()):
                raise ProfileError(f"profile staging directory is not empty: {target}")
        else:
            target.mkdir(parents=True)

        records = []
        for name, source in sources:
            staged = target / source.name
            shutil.copyfile(source, staged, follow_symlinks=True)
            records.append(
                ProfileCopy(
                    name=name,
                    source=source.absolute(),
                    destination=staged.absolute(),
                    sha256=_sha256(staged),
                )
            )
        return tuple(records)

    def _read_existing(self) -> dict[str, ProfileData]:
        profiles = {}
        self._require_distinct_paths()
        for name, path, units, required in self._specifications():
            if not path.is_file():
                if required:
                    raise ProfileError(f"{path}: required profile file is missing")
                continue
            profiles[name] = _read_profile(name, units, path)
        return profiles

    def _require_distinct_paths(self) -> None:
        paths = [path for _name, path, _units, _required in self._specifications()]
        duplicates = {path for path in paths if paths.count(path) > 1}
        if duplicates:
            names = ", ".join(str(path) for path in sorted(duplicates))
            raise ProfileError(f"{names}: one file cannot serve multiple profile roles")

    def _specifications(self) -> tuple[tuple[str, Path, str, bool], ...]:
        return (
            ("n", self.density, "1/cm^3", True),
            ("Te", self.electron_temperature, "eV", True),
            ("Ti", self.ion_temperature, "eV", True),
            ("q", self.safety_factor, "1", True),
            ("Er", self.radial_electric_field, "statV/cm", False),
            ("Vz", self.toroidal_velocity, "cm/s", False),
        )

    @staticmethod
    def _require_matching_grid(profile: ProfileData, reference: ProfileData) -> None:
        if not np.array_equal(profile.radius, reference.radius):
            raise ProfileError(
                f"{profile.path}: radial grid must exactly match {reference.path.name}; "
                "the current Fortran reader associates values by row"
            )

    @staticmethod
    def _require_coverage(
        profile: ProfileData, radial_minimum: float, plasma_radius: float
    ) -> None:
        if profile.radius[0] > radial_minimum or profile.radius[-1] < plasma_radius:
            raise ProfileError(
                f"{profile.path}: radius [cm] range "
                f"[{profile.radius[0]:g}, {profile.radius[-1]:g}] does not cover requested "
                f"domain [{radial_minimum:g}, {plasma_radius:g}] cm"
            )

    @staticmethod
    def _validate_positive(profile: ProfileData) -> None:
        invalid = np.flatnonzero(profile.values <= 0.0)
        if invalid.size:
            row = int(invalid[0]) + 1
            raise ProfileError(
                f"{profile.path}: row {row} {profile.name} [{profile.units}] must be positive"
            )

    @staticmethod
    def _validate_density_units(profile: ProfileData) -> None:
        maximum = float(np.max(profile.values))
        if maximum > 1.0e17:
            raise ProfileError(
                f"{profile.path}: density [1/cm^3] maximum {maximum:.6g} exceeds 1e17; "
                "values appear to use SI 1/m^3 units"
            )

    @staticmethod
    def _find_resonance(
        profile: ProfileData,
        radial_minimum: float,
        plasma_radius: float,
        m_mode: int,
        n_mode: int,
    ) -> float:
        target = abs(m_mode / n_mode)
        inside = (profile.radius > radial_minimum) & (profile.radius < plasma_radius)
        radius = np.concatenate(([radial_minimum], profile.radius[inside], [plasma_radius]))
        values = np.concatenate(
            (
                [np.interp(radial_minimum, profile.radius, profile.values)],
                profile.values[inside],
                [np.interp(plasma_radius, profile.radius, profile.values)],
            )
        )
        difference = np.abs(values) - target
        exact = np.flatnonzero(difference == 0.0)
        if exact.size:
            return float(radius[exact[0]])
        crossings = np.flatnonzero(difference[:-1] * difference[1:] < 0.0)
        if not crossings.size:
            raise ProfileError(
                f"{profile.path}: |q| has no |m/n| = {target:g} crossing inside "
                f"[{radial_minimum:g}, {plasma_radius:g}] cm"
            )
        index = int(crossings[0])
        fraction = -difference[index] / (difference[index + 1] - difference[index])
        return float(radius[index] + fraction * (radius[index + 1] - radius[index]))


def _read_profile(name: str, units: str, path: Path) -> ProfileData:
    rows: list[tuple[float, float]] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "!")):
            continue
        columns = stripped.split()
        if len(columns) != 2:
            raise ProfileError(
                f"{path}: row {line_number} must contain exactly two numeric columns "
                f"(r_eff [cm], {name} [{units}])"
            )
        try:
            radius, value = (
                float(column.replace("D", "E").replace("d", "e")) for column in columns
            )
        except ValueError as error:
            raise ProfileError(
                f"{path}: row {line_number} must contain exactly two numeric columns "
                f"(r_eff [cm], {name} [{units}])"
            ) from error
        if not np.isfinite(radius) or not np.isfinite(value):
            raise ProfileError(
                f"{path}: row {line_number} r_eff [cm] and {name} [{units}] must be finite"
            )
        if rows and radius <= rows[-1][0]:
            raise ProfileError(f"{path}: row {line_number} r_eff [cm] must be strictly increasing")
        rows.append((radius, value))

    if len(rows) < 2:
        raise ProfileError(f"{path}: profile must contain at least two rows")
    data = np.asarray(rows, dtype=np.float64)
    return ProfileData(name=name, units=units, path=path, radius=data[:, 0], values=data[:, 1])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
