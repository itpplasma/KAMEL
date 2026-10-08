"""Prepare KIM inputs from a MARS-F profile quartet and an equilibrium table.

Install from the repository root with ``pip install -e ./KIM/python``. Edit
the paths and source declarations below, then run ``python prepare_marsf.py``.
Relative paths are resolved from the directory where you run the script.
The JSON request supplies your plasma, mode, field, and numerical settings.
Preparation writes a new directory; it does not run the scientific solver.
"""

from pathlib import Path

from kim import MarsFMetadata, SimulationConfig, prepare_marsf_case, read_marsf_profiles

# Change these four paths. The output directory must not already exist.
MARSF_DIRECTORY = Path("./marsf-input")
EQUILIBRIUM_FILE = Path("./equil_r_q_psi.dat")
REQUEST_FILE = Path("./request.json")
DESTINATION = Path("./prepared-case")

# Check these declarations against your source data; nothing is auto-detected.
# The source directory must contain PROFDEN.IN, PROFTE.IN, PROFTI.IN, PROFROT.IN.
# Each has one header line followed by two columns: radial coordinate and value.
METADATA = MarsFMetadata(
    source="my-marsf-case",
    coordinate="sqrt_psiN",
    coordinate_unit="1",
    density_unit="1/m^3",
    electron_temperature_unit="eV",
    ion_temperature_unit="eV",
    toroidal_velocity_unit="m/s",  # Linear toroidal velocity, not angular rotation.
    equilibrium_provenance="my-equilibrium-calculation",
)


def main() -> None:
    config = SimulationConfig.model_validate_json(REQUEST_FILE.read_text(encoding="utf-8"))
    source = read_marsf_profiles(MARSF_DIRECTORY, METADATA)
    case = prepare_marsf_case(
        source,
        config,
        DESTINATION,
        equilibrium_file=EQUILIBRIUM_FILE,
        q_operation="preserve",  # Keep the equilibrium's q sign.
    )
    print(f"Prepared case: {case.directory}")
    print(f"KIM request: {case.request}")
    print(f"Conversion and provenance report: {case.report}")


if __name__ == "__main__":
    main()
