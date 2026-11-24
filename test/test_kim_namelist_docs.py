"""Keep the KIM namelist reference and example namelists in sync with the source.

KIM reads its namelist groups without error recovery, so an example that sets an
option the code no longer declares aborts at startup. These checks compare the
`namelist /GROUP/ ...` declarations in KIM/src with KIM/nmls/README.md and with
every KIM/nmls/*.nml file.
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
KIM_SOURCE_DIR = REPO_ROOT / "KIM" / "src"
NAMELIST_DIR = REPO_ROOT / "KIM" / "nmls"
REFERENCE = NAMELIST_DIR / "README.md"

# Groups that kim_read_config reads sequentially, without rewinding, in this order.
REQUIRED_GROUPS_IN_READ_ORDER = (
    "kim_config",
    "wkb_dispersion",
    "kim_io",
    "kim_setup",
    "kim_grid",
    "kim_profiles",
)

NAMELIST_DECLARATION = re.compile(
    r"^\s*namelist\s*/\s*(\w+)\s*/(.*)$", re.IGNORECASE | re.MULTILINE
)
DOCUMENTED_OPTION = re.compile(r"^\s*-\s*(\w+)\s*\.\.\.", re.MULTILINE)
GROUP_BLOCK = re.compile(r"^\s*[&$](\w+)(.*?)(?:^\s*/|\s/[ \t]*$)", re.MULTILINE | re.DOTALL)
ASSIGNED_KEY = re.compile(r"(?:^|[\s,])([A-Za-z_]\w*)\s*(?:\([^)]*\))?\s*=")


def strip_fortran_comments(text):
    """Remove `!` comments outside character literals."""
    stripped_lines = []
    for line in text.splitlines():
        quote = None
        for index, char in enumerate(line):
            if quote:
                if char == quote:
                    quote = None
            elif char in "'\"":
                quote = char
            elif char == "!":
                line = line[:index]
                break
        stripped_lines.append(line)
    return "\n".join(stripped_lines)


def join_continuation_lines(text):
    return re.sub(r"&[ \t]*\n[ \t]*&?", " ", text)


def parse_namelist_declarations(source_files):
    groups = {}
    for path in source_files:
        text = join_continuation_lines(strip_fortran_comments(path.read_text(errors="replace")))
        for match in NAMELIST_DECLARATION.finditer(text):
            variables = {name.strip().lower() for name in match.group(2).split(",") if name.strip()}
            groups.setdefault(match.group(1).lower(), set()).update(variables)
    return groups


def parse_reference_sections(markdown):
    sections = {}
    for chunk in re.split(r"^## ", markdown, flags=re.MULTILINE)[1:]:
        heading, _, body = chunk.partition("\n")
        options = {name.lower() for name in DOCUMENTED_OPTION.findall(body)}
        sections[heading.strip().lower()] = options
    return sections


def parse_namelist_file(path):
    """Return the groups of a namelist file in file order as (name, keys) pairs."""
    text = strip_fortran_comments(path.read_text(errors="replace"))
    return [
        (match.group(1).lower(), {key.lower() for key in ASSIGNED_KEY.findall(match.group(2))})
        for match in GROUP_BLOCK.finditer(text)
    ]


@pytest.fixture(scope="module")
def declared_groups():
    groups = parse_namelist_declarations(sorted(KIM_SOURCE_DIR.rglob("*.f90")))
    assert set(REQUIRED_GROUPS_IN_READ_ORDER) <= set(groups), "KIM namelist declarations not found"
    return groups


@pytest.fixture(scope="module")
def reference_sections():
    return parse_reference_sections(REFERENCE.read_text())


def test_every_declared_group_has_a_reference_section(declared_groups, reference_sections):
    missing = sorted(set(declared_groups) - set(reference_sections))
    assert not missing, f"KIM/nmls/README.md has no section for namelist groups {missing}"


def test_reference_documents_exactly_the_declared_options(declared_groups, reference_sections):
    problems = []
    for group, declared in sorted(declared_groups.items()):
        documented = reference_sections.get(group, set())
        undocumented = sorted(declared - documented)
        obsolete = sorted(documented - declared)
        if undocumented:
            problems.append(f"{group}: undocumented {undocumented}")
        if obsolete:
            problems.append(f"{group}: documented but not read by KIM {obsolete}")
    assert not problems, "KIM/nmls/README.md is out of sync:\n" + "\n".join(problems)


@pytest.mark.parametrize(
    "namelist_path", sorted(NAMELIST_DIR.glob("*.nml")), ids=lambda path: path.name
)
def test_example_namelist_matches_declarations(namelist_path, declared_groups):
    groups = parse_namelist_file(namelist_path)
    problems = []
    for group, keys in groups:
        if group not in declared_groups:
            problems.append(f"unknown group &{group}")
            continue
        unknown = sorted(keys - declared_groups[group])
        if unknown:
            problems.append(f"&{group}: options not read by KIM {unknown}")

    required_positions = [name for name, _ in groups if name in REQUIRED_GROUPS_IN_READ_ORDER]
    if required_positions != list(REQUIRED_GROUPS_IN_READ_ORDER):
        problems.append(
            f"required groups must appear once in the order {list(REQUIRED_GROUPS_IN_READ_ORDER)}; "
            f"found {required_positions}"
        )
    assert not problems, f"{namelist_path.name}:\n" + "\n".join(problems)
