# KAMEL Public Release Preparation

**Date:** 2026-03-04
**Goal:** Prepare the KAMEL repository for public visibility on GitHub.

## Scope

| Item | Action |
|------|--------|
| `CONTRIBUTING.md` | Create new (~50 lines) |
| `CODE_OF_CONDUCT.md` | Create new (Contributor Covenant v2.1, plasma.itp@tugraz.at) |
| `README.md` | Add Contributing + License sections |
| `LICENSE` | Update copyright year to 2026 |
| Tag `v1.0.0` | Annotated tag on main after commit |

## Decisions

- **Primary goal:** Academic visibility and citation, with contribution-ready structure.
- **Institutional approval:** Obtained from TU Graz ITP and all contributors.
- **AGENTS.md:** Keep in public repo as-is.
- **No secrets or sensitive data found** in the repository.
- **Third-party dependencies** all fetched at build time via CMake FetchContent, compatible with MIT.
