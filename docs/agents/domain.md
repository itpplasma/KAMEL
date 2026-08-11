# Domain Documentation

KAMEL uses a multi-context domain-document layout.

## Before exploring

1. Read `CONTEXT-MAP.md` at the repository root.
2. Read every component `CONTEXT.md` selected by the map for the task.
3. Read system-wide ADRs under `docs/adr/` and relevant component ADRs under `<component>/docs/adr/` when those directories exist.
4. For cross-code work, read every participating context rather than choosing only the initiating component.

Missing ADR directories are not errors. They are created when a durable architecture decision is first recorded.

## Use canonical vocabulary

Use the terms defined by the relevant `CONTEXT.md` in issue titles, acceptance criteria, tests, code, and documentation. Do not replace a defined term with an approximate synonym. If contexts use apparently conflicting terms, surface the conflict before changing an interface.

## Respect ownership boundaries

Each context identifies which component owns a calculation or policy. A cross-code feature should keep that responsibility in its owning component and pass a typed result across the boundary. If an implementation would move ownership, record and review an ADR first.

## Flag decision conflicts

If proposed work conflicts with an existing ADR or a settled statement in the active roadmap, state the conflict explicitly. Do not silently supersede the documented decision.
