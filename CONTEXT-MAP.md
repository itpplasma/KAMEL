# KAMEL Context Map

KAMEL is a multi-context scientific-computing repository. Read the context documents for every component touched by a task.

| Area | Context document | Read when working on |
|---|---|---|
| KIM | `KIM/CONTEXT.md` | Integral plasma response, periodic/global solvers, susceptibilities, fields, currents, and KIM transport coefficients |
| QL-Balance | `QL-Balance/CONTEXT.md` | Global radial transport, wave-code adapters, profile evolution, shielding normalization, and diffusion-coefficient assembly |
| KiLCA | `KiLCA/CONTEXT.md` | Cylindrical electromagnetic response, vacuum fields, and magnetic-drive provenance |
| Shared/common | Read every affected component context | Code under `common/` or cross-code data types and conventions |
| Python/KAMEL orchestration | Read every solver context used by the workflow | Configuration, preprocessing, execution, and HDF5 post-processing under `python/` |
| Cross-code features | Read all participating contexts | Interfaces between KIM, KiLCA, QL-Balance, shared code, or Python |

System-wide decisions belong under `docs/adr/`. Component-specific decisions belong under `<component>/docs/adr/`. These directories may be created when an architectural decision is first recorded.

The implementation plan for periodic KIM integral transport and QL-Balance coupling is `roadmap.md` on its feature branch. Treat it as the current design source until its decisions are replaced by an ADR or merged documentation.
