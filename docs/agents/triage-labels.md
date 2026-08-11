# Triage Labels

The engineering skills use five canonical workflow roles. Their GitHub label mapping is:

| Canonical role | GitHub label | Meaning |
|---|---|---|
| `needs-triage` | `needs-triage` | Maintainer evaluation required |
| `needs-info` | `needs-info` | Waiting for information from the reporter |
| `ready-for-agent` | `ready-for-agent` | Fully specified and ready for an autonomous agent |
| `ready-for-human` | `ready-for-human` | Ready for human implementation or review |
| `wontfix` | `wontfix` | Will not be actioned |

When a skill marks an issue AFK-ready, apply `ready-for-agent`. Physics or architecture decisions that still need a person use `ready-for-human` until resolved.
