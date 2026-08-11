# Issue tracker: GitHub

Issues and PRDs for this repository live in GitHub Issues under `itpplasma/KAMEL`. Use the `gh` CLI for all operations.

## Conventions

- Create: `gh issue create --repo itpplasma/KAMEL --title "..." --body-file <file>`.
- Read: `gh issue view <number> --repo itpplasma/KAMEL --comments`.
- List: `gh issue list --repo itpplasma/KAMEL --state open --json number,title,body,labels,comments`.
- Comment: `gh issue comment <number> --repo itpplasma/KAMEL --body "..."`.
- Label: `gh issue edit <number> --repo itpplasma/KAMEL --add-label "..."`.
- Close: `gh issue close <number> --repo itpplasma/KAMEL --comment "..."`.

Use `--body-file` or a quoted heredoc for multiline issue bodies. Never interpolate untrusted text through the shell.

## Skill mappings

- When a skill says “publish to the issue tracker,” create a GitHub issue in `itpplasma/KAMEL`.
- When a skill says “fetch the relevant ticket,” read its complete body, labels, and comments.
- Publish dependency blockers first so later issues can reference actual issue numbers.
- Do not close or rewrite a parent/tracking issue while creating implementation slices.
