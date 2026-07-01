#!/usr/bin/env python3
"""Up/down migration of golden-record input files across input-API drift.

Input decks drift relative to the binary they run against. Two directions:

  up    old historical deck  -> current binary  (rename/drop/add keys the new
        binary expects). Brings a stored deck to the repo's current schema.
  down  current deck         -> legacy binary   (strip/rename keys the old
        binary does not know). Used by the 3-version A/B/C legacy reruns.

Migrations live in <repo>/migrations/NNNN_slug.{yaml,yml,json}, ordered by the
numeric prefix. Each declares ops on Fortran-namelist keys. up applies the ops
in order; down applies them in reverse with inverted semantics, so one
definition serves both directions.

Ops (each names the namelist block it acts in; block match is case-insensitive):
  rename  {namelist, from, to}          up from->to    down to->from
  drop    {namelist, key, default?}     up remove key  down add key=default
  add     {namelist, key, value}        up add         down remove
  set     {namelist, key, old, new}     up value->new  down value->old

Namelist edits are line-based and formatting-preserving: only the touched key
line changes. Applying a direction is idempotent (a rename whose source key is
absent is a no-op), so reruns are safe.
"""
import argparse
import json
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

BLOCK_START = re.compile(r"^\s*&(\w+)", re.IGNORECASE)
BLOCK_END = re.compile(r"^\s*(/|&end)\s*$", re.IGNORECASE)


def _key_re(key):
    # match "  key = value   ! comment", capturing indent, key, sep, rest
    return re.compile(
        r"^(?P<indent>\s*)(?P<key>" + re.escape(key) + r")(?P<sep>\s*=\s*)(?P<rest>.*)$",
        re.IGNORECASE,
    )


def load_migrations(mig_dir):
    mig_dir = Path(mig_dir)
    out = []
    if not mig_dir.is_dir():
        return out
    for p in sorted(mig_dir.iterdir()):
        if p.suffix.lower() in (".yaml", ".yml"):
            if yaml is None:
                raise SystemExit("pyyaml not available; convert %s to .json" % p.name)
            data = yaml.safe_load(p.read_text())
        elif p.suffix.lower() == ".json":
            data = json.loads(p.read_text())
        else:
            continue
        if not data:
            continue
        data.setdefault("id", p.stem)
        out.append(data)
    return out


def invert_op(op):
    kind = op["op"]
    if kind == "rename":
        return {"op": "rename", "namelist": op["namelist"], "from": op["to"], "to": op["from"]}
    if kind == "drop":
        if "default" not in op:
            return None  # cannot re-add without a default; down is a no-op for this op
        return {"op": "add", "namelist": op["namelist"], "key": op["key"], "value": op["default"]}
    if kind == "add":
        return {"op": "drop", "namelist": op["namelist"], "key": op["key"], "default": op.get("value")}
    if kind == "set":
        return {"op": "set", "namelist": op["namelist"], "key": op["key"], "old": op["new"], "new": op["old"]}
    raise ValueError("unknown op: %r" % kind)


def _in_block(lines, namelist):
    """Yield (idx, is_end_idx) markers: returns list of (start, end) line ranges
    for the named block. Range is [start+1, end) i.e. body lines exclusive of
    the &name and the closing /."""
    ranges = []
    cur = None
    for i, ln in enumerate(lines):
        m = BLOCK_START.match(ln)
        if m:
            cur = (m.group(1).lower(), i)
            continue
        if cur and BLOCK_END.match(ln):
            if cur[0] == namelist.lower():
                ranges.append((cur[1], i))
            cur = None
    return ranges


def apply_op(lines, op, log):
    nl = op["namelist"]
    ranges = _in_block(lines, nl)
    if not ranges:
        log.append("  skip (%s): namelist &%s not found" % (op["op"], nl))
        return lines, False
    changed = False
    kind = op["op"]
    for start, end in ranges:
        if kind in ("rename", "drop", "set"):
            key = op.get("from") or op.get("key")
            kre = _key_re(key)
            for i in range(start + 1, end):
                m = kre.match(lines[i])
                if not m:
                    continue
                if kind == "rename":
                    lines[i] = m.group("indent") + op["to"] + m.group("sep") + m.group("rest")
                    log.append("  rename &%s %s -> %s" % (nl, key, op["to"]))
                elif kind == "drop":
                    log.append("  drop &%s %s" % (nl, key))
                    lines[i] = None
                elif kind == "set":
                    rest = m.group("rest")
                    rest2 = re.sub(re.escape(str(op["old"])), str(op["new"]), rest, count=1)
                    lines[i] = m.group("indent") + m.group("key") + m.group("sep") + rest2
                    log.append("  set &%s %s = %s" % (nl, key, op["new"]))
                changed = True
        elif kind == "add":
            key = op["key"]
            kre = _key_re(key)
            present = any(kre.match(lines[i]) for i in range(start + 1, end) if lines[i] is not None)
            if present:
                log.append("  skip add &%s %s (already present)" % (nl, key))
                continue
            indent = "    "
            for i in range(start + 1, end):
                if lines[i] is not None and lines[i].strip():
                    mi = re.match(r"^(\s*)", lines[i])
                    indent = mi.group(1) or indent
                    break
            lines.insert(end, "%s%s = %s" % (indent, key, op["value"]))
            log.append("  add &%s %s = %s" % (nl, key, op["value"]))
            changed = True
            # ranges computed before insert; recompute to stay correct
            ranges = _in_block(lines, nl)
    return [l for l in lines if l is not None], changed


def migrate_file(path, migrations, direction, log):
    path = Path(path)
    lines = path.read_text().splitlines()
    any_change = False
    chain = migrations if direction == "up" else list(reversed(migrations))
    for mig in chain:
        if mig.get("file") and Path(mig["file"]).name != path.name:
            continue
        ops = mig.get("ops", [])
        seq = ops if direction == "up" else [invert_op(o) for o in reversed(ops)]
        for op in seq:
            if op is None:
                continue
            lines, ch = apply_op(lines, op, log)
            any_change = any_change or ch
    if any_change:
        text = "\n".join(lines)
        if not text.endswith("\n"):
            text += "\n"
        path.write_text(text)
    return any_change


def cmd_apply(args):
    migrations = load_migrations(args.migrations)
    if args.id:
        migrations = [m for m in migrations if m["id"] in args.id]
    if not migrations:
        print("no migrations to apply", file=sys.stderr)
        return 0
    targets = []
    if args.file:
        targets = [Path(args.file)]
    else:
        d = Path(args.dir)
        names = {Path(m["file"]).name for m in migrations if m.get("file")}
        targets = [d / n for n in names if (d / n).is_file()]
    if not targets:
        print("no target input files found for the loaded migrations", file=sys.stderr)
        return 0
    rc = 0
    for t in targets:
        log = []
        try:
            ch = migrate_file(t, migrations, args.direction, log)
        except Exception as e:  # noqa: BLE001
            print("ERROR migrating %s: %s" % (t, e), file=sys.stderr)
            rc = 1
            continue
        tag = "changed" if ch else "no-op"
        print("[%s %s] %s" % (args.direction, tag, t))
        for line in log:
            print(line)
    return rc


def cmd_list(args):
    for m in load_migrations(args.migrations):
        print("%-28s %s" % (m["id"], m.get("description", "")))
        for op in m.get("ops", []):
            print("   - %s" % op)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("apply", help="apply migrations up or down")
    a.add_argument("--direction", choices=("up", "down"), required=True)
    a.add_argument("--migrations", default="migrations")
    g = a.add_mutually_exclusive_group(required=True)
    g.add_argument("--file", help="single input file to migrate")
    g.add_argument("--dir", help="case/run dir; migrate every file a migration names")
    a.add_argument("--id", nargs="*", help="restrict to these migration ids")
    a.set_defaults(func=cmd_apply)
    l = sub.add_parser("list", help="list migrations")
    l.add_argument("--migrations", default="migrations")
    l.set_defaults(func=cmd_list)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
