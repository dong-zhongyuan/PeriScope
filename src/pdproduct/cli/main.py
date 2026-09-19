"""pdproduct CLI — thin entrypoint; all logic lives in library modules."""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys

from pdproduct.contracts.loader import SCHEMA_FILES, validate_file
from pdproduct.provenance import new_run_manifest
from pdproduct.registry import methods as reg


def _json_default(value: object) -> str:
    """Serialize YAML date values without hiding unexpected objects."""
    if isinstance(value, (_dt.datetime, _dt.date)):
        return value.isoformat()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _json_dump(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, indent=2, default=_json_default)


def _cmd_validate_manifest(args: argparse.Namespace) -> int:
    ok = validate_file(args.path, args.schema)
    print(f"SCHEMA_OK {args.path} against {args.schema}" if ok else "FAILED")
    return 0


def _cmd_registry(args: argparse.Namespace) -> int:
    if args.reg_cmd == "summary":
        print(_json_dump({"readiness": reg.readiness_summary(), "meta": reg.meta()}))
        return 0
    if args.reg_cmd == "query":
        result = reg.query(
            category=args.category,
            mvp_role_contains=args.mvp_role,
            evidence=args.evidence,
            readiness=args.readiness,
        )
        for mid in sorted(result):
            m = result[mid]
            print(f"{mid}	{m.get('category')}	{m.get('readiness')}	{m.get('mvp_role')}")
        return 0
    if args.reg_cmd == "get":
        print(_json_dump(reg.get(args.method_id)))
        return 0
    return 2


def _cmd_run(args: argparse.Namespace) -> int:
    if args.run_cmd == "init":
        manifest = new_run_manifest(
            args.run_id,
            approved_scope=args.scope,
            command=args.command,
            operator=args.operator,
            seed=args.seed,
        )
        text = _json_dump(manifest)
        if args.out:
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text + "\n")
            print(f"RUN_MANIFEST_WRITTEN {args.out}")
        else:
            print(text)
        return 0
    return 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pdproduct",
        description="NeuroBridge-VT: evidence-aware virtual tissue inference and assay design",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("validate-manifest", help="validate a JSON file against a contract schema")
    p.add_argument("path")
    p.add_argument("--schema", default="dataset_manifest", choices=sorted(SCHEMA_FILES))
    p.set_defaults(func=_cmd_validate_manifest)

    p = sub.add_parser("registry", help="query the method capability registry")
    reg_sub = p.add_subparsers(dest="reg_cmd", required=True)
    q = reg_sub.add_parser("query")
    q.add_argument("--category")
    q.add_argument("--mvp-role")
    q.add_argument("--evidence")
    q.add_argument("--readiness")
    q = reg_sub.add_parser("get")
    q.add_argument("method_id")
    reg_sub.add_parser("summary")
    p.set_defaults(func=_cmd_registry)

    p = sub.add_parser("run", help="run lifecycle helpers")
    run_sub = p.add_subparsers(dest="run_cmd", required=True)
    r = run_sub.add_parser("init", help="create a schema-valid run manifest draft")
    r.add_argument("run_id")
    r.add_argument("--scope", required=True, help="approved scope, e.g. phaseB-W1-smoke (ADR-0001)")
    r.add_argument("--command", required=True)
    r.add_argument("--operator")
    r.add_argument("--seed", type=int)
    r.add_argument("--out", help="write manifest JSON to this path (run dir must exist)")
    p.set_defaults(func=_cmd_run)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
