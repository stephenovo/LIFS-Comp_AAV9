"""Command-line entry points for data auditing and candidate ranking."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from aav9_sma.data.audit import audit_csv
from aav9_sma.screening.score import rank_candidates


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aav9-sma")
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit_parser = subparsers.add_parser("audit-data", help="Audit a canonical CSV file")
    audit_parser.add_argument("input", type=Path)
    audit_parser.add_argument("--output", type=Path)

    rank_parser = subparsers.add_parser("rank-candidates", help="Rank model predictions")
    rank_parser.add_argument("input", type=Path)
    rank_parser.add_argument("--packaging-threshold", type=float, required=True)
    rank_parser.add_argument("--output", type=Path, required=True)
    return parser


def _write_json(payload: dict[str, object], output: Path | None) -> None:
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if output is None:
        print(rendered)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered + "\n", encoding="utf-8")


def main() -> None:
    args = _build_parser().parse_args()
    if args.command == "audit-data":
        _write_json(audit_csv(args.input).to_dict(), args.output)
        return
    if args.command == "rank-candidates":
        predictions = pd.read_csv(args.input)
        ranked = rank_candidates(predictions, packaging_threshold=args.packaging_threshold)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        ranked.to_csv(args.output, index=False)
        return
    raise RuntimeError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()

