#!/usr/bin/env bash
set -euo pipefail

python -m pip install -e ".[dev]"
aav9-sma demo --output-dir artifacts/demo
aav9-sma audit-data artifacts/demo/demo_canonical.csv --output artifacts/demo/audit_recheck.json
pytest -q tests/test_demo.py
