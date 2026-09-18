# Reproducibility and smoke test

The full Fit4Function reconstruction depends on external release data and SRA/ENA downloads. It is not a suitable first command for a new checkout.

## Five-minute local check

From a clean checkout:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
aav9-sma demo --output-dir artifacts/demo
pytest -q
```

The `demo` command generates deterministic synthetic data, runs the canonical audit, creates prediction columns, applies the packaging gate, and writes a ranked candidate table. It never downloads research data and its outputs must not be used for biological claims.

The same flow is available as:

```bash
bash scripts/reproduce_demo.sh
```

## Full-data run

The full run requires the external Fit4Function checkout, the pinned public source release, and SRA/ENA inputs described in [FIT4FUNCTION_RUNBOOK.md](FIT4FUNCTION_RUNBOOK.md). Before publishing a full result, record:

- this repository commit;
- the exact Fit4Function source commit;
- the input data manifest and SHA256 checksums;
- Python and dependency versions;
- command-line arguments and random seeds;
- model and output artifact checksums.

The Fit4Function checkout is pinned in [source_manifest.json](source_manifest.json)
and must be checked out at that commit before running the full-data commands.

The demo is an installation and interface check, not a substitute for the data gate, sequence-aware validation, or experimental validation.

## Checksum and full-data run

Full-data inputs are intentionally external and are not downloaded by the
repository. Copy [`data_manifest.example.json`](data_manifest.example.json),
replace each placeholder with a real SHA256 digest, and verify it before a run:

```bash
aav9-sma verify-manifest docs/data_manifest.json --root .
```

The non-destructive full-data wrapper checks the pinned Fit4Function commit,
verifies the manifest, runs the release audit and both Animal 4 development
benchmarks, then runs the virtual screen:

```bash
bash scripts/reproduce_full.sh \
  --fit4function-dir data/raw/fit4function_official \
  --manifest docs/data_manifest.json \
  --data-root . \
  --reconstructed data/processed/fit4function_multiorgan_reconstructed.csv.gz \
  --screen-csv data/raw/fit4function_official/data/fit4function_library_screens.csv
```

This wrapper never downloads raw data and never calls a result from Animal 4 a
final blind validation. The Animal 4 role policy is documented in
[`ANIMAL4_POLICY.md`](ANIMAL4_POLICY.md).
