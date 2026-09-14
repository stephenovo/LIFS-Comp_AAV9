"""NCBI SRA metadata retrieval for the Fit4Function BioProject."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from typing import Any

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _package_to_row(package: ET.Element) -> dict[str, Any]:
    experiment = package.find("EXPERIMENT")
    sample = package.find("SAMPLE")
    run = package.find("./RUN_SET/RUN")
    if experiment is None or sample is None or run is None:
        raise ValueError("Incomplete SRA experiment package")
    attributes = {
        node.findtext("TAG"): node.findtext("VALUE")
        for node in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
    }
    return {
        "experiment_accession": experiment.get("accession"),
        "experiment_alias": experiment.get("alias"),
        "experiment_title": experiment.findtext("TITLE"),
        "biosample_accession": sample.get("accession"),
        "run_accession": run.get("accession"),
        "bytes": int(run.get("size", 0)),
        "spots": int(run.get("total_spots", 0)),
        "library": attributes.get("library"),
        "assay_type": attributes.get("type"),
        "sample": attributes.get("sample"),
        "lab_host": attributes.get("lab_host"),
    }


def fetch_sra_manifest(bioproject: str = "PRJNA1131359") -> list[dict[str, Any]]:
    """Fetch all public SRA experiment metadata for a BioProject."""
    search = requests.get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "sra",
            "term": f"{bioproject}[bioproject]",
            "retmax": 10000,
            "retmode": "json",
        },
        timeout=60,
    )
    search.raise_for_status()
    identifiers = search.json()["esearchresult"]["idlist"]
    rows: list[dict[str, Any]] = []
    for start in range(0, len(identifiers), 100):
        response = requests.get(
            f"{EUTILS}/efetch.fcgi",
            params={
                "db": "sra",
                "id": ",".join(identifiers[start : start + 100]),
                "retmode": "xml",
            },
            timeout=120,
        )
        response.raise_for_status()
        root = ET.fromstring(response.content)
        rows.extend(_package_to_row(package) for package in root.findall("EXPERIMENT_PACKAGE"))
    return rows


def summarize_sra_manifest(rows: list[dict[str, Any]]) -> dict[str, object]:
    grouped: defaultdict[tuple[str | None, str | None], dict[str, int]] = defaultdict(
        lambda: {"runs": 0, "bytes": 0}
    )
    for row in rows:
        group = grouped[(row.get("library"), row.get("assay_type"))]
        group["runs"] += 1
        group["bytes"] += int(row.get("bytes", 0))
    return {
        "run_count": len(rows),
        "total_bytes": sum(int(row.get("bytes", 0)) for row in rows),
        "libraries": dict(Counter(str(row.get("library")) for row in rows)),
        "assay_types": dict(Counter(str(row.get("assay_type")) for row in rows)),
        "groups": {
            f"{library}__{assay_type}": values
            for (library, assay_type), values in sorted(
                grouped.items(), key=lambda item: (str(item[0][0]), str(item[0][1]))
            )
        },
    }
