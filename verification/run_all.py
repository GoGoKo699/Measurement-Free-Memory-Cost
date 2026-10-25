#!/usr/bin/env python3
"""Run the scientific finite checks and compare committed reference outputs."""
from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CHECKS = (
    ("pauli_algebra", "verify_pauli_algebra.py"),
    ("archive", "verify_archive.py"),
    ("five_qubit", "verify_five_qubit.py"),
)
METADATA = {"python_version", "numpy_version"}


def compare(actual, expected, path="result"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"{path}: JSON keys differ")
        for key, value in expected.items():
            if key not in METADATA:
                compare(actual[key], value, f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AssertionError(f"{path}: list shape differs")
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(expected, float):
        if not isinstance(actual, (int, float)) or not math.isclose(
            actual, expected, rel_tol=1e-9, abs_tol=1e-9
        ):
            raise AssertionError(f"{path}: {actual!r} != {expected!r}")
    elif actual != expected:
        raise AssertionError(f"{path}: {actual!r} != {expected!r}")


def main():
    summary = []
    for name, filename in CHECKS:
        completed = subprocess.run(
            [sys.executable, str(HERE / filename)],
            check=True, capture_output=True, text=True
        )
        actual = json.loads(completed.stdout)
        expected_path = HERE / "expected" / f"{name}.json"
        expected = json.loads(expected_path.read_text(encoding="utf-8"))
        compare(actual, expected, name)
        summary.append({
            "check": name,
            "status": "PASS",
            "reference_output": str(expected_path.relative_to(HERE.parent))
        })
    print(json.dumps({
        "status": "PASS",
        "checks": summary,
        "scope": "Finite identities and scalar certificates; see the written proofs."
    }, indent=2))


if __name__ == "__main__":
    main()

