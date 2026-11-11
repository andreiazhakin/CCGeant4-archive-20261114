"""Compare the Xe C++ table probe with the independent binary-table reference."""
import argparse
import csv
import math
from pathlib import Path

FIELDS = ("M0", "M1", "M2", "full", "eq", "B")
KEY = ("Z", "Zt", "N", "E", "cut")
REFERENCE = Path(__file__).resolve().parent / "reference" / "xe_si_kernel.csv"


def load(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = {}
        for row in csv.DictReader(stream):
            key = tuple(format(float(row[name]), ".12g") for name in KEY)
            if key in rows:
                raise ValueError(f"Duplicate probe configuration {key}")
            rows[key] = row
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe_csv", type=Path)
    args = parser.parse_args()
    expected, actual = load(REFERENCE), load(args.probe_csv)
    if not expected or expected.keys() != actual.keys():
        raise ValueError(f"Inconsistent probe grid ({len(expected)} reference, {len(actual)} actual)")
    max_relative = 0.0
    for key, ref in expected.items():
        for field in FIELDS:
            reference, value = float(ref[field]), float(actual[key][field])
            if not math.isclose(reference, value, rel_tol=3e-10, abs_tol=0):
                raise ValueError(f"{key} {field}: expected {reference}, found {value}")
            if reference != 0:
                max_relative = max(max_relative, abs(value / reference - 1))
    print(f"PASS: {len(expected)} Xe/Si configurations; maximum relative error {max_relative:.3g}")


if __name__ == "__main__":
    main()
