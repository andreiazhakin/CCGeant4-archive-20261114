"""Transport fully stripped Cu, Xe and U through 100 micrometres of silicon."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from prepare_tables import EXAMPLE, unpack


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--geant4-bin", type=Path, help="DLL directory of the matching Windows installation")
    parser.add_argument("--events", type=int, default=8, help="Events per projectile")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--tables", type=Path, help="Use existing inputs in the example directory layout")
    parser.add_argument("--output", type=Path, default=EXAMPLE / "results")
    args = parser.parse_args()
    if args.events < 1 or args.threads < 1:
        parser.error("--events and --threads must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tables = args.tables.resolve() if args.tables else unpack(output / "tables")
    species = ((29, 63), (54, 124), (92, 238))
    required = [tables / "global.csv", tables / "atima_low.csv"]
    for z, _ in species:
        required.extend(tables / name / f"Z{z:03d}.csv" for name in ("delta", "stopping"))
    for path in required:
        if not path.is_file():
            parser.error(f"missing example input: {path}")
    env = os.environ.copy()
    if args.geant4_bin:
        env["PATH"] = str(args.geant4_bin.resolve()) + os.pathsep + env.get("PATH", "")
    results = []
    summary = output / "summary.json"
    summary.unlink(missing_ok=True)
    for z, a in species:
        csv_path = output / f"Z{z:03d}.csv"
        log_path = csv_path.with_suffix(".log")
        command = [
            str(args.executable.resolve()), "--delta", str(tables / "delta"),
            "--stopping", str(tables / "stopping"), "--global", str(tables / "global.csv"),
            "--atima-low", str(tables / "atima_low.csv"), "--Z", str(z), "--A", str(a),
            "--N", "0", "--energy", "400", "--thickness", "100", "--cut", "25",
            "--seed", "6110101", "--events", str(args.events), "--threads", str(args.threads),
            "--kernel", "0", "--out", str(csv_path),
        ]
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        passed = result.returncode == 0 and "DIRAC_TRANSPORT_PASS" in log_path.read_text(
            encoding="utf-8", errors="replace"
        )
        results.append({"Z": z, "A": a, "events": args.events, "passed": passed})
        summary.write_text(json.dumps({
            "passed": len(results) == len(species) and all(run["passed"] for run in results),
            "threads": args.threads, "seed": 6110101, "runs": results,
        }, indent=2) + "\n", encoding="utf-8")
        if not passed:
            raise SystemExit(f"Z={z} failed; see {log_path}")
        print(f"Z={z}: PASS")


if __name__ == "__main__":
    main()
