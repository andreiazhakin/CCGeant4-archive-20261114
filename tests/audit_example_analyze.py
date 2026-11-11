"""Check event-level energy accounting and summarize a completed transport example."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics


def summarize(directory):
    setup = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    species = [(29, 63), (54, 124), (92, 238)]
    if not setup["passed"] or [(r["Z"], r["A"]) for r in setup["runs"]] != species:
        raise ValueError("The transport run did not complete for all three ions")
    result = []
    for run in setup["runs"]:
        path = directory / f"Z{run['Z']:03d}.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
        events = run["events"]
        if events <= 0 or len(rows) != events or {r["event_id"] for r in rows} != set(range(events)):
            raise ValueError(f"{path.name}: unexpected event IDs or count")
        total_balance = []
        for row in rows:
            if not all(math.isfinite(v) for v in row.values()):
                raise ValueError(f"{path.name}: non-finite event value")
            if (row["Z"], row["A"], row["N_initial"], row["E_u_MeV"]) != (run["Z"], run["A"], 0, 400):
                raise ValueError(f"{path.name}: wrong projectile or initial state")
            if (row["K_final_MeV"] < 0 or row["edep_all_MeV"] < 0
                    or not 0 <= row["q_final"] <= run["Z"]
                    or not 0 <= row["max_delta_T_over_Tmax"] <= 1 + 1e-10):
                raise ValueError(f"{path.name}: invalid event {int(row['event_id'])}")
            total_balance.append(abs(row["K_initial_MeV"] + row["injected_kinetic_MeV"]
                                     - row["edep_all_MeV"] - row["escaped_kinetic_MeV"]))
        primary_balance = max(abs(r["primary_balance_MeV"]) for r in rows)
        if max(max(total_balance), primary_balance) >= 1e-6:
            raise ValueError(f"{path.name}: energy balance exceeds 1e-6 MeV")
        result.append({"Z": run["Z"], "A": run["A"], "events": events,
                       "mean_primary_loss_MeV": statistics.fmean(r["primary_loss_MeV"] for r in rows),
                       "mean_slab_deposit_MeV": statistics.fmean(r["edep_slab_MeV"] for r in rows),
                       "total_explicit_delta_electrons": sum(int(r["delta_count"]) for r in rows),
                       "max_event_balance_MeV": max(total_balance), "max_primary_balance_MeV": primary_balance})
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(args.results), indent=2))
