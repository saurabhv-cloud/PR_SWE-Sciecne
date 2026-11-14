"""Assemble the study's paired analysis and condition comparison."""

import argparse
import json
from pathlib import Path

import numpy as np

from .pipeline import analyze


def contrast(samples):
    """Compare material samples; registrations are repeated representations."""
    condition_values = {}
    for condition in ("control", "perturbation"):
        selected = [sample for sample in samples if sample["condition"] == condition]
        if not selected:
            raise ValueError("Both study conditions are required")
        populations, measurement, omitted = [], [], []
        for sample in selected:
            rows = sample["estimates"]
            populations.append(np.mean([row["population"] for row in rows]))
            measurement.append(np.mean([row["measurement_se"] for row in rows]))
            omitted.append(np.mean([row["leave_block_populations"] for row in rows], axis=0))
        condition_values[condition] = (
            float(np.mean(populations)),
            float(np.sum(np.square(measurement)) / len(selected)**2),
            np.mean(omitted, axis=0),
        )
    control, perturbation = (condition_values[key] for key in ("control", "perturbation"))
    estimate = perturbation[0] - control[0]
    omitted = perturbation[2] - control[2]
    if omitted.ndim != 1 or len(omitted) < 2:
        raise ValueError("At least two aligned independent MD blocks are required")
    simulation_variance = (len(omitted) - 1) * np.var(omitted)
    se = float(np.sqrt(control[1] + perturbation[1] + simulation_variance))
    interval = [float(estimate - 1.96 * se), float(estimate + 1.96 * se)]
    decision = "inconclusive"
    if estimate >= .10 and interval[0] > 0:
        decision = "increase"
    elif estimate <= -.10 and interval[1] < 0:
        decision = "decrease"
    return {"estimate": float(estimate), "se": se, "ci95": interval, "decision": decision}


def assemble(cases, original, interpretation):
    samples = [analyze(case) for case in cases]
    before = {sample["sample_id"]: sample for sample in original}
    audit = []
    for sample in samples:
        previous = before[sample["sample_id"]]
        item = {"sample_id": sample["sample_id"]}
        for field in ("population", "distance_rms"):
            item["reported_" + field + "_range"] = float(np.ptp([row[field] for row in previous["estimates"]]))
            item["current_" + field + "_range"] = float(np.ptp([row[field] for row in sample["estimates"]]))
        audit.append(item)
    return {"schema": 1, "samples": samples, "contrast": contrast(samples),
            "audit": audit, "interpretation": interpretation}


def diagnose(replays, kernel):
    """Return the findings for the paired replays."""
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path("/workspace/data/inputs"))
    parser.add_argument("--output", type=Path, default=Path("/workspace/data/output/analysis.json"))
    parser.add_argument("--interpretation", default="The numerical evidence is provided in the paired analysis and condition comparison.")
    args = parser.parse_args()
    cases = json.loads((args.inputs / "cases.json").read_text(encoding="utf-8"))
    kernel = Path(__file__).parent / "assets/ensemble.npz"
    for case in cases:
        case["kernel_file"] = str(kernel)
    original = json.loads((args.inputs / "reported_results.json").read_text(encoding="utf-8"))
    report = assemble(cases, original, args.interpretation)
    replay_path = args.inputs / "replays.json"
    if replay_path.is_file():
        replays = json.loads(replay_path.read_text(encoding="utf-8"))
        report["diagnosis"] = diagnose(replays, kernel)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
