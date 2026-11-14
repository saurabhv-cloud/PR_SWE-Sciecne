"""Run the study analysis from raw spectra and registered crystal refinements."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import least_squares

from .model import CONTACT_CUTOFFS, GRID, PAIR_INDEX, ROTATIONS, oxygen_centers
from .structure import displacement_targets


def optical_targets(case):
    means, errors, reading_features = [], [], {}
    for channel in ("raman", "fluorescence"):
        calibration = case["calibration"][channel]
        readings, variances = [], []
        for record in case["spectra"]:
            if record["channel"] != channel:
                continue
            counts = np.asarray(record["counts"], dtype=float)
            blank = np.asarray(record["blank"], dtype=float)
            shape = np.asarray(calibration["template"]) * record["dose"] * calibration["gain"]
            weights = 1 / np.maximum(counts, 1.)
            precision = np.sum(weights * shape**2)
            value = np.sum(weights * shape * (counts - blank)) / precision
            reading_features[record["reading"]] = float(value)
            readings.append(value)
            variances.append(1 / precision)
        precision = 1 / np.asarray(variances)
        means.append(np.average(readings, weights=precision))
        errors.append(np.hypot(1 / np.sqrt(precision.sum()), calibration["feature_sd"]))
    return np.array(means), np.array(errors), reading_features


def parent_components(moment_matrices):
    """Pool symmetry replicas in the parent frame as one refinement observation."""
    transformed = ROTATIONS.transpose(0, 2, 1) @ moment_matrices @ ROTATIONS
    average = transformed.mean(axis=-3)
    return average[..., PAIR_INDEX[0], PAIR_INDEX[1]]


def load_kernel(path):
    with np.load(path, allow_pickle=False) as data:
        features = data["features"]
        contacts = data["contacts"]
    return features, contacts


def fit(features, contacts, observations, errors, mode):
    selected = {"joint": slice(None), "optics": slice(0, 2), "structure": slice(2, None)}[mode]
    bounds = np.array([(axis[0], axis[-1]) for axis in GRID]).T
    centers = bounds.mean(axis=0)
    scales = bounds[1] - bounds[0]
    observed = observations[selected]
    uncertainty = errors[selected]
    feature_model = RegularGridInterpolator(GRID, features[..., selected])
    contact_model = RegularGridInterpolator(GRID, contacts)

    def residuals(parameters):
        difference = (feature_model(parameters)[0] - observed) / uncertainty
        return np.concatenate([difference, .001 * (parameters - centers) / scales])

    result = least_squares(residuals, centers, bounds=bounds, x_scale=scales)
    if not result.success:
        raise RuntimeError(result.message)
    summary = contact_model(result.x)[0]
    summary[4] = np.sqrt(max(0., summary[4] - summary[3] ** 2))
    discrepancy = np.square(result.fun[:len(observed)]).sum()
    return result.x, summary, float(discrepancy)


def analyze(case, kernels=None, structure_reader=displacement_targets, mode="joint", uncertainty=True):
    """Return JSON-serializable populations and diagnostics for every representation."""
    if case["schema"] != 1:
        raise ValueError("Unsupported study case schema")
    if kernels is None:
        kernels = load_kernel(case["kernel_file"])
    feature_blocks, contact_blocks = kernels
    features = feature_blocks.mean(axis=-2)
    contacts = contact_blocks.mean(axis=-2)
    blocks = feature_blocks.shape[-2]
    geometry = np.asarray(case["geometry"])
    centers = oxygen_centers(geometry)
    optical, optical_error, readings = optical_targets(case)
    errors = np.r_[optical_error, np.full(6, case["structure_sd"])]
    estimates = []
    for representation in case["representations"]:
        moments = structure_reader(representation, centers)
        targets = np.r_[optical, parent_components(moments)]
        parameters, population, discrepancy = fit(features, contacts, targets, errors, mode)
        measurement_se = simulation_se = se = interval = [None, None]
        jackknife = []
        if uncertainty:
            perturbations = []
            for index, error in enumerate(errors):
                offset = np.zeros_like(errors)
                offset[index] = error
                _, lower, _ = fit(features, contacts, targets - offset, errors, mode)
                _, upper, _ = fit(features, contacts, targets + offset, errors, mode)
                perturbations.append((upper[[1, 4]] - lower[[1, 4]]) / 2)
            measurement_se = np.linalg.norm(perturbations, axis=0).tolist()
            for block in range(blocks):
                features_without = (blocks * features - feature_blocks[..., block, :]) / (blocks - 1)
                contacts_without = (blocks * contacts - contact_blocks[..., block, :]) / (blocks - 1)
                _, estimate, _ = fit(features_without, contacts_without, targets, errors, mode)
                jackknife.append(estimate[[1, 4]])
            simulation_se = np.sqrt((blocks - 1) * np.var(jackknife, axis=0)).tolist()
            se = np.hypot(measurement_se, simulation_se)
            ends = population[[1, 4], None] + 1.96 * se[:, None] * np.array([-1, 1])
            interval = [np.clip(ends[0], 0, 1).tolist(), np.maximum(ends[1], 0).tolist()]
        estimates.append({
            "representation": representation["id"],
            "population": float(population[1]),
            "se": se[0],
            "ci95": interval[0],
            "measurement_se": measurement_se[0],
            "simulation_se": simulation_se[0],
            "leave_block_populations": [float(value[0]) for value in jackknife],
            "distance_rms": float(population[4]),
            "distance_rms_se": se[1],
            "distance_rms_ci95": interval[1],
            "distance_rms_measurement_se": measurement_se[1],
            "distance_rms_simulation_se": simulation_se[1],
            "leave_block_distance_rms": [float(value[1]) for value in jackknife],
            "mean_distance": float(population[3]),
            "contact_sensitivity": dict(zip((str(value) for value in CONTACT_CUTOFFS), population[:3].tolist())),
            "parameters": parameters.tolist(),
            "standardized_residual_ss": discrepancy,
            "oxygen_displacement_trace": np.trace(moments, axis1=1, axis2=2).tolist(),
        })
    return {
        "sample_id": case["sample_id"],
        "condition": case["condition"],
        "mode": mode,
        "optical_features": dict(zip(("raman", "fluorescence"), optical.tolist())),
        "optical_se": dict(zip(("raman", "fluorescence"), optical_error.tolist())),
        "technical_reading_features": readings,
        "estimates": estimates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--mode", choices=("joint", "optics", "structure"), default="joint")
    args = parser.parse_args()
    case = json.loads(args.case.read_text(encoding="utf-8"))
    case["kernel_file"] = str(args.case.parent / case["kernel_file"])
    print(json.dumps(analyze(case, mode=args.mode), indent=2))


if __name__ == "__main__":
    main()
