"""Metrics for a predicted protein–ligand complex and its crystal counterpart."""

import os
import shutil
import tempfile
from pathlib import Path
from typing import ClassVar

import numpy as np
from biotite.structure.io import pdbx
from pydantic import Field

from blocks.metrics.chirality import ChiralityCheckBlock
from blocks.networks.prediction.helpers import _match_aligned_rmsd
from blocks.networks.prediction.validation.lddt import rf3_allatom_lddt
from blocks.util_blocks.remap_ligand import RemapLigandBlock
from blocks.utils.geometry import get_ligand_mask
from chaperone import InlineMetricBlock, JobArgs, PipelineDataHandler
from chaperone.data import make_row
from chaperone.io import parse_to_atom_array


def _reference_parts(array, ligand):
    protein = array.is_residue & ~array.hetero
    selected = get_ligand_mask(array)
    return array[protein], array[selected]


def _write_cif(array, path):
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, array, data_block="reference", include_bonds=True)
    cif.write(path)


def _prediction_path(prediction, work):
    # Present local files to the pipeline with relative paths.
    try:
        return os.path.relpath(prediction)
    except ValueError:
        local = work / Path(prediction).name
        shutil.copyfile(prediction, local)
        return os.path.relpath(local)


def _interface_score(ref_protein, ref_ligand, pred_protein, pred_ligand):
    pred_ligand = pred_ligand.copy()
    pred_ligand.set_annotation("chain_id", np.repeat(ref_ligand.chain_id[0], len(pred_ligand)))
    pred_ligand.res_id[:] = ref_ligand.res_id[0]
    reference = ref_protein + ref_ligand
    prediction = pred_protein + pred_ligand
    reference.set_annotation(
        "chain_iid", np.repeat(["protein", "ligand"], [len(ref_protein), len(ref_ligand)])
    )
    report = rf3_allatom_lddt(
        reference,
        prediction,
        interfaces_to_score=[("protein", "ligand", "protein-ligand")],
        pn_units_to_score=[],
    )
    return float(report["rf3_by_type_lddt"][0]["best_of_1_lddt"])


def audit_ligand_pose(prediction, reference, *, ligand, core_atoms=None) -> dict:
    """Audit one pair of mmCIF files for the specified CCD component.

    Pass reference atom names as ``core_atoms`` to also score that subset.
    """
    crystal = parse_to_atom_array(reference, assign_types=True)
    ref_protein, ref_ligand = _reference_parts(crystal, ligand)
    if not len(ref_protein) or not len(ref_ligand):
        raise ValueError("The reference must contain protein and the requested component")

    predicted = parse_to_atom_array(prediction, assign_types=True)
    pred_protein = predicted[predicted.is_residue & ~predicted.hetero]
    pred_ligand = predicted[predicted.is_ligand]
    if len(pred_ligand) != len(ref_ligand):
        raise ValueError("The ligand atom counts differ")
    ref_indices = {str(name): index for index, name in enumerate(ref_ligand.atom_name)}
    core = None
    if core_atoms is not None:
        core_atoms = list(core_atoms)
        if not core_atoms or not set(core_atoms) <= ref_indices.keys():
            raise ValueError("Core names must identify a nonempty subset of the component")
        core = [ref_indices[name] for name in dict.fromkeys(core_atoms)]

    temp_root = Path(tempfile.gettempdir())
    if temp_root.drive.lower() != Path.cwd().drive.lower():
        temp_root = Path.cwd()
    with tempfile.TemporaryDirectory(dir=temp_root) as directory:
        work = Path(directory)
        ref_cif = str(work / "reference.cif")
        _write_cif(ref_ligand, ref_cif)
        input_path = _prediction_path(prediction, work)
        handler = PipelineDataHandler([make_row(name="prediction", path=input_path)])
        remap = RemapLigandBlock(reference_cif=ref_cif, target_name=ligand, core_atoms=core_atoms)
        renamed = remap.execute(handler, work / "mapped", 0)
        loaded = next(iter(renamed.get_data(which="last", load=True)))
        mapped = loaded.atom_array[loaded.atom_array.res_name == ligand]
        if len(mapped) != len(pred_ligand) or not np.array_equal(mapped.coord, pred_ligand.coord):
            raise ValueError("Mapped ligand atoms must retain their coordinates")
        atom_map = dict(zip(map(str, pred_ligand.atom_name), map(str, mapped.atom_name)))
        if len(atom_map) != len(ref_indices) or set(atom_map.values()) != ref_indices.keys():
            raise ValueError("The component atom map must be complete")

        match = [ref_indices[name] for name in atom_map.values()]
        match_rmsd = loaded.metrics.get(remap.name, {}).get("match_rmsd")
        full_rmsd = match_rmsd if core is None else None
        if full_rmsd is None:
            full_rmsd = _match_aligned_rmsd(match, ref_ligand.coord, pred_ligand.coord)
        core_rmsd = None
        if core is not None:
            core_rmsd = match_rmsd
            if match_rmsd is None or len(core) < 3:
                core_rmsd = _match_aligned_rmsd(match, ref_ligand.coord, pred_ligand.coord, core)

        chirality = ChiralityCheckBlock(reference_cif=ref_cif)
        context = chirality.setup(renamed, work, 0)
        agreement = chirality.compute_metrics(loaded, context=context)["correct_chirality"]
        interface = _interface_score(ref_protein, ref_ligand, pred_protein, mapped)

    return {
        "ligand_rmsd": float(full_rmsd),
        "core_rmsd": None if core_rmsd is None else float(core_rmsd),
        "atom_map": atom_map,
        "correct_chirality": float(agreement),
        "interface_lddt": interface,
    }


class PoseAuditBlock(InlineMetricBlock):
    """Audit the last structure of each pipeline row against a crystal file."""

    _name = "pose_audit"
    _display_name = "Ligand Pose Audit"
    _description = "Compare predicted complexes with a crystal reference"
    _resources: ClassVar[JobArgs] = JobArgs(partition="cpu")

    reference_cif: str = Field(..., description="Path to the crystal mmCIF file.")
    ligand: str = Field(..., description="CCD component identifier.")
    core_atoms: list[str] | None = Field(None, description="Reference atom names for the core score.")

    def get_data_kwargs(self):
        return {"load": False, "which": "last"}

    def compute_metrics(self, loaded, *, context=None, index=0):
        return audit_ligand_pose(
            loaded.path,
            self.reference_cif,
            ligand=self.ligand,
            core_atoms=self.core_atoms,
        )
