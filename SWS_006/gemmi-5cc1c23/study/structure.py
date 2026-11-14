"""Prepare O2 displacement measurements in the MD crystal frame."""

import gemmi
import numpy as np

from .model import CELL


SIX = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def matrix(elements):
    result = np.zeros((3, 3))
    for (i, j), value in zip(SIX, elements, strict=True):
        result[i, j] = result[j, i] = value
    return result


def read_structure(representation):
    document = gemmi.cif.read_string(representation["mmcif"])
    return gemmi.make_structure_from_block(document.sole_block())


def ordered_moments(sites, cell, axes, expected_centers):
    frac = np.array(cell.frac.mat.tolist())
    cart = np.array(cell.orth.mat.tolist())
    unnormalize = axes @ cart @ np.diag(np.linalg.norm(frac, axis=1))
    moments = np.empty((len(expected_centers), 3, 3))
    matched = set()
    for site in sites:
        position = axes @ np.array(cell.orthogonalize(site.fract).tolist())
        residual = (position - expected_centers) / CELL.diagonal()
        residual = (residual - np.round(residual)) * CELL.diagonal()
        errors = np.linalg.norm(residual, axis=1)
        index = int(np.argmin(errors))
        if errors[index] > .001 or index in matched:
            raise ValueError("A crystallographic site does not match the resolved oxygen map")
        matched.add(index)
        moments[index] = unnormalize @ matrix(site.aniso.elements_pdb()) @ unnormalize.T
    if len(matched) != len(expected_centers):
        raise ValueError("The resolved oxygen map is incomplete")
    return moments


def displacement_targets(representation, expected_centers):
    model = read_structure(representation)
    model.setup_cell_images()
    small = gemmi.mx_to_sx_structure(model)
    # The downstream consumer reads a small-molecule CIF. Its serialized
    # components use reciprocal-axis normalization, including on image sites.
    expanded = small.get_all_unit_cell_sites()
    export = gemmi.SmallStructure()
    export.cell = small.cell
    export.spacegroup_hm = "P 1"
    for index, site in enumerate(expanded):
        site.label = f"O{index + 1}"
        export.add_site(site)
    decoded = gemmi.make_small_structure_from_block(export.make_cif_block())
    # The small-molecule CIF has the lattice geometry. Registration and the
    # explicit Cartesian orthogonalization are retained in the model metadata.
    return ordered_moments(decoded.sites, small.cell, np.asarray(representation["lab_to_crystal"]), expected_centers)

