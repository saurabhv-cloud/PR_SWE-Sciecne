"""Reduced rigid-body molecular model and optical forward observables."""

import numpy as np
from scipy.linalg import expm


TEMPERATURE = 290.0
KT = 0.8314463 * TEMPERATURE  # amu angstrom**2 / ps**2
DT = 0.20  # ps, exact harmonic propagator interval
CELL = np.diag([5.2521, 11.4126, 16.5403])
ROTATIONS = np.array([np.diag(x) for x in [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]]])
TRANSLATIONS = np.array([[0., 0., 0.], [.5, 0., .5], [.5, .5, 0.], [0., .5, .5]])
NEIGHBOR = 3
IMAGE_SHIFT = CELL @ [0, -1, 0]
Q = np.linspace(.05, .95, 5)
ANGLES = np.linspace(-.06, .06, 7)
GRID = (Q, ANGLES, ANGLES)
CONTACT_CUTOFFS = np.array([2.75, 2.80, 2.85])
PAIR_INDEX = np.triu_indices(3)
MASS = np.array([206., 206., 500., 700.])
HADAMARD = .5 * np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]])
TENSOR = np.array([[1., .3, .25], [.3, .5, .35], [.25, .35, .2]])


def unit(vector):
    return vector / np.linalg.norm(vector)


RAMAN_E = np.array([1., 0., 0.])
FLUORESCENCE_E = np.array([0., 0., 1.])
DIPOLE = unit(np.array([.2, .6, .77]))


def directions(geometry):
    nitrogen, _, oxygen = geometry
    oxygen_image = ROTATIONS[NEIGHBOR] @ oxygen + CELL @ TRANSLATIONS[NEIGHBOR] + IMAGE_SHIFT
    along = unit(oxygen_image - nitrogen)
    across = unit(np.cross(along, ROTATIONS[NEIGHBOR] @ along))
    return np.stack([along, across])


def simulate(q, blocks, frames, rng):
    """Draw stationary positions from exact coupled Langevin dynamics in ps."""
    variances = np.array([.064 * q, .064 * (1 - q), .035**2, .035**2])
    translational = np.array([1., 1.35, 1.35, 1.35])
    factors = np.column_stack([translational, translational, np.ones((4, 2))])
    variances = variances * np.mean(1. / factors, axis=0) ** -1 / factors
    velocity_variances = KT / MASS
    update = np.empty((4, 4, 2, 2))
    noise = np.empty_like(update)
    for molecule_mode in range(4):
        for coordinate in range(4):
            variance = variances[molecule_mode, coordinate]
            speed_variance = velocity_variances[coordinate]
            stationary = np.diag([variance, speed_variance])
            rate = np.array([[0., 1.], [-speed_variance / variance, -3.]])
            step = expm(DT * rate)
            update[molecule_mode, coordinate] = step
            noise[molecule_mode, coordinate] = np.linalg.cholesky(stationary - step @ stationary @ step.T)
    scales = np.sqrt(np.stack([variances, np.broadcast_to(velocity_variances, (4, 4))], axis=-1))
    state = rng.normal(size=(blocks, 4, 4, 2)) * scales
    series = np.empty((blocks, frames, 4, 4), dtype=np.float32)
    for index in range(frames):
        series[:, index] = state[..., 0]
        innovation = rng.normal(size=state.shape)
        state = np.einsum("ijab,kijb->kija", update, state) + np.einsum("ijab,kijb->kija", noise, innovation)
    return np.einsum("ij,btjk->btik", HADAMARD, series, dtype=np.float32, casting="unsafe")


def body_rotations(r, f):
    """Rotation Rz(r) Rx(f) in each molecule's local crystal frame."""
    cr, sr = np.cos(r), np.sin(r)
    cf, sf = np.cos(f), np.sin(f)
    result = np.empty((*r.shape, 3, 3))
    result[..., 0, :] = np.stack([cr, -sr * cf, sr * sf], axis=-1)
    result[..., 1, :] = np.stack([sr, cr * cf, -cr * sf], axis=-1)
    result[..., 2, :] = np.stack([np.zeros_like(sf), sf, cf], axis=-1)
    return result


def moments(coordinates, geometry):
    shifts = coordinates[..., :2] @ directions(geometry)
    image_shifts = np.einsum("mij,btmj->btmi", ROTATIONS, shifts)
    # A fixed reference mean at the nominal O2 position is used in both routes.
    return np.einsum("btmi,btmj->bmij", image_shifts, image_shifts) / shifts.shape[1]


def optical_and_contacts(coordinates, geometry, r, f):
    rotations = body_rotations(coordinates[..., 2] + r, coordinates[..., 3] + f)
    crystal_rotations = ROTATIONS @ rotations
    polarization = np.einsum("i,...ij->...j", RAMAN_E, crystal_rotations)
    raman = np.square(np.einsum("...i,ij,...j->...", polarization, TENSOR, polarization))
    fluorescence = np.square(FLUORESCENCE_E @ crystal_rotations @ DIPOLE)
    shifts = coordinates[..., :2] @ directions(geometry)
    nitrogen, hydrogen, oxygen = geometry
    asu_rotation = rotations[..., 0, :, :]
    asu_origin = oxygen + shifts[..., 0, :]
    n_pos = asu_origin + asu_rotation @ (nitrogen - oxygen)
    h_pos = asu_origin + asu_rotation @ (hydrogen - oxygen)
    o_pos = (oxygen + shifts[..., NEIGHBOR, :]) @ ROTATIONS[NEIGHBOR].T
    o_pos += CELL @ TRANSLATIONS[NEIGHBOR] + IMAGE_SHIFT
    distance = np.linalg.norm(n_pos - o_pos, axis=-1)
    hn = n_pos - h_pos
    ho = o_pos - h_pos
    cosine = np.sum(hn * ho, axis=-1) / (np.linalg.norm(hn, axis=-1) * np.linalg.norm(ho, axis=-1))
    contacts = (distance[..., None] <= CONTACT_CUTOFFS) & (cosine[..., None] <= np.cos(np.deg2rad(150)))
    optics = np.column_stack([raman.mean(axis=(1, 2)), fluorescence.mean(axis=(1, 2))])
    interactions = np.column_stack([contacts.mean(axis=1), distance.mean(axis=1), np.square(distance).mean(axis=1)])
    return optics, interactions


def oxygen_centers(geometry):
    return ROTATIONS @ geometry[2] + TRANSLATIONS @ CELL.T
