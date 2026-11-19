# Copyright (c) MyoSuite Authors. All rights reserved.
#
# This source code is licensed under the Apache 2 license found in the
# LICENSE file in the root directory of this source tree.
"""Muscle fatigue of a recorded work task, from its landmark recording.

:func:`task_fatigue` fits a landmark recording (:class:`.retarget.ClipTargets`) on a
model, derives the muscle commands of every frame and returns the fatigue state of
the muscles over repeated performances of the task::

    from myosuite.biomechanics.task_fatigue import task_fatigue

    out = task_fatigue(model, lmap, targets, repetitions=3, rest_s=5.0)
    out["qpos"]  # (T, nq) fitted poses
    out["MF"][-1, -1]  # fatigued fraction of every tracked muscle at the end
"""

from __future__ import annotations

from collections.abc import Sequence

import mujoco
import numpy as np

from myosuite.biomechanics.inverse_dynamics import clip_kinematics
from myosuite.biomechanics.retarget import (
    ClipTargets,
    LandmarkIK,
    LandmarkMap,
    RefitConfig,
)
from myosuite.biomechanics.static_optimization import (
    StaticOptimizationConfig,
    StaticOptimizer,
)
from myosuite.core.muscle_conditions import CumulativeFatigue


def task_fatigue(
    model: mujoco.MjModel,
    lmap: LandmarkMap,
    targets: ClipTargets,
    repetitions: int,
    rest_s: float,
    *,
    muscles: Sequence[int | str] | None = None,
    config: RefitConfig | None = None,
    lowpass_hz: float | None = None,
    so_config: StaticOptimizationConfig | None = None,
) -> dict[str, np.ndarray]:
    """Fatigue state of the muscles over repeated performances of a recorded task.

    Args:
        model: Compiled model.
        lmap: Landmark map of the model (names its bodies and sites).
        targets: The recording: ``frequency`` [Hz], landmark positions ``pos``
            ``(T, L, 3)`` [m], orientations ``quat`` ``(T, O, 4)`` (wxyz), the
            ``seed`` pose and optional ``posture`` values. ``T >= 4``.
        repetitions: Number of performances of the recording (at least 1).
        rest_s: Rest between consecutive performances [s] (not negative).
        muscles: Muscle actuators to track (ids or names); ``None``: every muscle
            in model order.
        config: Solver options of the landmark fit (``None``: library settings).
        lowpass_hz: Cut-off [Hz] of the kinematics filter (``None``: unfiltered).
        so_config: Static optimization settings (``None``: library settings).

    Returns:
        A dict with ``qpos`` ``(T, nq)`` fitted poses; ``muscle_ids`` ``(na,)``
        actuator ids of the tracked muscles, the order of the last axis below;
        ``command`` ``(repetitions, T, na)`` muscle activations of each frame;
        ``MA``, ``MR``, ``MF`` ``(repetitions, T, na)`` active, resting and
        fatigued fractions after each frame.

    Raises:
        ValueError: If ``repetitions``, ``rest_s`` or ``targets.frequency`` is out of
            range, or no muscle is tracked.
    """
    repetitions = int(repetitions)
    if repetitions < 1 or not rest_s >= 0.0 or not targets.frequency > 0.0:
        raise ValueError("need repetitions >= 1, rest_s >= 0 and frequency > 0")
    dt = 1.0 / float(targets.frequency)
    timestep = float(model.opt.timestep)
    step = timestep * max(int(round(dt / timestep)), 1)

    fatigue = CumulativeFatigue(model, muscles=muscles)
    ids = fatigue.muscle_ids
    if ids.size == 0:
        raise ValueError("no muscle actuators to track")

    ik = LandmarkIK(model, lmap, config)
    qpos, _ = ik.fit(targets.pos, targets.quat, targets.seed, targets.posture)
    kin = clip_kinematics(model, qpos, dt, lowpass_hz=lowpass_hz)
    optimizer = StaticOptimizer(model, so_config)

    n_frames = len(qpos)
    shape = (repetitions, n_frames, ids.size)
    command, MA, MR, MF = (np.empty(shape) for _ in range(4))
    for rep in range(repetitions):
        for t in range(n_frames):
            cap = np.full(model.nu, np.inf)
            cap[ids] = 1.0 - fatigue.MF
            frame = optimizer.solve(
                kin.qpos[t], kin.qvel[t], kin.qacc[t], max_activation=cap
            )
            command[rep, t] = frame.activation[ids]
            ma, mr, mf = fatigue.compute_act(command[rep, t], dt=step)
            MA[rep, t], MR[rep, t], MF[rep, t] = ma, mr, mf
        if rep < repetitions - 1 and rest_s > 0.0:
            fatigue.compute_act(np.zeros(ids.size), dt=float(rest_s))
    return {
        "qpos": qpos,
        "muscle_ids": ids.copy(),
        "command": command,
        "MA": MA,
        "MR": MR,
        "MF": MF,
    }
