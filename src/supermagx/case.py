"""One-way module wiring for a user-defined magnet calculation."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from .electromagnetics import AxisymmetricMagnetostatics
from .excitation import PiecewiseLinearCurrent
from .geometry import Magnet, RadialGroup, RectTurn, radial_stack
from .losses import (coupled_network_loss, eddy_energy,
                     eddy_energy_slab_modes, joule_energy)


def _geometry(spec: Mapping[str, Any]) -> Magnet:
    kind = spec.get("kind")
    if kind == "radial_stack":
        return radial_stack(float(spec["r_start_m"]),
                            (RadialGroup(**entry) for entry in spec["groups"]))
    if kind == "explicit_turns":
        return Magnet(tuple(RectTurn(**entry) for entry in spec["turns"]))
    raise ValueError("geometry.kind must be radial_stack or explicit_turns")


def run_case(case: Mapping[str, Any]) -> dict[str, Any]:
    """Run phase-one components and return JSON-serializable results.

    Thermal and mechanics remain separate future solvers. Their status is
    explicit in the output so an energy estimate cannot be read as a hotspot
    or a Lorentz-force proxy as a stress result.
    """
    magnet = _geometry(case["geometry"])
    pulse = case["pulse"]
    waveform = PiecewiseLinearCurrent(tuple(pulse["times_s"]),
                                      pulse["currents_a"])
    rho = float(case["resistivity_ohm_m"])
    em_options = case.get("electromagnetics", {})
    model = AxisymmetricMagnetostatics(
        magnet, rtol=float(em_options.get("rtol", 1e-6)),
        max_order=int(em_options.get("max_order", 64)))
    if "field_currents_a" in case:
        field_current = dict(case["field_currents_a"])
        field_time = None
    else:
        amps = np.array([waveform.values(z) for z in magnet.zones])
        index = int(np.argmax(np.linalg.norm(amps, axis=0)))
        field_current = {z: float(amps[i, index])
                         for i, z in enumerate(magnet.zones)}
        field_time = waveform.times[index]
    probes = []
    for r, z in case.get("probes_m", []):
        field = model.field(float(r), float(z), field_current)
        probes.append({"r_m": float(r), "z_m": float(z),
                       "br_t": field.br, "bz_t": field.bz,
                       "quadrature_order": field.order,
                       "quadrature_relative_change": field.relative_change})
    qj = joule_energy(magnet, waveform, resistivity=rho)
    losses: dict[str, Any] = {
        "joule_j": qj.total_j,
        "joule_by_zone_j": qj.by_zone,
        "resistance_by_zone_ohm": qj.resistances_ohm,
    }
    if "eddy" in case:
        eddy = case["eddy"]
        kind = eddy.get("model", "calibrated_one_pole")
        if kind == "coupled_network":
            settings = {"resistivity": rho,
                        "mutual_rtol": float(eddy.get("mutual_rtol", 1e-4)),
                        "mutual_max_order": int(eddy.get("mutual_max_order", 24))}
            qn = coupled_network_loss(
                magnet, waveform, radial_cells=eddy.get("radial_cells", 2),
                axial_cells=eddy.get("axial_cells", 2), **settings)
            losses.update(joule_j=qn.transport_j,
                          joule_by_zone_j=qn.transport_by_zone_j,
                          resistance_by_zone_ohm=qn.resistance_by_zone_ohm,
                          eddy_j=qn.eddy_j,
                          eddy_before_tail_j=qn.eddy_before_tail_j,
                          eddy_tail_j=qn.tail_j, eddy_model=qn.model,
                          total_j=qn.total_j,
                          network_cells=qn.cells, network_modes=qn.modes,
                          shortest_eddy_tau_s=qn.shortest_tau_s,
                          longest_eddy_tau_s=qn.longest_tau_s)
            if "refinement_meshes" in eddy:
                history = []
                for shape in eddy["refinement_meshes"]:
                    if len(shape) != 2:
                        raise ValueError("each refinement mesh needs [radial, axial]")
                    trial = coupled_network_loss(
                        magnet, waveform, radial_cells=shape[0],
                        axial_cells=shape[1], **settings)
                    previous = history[-1] if history else None
                    history.append({
                        "radial_cells": shape[0], "axial_cells": shape[1],
                        "eddy_j": trial.eddy_j, "total_j": trial.total_j,
                        "tail_j": trial.tail_j,
                        "relative_eddy_change": (abs(trial.eddy_j-previous["eddy_j"])
                                                  / max(abs(trial.eddy_j), 1e-30)
                                                  if previous else None),
                        "relative_total_change": (abs(trial.total_j-previous["total_j"])
                                                   / max(abs(trial.total_j), 1e-30)
                                                   if previous else None),
                    })
                losses["network_refinement"] = history
        elif kind == "slab_modes":
            qe = eddy_energy_slab_modes(
                magnet, waveform, resistivity=rho,
                modes=eddy.get("modes", 16),
                field_order=int(eddy.get("field_order", 6)))
        elif kind == "calibrated_one_pole":
            qe = eddy_energy(magnet, waveform, resistivity=rho,
                             gamma=float(eddy["gamma"]),
                             field_order=int(eddy.get("field_order", 6)))
        elif kind != "coupled_network":
            raise ValueError("eddy.model must be coupled_network, slab_modes "
                             "or calibrated_one_pole")
        if kind != "coupled_network":
            losses.update(eddy_j=qe.total_j, eddy_tail_j=qe.tail_j,
                          eddy_by_zone_j=qe.by_zone, eddy_model=qe.model,
                          total_j=qj.total_j+qe.total_j)
    if case.get("compute_mutual", False):
        m = model.mutual_matrix(
            rtol=float(em_options.get("mutual_rtol", 1e-4)),
            max_order=int(em_options.get("mutual_max_order", 24)))
        mutual = m.tolist()
    else:
        mutual = None
    return {
        "geometry": {"turns": len(magnet.turns),
                     "zones": list(magnet.zones)},
        "field_current_a": field_current,
        "field_time_s": field_time,
        "field_probes": probes,
        "mutual_matrix_h": mutual,
        "losses": losses,
        "thermal_status": "not solved",
        "mechanics_status": "not solved",
    }
