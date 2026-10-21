"""Screen 5--8 mm radial CuNb conductors against the op4 winding family.

Usage: python scan_narrow.py
Outputs one auditable row per geometry; this is a semi-analytic screen, not a
structural or thermal FEM result. The op4 reference is recomputed on startup.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from supermagx.electromagnetics import AxisymmetricMagnetostatics  # noqa: E402
from supermagx.electromagnetics.axisymmetric import (  # noqa: E402
    _self_mutual,
    loop_mutual,
)
from supermagx.geometry import RadialGroup, radial_stack  # noqa: E402
from supermagx.losses.eddy import _filtered_slope_integral  # noqa: E402

RHO = 6e-9
DENSITY = 8900.0
CP = 200.0
MU0 = 4e-7 * math.pi
GAMMA = 0.75216163  # op3 single-pulse calibrated one-pole time factor
B_TARGET = 20.0
FIELD_CURRENT_FACTOR = 1.0154074  # op3-to-FEM centre-field correction
FLAT = .001
R_START = .125
ZONE_LAYERS = (2, 2, 1, 1, 1, 1, 1, 1)
ZYLON_MM = (18, 16, 14, 12, 10, 8, 8, 25)
LAYOUTS = {
    "N340": (33, 33, 34, 34, 34, 35, 35, 36),
    "N313": (30, 30, 31, 31, 32, 32, 33, 34),
}


def make_geometry(width_mm: float, height_mm: float, layout: str):
    if width_mm < 5:
        raise ValueError("radial CuNb width must be at least 5 mm")
    if layout not in LAYOUTS:
        raise ValueError(f"unknown layout {layout}")
    width = width_mm / 1000
    height = height_mm / 1000
    pitch_r = width + .001
    pitch_z = height + .001
    counts = LAYOUTS[layout]
    groups = [
        RadialGroup(
            zone=k + 1,
            layers=ZONE_LAYERS[k],
            turns_per_layer=counts[k],
            metal_width=width,
            metal_height=height,
            radial_pitch=pitch_r,
            axial_pitch=pitch_z,
            gap_after=(ZYLON_MM[k] + (3 if k < 7 else 0)) / 1000,
        )
        for k in range(8)
    ]
    magnet = radial_stack(R_START, groups)
    r_outer = R_START + 10 * pitch_r + (sum(ZYLON_MM) + 7 * 3) / 1000
    axial = max(counts) * pitch_z + .010
    return magnet, r_outer, axial


def zone_mutual_line(magnet) -> np.ndarray:
    """Filament pair mutuals with finite rectangular self terms.

    The op4 7.5 x 12.5 mm reference row sums agree within 0.04%; near-pair
    finite-cross-section corrections are still a model uncertainty.
    """
    turns = magnet.turns
    radius = np.array([turn.r_center for turn in turns])
    z = np.array([turn.z_center for turn in turns])
    zone = np.array([turn.zone - 1 for turn in turns])
    dz = z[:, None] - z[None, :]
    np.fill_diagonal(dz, 1e-6)
    loops = loop_mutual(radius[:, None], radius[None, :], dz)
    np.fill_diagonal(loops, [_self_mutual(turn, 12) for turn in turns])
    winding = np.eye(8)[zone]
    matrix = winding.T @ loops @ winding
    return (matrix + matrix.T) / 2


def pulse_energy(
    ramp: float,
    current: float,
    resistance: float,
    kr: np.ndarray,
    kz: np.ndarray,
    tau_r: float,
    tau_z: float,
) -> tuple[float, float, float, float, np.ndarray]:
    durations = np.array([ramp, FLAT, ramp])
    slopes = np.array([1 / ramp, 0, -1 / ramp])
    ir, tail_r = _filtered_slope_integral(slopes, durations, tau_r)
    iz, tail_z = _filtered_slope_integral(slopes, durations, tau_z)
    qj = current**2 * resistance * (FLAT + 2 * ramp / 3)
    qeddy_turn = (kr * (ir + tail_r) + kz * (iz + tail_z)) / RHO
    qeddy = float(qeddy_turn.sum())
    tail = float((kr.sum() * tail_r + kz.sum() * tail_z) / RHO)
    return qj + qeddy, qj, qeddy, tail, qeddy_turn


def evaluate(item: tuple[float, float, str]) -> dict:
    width_mm, height_mm, layout = item
    magnet, r_outer, axial = make_geometry(width_mm, height_mm, layout)
    turns = magnet.turns
    model = AxisymmetricMagnetostatics(magnet, rtol=1e-7)
    b_per_amp = model.field(0, 0, {zone: 1.0 for zone in magnet.zones}).bz
    current = FIELD_CURRENT_FACTOR * B_TARGET / b_per_amp
    resistance = np.array([turn.resistance(RHO) for turn in turns])
    volume = np.array([turn.metal_volume for turn in turns])
    radius = np.array([turn.r_center for turn in turns])
    zone = np.array([turn.zone - 1 for turn in turns])
    kr = np.zeros(len(turns))
    kz = np.zeros(len(turns))
    fr = np.zeros(len(turns))
    fz = np.zeros(len(turns))
    for k, turn in enumerate(turns):
        br_by_zone, bz_by_zone = model.screening_coefficients(
            turn.r_center, turn.z_center, order=4,
        )
        br = current * float(br_by_zone.sum())
        bz = current * float(bz_by_zone.sum())
        kr[k] = volume[k] * turn.width**2 * bz**2 / 12
        kz[k] = volume[k] * turn.height**2 * br**2 / 12
        j = current / (turn.width * turn.height)
        fr[k] = j * bz
        fz[k] = -j * br
    tau_r = GAMMA * MU0 * (width_mm / 1000)**2 / (math.pi**2 * RHO)
    tau_z = GAMMA * MU0 * (height_mm / 1000)**2 / (math.pi**2 * RHO)
    resistance_total = float(resistance.sum())

    def objective(ramp: float) -> float:
        return pulse_energy(ramp, current, resistance_total, kr, kz,
                            tau_r, tau_z)[0]

    optimum = minimize_scalar(objective, bounds=(.0005, .08),
                              method="bounded", options={"xatol": 1e-8})
    candidates = [(float(optimum.x), float(optimum.fun)),
                  (.0005, objective(.0005)), (.08, objective(.08))]
    ramp, _ = min(candidates, key=lambda pair: pair[1])
    total, qj, qeddy, tail, qeddy_turn = pulse_energy(
        ramp, current, resistance_total, kr, kz, tau_r, tau_z,
    )
    zone_resistance = np.bincount(zone, weights=resistance, minlength=8)
    zone_qj = current**2 * zone_resistance * (FLAT + 2 * ramp / 3)
    zone_qeddy = np.bincount(zone, weights=qeddy_turn, minlength=8)
    mutual = zone_mutual_line(magnet)
    l_common = mutual.sum(axis=1)
    voltage = current * l_common / ramp + current * zone_resistance
    qturn = current**2 * resistance * (FLAT + 2 * ramp / 3) + qeddy_turn
    turn_mass = volume * DENSITY
    counts = LAYOUTS[layout]
    pitch_z_mm = height_mm + 1
    axial_winding_mm = np.array(counts) * pitch_z_mm
    zylon_area = float(sum(
        ZYLON_MM[k] * (axial_winding_mm[k] + 10) for k in range(8)))
    gfrp_area = float(sum(
        3 * (max(axial_winding_mm[k], axial_winding_mm[k + 1]) + 10)
        for k in range(7)))
    pitch_area = (width_mm + 1) * (height_mm + 1)
    metal_area = len(turns) * width_mm * height_mm
    insulation_area = len(turns) * (pitch_area - width_mm * height_mm)
    material_area = metal_area + insulation_area + zylon_area + gfrp_area
    return {
        "name": f"{layout}_W{width_mm:g}_H{height_mm:g}",
        "layout": layout, "width_mm": width_mm, "height_mm": height_mm,
        "turns": len(turns),
        "zone_turns": "/".join(str(int(v)) for v in np.bincount(zone, minlength=8)),
        "current_kA": current / 1000,
        "b_center_calibrated_t": B_TARGET,
        "outer_diameter_mm": 2000 * r_outer,
        "axial_length_mm": 1000 * axial,
        "ramp_ms": 1000 * ramp,
        "voltage_max_kV": float(max(abs(voltage))) / 1000,
        "voltage_zone": int(np.argmax(abs(voltage))) + 1,
        "energy_MJ": float(.5 * current**2 * mutual.sum()) / 1e6,
        "qtotal_kJ": total / 1000,
        "qj_kJ": qj / 1000,
        "qeddy_kJ": qeddy / 1000,
        "qtail_kJ": tail / 1000,
        "qzone_max_kJ": float(max(zone_qj + zone_qeddy)) / 1000,
        "fr_max_GNm3": float(max(abs(fr))) / 1e9,
        "fz_max_GNm3": float(max(abs(fz))) / 1e9,
        "hoop_proxy_MPa": float(max(abs(radius * fr))) / 1e6,
        "fr_integrated_total_MN": float(sum(fr * volume)) / 1e6,
        "fz_abs_integrated_total_MN": float(sum(abs(fz) * volume)) / 1e6,
        "metal_J_Amm2": current / (width_mm * height_mm),
        "cell_J_Amm2": current / pitch_area,
        "all_material_J_Amm2": len(turns) * current / material_area,
        "envelope_J_Amm2": len(turns) * current
        / ((r_outer - R_START) * 1000 * axial * 1000),
        "metal_area_mm2": metal_area,
        "insulation_area_mm2": insulation_area,
        "zylon_area_mm2": zylon_area,
        "gfrp_area_mm2": gfrp_area,
        "all_material_area_mm2": material_area,
        "conductor_length_m": float(sum(2 * math.pi * radius)),
        "metal_mass_kg": float(turn_mass.sum()),
        "dTavg_K": total / (float(turn_mass.sum()) * CP),
        "dTturn_max_K": float(max(qturn / (turn_mass * CP))),
        "resistance_sum_mohm": resistance_total * 1000,
        "zone_resistance_mohm": "/".join(f"{x*1000:.9g}" for x in zone_resistance),
        "zone_mutual_rowsum_mH": "/".join(f"{x*1000:.9g}" for x in l_common),
        "zone_voltage_kV": "/".join(f"{x/1000:.9g}" for x in voltage),
        "zone_qj_kJ": "/".join(f"{x/1000:.9g}" for x in zone_qj),
        "zone_qeddy_kJ": "/".join(f"{x/1000:.9g}" for x in zone_qeddy),
        "kr": float(kr.sum()), "kz": float(kz.sum()),
        "tau_r_ms": 1000 * tau_r, "tau_z_ms": 1000 * tau_z,
    }


def validate_reference(row: dict) -> dict:
    expected = {
        "current_kA": 29.78411723134254,
        "qtotal_kJ": 703.478662752696,
        "voltage_max_kV": 5.830411290942461,
        "fr_max_GNm3": 6.344346612550876,
        "hoop_proxy_MPa": 820.0067996722007,
        "energy_MJ": 9.923719917176491,
    }
    errors = {key: 100 * (row[key] / value - 1)
              for key, value in expected.items()}
    # The old source is missing. This check protects the reconstructed
    # implementation against accidentally changing the op4 reference case.
    if max(abs(value) for value in errors.values()) > .5:
        raise ArithmeticError(f"op4 reference reconstruction drifted: {errors}")
    return errors


def pareto(rows: list[dict]) -> list[dict]:
    metrics = ("qtotal_kJ", "fr_max_GNm3", "voltage_max_kV",
               "dTturn_max_K", "metal_mass_kg")
    feasible = [row for row in rows if row["current_kA"] < 40
                and row["voltage_max_kV"] < 40]
    return [row for row in feasible if not any(
        other is not row
        and all(other[k] <= row[k] for k in metrics)
        and any(other[k] < row[k] for k in metrics)
        for other in feasible
    )]


def force_nonworsening_shortlist(rows: list[dict]) -> list[dict]:
    """Compare new widths with the existing op4 screened force/heat proxies.

    This is a comparison rule, not an assertion of <700 MPa von Mises stress.
    """
    return sorted((row for row in rows
                   if row["fr_max_GNm3"] <= 6.344346612550876
                   and row["hoop_proxy_MPa"] <= 820.0067996722007
                   and row["dTturn_max_K"] <= 22.741529989842892
                   and row["current_kA"] < 40
                   and row["voltage_max_kV"] < 40),
                  key=lambda row: row["qtotal_kJ"])


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    here = Path(__file__).resolve().parent
    start = perf_counter()
    reference = evaluate((7.5, 12.5, "N340"))
    errors = validate_reference(reference)
    (here / "reference_check.json").write_text(
        json.dumps({"reference": reference, "relative_errors_percent": errors},
                   indent=2) + "\n", encoding="utf-8")
    print("op4 reference reconstruction, max error %:",
          max(abs(value) for value in errors.values()), flush=True)
    items = [(w, h, layout)
             for layout in LAYOUTS
             for w in (5, 5.5, 6, 6.5, 7, 7.5, 8)
             for h in (12, 12.5, 13, 13.5, 14)]
    rows = [reference]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(evaluate, item): item for item in items
                   if item != (7.5, 12.5, "N340")}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(row["name"], "Q(kJ)", round(row["qtotal_kJ"], 2),
                  "Fr(GN/m3)", round(row["fr_max_GNm3"], 2),
                  flush=True)
    rows.sort(key=lambda row: (row["layout"], row["width_mm"],
                               row["height_mm"]))
    front = pareto(rows)
    shortlist = force_nonworsening_shortlist(rows)
    write_csv(here / "scan_results.csv", rows)
    write_csv(here / "pareto.csv", front)
    write_csv(here / "force_nonworsening_shortlist.csv", shortlist)
    print("cases", len(rows), "Pareto", len(front),
          "force-nonworsening", len(shortlist),
          "best", shortlist[0]["name"],
          "elapsed_s", round(perf_counter() - start, 1), flush=True)


if __name__ == "__main__":
    main()
