"""Match the exact ring geometry and pulse between COMSOL and SuperMagX."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn
from supermagx.losses import compile_coupled_loss

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"


def compare(name: str, meshes: list[int], distribution: str) -> dict:
    fem = json.loads((OUTPUT / f"{name}_summary.json").read_text(encoding="utf-8"))
    fields = fem["spec"].split(",")
    rho, width, height, ramp = map(float, fields[1:5])
    flat = float(fields[8]) if len(fields) in (9, 10, 13) else .001
    if flat != .001:
        raise ValueError("comparison assumes a 1 ms flat top")
    turns = [RectTurn(.1, .1+width, -height/2, height/2, "zone")]
    times = [0, ramp, ramp+flat, 2*ramp+flat]
    peaks = {"zone": 1000.0}
    delays = {"zone": 0.0}
    if len(fields) == 13:
        r2, z2, scale, delay = map(float, fields[9:13])
        turns.append(RectTurn(r2, r2+width, z2-height/2, z2+height/2,
                              "second"))
        peaks["second"] = scale*1000
        delays["second"] = delay
        times.extend([delay, delay+ramp, delay+ramp+flat,
                      delay+2*ramp+flat])
    magnet = Magnet(tuple(turns))
    times = sorted(set(times))

    def current(t: float, peak: float, delay: float) -> float:
        return peak*(min(max((t-delay)/ramp, 0), 1)
                     - min(max((t-delay-ramp-flat)/ramp, 0), 1))

    waveform = PiecewiseLinearCurrent(
        tuple(times), {zone: tuple(current(t, peak, delays[zone])
                                   for t in times)
                       for zone, peak in peaks.items()})
    transport = sum(turn.resistance(rho)*peaks[turn.zone]**2
                    *(flat+2*ramp/3) for turn in turns)
    fem_eddy = fem["total_heat_to_end_j"]-transport
    rows = []
    for n in meshes:
        operator = compile_coupled_loss(magnet, radial_cells=n,
                                        axial_cells=n,
                                        cell_distribution=distribution,
                                        mutual_rtol=1e-4,
                                        mutual_max_order=32)
        loss = operator.evaluate(waveform, resistivity=rho)
        coefficients = operator.low_frequency_coefficients()
        if len(turns) == 1:
            slow_eddy = float(coefficients.eddy_kernel_times_rho[0, 0]
                              /rho*2*1000**2/ramp)
        else:
            slow_eddy = None
        rows.append({
            "cells_per_axis": n,
            "cells": loss.cells,
            "transport_j": loss.transport_j,
            "eddy_j": loss.eddy_j,
            "tail_j": loss.tail_j,
            "total_j": loss.total_j,
            "slow_ramp_eddy_j": slow_eddy,
            "relative_total_error_vs_fem":
                loss.total_j/fem["total_heat_to_end_j"]-1,
            "relative_eddy_error_vs_fem": loss.eddy_j/fem_eddy-1,
            "relative_tail_error_vs_fem":
                loss.tail_j/fem["tail_heat_to_end_j"]-1,
        })
    return {
        "case": name, "distribution": distribution, "rho_ohm_m": rho,
        "width_m": width, "height_m": height, "ramp_s": ramp,
        "fem": {
            "transport_j_analytic": transport,
            "eddy_j_by_subtraction": fem_eddy,
            "tail_j": fem["tail_heat_to_end_j"],
            "total_j": fem["total_heat_to_end_j"],
            "last_power_w": fem["end_power_w"],
            "samples": fem["samples"],
        },
        "network": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+")
    parser.add_argument("--meshes", nargs="+", type=int, default=[8, 16, 24])
    parser.add_argument("--distribution", choices=["uniform", "edge_cosine"],
                        default="uniform")
    parser.add_argument("--output", type=Path,
                        help="optional output path; keeps earlier mesh sweeps intact")
    args = parser.parse_args()
    all_results = [compare(name, args.meshes, args.distribution)
                   for name in args.cases]
    path = args.output or OUTPUT / f"comparisons_{args.distribution}.json"
    path.write_text(json.dumps(all_results, indent=2) + "\n", encoding="utf-8")
    for result in all_results:
        print(result["case"], result["fem"])
        for row in result["network"]:
            print("  ", row)
    print("Saved", path)
