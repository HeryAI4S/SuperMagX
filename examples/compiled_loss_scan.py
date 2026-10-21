"""Reuse one magnetic-diffusion operator across ramp and resistivity scans.

Run: python examples/compiled_loss_scan.py
The numbers describe the idealized ring below, not an APEX design candidate.
"""

from __future__ import annotations

from time import perf_counter

import numpy as np

from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn
from supermagx.losses import compile_coupled_loss


def trapezoid(ramp: float) -> PiecewiseLinearCurrent:
    return PiecewiseLinearCurrent(
        (0.0, ramp, ramp + .001, 2*ramp + .001),
        {"supply": (0.0, 1000.0, 1000.0, 0.0)},
    )


magnet = Magnet((RectTurn(.1, .12, -.01, .01, "supply"),))
start = perf_counter()
operator = compile_coupled_loss(
    magnet, radial_cells=16, axial_cells=16,
    cell_distribution="edge_cosine",
)
build_s = perf_counter() - start

for rho in (3e-9, 6e-9, 12e-9):
    start = perf_counter()
    candidates = [
        (ramp, operator.evaluate(trapezoid(float(ramp)), resistivity=rho))
        for ramp in np.geomspace(.001, .1, 100)
    ]
    evaluate_s = perf_counter() - start
    ramp, result = min(candidates, key=lambda row: row[1].total_j)
    boundary = ramp == candidates[0][0] or ramp == candidates[-1][0]
    print(
        f"rho={rho:.3g} ohm m; ramp={ramp:.6g} s; "
        f"Q={result.total_j:.9g} J; transport={result.transport_j:.9g} J; "
        f"eddy={result.eddy_j:.9g} J; tail={result.tail_j:.9g} J; "
        f"at_scan_boundary={boundary}; "
        f"operator build={build_s:.3g} s; 100 evaluations={evaluate_s:.3g} s"
    )
