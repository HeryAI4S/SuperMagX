"""Create a compact, auditable record from local COMSOL and cell-model runs.

Run after run_validation.py and compare.py. Native MPH and raw traces remain
local under output/; this JSON records their checksums and the measured data.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
CASES = {
    "single_a": ("single_a_ultra", "single_a_fine", "single_a_coarse"),
    "narrow": ("narrow_ultra", "narrow_fine", "narrow_coarse"),
    "highrho": ("highrho_ultra", "highrho_fine", "highrho_coarse"),
    "two": ("two_ultra", "two_fine", "two_coarse"),
}


def read_json(path: Path) -> dict | list:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_sha256(path: Path) -> str:
    """Hash tracked UTF-8 source with Git's LF line-ending convention."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def observed_order(values: dict[int, float]) -> float | None:
    # Positive monotone differences are needed for the asymptotic ratio.
    ratio = (values[16] - values[20]) / (values[20] - values[24])
    if ratio <= 0:
        return None

    def residual(order: float) -> float:
        return ((16**-order - 20**-order)
                / (20**-order - 24**-order) - ratio)

    if residual(.01) * residual(6) > 0:
        return None
    return float(brentq(residual, .01, 6))


def extrapolate(values: dict[int, float]) -> dict:
    # The predicted O(h^2) limit uses only the two finest edge-graded meshes.
    limit = float(values[24] + (values[24] - values[20])
                  / ((24 / 20) ** 2 - 1))
    return {"order_from_16_20_24": observed_order(values),
            "n24": float(values[24]), "n20": float(values[20]),
            "limit_p2": limit,
            "p2_correction": limit - float(values[24])}


def main() -> None:
    n16 = {row["case"]: row for row in
           read_json(OUTPUT / "edge_cosine_16.json")}
    n20_24 = {row["case"]: row for row in
              read_json(OUTPUT / "comparisons_edge_cosine.json")}
    record = {
        "description": "Axisymmetric solid-ring magnetic diffusion; conditional p=2 mesh extrapolation",
        "comsol_version": "6.4.0.429",
        "power_expression": "IntSurface(2*pi*r*mf.Qh)",
        "source_hash_convention": "SHA-256 of Git-normalized LF source bytes",
        "source_sha256": {
            str(path.relative_to(HERE.parent.parent)): source_sha256(path)
            for path in (
                HERE / "CoupledLossBenchmark.java",
                HERE / "run_validation.py",
                HERE / "compare.py",
                HERE / "summarize.py",
                HERE.parent.parent / "src/supermagx/geometry.py",
                HERE.parent.parent / "src/supermagx/electromagnetics/axisymmetric.py",
                HERE.parent.parent / "src/supermagx/losses/network.py",
                HERE.parent.parent / "src/supermagx/losses/eddy.py",
            )},
        "cases": {},
    }
    for short, (reference, fine, coarse) in CASES.items():
        fem = {level: read_json(OUTPUT / f"{name}_summary.json")
               for level, name in (("reference", reference),
                                   ("fine", fine), ("coarse", coarse))}
        first_comparison = n16.get(reference) or n16[fine]
        rows = [first_comparison["network"][0]]
        comparison = n20_24.get(reference)
        if comparison is None:
            # compare.py may have been run against the previous FEM refinement.
            comparison = n20_24[fine]
        rows.extend(comparison["network"])
        by_n = {row["cells_per_axis"]: row for row in rows}
        if set(by_n) != {16, 20, 24}:
            raise ValueError(f"missing mesh level for {short}")
        transport = float(first_comparison["fem"]["transport_j_analytic"])
        fem_total = float(fem["reference"]["total_heat_to_end_j"])
        fem_tail = float(fem["reference"]["tail_heat_to_end_j"])
        fem_eddy = fem_total - transport
        eddy = extrapolate({n: row["eddy_j"] for n, row in by_n.items()})
        tail = extrapolate({n: row["tail_j"] for n, row in by_n.items()})
        total = transport + eddy["limit_p2"]
        ref_path = HERE.parent.parent / fem["reference"]["mph"]
        record["cases"][short] = {
            "reference_spec": fem["reference"]["spec"],
            "network_source_cases": [first_comparison["case"],
                                     comparison["case"]],
            "native_mph_sha256": sha256(ref_path),
            "fem": {"transport_analytic_j": transport,
                    "eddy_by_subtraction_j": fem_eddy,
                    "tail_to_end_j": fem_tail,
                    "total_to_end_j": fem_total,
                    "end_power_w": fem["reference"]["end_power_w"],
                    "mesh_time_change_total_percent":
                        100 * (fem_total / fem["fine"]["total_heat_to_end_j"] - 1),
                    "mesh_time_change_tail_percent":
                        100 * (fem_tail / fem["fine"]["tail_heat_to_end_j"] - 1),
                    "coarse_to_fine_total_percent":
                        100 * (fem["fine"]["total_heat_to_end_j"]
                               / fem["coarse"]["total_heat_to_end_j"] - 1)},
            "edge_mesh": {"eddy": eddy, "tail": tail},
            "extrapolated": {
                "transport_j": transport,
                "eddy_j": eddy["limit_p2"],
                "tail_j": tail["limit_p2"],
                "total_j": total,
                "relative_total_error_percent": 100 * (total / fem_total - 1),
                "relative_eddy_error_percent":
                    100 * (eddy["limit_p2"] / fem_eddy - 1),
                "relative_tail_error_percent":
                    100 * (tail["limit_p2"] / fem_tail - 1),
            },
        }
    dest = HERE / "benchmark_results.json"
    dest.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    for name, case in record["cases"].items():
        ext = case["extrapolated"]
        print(name, "error %: total", round(ext["relative_total_error_percent"], 3),
              "eddy", round(ext["relative_eddy_error_percent"], 3),
              "tail", round(ext["relative_tail_error_percent"], 3))
    print(dest)


if __name__ == "__main__":
    main()
