"""Compile and run reproducible COMSOL 6.4 loss benchmarks.

COMSOL batch does not pass Java main arguments and its Java security policy
rejects ordinary file/env reads. The case specification is therefore baked
into a temporary compiled copy of the Java source for each run. Raw stdout
and native MPH files stay under ignored ``validation/comsol/output``.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
COMSOL = Path(r"C:\software\comsol6.4\COMSOL64\Multiphysics\bin\win64")
OUTPUT = HERE / "output"
GENERATED = HERE / "generated"
CASES = {
    "single_a_coarse": "single_a_coarse,6e-9,0.02,0.02,0.01,0.001,0.06,0.00025",
    "single_a_fine": "single_a_fine,6e-9,0.02,0.02,0.01,0.0005,0.06,0.000125",
    "single_a_ultra": "single_a_ultra,6e-9,0.02,0.02,0.01,0.00025,0.06,0.0000625",
    "single_a_timeonly": "single_a_timeonly,6e-9,0.02,0.02,0.01,0.0005,0.06,0.0000625",
    "single_a_meshonly": "single_a_meshonly,6e-9,0.02,0.02,0.01,0.00025,0.06,0.000125",
    "single_a_farfield": "single_a_farfield,6e-9,0.02,0.02,0.01,0.0005,0.06,0.000125,0.001,2.0",
    "narrow_coarse": "narrow_coarse,6e-9,0.008,0.016,0.003,0.0004,0.03,0.0001",
    "narrow_fine": "narrow_fine,6e-9,0.008,0.016,0.003,0.0002,0.03,0.00005",
    "narrow_ultra": "narrow_ultra,6e-9,0.008,0.016,0.003,0.0001,0.03,0.000025",
    "highrho_coarse": "highrho_coarse,1.2e-8,0.02,0.02,0.005,0.00075,0.04,0.0001",
    "highrho_fine": "highrho_fine,1.2e-8,0.02,0.02,0.005,0.000375,0.04,0.00005",
    "highrho_ultra": "highrho_ultra,1.2e-8,0.02,0.02,0.005,0.0001875,0.04,0.000025",
    "two_coarse": "two_coarse,6e-9,0.02,0.02,0.01,0.001,0.07,0.00025,0.001,0.15,0.015,0.7,0.003",
    "two_fine": "two_fine,6e-9,0.02,0.02,0.01,0.0005,0.07,0.000125,0.001,0.15,0.015,0.7,0.003",
    "two_ultra": "two_ultra,6e-9,0.02,0.02,0.01,0.00025,0.07,0.0000625,0.001,0.15,0.015,0.7,0.003",
}


def run_case(case_name: str) -> dict:
    spec = CASES[case_name]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    GENERATED.mkdir(parents=True, exist_ok=True)
    template = (HERE / "CoupledLossBenchmark.java").read_text(encoding="utf-8")
    source = re.sub(r'static final String DEFAULT_SPEC = "[^"]+";',
                    f'static final String DEFAULT_SPEC = "{spec}";', template)
    if source == template and spec not in template:
        raise RuntimeError("could not insert case specification")
    java = GENERATED / "CoupledLossBenchmark.java"
    java.write_text(source, encoding="utf-8")
    subprocess.run([COMSOL / "comsolcompile.exe", java], cwd=ROOT,
                   check=True, timeout=120)
    stdout = OUTPUT / f"{case_name}_stdout.log"
    with stdout.open("w", encoding="utf-8") as log:
        subprocess.run([
            COMSOL / "comsolbatch.exe", "-inputfile",
            GENERATED / "CoupledLossBenchmark.class", "-outputfile",
            OUTPUT / f"{case_name}_batch.mph", "-batchlog",
            OUTPUT / f"{case_name}.log", "-batchlogout",
        ], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
            check=True, timeout=600)
    raw = stdout.read_text(encoding="utf-8", errors="replace")
    if f"OUTPUT={case_name}" not in raw or "POWER_SHAPE=" not in raw:
        raise RuntimeError(f"COMSOL result missing; inspect {stdout}")
    values = np.array([[float(a), float(b)] for a, b in
                       re.findall(r"^DATA,([^,]+),([^\r\n]+)", raw, re.M)])
    if len(values) < 3 or not np.all(np.diff(values[:, 0]) > 0):
        raise RuntimeError("invalid time-power trace")
    fields = spec.split(",")
    end = float(fields[6])
    if abs(values[-1, 0]-end) > 1e-12:
        raise RuntimeError("transient ended before requested time")
    pulse_end = 2*float(fields[4]) + .001
    if len(fields) == 13:
        pulse_end += float(fields[12])
    before = values[:, 0] <= pulse_end+1e-12
    after = values[:, 0] >= pulse_end-1e-12
    def integrate(rows: np.ndarray) -> float:
        return float(np.sum((rows[1:, 1]+rows[:-1, 1])
                            *np.diff(rows[:, 0])/2))

    result = {
        "case": case_name, "spec": spec, "samples": len(values),
        "comsol_version": "6.4.0.429",
        "power_expression": "IntSurface(2*pi*r*mf.Qh)",
        "pulse_heat_j": integrate(values[before]),
        "tail_heat_to_end_j": integrate(values[after]),
        "total_heat_to_end_j": integrate(values),
        "end_power_w": float(values[-1, 1]),
        "stdout": str(stdout.relative_to(ROOT)),
        "mph": str((OUTPUT / f"{case_name}.mph").relative_to(ROOT)),
    }
    (OUTPUT / f"{case_name}_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="+", choices=CASES)
    selected = parser.parse_args().cases
    for name in selected:
        print(json.dumps(run_case(name)), flush=True)
