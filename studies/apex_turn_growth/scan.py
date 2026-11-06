"""M31: increasing turns using existing exact finite-rectangle on-axis kernels.

Run: python studies/apex_turn_growth/scan.py
Pure magnetostatics: NO stress, terminal-voltage, or thermal certification.
"""
from __future__ import annotations
import csv
from pathlib import Path
from math import pi
from supermagx.geometry import RectTurn
from supermagx.design.turn_axis import exact_turn_axis_coefficient

INITIAL_RADIUS_M = .130
RADIAL_PITCH_M = .0043
AXIAL_PITCH_M = .0083
METAL_RADIAL_M = .004
METAL_AXIAL_M = .008
MODULE_GAPS_M = (.010, .012, .008)
EXTERNAL_ZYLON_M = .016
TARGET_FIELD_T = 20.
CORNER_RADIUS_M = .0005
METAL_AREA_ROUNDED_M2 = METAL_RADIAL_M*METAL_AXIAL_M - (4-pi)*CORNER_RADIUS_M**2

def study_case(radial_layers: tuple[int,int,int,int], turns_per_layer: int) -> dict:
    """Four physical modules; PS3 drives modules 3 and 4 in series.

    All THREE power circuits are assigned equal currents solely to determine
    the minimum possible peak current (hence max J) for this specific fixed
    all-positive-coefficient geometry with identical conductor sections.

    Fields are for EXACT rectangular metal turns (uniform J per rectangle).
    Rounded-corner area is used only as a current-density SCREEN and does
    not claim to reproduce the actual rounded-wire magnetic source.
    """
    if len(radial_layers)!=4 or any(type(n) is not int or n<1 for n in radial_layers):
        raise ValueError("radial_layers must have four positive integers")
    if type(turns_per_layer) is not int or turns_per_layer<1:
        raise ValueError("turns_per_layer must be a positive integer")
    ks={"PS1":0., "PS2":0., "PS3":0.}
    r=INITIAL_RADIUS_M
    for module,n_layers in enumerate(radial_layers):
        zone=f"PS{min(module+1,3)}"
        for k in range(n_layers):
            rinner=r+k*RADIAL_PITCH_M
            for j in range(turns_per_layer):
                z=(j-(turns_per_layer-1)/2)*AXIAL_PITCH_M
                t=RectTurn(rinner,rinner+METAL_RADIAL_M,
                           z-METAL_AXIAL_M/2,z+METAL_AXIAL_M/2,
                           zone,source_profile="uniform")
                ks[zone]+=exact_turn_axis_coefficient(t)
        r+=n_layers*RADIAL_PITCH_M
        if module<3:
            r+=MODULE_GAPS_M[module]
    common_i=TARGET_FIELD_T/sum(ks.values())
    return {
        "layers": "+".join(str(n) for n in radial_layers),
        "radial_layers_total": sum(radial_layers),
        "turns_per_layer": turns_per_layer,
        "turns_total": sum(radial_layers)*turns_per_layer,
        "axial_envelope_mm": 1000*((turns_per_layer-1)*AXIAL_PITCH_M+METAL_AXIAL_M),
        "outer_diameter_mm": 2000*(r+EXTERNAL_ZYLON_M),
        "required_equal_supply_current_kA": common_i/1000,
        "current_density_rectangular_A_mm2": common_i/(METAL_RADIAL_M*METAL_AXIAL_M)/1e6,
        "rounded_area_screen_A_mm2": common_i/METAL_AREA_ROUNDED_M2/1e6,
        "under_600_A_mm2_area_screen": common_i<=600e6*METAL_AREA_ROUNDED_M2,
        "under_25kA": common_i<25000,
        "PS1_field_coefficient_T_per_kA": ks["PS1"]*1000,
        "PS2_field_coefficient_T_per_kA": ks["PS2"]*1000,
        "PS3_field_coefficient_T_per_kA": ks["PS3"]*1000,
    }

def cases() -> list[dict]:
    candidates=((2,2,2,2),(3,3,2,2),(3,3,3,3),(3,3,4,4),
                (4,4,4,4),(4,4,5,5),(5,5,5,5),(6,6,6,6))
    return [study_case(layers,axial) for layers in candidates
            for axial in (60,64,70,80,100,120,160)]

if __name__=="__main__":
    results=cases()
    path=Path(__file__).with_name("results.csv")
    with path.open("w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    print(f"{len(results)} deterministic magnetostatic candidates -> {path}")
