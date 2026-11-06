"""Exact M26-A fixed-topology 20 T field feasibility, no FEM."""
from supermagx.geometry import Magnet, RectTurn
from supermagx.design import axis_field_envelope, zone_axis_coefficients

radii=((132.,136.3),(150.6,154.9),(171.2,175.5),(187.8,192.1))
turns=[]
for mod,rs in enumerate(radii):
    for r in rs:
        for j in range(60):
            z=-244.85+8.3*j
            turns.append(RectTurn((r-2)*1e-3,(r+2)*1e-3,
                                  (z-4)*1e-3,(z+4)*1e-3,
                                  f"PS{min(mod+1,3)}",source_profile="uniform"))
magnet=Magnet(tuple(turns))
currents={"PS1":14600.785369366342,"PS2":17402.973451509286,
          "PS3":23740.635279040147}
for z in (0.,.0025):
    k=zone_axis_coefficients(magnet,z)
    print(f"Bz(0,{z*1000:g} mm): {sum(k[c]*I for c,I in currents.items()):.9f} T")
    if z==0:
        print("Field coefficients [T/kA]:",{c:1000*x for c,x in k.items()})
e=axis_field_envelope(magnet,20.,600e6,supply_current_ceiling_a=25000.,
                       expected_supply_count=3)
print("Current caps [kA]:",{z:x/1000 for z,x in e.current_ceilings_a.items()})
print(f"20 T fixed-topology B supremum [T]: {e.field_supremum_t:.9f}")
print(f"Required fraction of caps: {e.required_fraction_of_ceilings:.6f}")
print("WARNING: exact for modeled RECTANGLES, not physical R0.5 rounded corners.")
