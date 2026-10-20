"""Finite-cross-section, axisymmetric free-space magnetostatics.

The loops are integrated over each conductor's actual r-z rectangle. The
prescribed transport distribution is selected per turn; it is *not* an FEM
solution for current redistribution, ferromagnets, or transient eddy currents.
Numerical quadrature is refined and nonconvergence is reported explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import pi
from typing import Hashable, Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ellipk, ellipkm1, ellipe

from supermagx.geometry import Magnet, RectTurn

MU0 = 4e-7 * pi


class ConvergenceError(ArithmeticError):
    """Requested quadrature tolerance was not reached."""


@dataclass(frozen=True)
class FieldResult:
    br: float
    bz: float
    order: int
    relative_change: float


def _ellipk(m):
    m = np.asarray(m)
    return np.where(m > .9, ellipkm1(1 - m), ellipk(m))


def loop_field(source_r, dz, target_r):
    """(Br, Bz) in tesla from a 1 A circular filament, including the axis."""
    a, dz, r = np.broadcast_arrays(np.asarray(source_r, float),
                                    np.asarray(dz, float),
                                    np.asarray(target_r, float))
    if np.any(a <= 0) or np.any(r < 0):
        raise ValueError("radii must be nonnegative, and source radius positive")
    axis = r == 0
    rsafe = np.where(axis, 0.0, r)
    d2 = (a - rsafe)**2 + dz**2
    if np.any(~axis & (d2 == 0)):
        raise ValueError("filament field is singular at its source")
    q = np.sqrt((a + rsafe)**2 + dz**2)
    m = 4 * a * rsafe / q**2
    k, e = _ellipk(m), ellipe(m)
    common = MU0 / (2 * pi * q)
    br = common * dz / np.where(axis, 1.0, rsafe) * (
        -k + (a*a + rsafe*rsafe + dz*dz) / d2 * e)
    bz = common * (k + (a*a - rsafe*rsafe - dz*dz) / d2 * e)
    bz_axis = MU0 * a*a / (2 * (a*a + dz*dz)**1.5)
    return np.where(axis, 0.0, br), np.where(axis, bz_axis, bz)


def loop_mutual(source_r, target_r, dz):
    """Mutual inductance of two circular filaments in henries.

    Coincident filaments are singular; use finite-turn integration for self M.
    """
    a, b, dz = np.broadcast_arrays(np.asarray(source_r, float),
                                    np.asarray(target_r, float),
                                    np.asarray(dz, float))
    if np.any(a <= 0) or np.any(b <= 0):
        raise ValueError("loop radii must be positive")
    m = 4*a*b / ((a+b)**2 + dz*dz)
    if np.any(m >= 1):
        raise ValueError("coincident filaments have divergent mutual inductance")
    # A direct K-E subtraction loses significant digits when m is tiny.
    small = m < 1e-6
    safe_m = np.where(small, .5, m)
    numerator = (2-safe_m)*_ellipk(safe_m) - 2*ellipe(safe_m)
    usual = MU0*np.sqrt(a*b)*numerator/np.sqrt(safe_m)
    series = MU0*np.sqrt(a*b)*(pi/16*m**1.5 + 3*pi/64*m**2.5)
    return np.where(small, series, usual)


@lru_cache(maxsize=32)
def _unit_gauss(order: int):
    x, w = leggauss(order)
    return (x + 1)/2, w/2


@lru_cache(maxsize=4096)
def _turn_nodes(turn: RectTurn, order: int):
    x, w = _unit_gauss(order)
    r = turn.r_inner + turn.width*x[:, None]
    z = turn.z_lower + turn.height*x[None, :]
    rr = np.broadcast_to(r, (order, order)).ravel()
    zz = np.broadcast_to(z, (order, order)).ravel()
    if turn.source_profile == "uniform":
        weights = np.broadcast_to(w[:, None]*w[None, :],
                                  (order, order)).ravel()
    else:
        weights = (turn.width*turn.height * w[:, None]*w[None, :]
                   / np.broadcast_to(r, (order, order))
                   / turn.h_integral).ravel()
    return rr, zz, weights


def _turn_field(turn: RectTurn, r: float, z: float, order: int):
    sr, sz, weights = _turn_nodes(turn, order)
    br, bz = loop_field(sr, z-sz, r)
    return float(weights@br), float(weights@bz)


def _pair_mutual(a: RectTurn, b: RectTurn, order: int):
    # Shift target/source grids to avoid coincident quadrature nodes.
    ar, az, aw = _turn_nodes(a, order)
    br, bz, bw = _turn_nodes(b, order+1)
    kernel = loop_mutual(ar[:, None], br[None, :], az[:, None]-bz[None, :])
    return float(aw @ kernel @ bw)


def _self_mutual(turn: RectTurn, order: int):
    """Duffy quadrature of the integrable logarithmic self singularity.

    Translation invariance in z collapses one integral. Radial pair symmetry
    leaves the triangle r2>=r1. Two Duffy triangles resolve (dr,dz)=(0,0).
    """
    x, w = _unit_gauss(order)
    r1 = turn.r_inner + turn.width*x[:, None, None]
    umax = turn.r_outer-r1
    s = x[None, :, None]
    t = x[None, None, :]
    factor = (turn.width*w[:, None, None] * w[None, :, None]
              * w[None, None, :] * umax*turn.height*s)
    out = 0.0
    for u, v in ((umax*s, turn.height*s*t),
                 (umax*s*t, turn.height*s)):
        r2 = r1 + u
        kernel = loop_mutual(r1, r2, v)
        profile = (1.0 if turn.source_profile == "uniform" else 1/(r1*r2))
        out += float(np.sum(factor*(turn.height-v)*kernel*profile))
    norm = (turn.width*turn.height if turn.source_profile == "uniform"
            else turn.h_integral)
    return 4*out/norm**2


class AxisymmetricMagnetostatics:
    """Independent source zones; no fixed layer count or winding topology."""

    def __init__(self, magnet: Magnet, *, rtol: float = 1e-6,
                 max_order: int = 64):
        if not (0 < rtol < 1) or max_order < 8:
            raise ValueError("require 0 < rtol < 1 and max_order >= 8")
        self.magnet = magnet
        self.rtol = rtol
        self.max_order = max_order

    def _zone_field(self, r: float, z: float, order: int):
        br = {zone: 0.0 for zone in self.magnet.zones}
        bz = br.copy()
        for turn in self.magnet.turns:
            b_r, b_z = _turn_field(turn, r, z, order)
            br[turn.zone] += turn.direction*b_r
            bz[turn.zone] += turn.direction*b_z
        return br, bz

    def field(self, r: float, z: float,
              currents: Mapping[Hashable, float]) -> FieldResult:
        """Field in free space; raises if quadrature has not converged.

        Points inside metal need a dedicated singular-cell integration and are
        deliberately rejected by this precision-controlled public method.
        """
        if not np.isfinite(r) or not np.isfinite(z) or r < 0:
            raise ValueError("observation coordinates must be finite, r >= 0")
        unknown = set(currents)-set(self.magnet.zones)
        if unknown:
            raise ValueError(f"unknown zone: {unknown}")
        if any(not np.isfinite(v) for v in currents.values()):
            raise ValueError("currents must be finite")
        if any(t.r_inner <= r <= t.r_outer and
               t.z_lower <= z <= t.z_upper for t in self.magnet.turns):
            raise ValueError("field point is inside a conductor")
        if not any(currents.values()):
            return FieldResult(0.0, 0.0, 0, 0.0)
        previous = None
        for order in (4, 8, 16, 32, 64, 128, 256):
            if order > self.max_order:
                break
            br_z, bz_z = self._zone_field(r, z, order)
            value = np.array((sum(br_z[k]*currents.get(k, 0.0)
                                  for k in self.magnet.zones),
                              sum(bz_z[k]*currents.get(k, 0.0)
                                  for k in self.magnet.zones)))
            if previous is not None:
                change = float(np.linalg.norm(value-previous)
                               / max(np.linalg.norm(value), 1e-15))
                if change <= self.rtol:
                    return FieldResult(float(value[0]), float(value[1]),
                                       order, change)
            previous = value
        raise ConvergenceError(f"field quadrature did not meet rtol={self.rtol}")

    def screening_coefficients(self, r: float, z: float, *, order: int = 6):
        """Turn-centre B per ampere by zone; for loss/force screening only.

        This evaluates within the conductor and therefore has no certified
        quadrature error; it is intentionally separate from field().
        """
        br, bz = self._zone_field(r, z, order)
        return np.array([br[k] for k in self.magnet.zones]), np.array(
            [bz[k] for k in self.magnet.zones])

    def mutual_matrix(self, *, rtol: float = 1e-4,
                      max_order: int = 24) -> np.ndarray:
        """Zone mutual matrix, H, by adaptive finite-turn quadrature.

        This O(N_turns^2) reference calculation can be slow for large coils.
        Every turn pair is checked; an unresolved pair raises ConvergenceError.
        """
        if not (0 < rtol < 1) or max_order < 8:
            raise ValueError("invalid mutual quadrature settings")
        turns = self.magnet.turns
        zones = self.magnet.zones
        zi = {zone: i for i, zone in enumerate(zones)}
        matrix = np.zeros((len(zones), len(zones)))
        orders = (2, 4, 6, 8, 12, 16, 24, 32, 48)
        for i, a in enumerate(turns):
            for j in range(i, len(turns)):
                b = turns[j]
                last = None
                converged = False
                for n in orders:
                    if n > max_order:
                        break
                    value = (_self_mutual(a, n) if i == j
                             else _pair_mutual(a, b, n))
                    if last is not None and abs(value-last) <= rtol*abs(value):
                        converged = True
                        break
                    last = value
                if not converged:
                    raise ConvergenceError(
                        f"mutual quadrature did not converge for turns {i},{j}")
                value *= a.direction*b.direction
                ia, ib = zi[a.zone], zi[b.zone]
                matrix[ia, ib] += value
                if i != j:
                    matrix[ib, ia] += value
        return matrix
