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


def _ellipk(m, complement):
    m = np.asarray(m)
    # 1-m loses all bits when the two loops nearly coincide. Form the
    # complementary parameter from squared physical separation instead.
    return np.where(m > .9, ellipkm1(complement), ellipk(m))


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
    complement = d2/q**2
    m = np.where(complement < .01, 1-complement, 4*a*rsafe/q**2)
    k, e = _ellipk(m, complement), ellipe(m)
    common = MU0 / (2 * pi * q)
    br = common * dz / np.where(axis, 1.0, rsafe) * (
        -k + (a*a + rsafe*rsafe + dz*dz) / d2 * e)
    bz = common * (k + (a*a - rsafe*rsafe - dz*dz) / d2 * e)
    axis_d2 = a*a + dz*dz
    bz_axis = MU0 * a*a / (2 * axis_d2**1.5)
    # Maxwell's source-free on-axis expansion avoids cancellation of two
    # nearly equal elliptic-integral terms in Br as r -> 0.
    near_axis = rsafe < 1e-4*a
    br_axis = 3*MU0*a*a*dz*rsafe/(4*axis_d2**2.5)
    bz_axis += 3*MU0*a*a*rsafe**2*(a*a-4*dz*dz)/(8*axis_d2**3.5)
    return np.where(near_axis, br_axis, br), np.where(near_axis, bz_axis, bz)


def loop_mutual(source_r, target_r, dz):
    """Mutual inductance of two circular filaments in henries.

    Coincident filaments are singular; use finite-turn integration for self M.
    """
    a, b, dz = np.broadcast_arrays(np.asarray(source_r, float),
                                    np.asarray(target_r, float),
                                    np.asarray(dz, float))
    if np.any(a <= 0) or np.any(b <= 0):
        raise ValueError("loop radii must be positive")
    q2 = (a+b)**2 + dz*dz
    complement = ((a-b)**2 + dz*dz)/q2
    if np.any(complement == 0):
        raise ValueError("coincident filaments have divergent mutual inductance")
    m = np.where(complement < .01, 1-complement, 4*a*b/q2)
    # A direct K-E subtraction loses significant digits when m is tiny.
    small = m < 1e-6
    safe_m = np.where(small, .5, m)
    numerator = ((2-safe_m)*_ellipk(safe_m, complement)
                 - 2*ellipe(safe_m))
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


@lru_cache(maxsize=2048)
def _radial_nodes(r_inner: float, r_outer: float, profile: str, order: int):
    x, w = _unit_gauss(order)
    width = r_outer-r_inner
    radii = r_inner + width*x
    if profile == "dc_equilibrium":
        w = width*w/(radii*np.log1p(width/r_inner))
    return radii, w


def _axial_average(source_r, target_r, centre_dz, source_h, target_h,
                   order: int):
    """Exact axial convolution weight, with one-dimensional Gauss quadrature.

    For two uniform axial intervals the difference-coordinate density is a
    trapezoid. Its three linear pieces replace a two-dimensional z integral.
    Only the circular-loop kernel is numerically integrated.
    """
    half_sum = (source_h+target_h)/2
    half_difference = abs(source_h-target_h)/2
    breaks = (-half_sum, -half_difference, half_difference, half_sum)
    x, w = _unit_gauss(order)
    result = np.zeros(np.broadcast_shapes(np.shape(source_r),
                                          np.shape(target_r)))
    for lower, upper in zip(breaks[:-1], breaks[1:]):
        if upper == lower:
            continue
        difference = lower+(upper-lower)*x
        overlap = np.minimum(source_h, target_h)
        overlap = np.minimum(overlap, half_sum-np.abs(difference))
        weights = (upper-lower)*w*overlap/(source_h*target_h)
        kernel = loop_mutual(np.asarray(source_r)[..., None],
                             np.asarray(target_r)[..., None],
                             centre_dz+difference)
        result += np.sum(kernel*weights, axis=-1)
    return result


@lru_cache(maxsize=8192)
def _pair_mutual_shape(shape_a: tuple, shape_b: tuple,
                       centre_dz: float, order: int):
    ar, aw = _radial_nodes(shape_a[0], shape_a[1], shape_a[3], order)
    br, bw = _radial_nodes(shape_b[0], shape_b[1], shape_b[3], order)
    axial = _axial_average(ar[:, None], br[None, :],
                            centre_dz, shape_a[2], shape_b[2], order)
    return float(np.sum(aw[:, None]*bw[None, :]*axial))


def _pair_mutual(a: RectTurn, b: RectTurn, order: int):
    shape_a = (a.r_inner, a.r_outer, a.height, a.source_profile)
    shape_b = (b.r_inner, b.r_outer, b.height, b.source_profile)
    if shape_a > shape_b:
        shape_a, shape_b = shape_b, shape_a
    return _pair_mutual_shape(shape_a, shape_b,
                              abs(a.z_center-b.z_center), order)


@lru_cache(maxsize=2048)
def _self_mutual_shape(r_inner: float, r_outer: float, height: float,
                       profile: str, order: int):
    """Duffy quadrature of the integrable logarithmic self singularity.

    Translation invariance in z collapses one integral. Radial pair symmetry
    leaves the triangle r2>=r1. Two Duffy triangles resolve (dr,dz)=(0,0).
    """
    turn = RectTurn(r_inner, r_outer, 0.0, height, "cache",
                    source_profile=profile)
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


def _self_mutual(turn: RectTurn, order: int):
    return _self_mutual_shape(turn.r_inner, turn.r_outer, turn.height,
                              turn.source_profile, order)


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
