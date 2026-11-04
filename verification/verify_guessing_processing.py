#!/usr/bin/env python3
"""Finite checks for BB84 score-to-envelope-rank transfer.

The general proof is analytic. These checks test the reduction's finite scalar
inequalities, averaging, and transfer to an independently solved injection LP.
They do not certify the physical noise or controller assumptions.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def branch_bound(k: int, s: int, c: int) -> float:
    return min(1.0, 2.0 ** ((s - k) / 6), 2.0 ** (s + c - k))


def h(k: int, r: float) -> float:
    if r <= 5 * k / 6:
        return 2.0 ** (r - k)
    if r <= 2 * k:
        return 2.0 ** ((r - 2 * k) / 7)
    return 1.0


def concave_majorant(k: int, r: float) -> float:
    b = 2.0 ** (-k)
    H = 2.0 ** (-k / 6)
    t = 5 * k / 6
    r = max(0.0, r)
    if r >= 2 * k:
        return 1.0
    if k >= 8:
        return b + (1 - b) * r / (2 * k)
    if r <= t:
        return b + (H - b) * r / t
    return H + (1 - H) * (r - t) / (2 * k - t)


def polygon_floor(k: int, p: float) -> float:
    b = 2.0 ** (-k)
    H = 2.0 ** (-k / 6)
    t = 5 * k / 6
    if p <= b:
        return 0.0
    if k >= 8:
        return 2 * k * (p - b) / (1 - b)
    if p <= H:
        return t * (p - b) / (H - b)
    return t + (2 * k - t) * (p - H) / (1 - H)


def affine_floor(k: int, p: float) -> float:
    b = 2.0 ** (-k)
    H = 2.0 ** (-k / 6)
    t = 5 * k / 6
    slope = max((H - b) / t, (1 - b) / (2 * k))
    return max(0.0, p - b) / slope


def injection_formula(D: int, M: float, a: float, q: float) -> float:
    if D == 0:
        return 0.0
    return max(0.0, max(a / q + (D - ell) * (1 - q) * a / q
                        - M * q ** (ell - 1) for ell in range(1, D + 1)))


def injection_lp(D: int, M: float, a: float, q: float) -> float:
    if D == 0:
        return 0.0
    matrix = np.zeros((D, D))
    rhs = np.zeros(D)
    for j in range(1, D + 1):
        rhs[j - 1] = a - M * q**j
        for i in range(1, j + 1):
            matrix[j - 1, i - 1] = q ** (j - i + 1)
    result = linprog(np.ones(D), A_ub=-matrix, b_ub=-rhs,
                     bounds=(0, None), method="highs")
    if not result.success:
        raise AssertionError(result.message)
    return float(result.fun)


def upper_concave_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Least concave polygon majorizing a finite x-sorted set."""
    hull: list[tuple[float, float]] = []
    for point in points:
        while len(hull) > 1:
            left = (hull[-1][1] - hull[-2][1]) * (point[0] - hull[-1][0])
            right = (point[1] - hull[-1][1]) * (hull[-1][0] - hull[-2][0])
            if left > right + 1e-14:
                break
            hull.pop()
        hull.append(point)
    return hull


def inverse_hull(hull: list[tuple[float, float]], p: float) -> float:
    if p <= hull[0][1]:
        return hull[0][0]
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        if p <= y1:
            return x0 + (x1 - x0) * (p - y0) / (y1 - y0)
    return hull[-1][0]


def inherited_recovery_bound(k: int, r: int, s: int | None = None,
                             include_bcw: bool = True) -> float:
    """Primary-literature reductions supplied by the recovery review.

    This numerical comparison checks the consequences of the stated formulas.
    Attribution and the analytic validity of those imported formulas require
    the separate primary-source review; this routine is not a source audit.
    """
    d = 2**k
    f = min(1.0, 2.0 ** (r - 2 * k))
    dimension_bound = 2.0 ** (r - k)
    if s is not None:
        f = min(f, 2.0 ** (s - k))
        dimension_bound = 2.0 ** (r - s - k)
    u = min((1 + f) / 2,
            min(sum(math.comb(k, j) for j in range(ell + 1)) * f
                + (2.0 ** (-ell - 1) if ell < k else 0.0)
                for ell in range(k + 1)))
    if include_bcw:
        u = min(u, (d * f + 1 / d) / 2)
    inverse_g = (1 + math.sqrt(max(0.0, (d - 1) * (d * u - 1)))) / d
    return min(1.0, dimension_bound, inverse_g)


def pow2_fraction(exponent: int) -> Fraction:
    return Fraction(2**exponent, 1) if exponent >= 0 else Fraction(1, 2**(-exponent))


def exact_recovery_chord_certificates() -> list[dict]:
    """Exact small-k boundary certificates for the inherited reduction.

    For each positive integer rank, the inherited bounds imply p^2 <= a
    listed rational candidate. A candidate <= (r/(2k))^2 certifies p<=r/(2k),
    which is stronger than the desired baseline-to-one endpoint chord.
    Rank zero is treated separately by its exact random-guess bound 2^-k.
    """
    certificates = []
    for k in (9, 10):
        witnesses = []
        for r in range(1, 2 * k + 1):
            f = pow2_fraction(r - 2 * k)
            candidates = {
                "dimension_squared": pow2_fraction(2 * (r - k)),
                "pair_MUB_PGM": (1 + f) / 2,
                "BCW_R3": (2**k * f + Fraction(1, 2**k)) / 2,
            }
            for ell in range(k + 1):
                tail = pow2_fraction(-ell - 1) if ell < k else Fraction(0)
                candidates[f"DFW_ell_{ell}"] = sum(math.comb(k, j)
                                                   for j in range(ell + 1)) * f + tail
            name, value = min(candidates.items(), key=lambda item: item[1])
            target = Fraction(r, 2 * k) ** 2
            assert value <= target
            witnesses.append({"r": r, "bound": name,
                              "value": str(value), "target": str(target)})
        certificates.append({"k": k, "arithmetic": "exact Fraction",
                             "all_positive_integer_ranks_certified": True,
                             "rank_zero_score": str(pow2_fraction(-k)),
                             "witnesses": witnesses})
    return certificates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        help='Save this run without replacing the reference output.')
    args = parser.parse_args()
    rng = np.random.default_rng(20261104)
    finite_count = 0
    mixture_count = 0
    worst_branch_margin = 0.0
    worst_mixture_margin = 0.0
    for k in range(1, 41):
        b = 2.0 ** (-k)
        H = 2.0 ** (-k / 6)
        t = 5 * k / 6
        hinge_minus_chord = H - (b + (1 - b) * t / (2 * k))
        assert (hinge_minus_chord > 0) == (k <= 7)
        for r in np.linspace(0, 2 * k, 2001):
            assert h(k, r) <= concave_majorant(k, r) + 2e-14
        for n in sorted({k, k + 2, 2 * k}):
            entries = []
            for s in range(n + 1):
                for c in range(n - s + 1):
                    r = 2 * s + c
                    p = branch_bound(k, s, c)
                    assert p <= h(k, r) + 2e-14
                    assert p <= concave_majorant(k, r) + 2e-14
                    assert polygon_floor(k, p) <= r + 2e-12
                    assert affine_floor(k, p) <= r + 2e-12
                    worst_branch_margin = min(worst_branch_margin,
                                              r - polygon_floor(k, p))
                    entries.append((r, p))
                    finite_count += 1
            entries = np.asarray(entries)
            for _ in range(50):
                size = int(rng.integers(2, 15))
                chosen = entries[rng.integers(0, len(entries), size)]
                weights = rng.dirichlet(np.ones(size))
                R = float(weights @ chosen[:, 0])
                # Scores may be below the upper branch bounds; all are
                # unconditional weighted averages with fixed probabilities.
                p = float(weights @ (chosen[:, 1] * rng.random(size)))
                margin = R - polygon_floor(k, p)
                assert margin >= -2e-12
                assert affine_floor(k, p) <= R + 2e-12
                worst_mixture_margin = min(worst_mixture_margin, margin)
                mixture_count += 1

    # Endpoint mixtures expose the direction of the required averaging step.
    # A convex exponential cannot be evaluated at expected rank as an upper
    # bound. The concave majorant is tight on this scalar relaxation for k>=8.
    endpoint_mixtures = []
    for k in (8, 12, 40):
        b = 2.0 ** (-k)
        for weight in (0.01, 0.2, 0.4, 0.9):
            R = weight * 2 * k
            p = b * (1 - weight) + weight
            assert abs(polygon_floor(k, p) - R) < 2e-12
            endpoint_mixtures.append({"k": k, "full_rank_weight": weight,
                                      "expected_rank": R, "score": p})

    lp_count = 0
    lp_worst_absolute_error = 0.0
    for D in (1, 2, 3, 5, 10, 20):
        for q in (0.2, 0.5, 0.9, 0.99):
            for a in (0.2, 1.0, 3.5):
                for reserve in (1.0, 2.0, 10.0):
                    M = a * reserve
                    formula = injection_formula(D, M, a, q)
                    lp = injection_lp(D, M, a, q)
                    error = abs(formula - lp)
                    assert error < 1e-7 * max(1.0, formula)
                    lp_worst_absolute_error = max(lp_worst_absolute_error, error)
                    lp_count += 1

    k, p, n, gamma, D = 12, 0.4, 120, 0.01, 2000
    q = 1 - gamma
    a = polygon_floor(k, p)
    gap = max(0.0, D * -math.log(q) - math.log(2 * n / a))
    example = {
        "k": k, "p": p, "n": n, "gamma": gamma, "D": D,
        "random_guess_baseline": 2.0 ** (-k),
        "lattice_only_baseline": 2.0 ** (-k / 6),
        "lattice_only_floor": max(0.0, 2 * k * (p - 2.0 ** (-k / 6))
                                    / (1 - 2.0 ** (-k / 6))),
        "combined_rank_floor": a,
        "logarithmic_charge_bound": a * gap,
        "integer_gate_logarithmic_real_bound": gap / math.log1p(1 / a),
        "integer_gate_logarithmic_ceiling": math.ceil(gap / math.log1p(1 / a)),
        "discrete_injection_real_bound": injection_formula(D, 2 * n, a, q),
        "discrete_injection_ceiling": math.ceil(injection_formula(D, 2 * n, a, q)),
        "interval_example_D": 700,
        "interval_discrete_real_bound": injection_formula(700, 2 * n, a, q),
        "width_feasibility_check_a_le_2nq": a <= q * 2 * n,
    }
    recovery_comparison = []
    for k_cmp in range(1, 41):
        direct_points = [
            (r, max(branch_bound(k_cmp, s, r - 2 * s)
                    for s in range(r // 2 + 1)))
            for r in range(2 * k_cmp + 1)
        ]
        recovery_points = [(r, inherited_recovery_bound(k_cmp, r))
                           for r in range(2 * k_cmp + 1)]
        algebra_aware_points = [
            (r, max(inherited_recovery_bound(k_cmp, r, s)
                    for s in range(r // 2 + 1)))
            for r in range(2 * k_cmp + 1)
        ]
        without_bcw_points = [(r, inherited_recovery_bound(k_cmp, r, include_bcw=False))
                              for r in range(2 * k_cmp + 1)]
        assert all(abs(a[1] - b[1]) < 1e-13
                   for a, b in zip(recovery_points, without_bcw_points))
        assert all(abs(a[1] - b[1]) < 1e-13
                   for a, b in zip(algebra_aware_points, recovery_points))
        direct_hull = upper_concave_hull(direct_points)
        recovery_hull = upper_concave_hull(recovery_points)
        algebra_aware_hull = upper_concave_hull(algebra_aware_points)
        combined_hull = upper_concave_hull([
            (r, min(direct_points[r][1], recovery_points[r][1]))
            for r in range(2 * k_cmp + 1)
        ])
        p_cmp = 0.4
        recovery_comparison.append({
            "k": k_cmp,
            "direct_integer_rank_hull": direct_hull,
            "inherited_recovery_integer_rank_hull": recovery_hull,
            "algebra_aware_inherited_recovery_integer_rank_hull": algebra_aware_hull,
            "combined_integer_rank_hull": combined_hull,
            "example_score": p_cmp,
            "direct_rank_floor_at_example": inverse_hull(direct_hull, p_cmp),
            "inherited_rank_floor_at_example": inverse_hull(recovery_hull, p_cmp),
            "combined_rank_floor_at_example": inverse_hull(combined_hull, p_cmp),
            "inherited_hull_is_endpoint_chord": len(recovery_hull) == 2,
            "algebra_aware_family_equals_rank_only_family": True,
            "BCW_R3_addition_leaves_maximized_curve_unchanged": True,
        })
    output = {
        "status": "passed",
        "seed": 20261104,
        "scope": "Finite scalar checks; no physical-model certification or universal proof by computation.",
        "k_range": [1, 40],
        "discrete_algebra_cases": finite_count,
        "randomized_unconditional_mixtures": mixture_count,
        "sampled_continuous_majorant_cases": 40 * 2001,
        "worst_branch_floor_margin": worst_branch_margin,
        "worst_mixture_floor_margin": worst_mixture_margin,
        "endpoint_mixtures": endpoint_mixtures,
        "independent_injection_LP_cases": lp_count,
        "LP_relative_tolerance": 1e-7,
        "LP_checks_passed": True,
        "low_score_example": example,
        "recovery_comparison": recovery_comparison,
        "recovery_comparison_limit": "Only the supplied primary-source formula family was compared; finite k=1..40 does not prove the all-k simplification or absolute optimality.",
        "recovery_families_compared": {
            "rank_only": "f=min(1,2^(r-2k)); dimension score<=min(1,2^(r-k)); refined discrimination inverse applied to min(pair MUB PGM (1+f)/2, BCW R3 (d*f+1/d)/2, DFW binomial family)",
            "algebra_aware": "f=min(1,2^(r-2k),2^(s-k)); dimension score<=min(1,2^(r-s-k)); same refined discrimination and recovery family including BCW R3; maximize over integer s with c=r-2s>=0",
            "algebra_aware_equality_proof": "For r<=k choose s=0,c=r. For k<=r<=2k choose s=r-k,c=2k-r. Both choices attain the rank-only family. They fit every width n>=k when r<=2n. Consequently adding the existing coherent-block bottleneck does not alter the rank-only upper curve in this width regime.",
            "finite_width_limit": "The stored unconstrained-rank comparison allows every s with c>=0 and fits n>=k. If n<k, enforce s+c<=n; this can strengthen either family.",
        },
        "exact_recovery_chord_certificates": exact_recovery_chord_certificates(),
    }
    report = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding='utf-8')
    print(report, end='')


if __name__ == "__main__":
    main()
