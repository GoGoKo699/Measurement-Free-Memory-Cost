#!/usr/bin/env python3
"""Finite checks for achievability.md; no simulation is used as a theorem.

Only the Python standard library is required. Exact rational tree enumeration,
exhaustive five-qubit Pauli histories, explicit register allocation, and numeric
finite-parameter inequalities exercise different parts of the proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path


def pauli(text: str) -> tuple[int, int]:
    x = sum((c in "XY") << i for i, c in enumerate(text))
    z = sum((c in "ZY") << i for i, c in enumerate(text))
    return x, z


def mul(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    """Phaseless Pauli multiplication."""
    return a[0] ^ b[0], a[1] ^ b[1]


def symp(a: tuple[int, int], b: tuple[int, int]) -> int:
    return ((a[0] & b[1]).bit_count() + (a[1] & b[0]).bit_count()) & 1


def five_qubit_histories() -> dict:
    generators = [pauli(s) for s in ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")]
    stabilizer = {(0, 0)}
    for g in generators:
        stabilizer |= {mul(v, g) for v in list(stabilizer)}
    assert len(stabilizer) == 16
    assert all(symp(a, b) == 0 for a, b in product(generators, repeat=2))

    def syndrome(error: tuple[int, int]) -> int:
        return sum(symp(error, g) << j for j, g in enumerate(generators))

    def weight(error: tuple[int, int]) -> int:
        return (error[0] | error[1]).bit_count()

    logical_x, logical_z = pauli("XXXXX"), pauli("ZZZZZ")
    assert symp(logical_x, logical_z) == 1
    assert syndrome(logical_x) == syndrome(logical_z) == 0
    distance = min(
        weight((x, z))
        for x, z in product(range(32), repeat=2)
        if syndrome((x, z)) == 0 and (x, z) not in stabilizer
    )
    assert distance == 3

    correctable = [(0, 0)]
    for i in range(5):
        correctable.extend(((1 << i, 0), (0, 1 << i), (1 << i, 1 << i)))
    table = {syndrome(e): e for e in correctable}
    assert len(table) == 16
    histories = accumulated_over_one = final_syndrome_only_failures = 0
    for history in product(correctable, repeat=3):
        actual = estimate = (0, 0)
        previous_syndrome = 0
        for error in history:
            actual = mul(error, actual)
            current = syndrome(actual)
            increment = current ^ previous_syndrome
            estimate = mul(table[increment], estimate)
            previous_syndrome = current
        assert mul(actual, estimate) in stabilizer
        histories += 1
        accumulated_over_one += weight(actual) > 1
        final_syndrome_only_failures += mul(actual, table[syndrome(actual)]) not in stabilizer

    # Independent check of root's complementary coherent-recovery ANFs.
    # This is not required for the zero-T theorem; the symplectic convention
    # fixes the decoded logical frame before any explicit encoder is chosen.
    for syn, error in table.items():
        a, b, c, d = [(syn >> j) & 1 for j in range(4)]
        fx = a ^ d ^ (a & c) ^ (b & c) ^ (b & d)
        fz = b ^ c ^ (a & b) ^ (a & d) ^ (c & d)
        assert fx == symp(error, logical_z)
        assert fz == symp(error, logical_x)
    return {
        "distance": distance,
        "correctable_increment_count": len(correctable),
        "three_increment_histories_checked": histories,
        "histories_with_accumulated_weight_above_one": accumulated_over_one,
        "failures_if_only_final_syndrome_used": final_syndrome_only_failures,
        "deferred_history_decoder_failures": 0,
        "complementary_coherent_recovery_anf_rows_checked": len(table),
    }


def majority(bits: list[int]) -> int:
    return int(sum(bits) > len(bits) // 2)


def f3(u: Fraction) -> Fraction:
    return 3 * u * u - 2 * u * u * u


def exact_tree_error(height: int, initial_noise: bool) -> Fraction:
    """Enumerate all edge flips and perform recursive majority on true leaves."""
    edges = (3 ** (height + 1) - 3) // 2
    count = edges + int(initial_noise)
    wrong_by_weight = [0] * (count + 1)
    for mask in range(1 << count):
        offset = int(initial_noise)
        root = (mask & 1) if initial_noise else 0

        def decode(bit: int, depth: int) -> int:
            nonlocal offset
            if depth == 0:
                return bit
            decoded_children = []
            for _ in range(3):
                edge_flip = (mask >> offset) & 1
                offset += 1
                decoded_children.append(decode(bit ^ edge_flip, depth - 1))
            return majority(decoded_children)

        output = decode(root, height)
        assert offset == count
        if output:
            wrong_by_weight[mask.bit_count()] += 1
    e = Fraction(1, 32)
    return sum(
        (multiplicity * e**w * (1 - e) ** (count - w)
         for w, multiplicity in enumerate(wrong_by_weight)),
        Fraction(0),
    )


def broadcast_checks() -> dict:
    e, q = Fraction(1, 32), Fraction(15, 16)
    invariant = f3(e + q / 4)
    assert invariant == Fraction(22831, 131072) < Fraction(1, 4)
    eta = Fraction(0)
    exact = []
    for height in range(3):
        if height:
            eta = f3(e + q * eta)
        theta = e + q * eta
        assert exact_tree_error(height, False) == eta
        assert exact_tree_error(height, True) == theta
        assert theta <= Fraction(1, 2) - q / 4
        exact.append({"height": height, "root_error": str(eta), "true_syndrome_error": str(theta)})
    # A deep numerical recurrence is only a check of implementation stability.
    eta_float = 0.0
    for _ in range(10000):
        u = float(e) + float(q) * eta_float
        eta_float = 3 * u * u - 2 * u * u * u
        assert eta_float <= 0.25
    return {
        "endpoint_gamma": "1/16",
        "invariant_image": str(invariant),
        "enumerated_heights": exact,
        "noise_patterns_checked": sum(3 * (1 << ((3 ** (h + 1) - 3) // 2)) for h in range(3)),
        "root_error_after_10000_generations": eta_float,
        "true_syndrome_error_after_10000_generations": float(e) + float(q) * eta_float,
    }


def allocation_checks() -> dict:
    b, s = 3, 2
    cases = []
    for blocks in range(2, 10):
        capacity = b * s * (2 * 3 ** (blocks - 2) - 1)
        pool = set(range(capacity))
        live: set[int] = set()
        peak = preparations = 0
        for j in range(2, blocks + 1):
            next_size = 3 * len(live) + b * s
            targets = set(sorted(pool - live)[:next_size])
            assert len(targets) == next_size and not (targets & live)
            peak = max(peak, len(live | targets))
            preparations += next_size
            live = targets
            assert len(live) == b * s * (3 ** (j - 1) - 1) // 2
        assert peak == capacity
        assert preparations == b * s * (3**blocks - 2 * blocks - 1) // 4
        cases.append({"blocks": blocks, "archive_peak": peak, "target_preparations": preparations})
    return {"replicas": b, "syndrome_bits": s, "explicit_register_allocations": cases}


def h2(x: float) -> float:
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


R0 = 1 - h2(1 / 8) - math.log2(3) / 8


def logsumexp(values: list[float]) -> float:
    top = max(values)
    return top + math.log(sum(math.exp(v - top) for v in values))


def log_binomial_tail(n: int, t: int, p: float) -> float:
    terms = [
        math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
        + j * math.log(p) + (n - j) * math.log1p(-p)
        for j in range(t + 1, n + 1)
    ]
    return logsumexp(terms)


def parameter_checks() -> dict:
    maximum_entropy_excess = -math.inf
    for n in range(1, 257):
        ball = sum(3**j * math.comb(n, j) for j in range(2 * (n // 16) + 1))
        excess = math.log2(ball) - n * (1 - R0)
        maximum_entropy_excess = max(maximum_entropy_excess, excess)
        assert excess <= 1e-10
    examples = []
    for gamma_fraction, k, rounds, epsilon in (
        (Fraction(1, 16), 1, 1, 0.05),
        (Fraction(1, 16), 1, 10, 0.01),
        (Fraction(1, 64), 10, 100, 0.001),
        (Fraction(1, 160), 100, 1000, 0.01),
        (Fraction(1, 1000), 1, 1, 0.9),
    ):
        gamma = float(gamma_fraction)
        tau = (1 / (16 * gamma_fraction)).numerator // (1 / (16 * gamma_fraction)).denominator
        blocks = (rounds + tau - 1) // tau
        effective_gamma = -math.expm1(tau * math.log1p(-gamma))
        q, p = 1 - effective_gamma, 3 * effective_gamma / 4
        alpha = 1 / 16
        kl = alpha * math.log(alpha / p) + (1 - alpha) * math.log((1 - alpha) / (1 - p))
        n = math.ceil(max((k + 1) / R0, math.log(2 * blocks / epsilon) / kl))
        s, t = n - k, n // 16
        assert s > 0 and effective_gamma <= 1 / 16 + 1e-15
        assert blocks <= 32 * gamma * rounds + 1 + 1e-12
        actual_lengths = [min(j * tau, rounds) - min((j - 1) * tau, rounds) for j in range(1, blocks + 1)]
        assert sum(actual_lengths) == rounds and min(actual_lengths) >= 1
        assert max(actual_lengths) <= tau
        code_union_log_bound = math.log(4 / 3) + (k - R0 * n) * math.log(2)
        assert code_union_log_bound <= math.log(2 / 3) + 1e-12
        tail = log_binomial_tail(n, t, p)
        assert tail <= -n * kl + 1e-8
        data_union = math.exp(math.log(blocks) + tail)
        assert data_union <= epsilon / 2 + 1e-12
        if blocks > 1:
            count = s * (blocks - 1)
            replicas = math.ceil(8 / q**2 * math.log(2 * count / epsilon))
            if replicas % 2 == 0:
                replicas += 1
            archive_union = math.exp(math.log(count) - replicas * q**2 / 8)
            width = n + replicas * s * (2 * 3 ** (blocks - 2) - 1)
            assert archive_union <= epsilon / 2 + 1e-12
        else:
            replicas, archive_union, width = 0, 0.0, n
        examples.append({
            "gamma": str(gamma_fraction), "k": k, "physical_rounds": rounds,
            "epsilon": epsilon, "tau": tau, "blocks": blocks,
            "data_width": n, "correctable_weight": t, "root_replicas": replicas,
            "width_decimal": str(width), "log2_width": math.log2(width),
            "exact_data_tail_union_bound": data_union,
            "archive_hoeffding_union_bound": archive_union,
        })
    return {"R0": R0, "entropy_ball_lengths_checked": 256,
            "maximum_entropy_bound_excess": maximum_entropy_excess, "examples": examples}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = {
        "status": "all_checks_passed",
        "scope": "Finite verification; universal statements are proved in science/clifford_lifetime.md.",
        "five_qubit_pauli_history": five_qubit_histories(),
        "broadcast": broadcast_checks(),
        "workspace": allocation_checks(),
        "finite_parameters": parameter_checks(),
    }
    encoded = json.dumps(results, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
