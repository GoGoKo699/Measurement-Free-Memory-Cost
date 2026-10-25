#!/usr/bin/env python3
"""Finite checks for the explicit five-qubit coherent recovery.

Run from the repository root:
    python3 verification/verify_five_qubit.py
Only deterministic small-system identities are tested here. General lifetime
and concatenation guarantees are proved in science/coherent_recovery.md.
"""
from __future__ import annotations

import itertools
import json
import math
from collections import Counter

import numpy as np

GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
LOGICAL_X, LOGICAL_Z = "XXXXX", "ZZZZZ"
SINGLE = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.diag([1, -1]).astype(complex),
}


def sp(a: str, b: str) -> int:
    return sum(x != "I" and y != "I" and x != y
               for x, y in zip(a, b)) % 2


def matrix(label: str) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for ch in label:
        out = np.kron(out, SINGLE[ch])
    return out


def syndrome(label: str) -> int:
    return sum(sp(label, g) << j for j, g in enumerate(GENERATORS))


def correction(s: int) -> tuple[int, int]:
    a, b, c, d = ((s >> j) & 1 for j in range(4))
    return a ^ d ^ (a & c) ^ (b & c) ^ (b & d), \
        b ^ c ^ (a & b) ^ (a & d) ^ (c & d)


def pauli_bits(label: str) -> int:
    x = sum((ch in "XY") << j for j, ch in enumerate(label))
    z = sum((ch in "YZ") << j for j, ch in enumerate(label))
    return x | (z << len(label))


def binary_rank(rows: list[int]) -> int:
    basis = {}
    for row in rows:
        while row:
            j = row.bit_length() - 1
            if j in basis:
                row ^= basis[j]
            else:
                basis[j] = row
                break
    return len(basis)


def run() -> dict:
    labels = ["".join(p) for p in itertools.product("IXYZ", repeat=5)]
    assert binary_rank([pauli_bits(g) for g in GENERATORS]) == 4
    assert all(sp(a, b) == 0 for a in GENERATORS for b in GENERATORS)
    assert all(sp(g, l) == 0 for g in GENERATORS
               for l in (LOGICAL_X, LOGICAL_Z))
    assert sp(LOGICAL_X, LOGICAL_Z) == 1

    correctable = [label for label in labels
                   if sum(ch != "I" for ch in label) <= 1]
    assert len(correctable) == 16
    assert len({syndrome(label) for label in correctable}) == 16
    for label in correctable:
        assert correction(syndrome(label)) == (
            sp(label, LOGICAL_Z), sp(label, LOGICAL_X)
        )

    # Construct a full symplectic basis independently of the correction ANFs.
    destabilizers = []
    for j in range(4):
        candidates = [
            label for label in labels
            if all(sp(label, g) == int(i == j)
                   for i, g in enumerate(GENERATORS))
            and sp(label, LOGICAL_X) == sp(label, LOGICAL_Z) == 0
            and all(sp(label, old) == 0 for old in destabilizers)
        ]
        assert candidates
        destabilizers.append(candidates[0])

    eye = np.eye(32, dtype=complex)
    projector = eye.copy()
    for g in GENERATORS:
        projector = projector @ (eye + matrix(g)) / 2
    assert np.linalg.matrix_rank(projector, tol=1e-10) == 2
    zero = projector @ (eye + matrix(LOGICAL_Z))[:, 0] / 2
    zero /= np.linalg.norm(zero)
    one = matrix(LOGICAL_X) @ zero
    encoder = np.column_stack([zero, one])
    assert np.max(np.abs(encoder.conj().T @ encoder - np.eye(2))) < 1e-12

    # Columns are D(s)|b_L>, with index 16*b+s. This fixes a Clifford
    # decoder mapping the stated stabilizers/logical pair to standard wires.
    decoder = np.zeros((32, 32), dtype=complex)
    d_mats = [matrix(label) for label in destabilizers]
    for s in range(16):
        shift = eye.copy()
        for j in range(4):
            if (s >> j) & 1:
                shift = d_mats[j] @ shift
        for b in range(2):
            decoder[16*b+s, :] = (shift @ encoder[:, b]).conj()
    unitary_residual = float(np.max(np.abs(
        decoder @ decoder.conj().T - eye
    )))
    assert unitary_residual < 1e-12

    basis_residual = 0.0
    for j, (g, d) in enumerate(zip(GENERATORS, destabilizers)):
        expected_z = np.diag([(-1)**((i >> j) & 1) for i in range(32)])
        expected_x = np.zeros((32, 32), dtype=complex)
        for i in range(32):
            expected_x[i ^ (1 << j), i] = 1
        for label, expected in ((g, expected_z), (d, expected_x)):
            residual = np.max(np.abs(
                decoder @ matrix(label) @ decoder.conj().T - expected
            ))
            basis_residual = max(basis_residual, float(residual))
    for label, gate in ((LOGICAL_X, SINGLE["X"]),
                        (LOGICAL_Z, SINGLE["Z"])):
        expected = np.kron(gate, np.eye(16))
        residual = np.max(np.abs(
            decoder @ matrix(label) @ decoder.conj().T - expected
        ))
        basis_residual = max(basis_residual, float(residual))
    assert basis_residual < 1e-12

    bad_by_weight = Counter()
    total_by_weight = Counter()
    pauli_residual = 0.0
    single_identity_residual = 0.0
    for label in labels:
        s = syndrome(label)
        x, z = correction(s)
        decoded = decoder @ matrix(label) @ encoder
        projected = decoded[[s, 16+s], :]
        other = decoded.copy()
        other[[s, 16+s], :] = 0
        assert np.max(np.abs(other)) < 1e-12
        corrected = (np.linalg.matrix_power(SINGLE["Z"], z)
                     @ np.linalg.matrix_power(SINGLE["X"], x)
                     @ projected)
        rx, rz = x ^ sp(label, LOGICAL_Z), z ^ sp(label, LOGICAL_X)
        expected = (np.linalg.matrix_power(SINGLE["X"], rx)
                    @ np.linalg.matrix_power(SINGLE["Z"], rz))
        phase = np.trace(expected.conj().T @ corrected) / 2
        residual = float(np.max(np.abs(corrected - phase * expected)))
        pauli_residual = max(pauli_residual, residual)
        assert abs(abs(phase) - 1) < 1e-12
        assert residual < 1e-12
        w = sum(ch != "I" for ch in label)
        total_by_weight[w] += 1
        if rx or rz:
            bad_by_weight[w] += 1
        if w <= 1:
            single_identity_residual = max(
                single_identity_residual,
                float(np.max(np.abs(corrected - phase * np.eye(2))))
            )
    assert single_identity_residual < 1e-12
    assert [total_by_weight[w] for w in range(6)] == [
        math.comb(5, w) * 3**w for w in range(6)
    ]
    assert bad_by_weight[0] == bad_by_weight[1] == 0

    def exact_failure(p: float) -> float:
        return sum(bad_by_weight[w] * (p/3)**w * (1-p)**(5-w)
                   for w in range(6))

    for p in (0, .0001, .001, .01, .0375, .1, .3, .7, 1):
        assert exact_failure(p) <= 10*p*p + 1e-12

    parameter_cases = 0
    max_failure_ratio = 0.0
    for k, duration, gamma, eps in itertools.product(
        (1, 2, 7), (1, 2, 21, 101, 10001),
        (1e-4, .001, .01, .05), (.001, .01, .1, .25)
    ):
        tau = math.floor(1/(20*gamma))
        intervals = (duration + tau - 1) // tau
        assert tau >= 1/(40*gamma) - 1e-12
        assert intervals <= 40*gamma*duration + 1 + 1e-12
        a = max(1.0, math.log(k*intervals/(10*eps))/math.log(8/3))
        level = math.ceil(math.log2(a))
        fail_bound = k*intervals*math.exp((2**level)*math.log(3/8))/10
        assert fail_bound <= eps*(1+1e-12)
        assert 5**level <= 5*a**math.log2(5)*(1+1e-12)
        max_failure_ratio = max(max_failure_ratio, fail_bound/eps)
        parameter_cases += 1

    return {
        "status": "PASS",
        "stabilizer_rank": 4,
        "correctable_syndromes": 16,
        "all_pauli_errors_checked": len(labels),
        "destabilizers": destabilizers,
        "bad_paulis_by_weight": [bad_by_weight[w] for w in range(6)],
        "quadratic_terms": {"logical_X": 3, "logical_Z": 3},
        "T_count_per_recovery_node_upper": 42,
        "additional_work_qubits": 0,
        "parameter_cases": parameter_cases,
        "max_failure_bound_over_target": max_failure_ratio,
        "max_decoder_unitarity_residual": unitary_residual,
        "max_clifford_basis_residual": basis_residual,
        "max_residual_pauli_check_error": pauli_residual,
        "max_single_error_recovery_residual": single_identity_residual,
        "scope": "Finite algebra and dense-matrix identities, not an asymptotic proof."
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))

