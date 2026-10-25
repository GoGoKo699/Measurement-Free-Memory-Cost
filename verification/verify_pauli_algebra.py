#!/usr/bin/env python3
"""Independent finite verification of the Pauli-algebra memory bounds.

The dense-matrix calculations do not import the binary circuit simulator.
They verify all two-qubit subspace algebras, non-Clifford envelope updates,
mixed sampled channels, a joint nonstabilizer replacement reset, and a decoder
whose terminal support rank exceeds twice its output width. A separate exact
rational calculation verifies the finite-depth injection LP by matching primal
and dual certificates. These checks supplement, rather than replace, the proofs.

Run from any directory:
    python3 verify_pauli_algebra.py --output pauli_algebra_results.json

Dependency: NumPy. The tested environment is recorded in the JSON result.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import platform
from typing import Iterable

import numpy as np


N_QUBITS = 2
DIMENSION = 1 << N_QUBITS
TOLERANCE = 1e-9
SUPPORT_TOLERANCE = 1e-10
RANDOM_SEED = 72461
RANDOM_CIRCUITS = 200
STEPS_PER_CIRCUIT = 24

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
SINGLE_PAULIS = {(0, 0): I2, (1, 0): X, (1, 1): Y, (0, 1): Z}
H = (X + Z) / np.sqrt(2)
S = np.diag([1, 1j]).astype(complex)
T = np.diag([1, np.exp(1j * np.pi / 4)])
CNOT = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
    dtype=complex,
)
MAGIC_STATE = (I2 + (X + Y + Z) / np.sqrt(3)) / 2

Subspace = frozenset[int]


class VerificationError(RuntimeError):
    """An exact or numerical verification condition failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def residual(matrix: np.ndarray) -> float:
    return float(np.max(np.abs(matrix)))


def span(generators: Iterable[int]) -> Subspace:
    labels = {0}
    for generator in generators:
        labels.update(label ^ generator for label in tuple(labels))
    return frozenset(labels)


def dimension(subspace: Subspace) -> int:
    size = len(subspace)
    require(size > 0 and size & (size - 1) == 0, "Invalid subspace size")
    return size.bit_length() - 1


def binary_rank(rows: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def symplectic(v: int, w: int, n: int = N_QUBITS) -> int:
    mask = (1 << n) - 1
    return (((v & mask) & (w >> n)).bit_count()
            + ((v >> n) & (w & mask)).bit_count()) % 2


def quantum_rank(subspace: Subspace) -> int:
    """Return half the rank of the restricted symplectic form."""
    basis: list[int] = []
    generated = frozenset({0})
    for label in sorted(subspace):
        if label not in generated:
            basis.append(label)
            generated = span(basis)
    gram_rows = [sum(symplectic(v, w) << j for j, w in enumerate(basis))
                 for v in basis]
    rank = binary_rank(gram_rows)
    require(rank % 2 == 0, "Alternating form has odd rank")
    return rank // 2


def all_subspaces() -> list[Subspace]:
    spaces = {frozenset({0})}
    pending = [frozenset({0})]
    while pending:
        current = pending.pop()
        for label in range(1, 1 << (2 * N_QUBITS)):
            enlarged = span(itertools.chain(current, [label]))
            if enlarged not in spaces:
                spaces.add(enlarged)
                pending.append(enlarged)
    return sorted(spaces, key=lambda space: (len(space), tuple(sorted(space))))


def pauli(label: int, n: int = N_QUBITS) -> np.ndarray:
    operator = np.array([[1]], dtype=complex)
    for j in range(n):
        factor = SINGLE_PAULIS[((label >> j) & 1, (label >> (n + j)) & 1)]
        operator = np.kron(operator, factor)
    return operator


PAULIS = np.stack([pauli(v) for v in range(1 << (2 * N_QUBITS))])


def coefficients(operator: np.ndarray) -> np.ndarray:
    return np.einsum("vij,ij->v", PAULIS.conj(), operator) / DIMENSION


def support(operator: np.ndarray) -> set[int]:
    return set(np.flatnonzero(np.abs(coefficients(operator)) > SUPPORT_TOLERANCE))


def embed(operator: np.ndarray, qubit: int) -> np.ndarray:
    return np.kron(operator, I2) if qubit == 0 else np.kron(I2, operator)


def reset_adjoint(operator: np.ndarray, qubit: int, state: np.ndarray) -> np.ndarray:
    weighted = (embed(state, qubit) @ operator).reshape(2, 2, 2, 2)
    reduced = np.trace(weighted, axis1=qubit, axis2=qubit + 2)
    return np.kron(I2, reduced) if qubit == 0 else np.kron(reduced, I2)


def verify_algebras(spaces: list[Subspace]) -> dict[str, int | float]:
    min_choi_eigenvalue = 0.0
    max_trace_residual = 0.0
    max_twirl_residual = 0.0
    max_partial_transpose_residual = 0.0
    for space in spaces:
        choi = sum(np.kron(PAULIS[v], PAULIS[v].T) for v in space) / DIMENSION
        eigenvalues = np.linalg.eigvalsh(choi)
        min_choi_eigenvalue = min(min_choi_eigenvalue, float(eigenvalues.min()))
        require(eigenvalues.min() >= -TOLERANCE, "Conditional expectation is not CP")
        trace_out = np.trace(
            choi.reshape(DIMENSION, DIMENSION, DIMENSION, DIMENSION),
            axis1=0, axis2=2,
        )
        tp_error = residual(trace_out - np.eye(DIMENSION))
        max_trace_residual = max(max_trace_residual, tp_error)
        require(tp_error < TOLERANCE, "Conditional expectation is not TP")

        orthogonal = [w for w in range(len(PAULIS))
                      if all(symplectic(v, w) == 0 for v in space)]
        for v, operator in enumerate(PAULIS):
            twirled = sum(PAULIS[w] @ operator @ PAULIS[w] for w in orthogonal)
            twirled /= len(orthogonal)
            projected = operator if v in space else np.zeros_like(operator)
            error = residual(twirled - projected)
            max_twirl_residual = max(max_twirl_residual, error)
            require(error < TOLERANCE, "Pauli twirl differs from algebra projection")

        # Independent spectral check of the coherent factor in the algebra's
        # normal form: normalized Choi partial-transpose trace norm is 2**s.
        partial_transpose = (
            (choi / DIMENSION)
            .reshape(DIMENSION, DIMENSION, DIMENSION, DIMENSION)
            .transpose(0, 3, 2, 1)
            .reshape(DIMENSION**2, DIMENSION**2)
        )
        trace_norm = float(np.abs(np.linalg.eigvalsh(partial_transpose)).sum())
        error = abs(trace_norm - 2**quantum_rank(space))
        max_partial_transpose_residual = max(max_partial_transpose_residual, error)
        require(error < TOLERANCE, "Algebra coherent-factor dimension mismatch")

    return {
        "subspaces_checked": len(spaces),
        "minimum_conditional_expectation_choi_eigenvalue": min_choi_eigenvalue,
        "maximum_trace_preservation_residual": max_trace_residual,
        "maximum_twirl_projection_residual": max_twirl_residual,
        "maximum_choi_partial_transpose_trace_norm_error": max_partial_transpose_residual,
    }


def verify_gate_envelopes(spaces: list[Subspace]) -> dict[str, int]:
    t_cases = 0
    axis_cases = 0
    local_cases = 0
    axis_label = 1 | (1 << (N_QUBITS + 1))  # X on qubit 0, Z on qubit 1.
    angle = np.pi / 7
    axis_unitary = (np.cos(angle) * np.eye(DIMENSION)
                    + 1j * np.sin(angle) * PAULIS[axis_label])
    generic_axis = (X + 2 * Y + 3 * Z) / np.sqrt(14)
    generic_unitary = np.cos(angle) * I2 + 1j * np.sin(angle) * generic_axis
    for space in spaces:
        old_rank = dimension(space)
        for qubit in range(N_QUBITS):
            z_label = 1 << (N_QUBITS + qubit)
            envelope = span(itertools.chain(space, [z_label]))
            require(dimension(envelope) <= old_rank + 1, "T rank increment exceeds one")
            for single_gate in (T, T.conj().T):
                unitary = embed(single_gate, qubit)
                for label in space:
                    require(support(unitary.conj().T @ PAULIS[label] @ unitary) <= envelope,
                            "T adjoint escaped its envelope")
                t_cases += 1

            qubit_mask = (1 << qubit) | z_label
            local_envelope = span(itertools.chain(
                (v & ~qubit_mask for v in space), [1 << qubit, z_label],
            ))
            require(dimension(local_envelope) <= old_rank + 2,
                    "One-qubit unitary rank increment exceeds two")
            unitary = embed(generic_unitary, qubit)
            for label in space:
                require(support(unitary.conj().T @ PAULIS[label] @ unitary) <= local_envelope,
                        "Generic one-qubit adjoint escaped its envelope")
            local_cases += 1

        axis_envelope = span(itertools.chain(space, [axis_label]))
        require(dimension(axis_envelope) <= old_rank + 1,
                "Pauli-axis rotation rank increment exceeds one")
        for label in space:
            require(support(axis_unitary.conj().T @ PAULIS[label] @ axis_unitary) <= axis_envelope,
                    "Nonlocal Pauli-axis rotation escaped its envelope")
        axis_cases += 1
    return {"T_and_T_dagger_cases": t_cases,
            "nonlocal_Pauli_rotation_cases": axis_cases,
            "generic_one_qubit_unitary_cases": local_cases}


def verify_joint_reset() -> dict[str, int | float]:
    # A nonstabilizer entangled state on two reset qubits, with a third retained.
    vector = np.array([1, 0, np.exp(1j * np.pi / 7), 1], dtype=complex) / np.sqrt(3)
    state = np.outer(vector, vector.conj())
    kraus = [np.kron(np.outer(vector, np.eye(4)[j]), I2) for j in range(4)]
    largest_error = 0.0
    cases = 0
    for reset_pauli in PAULIS:
        for retained_pauli in (I2, X, Y, Z):
            observable = np.kron(reset_pauli, retained_pauli)
            pulled_back = sum(k.conj().T @ observable @ k for k in kraus)
            expected = np.trace(state @ reset_pauli) * np.kron(np.eye(4), retained_pauli)
            error = residual(pulled_back - expected)
            largest_error = max(largest_error, error)
            require(error < TOLERANCE, "Joint reset adjoint formula failed")
            cases += 1
    return {"joint_reset_Pauli_cases": cases,
            "maximum_joint_reset_adjoint_residual": largest_error}


def verify_decoder_support() -> dict[str, int | float]:
    unitary = CNOT @ embed(T @ H @ T, 1)
    decoded = [unitary.conj().T @ embed(p, 0) @ unitary for p in (I2, X, Y, Z)]
    envelope = span(label for operator in decoded for label in support(operator))
    expected_x = np.kron(X, Z / np.sqrt(2) + (X + Y) / 2)
    formula_error = max(residual(decoded[1] - expected_x), residual(decoded[3] - embed(Z, 0)))
    require(formula_error < TOLERANCE, "Explicit decoder pullback formula failed")
    require(dimension(envelope) == 4, "One-output-qubit decoder should have rank-four support")
    clifford_envelope = span(
        label for p in (I2, X, Y, Z)
        for label in support(CNOT.conj().T @ embed(p, 0) @ CNOT)
    )
    require(dimension(clifford_envelope) == 2, "Clifford decoder should have rank-two support")
    return {"output_qubits": 1,
            "non_Clifford_decoder_support_rank": dimension(envelope),
            "Clifford_decoder_support_rank": dimension(clifford_envelope),
            "maximum_decoder_formula_residual": formula_error}


def random_unitary(rng: np.random.Generator, size: int) -> np.ndarray:
    sample = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    unitary, triangular = np.linalg.qr(sample)
    phases = np.diag(triangular)
    return unitary * (phases / np.abs(phases))


def verify_logical_channel_rank(spaces: list[Subspace]) -> dict[str, int | float]:
    """Use Kraus channels with an averaged prefix, algebra cut, and decoder."""
    rng = np.random.default_rng(RANDOM_SEED + 1)
    max_identity_error = 0.0
    max_tp_error = 0.0
    max_radius = 0.0
    max_bound_violation = 0.0
    cases = 0
    for logical_dimension in (2, 4):
        for space in spaces:
            # Two independent isometries with fixed probabilities form a CPTP
            # prefix; this deliberately need not be an isometric encoder.
            first = random_unitary(rng, DIMENSION)[:, :logical_dimension]
            second = random_unitary(rng, DIMENSION)[:, :logical_dimension]
            prefix = [np.sqrt(0.3) * first, np.sqrt(0.7) * second]
            decoder_unitary = random_unitary(rng, DIMENSION)
            if logical_dimension == DIMENSION:
                decoder = [decoder_unitary]
            else:
                # Discard the second physical qubit after a generic unitary.
                decoder = [decoder_unitary[[j, 2 + j], :] for j in range(2)]
            orthogonal = [w for w in range(len(PAULIS))
                          if all(symplectic(v, w) == 0 for v in space)]
            kraus = [d @ PAULIS[w] @ p / np.sqrt(len(orthogonal))
                     for d in decoder for w in orthogonal for p in prefix]
            completeness = sum(operator.conj().T @ operator for operator in kraus)
            tp_error = residual(completeness - np.eye(logical_dimension))
            max_tp_error = max(max_tp_error, tp_error)
            require(tp_error < TOLERANCE, "Logical fixture is not trace preserving")

            superoperator = sum(np.kron(operator.conj(), operator) for operator in kraus)
            channel_rank = int(np.linalg.matrix_rank(superoperator, tol=TOLERANCE))
            require(channel_rank <= len(space), "Logical channel rank exceeds algebra dimension")
            radius = float(np.max(np.abs(np.linalg.eigvals(superoperator))))
            max_radius = max(max_radius, radius)
            require(radius <= 1 + TOLERANCE, "CPTP superoperator has an expansive eigenvalue")

            # Choi overlap is calculated from vectorized Kraus matrices, apart
            # from the superoperator representation used above.
            choi = sum(np.outer(operator.reshape(-1, order="F"),
                                operator.reshape(-1, order="F").conj())
                       for operator in kraus)
            bell = np.eye(logical_dimension).reshape(-1, order="F") / np.sqrt(logical_dimension)
            fidelity = float(np.real(bell.conj() @ choi @ bell) / logical_dimension)
            trace_fidelity = np.trace(superoperator) / logical_dimension**2
            identity_error = abs(fidelity - trace_fidelity)
            max_identity_error = max(max_identity_error, float(identity_error))
            require(identity_error < TOLERANCE, "Bell fidelity differs from superoperator trace")
            require(fidelity <= channel_rank / logical_dimension**2 + TOLERANCE,
                    "Superoperator rank fidelity bound failed")
            bound = min(1.0, len(space) / logical_dimension**2,
                        2**quantum_rank(space) / logical_dimension)
            violation = max(0.0, fidelity - bound)
            max_bound_violation = max(max_bound_violation, violation)
            require(violation < TOLERANCE, "Algebra fidelity bottleneck failed")
            cases += 1
    return {"logical_CPTP_channels": cases,
            "maximum_trace_preservation_residual": max_tp_error,
            "maximum_fidelity_superoperator_trace_residual": max_identity_error,
            "maximum_channel_eigenvalue_modulus": max_radius,
            "maximum_algebra_fidelity_bound_violation": max_bound_violation}


def verify_mixed_circuits() -> dict[str, int | float]:
    rng = np.random.default_rng(RANDOM_SEED)
    cliffords = [embed(gate, j) for gate in (H, S) for j in range(N_QUBITS)] + [CNOT]
    max_factor_error = 0.0
    max_replacement_error = 0.0
    growth_events = 0
    zero_envelope_cuts = 0
    for _ in range(RANDOM_CIRCUITS):
        space = frozenset(range(len(PAULIS)))
        images = PAULIS.copy()
        for _ in range(STEPS_PER_CIRCUIT):
            kind = int(rng.integers(5))
            qubit = int(rng.integers(N_QUBITS))
            old_rank = dimension(space)
            if kind == 0:
                unitary = cliffords[int(rng.integers(len(cliffords)))]
                mapped_labels = []
                for label in space:
                    mapped = support(unitary.conj().T @ PAULIS[label] @ unitary)
                    require(len(mapped) == 1, "Clifford image has multiple Pauli terms")
                    mapped_labels.extend(mapped)
                space = span(mapped_labels)
                images = unitary.conj().T @ images @ unitary
            elif kind == 1:
                gate = T if rng.random() < 0.5 else T.conj().T
                unitary = embed(gate, qubit)
                space = span(itertools.chain(space, [1 << (N_QUBITS + qubit)]))
                images = unitary.conj().T @ images @ unitary
                growth_events += int(dimension(space) > old_rank)
            elif kind in (2, 3):
                mask = (1 << qubit) | (1 << (N_QUBITS + qubit))
                if kind == 2:
                    space = span(v & ~mask for v in space)
                    state = MAGIC_STATE
                else:
                    space = frozenset(v for v in space if not v & mask)
                    state = I2 / 2
                images = np.stack([reset_adjoint(operator, qubit, state) for operator in images])
            else:
                measured_label = int(rng.integers(1, len(PAULIS)))
                measured = PAULIS[measured_label]
                space = frozenset(v for v in space if symplectic(v, measured_label) == 0)
                images = (images + measured @ images @ measured) / 2

            expansion = np.einsum("vij,aij->av", PAULIS.conj(), images) / DIMENSION
            forbidden = [v for v in range(len(PAULIS)) if v not in space]
            error = residual(expansion[:, forbidden]) if forbidden else 0.0
            max_factor_error = max(max_factor_error, error)
            require(error < TOLERANCE, "Sampled tail adjoint escaped its Pauli algebra")
            if dimension(space) == 0:
                scalar_images = expansion[:, 0, None, None] * np.eye(DIMENSION)
                error = residual(images - scalar_images)
                max_replacement_error = max(max_replacement_error, error)
                require(error < TOLERANCE, "Zero envelope failed replacement-channel criterion")
                zero_envelope_cuts += 1
    require(growth_events > 0 and zero_envelope_cuts > 0,
            "Random fixtures did not exercise both rank growth and replacement")
    return {"circuits": RANDOM_CIRCUITS,
            "steps_per_circuit": STEPS_PER_CIRCUIT,
            "seed": RANDOM_SEED,
            "observed_T_envelope_rank_growth_events": growth_events,
            "zero_envelope_cuts": zero_envelope_cuts,
            "maximum_factorization_residual": max_factor_error,
            "maximum_replacement_criterion_residual": max_replacement_error}


def verify_injection_case(depth: int, initial: Fraction, floor: Fraction,
                          q: Fraction, cap: Fraction) -> bool:
    """Check exact primal/dual certificates; return whether cap is feasible."""
    zero = Fraction(0)
    require(initial >= floor > 0 and cap >= initial and 0 < q < 1,
            "Invalid scalar verification instance")
    rank = initial
    injections: list[Fraction] = []
    pre_noise: list[Fraction] = []
    for _ in range(depth):
        injection = max(zero, floor / q - rank)
        injections.append(injection)
        pre_noise.append(rank + injection)
        rank = q * (rank + injection)
        require(rank >= floor, "Greedy scalar trajectory violated floor")

    objectives = [zero]
    for start in range(1, depth + 1):
        objectives.append(floor / q + (depth - start) * (1 - q) * floor / q
                          - initial * q**(start - 1))
    optimum = max(objectives)
    require(sum(injections, zero) == optimum, "Primal and dual scalar costs differ")

    if optimum > 0:
        start = objectives.index(optimum)
        dual = [zero] * depth
        for j in range(start - 1, depth - 1):
            dual[j] = (1 - q) / q
        dual[-1] = 1 / q
    else:
        dual = [zero] * depth

    dual_value = sum((dual[j - 1] * (floor - initial * q**j)
                      for j in range(1, depth + 1)), zero)
    require(dual_value == optimum, "Dual certificate has wrong objective")
    for i in range(1, depth + 1):
        coefficient = sum((dual[j - 1] * q**(j - i + 1)
                           for j in range(i, depth + 1)), zero)
        require(coefficient <= 1, "Dual certificate is infeasible")
        require(injections[i - 1] * (1 - coefficient) == 0,
                "Variable complementary slackness failed")
    for j in range(1, depth + 1):
        lhs = sum((q**(j - i + 1) * injections[i - 1]
                   for i in range(1, j + 1)), zero)
        slack = lhs - (floor - initial * q**j)
        require(slack >= 0, "Primal certificate is infeasible")
        require(dual[j - 1] * slack == 0, "Constraint complementary slackness failed")

    coast_rank = initial
    coast_steps = 0
    while coast_steps < depth and q * coast_rank >= floor:
        coast_rank *= q
        coast_steps += 1
    closed_form = (zero if coast_steps == depth else
                   floor / q - coast_rank
                   + (depth - coast_steps - 1) * (1 - q) * floor / q)
    require(closed_form == optimum, "Coasting closed form differs from LP optimum")
    additive_bound = max(zero, (1 - q) * floor * depth / q - (initial - floor))
    require(optimum >= additive_bound, "Additive corollary exceeded exact optimum")

    cap_feasible = depth == 0 or floor <= q * cap
    if cap_feasible:
        require(all(value <= cap for value in pre_noise), "Feasible cap exceeded")
    else:
        require(pre_noise[0] > cap, "Infeasible cap was not detected")
    return cap_feasible


def verify_injection_optimization() -> dict[str, int | str]:
    rates = (Fraction(1, 5), Fraction(1, 2), Fraction(3, 4), Fraction(9, 10))
    cases = 0
    cap_feasible_cases = 0
    equality_cases = 0
    for q in rates:
        for cap in (Fraction(2), Fraction(4), Fraction(8)):
            floors = {Fraction(1, 2), Fraction(1), Fraction(3, 2), cap / 2, q * cap, cap}
            for floor in sorted(floors):
                if floor > cap:
                    continue
                initials = {floor, (floor + cap) / 2, cap}
                for initial in sorted(initials):
                    for depth in range(11):
                        cap_feasible_cases += int(
                            verify_injection_case(depth, initial, floor, q, cap)
                        )
                        cases += 1
        # Boundary cases M q^t = a are handled without floating-point logarithms.
        for floor in (Fraction(1), Fraction(3, 2)):
            for coast_steps in range(4):
                initial = floor / q**coast_steps
                half = initial / 2
                cap = Fraction(2 * ((half.numerator + half.denominator - 1) // half.denominator))
                for depth in (coast_steps, coast_steps + 1, coast_steps + 4):
                    cap_feasible_cases += int(
                        verify_injection_case(depth, initial, floor, q, cap)
                    )
                    cases += 1
                    equality_cases += 1
    return {"arithmetic": "exact rational (fractions.Fraction)",
            "primal_dual_certificate_cases": cases,
            "cap_feasible_cases": cap_feasible_cases,
            "cap_infeasible_cases": cases - cap_feasible_cases,
            "exact_coasting_boundary_cases": equality_cases}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Also save the JSON report at this path")
    args = parser.parse_args()
    spaces = all_subspaces()
    require(len(spaces) == 67, "Enumeration did not produce all two-qubit subspaces")
    report = {
        "status": "PASS",
        "scope": "finite operator identities and exact scalar LP certificates",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "numerical_tolerance": TOLERANCE,
        "algebras": verify_algebras(spaces),
        "gate_envelopes": verify_gate_envelopes(spaces),
        "joint_reset": verify_joint_reset(),
        "terminal_decoder": verify_decoder_support(),
        "logical_channel_rank": verify_logical_channel_rank(spaces),
        "mixed_circuits": verify_mixed_circuits(),
        "injection_optimization": verify_injection_optimization(),
    }
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
