#!/usr/bin/env python3
"""Finite structural checks for the coherent Tanner compilation theorem.

This does not instantiate a large Tanner code or certify its decoder theorem.
The code-family/decoder premises are inherited; the universal argument is in
science/tanner_memory.md. Here we check a nonlinear lookup's full Choi channel,
explicit compilation counts, bounded-gap length formulas and scalar margins.
"""
from __future__ import annotations

from fractions import Fraction as F
import itertools
import json
import math
import numpy as np

from verify_five_qubit import correction, syndrome, sp, LOGICAL_X, LOGICAL_Z


def gates_for_table(table):
    """Two syndrome inputs4,5; minterm work6..9; output10,11."""
    gates = []
    for s in range(4):
        negatives = [4+j for j in range(2) if not (s >> j) & 1]
        for j in negatives:
            gates.append(('X', j))
        gates.append(('CCX', 4, 5, 6+s))
        for j in reversed(negatives):
            gates.append(('X', j))
        for j in range(2):
            if (table[s] >> j) & 1:
                gates.append(('CX', 6+s, 10+j))
    return gates


def gate_on_basis(x, gate):
    name, *bits = gate
    if name == 'X':
        return x ^ (1 << bits[0])
    if name == 'CX':
        return x ^ (((x >> bits[0]) & 1) << bits[1])
    if name == 'CCX':
        return x ^ ((((x >> bits[0]) & 1) & ((x >> bits[1]) & 1)) << bits[2])
    raise ValueError(name)


def apply_gate(psi, gate):
    if gate[0] in ('X', 'CX', 'CCX'):
        out = np.empty_like(psi)
        for j in range(len(psi)):
            out[gate_on_basis(j, gate)] = psi[j]
        return out
    if gate[0] == 'CZ':
        _, a, b = gate
        return np.array([(-v if ((j>>a)&1 and (j>>b)&1) else v)
                         for j, v in enumerate(psi)])
    if gate[0] == 'H':
        _, b = gate
        out = psi.copy()
        for j in range(len(psi)):
            if not (j >> b) & 1:
                k = j ^ (1 << b)
                out[j] = (psi[j]+psi[k])/math.sqrt(2)
                out[k] = (psi[j]-psi[k])/math.sqrt(2)
        return out
    raise ValueError(gate[0])


def exhaustive_lookup():
    cases = 0
    for table in itertools.product(range(4), repeat=4):
        gates = gates_for_table(table)
        assert sum(g[0] == 'CCX' for g in gates) == 4
        for s in range(4):
            x = s << 4
            y = x
            for gate in gates:
                y = gate_on_basis(y, gate)
            assert (y >> 10) == table[s]
            assert ((y >> 4) & 3) == s
            for gate in reversed(gates):
                y = gate_on_basis(y, gate)
            assert y == x
            cases += 1
    return cases


def choi_channel():
    # Nonlinear correction: fX=s0*s1, fZ=s0 XOR s1.
    table = tuple((((s&1)*((s>>1)&1)) | (((s&1)^((s>>1)&1)) << 1))
                  for s in range(4))
    gates = gates_for_table(table)
    psi = np.zeros(2**12, dtype=complex)
    initial = np.zeros(16, dtype=complex)
    for r in range(4):
        psi[r | (r << 2)] = .5
        initial[r | (r << 2)] = .5
    # Data2,3; references0,1. Extract commuting ZZ and XX syndromes4,5.
    extraction = [('CX', 2, 4), ('CX', 3, 4),
                  ('H', 2), ('H', 3), ('CX', 2, 5), ('CX', 3, 5),
                  ('H', 2), ('H', 3)]
    for gate in extraction + gates + [('CX', 10, 2), ('CZ', 11, 2)] + list(reversed(gates)):
        psi = apply_gate(psi, gate)
    dirty_work = float(sum(abs(v)**2 for j, v in enumerate(psi) if j >> 6))
    assert dirty_work < 1e-25
    reshaped = psi.reshape(-1, 16)
    actual = reshaped.T @ reshaped.conj()
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.diag([1, -1]).astype(complex)
    S0, S1 = np.kron(Z, Z), np.kron(X, X)
    X0, Z0 = np.kron(I2, X), np.kron(I2, Z)
    expected = np.zeros((16, 16), dtype=complex)
    for s in range(4):
        projector = ((np.eye(4)+(-1)**(s&1)*S0) @
                     (np.eye(4)+(-1)**((s>>1)&1)*S1))/4
        P = (Z0 if (table[s]>>1)&1 else np.eye(4)) @ (X0 if table[s]&1 else np.eye(4))
        vec = np.kron(P @ projector, np.eye(4)) @ initial
        expected += np.outer(vec, vec.conj())
    error = float(np.max(np.abs(actual-expected)))
    assert error < 1e-13
    return {'input_dimension': 4, 'reference_dimension': 4,
            'nonlinear_table': list(table), 'compute_uncompute_toffolis': 8,
            'prescribed_T_slots': 56, 'work_leak_probability': dirty_work,
            'choi_max_entry_error': error}


def scalar_and_resource_checks():
    checks = 0
    # These are algebraic-parameter checks, not claims these small local codes exist.
    for b in range(1, 9):
        for t in range(1, 17):
            L_per_N = F(4+8*t, b)
            S_per_N = L_per_N*(4*b*2**b+4*b)+32*(1+t)
            assert 7+S_per_N <= 64*(1+t)*b*2**b
            assert 14*L_per_N*(b-1)*2**b <= 112*(1+t)*b*2**b
            checks += 1
    for a, A, theta, d0 in itertools.product(
            [F(1,2), F(1,3), F(1,100)], [F(1,10), F(3), F(100)],
            [F(1,10), F(1,1000)], [F(1,10), F(1,1000)]):
        u=min(theta/(4*A), d0/8, F(1,8))
        v=a*u/(1-a)
        assert v<=u and A*(v+u)<=theta
        assert a*(v+u)==v and 2*v<d0/2
        for div in [F(2), F(5,2), F(7), F(100)]:
            gamma=u/div
            h=u//gamma
            assert h>=u/(2*gamma)
            assert 1-(1-gamma)**h <= h*gamma <= u
        checks += 1
    ratios=0
    for Q in [3,5,7,9]:
        for i in range(3,11):
            N0=Q**i*(Q**(2*i)-1)
            N1=Q**(i+1)*(Q**(2*i+2)-1)
            assert F(N1,N0)<2*Q**3
            ratios+=1
    return {'resource_and_invariant_cases': checks, 'exact_length_ratio_cases': ratios}


def sufficient_size_example():
    # A scalar resource fixture, NOT a claim that these are admissible Tanner
    # code parameters or an experimentally attained noise threshold.
    u, gamma = F(1, 32), F(1, 1024)
    D, epsilon, k, rate = 1_000_000, F(1, 2**20), 128, F(1, 4)
    h = u // gamma
    J = (D+h-1)//h
    # Hoeffding gives c >= u^2/2. Integer L has 2^-L <= eps/(2J),
    # hence also exp(-L) <= eps/(2J); this is an exact sufficient certificate.
    ratio = F(2*J, 1)/epsilon
    L = ratio.numerator.bit_length() - ratio.denominator.bit_length() + 1
    while F(1,2**L) > 1/ratio:
        L += 1
    N = max(math.ceil(k/rate), math.ceil(F(2*L,1)/(u*u)))
    assert u*u*N/2 >= L
    certified_failure = F(2*J, 2**L)
    assert certified_failure <= epsilon
    assert J <= 2*gamma*D/u+1
    return {'interpretation': 'Conditional scalar fixture only; no Tanner block simulated.',
            'u': str(u), 'gamma': str(gamma), 'D': D, 'epoch_length': h,
            'epochs': J, 'logical_qubits': k, 'rate_lower_bound': str(rate),
            'target_error': str(epsilon), 'sufficient_data_length': N,
            'failure_bound_exact': str(certified_failure),
            'certificate': 'Hoeffding c>=u^2/2; 2J exp(-u^2 N/2)<=2J*2^-L<=epsilon.'}


def hybrid_inner_checks():
    """Exact Pauli covariance and biased-noise checks for the inner recovery.

    The independently checked dense recovery identity lives in verify_five_qubit.
    Here the extra issue is that a nonzero outer residual already exists when
    fresh physical noise arrives. No independence of its X/Z bits is assumed.
    """
    labels = list(itertools.product('IXYZ', repeat=5))
    bit = {'I': 0, 'X': 1, 'Z': 2, 'Y': 3}
    char = {value: key for key, value in bit.items()}

    def residual(label):
        fx, fz = correction(syndrome(label))
        return (fx ^ sp(label, LOGICAL_Z)) | ((fz ^ sp(label, LOGICAL_X)) << 1)

    covariance_cases = 0
    for label in labels:
        original_syndrome = syndrome(label)
        original_residual = residual(label)
        for prior in 'IXYZ':
            # XXXXX and ZZZZZ represent logical X and Z; YYYYY differs from
            # their logical-Y product only by a physically irrelevant phase.
            shifted = tuple(char[bit[p] ^ bit[prior]] for p in label)
            assert syndrome(shifted) == original_syndrome
            assert residual(shifted) == original_residual ^ bit[prior]
            covariance_cases += 1

    depolarizing = (F(77,80), F(1,80), F(1,80), F(1,80))
    pure_x = (F(77,80), F(3,80), F(0), F(0))
    biased = (F(313,320), F(1,80), F(1,160), F(1,320))
    no_error = (F(1), F(0), F(0), F(0))
    fixtures = [
        ('depolarizing', [depolarizing]*5),
        ('pure_x', [pure_x]*5),
        ('biased', [biased]*5),
        ('different_positions', [depolarizing, pure_x, biased, no_error, biased]),
    ]
    results = []
    for name, distributions in fixtures:
        assert all(sum(p) == 1 and min(p) >= 0 for p in distributions)
        output = [F(0)]*4  # I,X,Z,Y order follows the binary residual label.
        for label in labels:
            probability = math.prod(distributions[j]['IXYZ'.index(p)]
                                    for j, p in enumerate(label))
            output[residual(label)] += probability
        assert sum(output) == 1
        physical_failure = [1-p[0] for p in distributions]
        pair_bound = sum(physical_failure[i]*physical_failure[j]
                         for i in range(5) for j in range(i+1,5))
        nonidentity = 1-output[0]
        assert nonidentity <= pair_bound <= 10*max(physical_failure)**2
        x_marginal, z_marginal = output[1]+output[3], output[2]+output[3]
        assert x_marginal <= nonidentity and z_marginal <= nonidentity
        correlation = output[3]-x_marginal*z_marginal
        if name == 'depolarizing':
            assert correlation != 0  # A concrete check that X/Z independence is false.
        results.append({'fixture': name, 'logical_nonidentity_probability': str(nonidentity),
                        'pair_union_bound': str(pair_bound),
                        'X_Z_joint_minus_product': str(correlation)})
    return {'prior_logical_pauli_covariance_cases': covariance_cases,
            'exact_biased_distribution_cases': len(fixtures),
            'distribution_results': results}


def hybrid_schedule_checks():
    parameter_cases = 0
    choices = []
    for u in [F(1,8), F(1,32), F(1,1000), F(1,2**20), F(1,10**12)]:
        level = 0
        while F(3,8)**(2**level)/10 > u/2:
            level += 1
        eta = F(3,8)**(2**level)/10
        W = 5**level
        assert eta <= u/2
        if level:
            assert F(3,8)**(2**(level-1))/10 > u/2
        node_count = (W-1)//4
        assert node_count == sum(5**j for j in range(level))
        choices.append({'outer_u': str(u), 'constant_inner_levels': level,
                        'inner_width_per_outer_position': W,
                        'inner_T_slots_per_outer_position': 42*node_count,
                        'effective_nonidentity_bound': str(eta)})
        for gamma in [F(1,20), F(1,21), F(1,100), F(1,1000)]:
            h = 1 // (20*gamma)
            assert h >= 1/(40*gamma) and 0 < h*gamma <= F(1,20)
            for batch_length in {1, h}:
                p0 = F(3,4)*(1-(1-gamma)**batch_length)
                assert p0 <= F(3,80)
                assert (10*p0)**(2**level)/10 <= eta <= u/2
            for D in [0, 1, h, h+1, 1_000_000]:
                J = (D+h-1)//h
                assert J <= 40*gamma*D+1
                assert (J == 0) == (D == 0)
                if J:
                    assert 1 <= D-(J-1)*h <= h
                # The combined count includes both inner and outer allocations.
                # These illustrative outer constants test accounting, not code existence.
                N, CW, CT = 120, 1024, 2048
                width = (W+CW)*N
                slots = (42*node_count+CT)*N*J
                assert width >= W*N and width >= CW*N
                assert slots >= 42*node_count*N*J and slots >= CT*N*J
                parameter_cases += 1
    return {'exact_schedule_and_resource_cases': parameter_cases,
            'constant_depth_choices': choices}


def hybrid_sufficient_size_example():
    # Conditional algebraic fixture only: these are not certified Tanner local codes.
    u, gamma = F(1,32), F(1,20)
    D, epsilon, k, rate = 1_000_000, F(1,2**20), 128, F(1,4)
    level = 0
    while F(3,8)**(2**level)/10 > u/2:
        level += 1
    h = 1 // (20*gamma)
    J = (D+h-1)//h
    exponent = 0
    while F(2*J, 2**exponent) > epsilon:
        exponent += 1
    N = max(math.ceil(k/rate), math.ceil(F(2*exponent)/(u*u)))
    failure = F(2*J, 2**exponent)
    assert u*u*N/2 >= exponent and failure <= epsilon
    assert gamma > u/2  # This fixture lies outside the direct Tanner noise range.
    return {'interpretation': 'Conditional hybrid scalar fixture only; no Tanner block simulated.',
            'u': str(u), 'gamma': str(gamma), 'D': D, 'epoch_length': h,
            'epochs': J, 'logical_qubits': k, 'rate_lower_bound': str(rate),
            'constant_inner_levels': level, 'inner_width_per_outer_position': 5**level,
            'effective_nonidentity_bound': str(F(3,8)**(2**level)/10),
            'target_error': str(epsilon), 'sufficient_outer_data_length': N,
            'failure_bound_exact': str(failure),
            'certificate': 'Inner residual <=u/2; 2J exp(-u^2 N/2)<=2J*2^-exponent<=epsilon.'}


def run():
    return {'status': 'PASS', 'exhaustive_two_input_two_output_lookup_cases': exhaustive_lookup(),
            'reference_channel': choi_channel(), 'scalar_checks': scalar_and_resource_checks(), 'sufficient_size_example': sufficient_size_example(),
            'hybrid_inner': hybrid_inner_checks(), 'hybrid_schedule': hybrid_schedule_checks(),
            'hybrid_sufficient_size_example': hybrid_sufficient_size_example(),
            'scope': 'Finite compilation/Choi identities and exact scalar checks; inherited Tanner family and decoder theorems are not numerically re-proved.'}


if __name__=='__main__':
    print(json.dumps(run(), indent=2))
