#!/usr/bin/env python3
"""Finite exact checks for the absorbing-envelope reliability ledger.

All binary algebra and noise probabilities are exact. Natural logarithms are
enclosed by rational atanh series with an explicit remainder. These finite
checks supplement the written dimension-independent argument; they neither
certify a physical memory implementation nor establish an optimal coefficient.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import random

from certify_reliability_examples import one_minus_root_interval


def require(condition, message):
    if not condition:
        raise ValueError(message)


def span(generators):
    """Return every vector, independently of a generator-matrix representation."""
    values = {0}
    for v in generators:
        values |= {x ^ v for x in values}
    return tuple(sorted(values))


def rank(space):
    require(space and space[0] == 0, 'invalid enumerated linear space')
    dimension = len(space).bit_length() - 1
    require(len(space) == 1 << dimension, 'space cardinality is not a power of two')
    return dimension


def all_spaces(bits):
    spaces = {(0,)}
    frontier = [(0,)]
    while frontier:
        current = frontier.pop()
        for v in range(1 << bits):
            if v not in current:
                enlarged = span((*current, v))
                if enlarged not in spaces:
                    spaces.add(enlarged)
                    frontier.append(enlarged)
    return sorted(spaces, key=lambda s: (len(s), s))


def coordinate_mask(n, qubits):
    return qubits | (qubits << n)


def erased(space, n, qubits):
    mask = coordinate_mask(n, qubits)
    return tuple(v for v in space if not (v & mask))


def reset(space, n, qubits):
    mask = coordinate_mask(n, qubits)
    return tuple(sorted({v & ~mask for v in space}))


def symplectic(v, w, n):
    mask = (1 << n) - 1
    return (((v & mask) & (w >> n)).bit_count()
            + ((v >> n) & (w & mask)).bit_count()) % 2


@lru_cache(maxsize=None)
def symplectic_half_rank(space, n):
    radical = tuple(v for v in space
                    if all(symplectic(v, w, n) == 0 for w in space))
    form_rank = rank(space)-rank(radical)
    require(form_rank % 2 == 0, 'alternating form has odd rank')
    return form_rank//2


def rotation(space, n, axis):
    # The refined update is essential: adding axis unconditionally would
    # resurrect the zero envelope even though the Heisenberg channel is unital.
    if all(symplectic(v, axis, n) == 0 for v in space):
        return space
    return span((*space, axis))


def hadamard(v, n, j):
    if ((v >> j) ^ (v >> (n+j))) & 1:
        v ^= (1 << j) | (1 << (n+j))
    return v


def phase(v, n, j):
    return v ^ (((v >> j) & 1) << (n+j))


def cnot(v, n, control, target):
    return (v ^ (((v >> control) & 1) << target)
            ^ (((v >> (n+target)) & 1) << (n+control)))


def clifford(space, n, op):
    name, *positions = op
    fn = {'H': hadamard, 'S': phase, 'CX': cnot}[name]
    return tuple(sorted(fn(v, n, *positions) for v in space))


def measurement(space, n, axis):
    return tuple(v for v in space if not symplectic(v, axis, n))


def pivot_qubits(space, n):
    """Find an injective physical-coordinate restriction, by subset search."""
    r = rank(space)
    for number in range(min(r, n)+1):
        for chosen in combinations(range(n), number):
            qmask = sum(1 << j for j in chosen)
            mask = coordinate_mask(n, qmask)
            if len({v & mask for v in space}) == len(space):
                return qmask
    raise ValueError('no at-most-r pivot-qubit set')


def mask_probability(mask, n, gamma):
    weight = mask.bit_count()
    return gamma**weight * (1-gamma)**(n-weight)


@lru_cache(maxsize=None)
def ln_interval(x, terms=72):
    """Exact rational lower/upper bounds, including x<1 and negative logs."""
    x = F(x)
    require(x > 0, 'logarithm requires positive input')
    exponent = 0
    while x >= 2:
        x /= 2
        exponent += 1
    while x < 1:
        x *= 2
        exponent -= 1

    def unit(y):
        t = (y-1)/(y+1)
        power, total = t, F(0)
        for j in range(terms):
            total += 2*power/(2*j+1)
            power *= t*t
        tail = 2*power/((2*terms+1)*(1-t*t))
        return total, total+tail

    low, high = unit(x)
    two_low, two_high = unit(F(2))
    if exponent >= 0:
        return low+exponent*two_low, high+exponent*two_high
    return low+exponent*two_high, high+exponent*two_low


def positive_log_lower(x):
    return max(F(0), ln_interval(x)[0])


def reliability_lower(D, M, gamma, delta):
    require(D >= 1 and 0 < gamma < 1 and 0 < delta < 1, 'bound domain')
    A = D*(1-delta)
    return gamma*A*positive_log_lower(A/delta)/ln_interval(1/gamma)[1]-M


def check_paired_pivots():
    """Independently search physical restrictions, without using the proof."""
    counts = defaultdict(int)
    families = {n: all_spaces(2*n) for n in (1, 2, 3)}
    for n, expected in ((1, 5), (2, 67), (3, 2825)):
        require(len(families[n]) == expected, 'incomplete exhaustive pivot family')
    rng = random.Random(6112026)
    samples = {(0,), tuple(range(256))}
    while len(samples) < 300:
        samples.add(span(rng.randrange(256) for _ in range(rng.randrange(9))))
    families[4] = sorted(samples, key=lambda V: (len(V), V))
    for n, spaces in families.items():
        for V in spaces:
            r, s = rank(V), symplectic_half_rank(V, n)
            require(pivot_qubits(V, n).bit_count() <= r-s,
                    'paired-pivot support exceeds r-s')
            # The proof first finds s full X/Z pairs carrying 2s independent
            # columns. Test the existence of such a physical projection too.
            require(any(len({v & coordinate_mask(n, sum(1 << j for j in chosen))
                             for v in V}) >= 1 << (2*s)
                        for chosen in combinations(range(n), s)),
                    'no independent collection of s physical coordinate pairs')
            counts[f'{n}_qubit_subspaces'] += 1
    return dict(sorted(counts.items()))


def check_subspaces():
    counts = defaultdict(int)
    expected = {1: 5, 2: 67}
    gammas = (F(1,100), F(1,20), F(1,3), F(2,3))
    rng = random.Random(6192026)
    families = {n: all_spaces(2*n) for n in (1, 2)}
    for n in (1, 2):
        require(len(families[n]) == expected[n], 'incomplete exhaustive subspace list')
    # The three-qubit samples include all rank classes, but are not exhaustive.
    sampled = {(0,), tuple(range(64))}
    for _ in range(180):
        sampled.add(span(rng.randrange(64) for __ in range(rng.randrange(7))))
    families[3] = sorted(sampled, key=lambda s: (len(s), s))

    for n, spaces in families.items():
        for V in spaces:
            r = rank(V)
            s = symplectic_half_rank(V, n)
            require(all((a ^ b) in V for a in V for b in V), 'not a subspace')
            pivots = pivot_qubits(V, n)
            require(pivots.bit_count() <= r-s, 'paired-pivot cardinality')
            mask = coordinate_mask(n, pivots)
            require(len({v & mask for v in V}) == len(V), 'pivot restriction not injective')
            for E in range(1 << n):
                survivor = erased(V, n, E)
                if E & pivots == pivots:
                    require(survivor == (0,), 'pivot erasure did not annihilate space')
                counts['erasure_masks'] += 1
            for gamma in gammas:
                expected_rank, extinction = F(0), F(0)
                for E in range(1 << n):
                    probability = mask_probability(E, n, gamma)
                    rr = rank(erased(V, n, E))
                    expected_rank += probability*rr
                    extinction += probability*(rr == 0)
                require(expected_rank <= (1-gamma)*r, 'mean contraction failed')
                if r:
                    require(extinction >= gamma**r, 'rank-dependent extinction hazard failed')
                    require(extinction >= gamma**(r-s), 'paired-pivot extinction hazard failed')
                    counts['positive_rank_hazards'] += 1
                counts['contraction_cases'] += 1

            # Every Pauli-axis rotation has the exact commuting support or two
            # anticommuting support labels; check enclosure directly labelwise.
            for axis in range(1, 1 << (2*n)):
                out = rotation(V, n, axis)
                require(rank(out) <= r+1, 'rotation rank charge failed')
                for v in V:
                    required = (v,) if not symplectic(v, axis, n) else (v, v ^ axis)
                    require(all(x in out for x in required), 'rotation support not enclosed')
                if not r:
                    require(out == (0,), 'zero envelope resurrected under rotation')
                counts['rotation_enclosures'] += 1
            for A in range(1 << n):
                out = reset(V, n, A)
                require(rank(out) <= r, 'reset increased rank')
                if not r:
                    require(out == (0,), 'zero envelope resurrected under reset')
                counts['reset_cases'] += 1
            counts[f'{n}_qubit_subspaces'] += 1
    return dict(sorted(counts.items()))


def apply_free(dist, n, op):
    out = defaultdict(F)
    for V, probability in dist.items():
        if op[0] in ('H', 'S', 'CX'):
            VV = clifford(V, n, op)
            require(rank(VV) == rank(V), 'Clifford did not preserve rank')
        elif op[0] == 'R':
            VV = reset(V, n, op[1])
        elif op[0] == 'M':
            VV = measurement(V, n, op[1])
        else:
            raise ValueError('unknown free operation')
        require(rank(VV) <= rank(V), 'free primitive increased rank')
        require(V != (0,) or VV == (0,), 'free primitive resurrected zero')
        out[VV] += probability
    return dict(out)


def rank_expectation(dist):
    return sum((p*rank(V) for V, p in dist.items()), F(0))


def check_schedule(n, gamma, start, schedule):
    dist = {start: F(1)}
    M = rank(start)
    D, T = 0, 0
    A, B, H = F(0), F(0), F(0)
    U, Q, paired_H = F(0), F(0), F(0)
    survival_product = F(1)
    conditional_jensen_cases = 0
    layer_records = []
    extinction_in_noise, extinction_in_free = F(0), F(0)
    total_noise_drop, total_gate_gain, total_free_drop = F(0), F(0), F(0)
    positive_ranks = set()
    for op in schedule:
        old_rank = rank_expectation(dist)
        old_zero = dist.get((0,), F(0))
        if op[0] == 'N':
            D += 1
            alive = 1-old_zero
            A += alive
            B += old_rank
            H += sum((p*gamma**rank(V) for V,p in dist.items() if rank(V)), F(0))
            quantum = sum((p*symplectic_half_rank(V, n) for V,p in dist.items()), F(0))
            uu = old_rank-quantum
            hazard = sum((p*gamma**(rank(V)-symplectic_half_rank(V, n))
                          for V,p in dist.items() if rank(V)), F(0))
            U += uu
            Q += quantum
            paired_H += hazard
            if alive:
                require(0 < hazard < alive, 'conditional hazard outside (0,1)')
                survival_product *= 1-hazard/alive
                positive_u = {rank(V)-symplectic_half_rank(V, n)
                              for V,p in dist.items() if p and rank(V)}
                if len(positive_u) == 1:
                    u = next(iter(positive_u))
                    require(uu == alive*u and hazard == alive*gamma**u,
                            'constant-u conditional Jensen identity')
                else:
                    require(uu*ln_interval(1/gamma)[0]
                            >= alive*ln_interval(alive/hazard)[1],
                            'conditional paired-hazard Jensen certificate')
                conditional_jensen_cases += 1
            positive_ranks.update(rank(V) for V,p in dist.items() if p and rank(V))
            out = defaultdict(F)
            for V, probability in dist.items():
                for E in range(1 << n):
                    VV = erased(V, n, E)
                    out[VV] += probability*mask_probability(E, n, gamma)
            dist = dict(out)
            new_rank = rank_expectation(dist)
            require(new_rank <= (1-gamma)*old_rank, 'schedule contraction')
            total_noise_drop += old_rank-new_rank
            extinction_in_noise += dist.get((0,), F(0))-old_zero
            layer_records.append({
                'enter_rank':old_rank, 'exit_rank':new_rank,
                'enter_survival':alive, 'exit_survival':1-dist.get((0,), F(0)),
                'T_before':T, 'hazard':hazard,
                'conditional_hazard':hazard/alive if alive else F(0)})
        elif op[0] == 'T':
            T += 1
            out = defaultdict(F)
            for V, probability in dist.items():
                VV = rotation(V, n, op[1])
                out[VV] += probability
            dist = dict(out)
            gain = rank_expectation(dist)-old_rank
            require(0 <= gain <= 1, 'schedule T charge')
            total_gate_gain += gain
        else:
            dist = apply_free(dist, n, op)
            total_free_drop += old_rank-rank_expectation(dist)
            extinction_in_free += dist.get((0,), F(0))-old_zero
        require(sum(dist.values(), F(0)) == 1, 'noise probabilities do not normalize')
        require(dist.get((0,), F(0)) >= old_zero, 'zero envelope not absorbing')
        require(1-dist.get((0,), F(0)) <= survival_product,
                'survival exceeds product of conditional hazard upper factors')

    final_rank = rank_expectation(dist)
    delta = dist.get((0,), F(0))
    require(extinction_in_noise+extinction_in_free == delta-(M == 0), 'first-extinction disjointness')
    require(H <= extinction_in_noise <= delta, 'sum of conditional hazards exceeds extinction')
    require(H <= paired_H <= extinction_in_noise, 'paired hazard sum exceeds extinction')
    require(B == U+Q, 'rank did not split into r-s and s cut times')
    require(survival_product >= 1-delta, 'exact survival product too small')
    require(A >= D*(1-delta), 'alive cut-time lower bound')
    require(M+total_gate_gain-total_free_drop-total_noise_drop == final_rank, 'rank ledger telescope')
    require(total_noise_drop >= gamma*B, 'rank contraction sum')
    require(T >= gamma*B-M+final_rank, 'T-count telescope')

    # These are actual interior cuts of the global sampled process: the
    # entering envelope distribution is retained, not reinitialized.
    intervals = ([(first, last) for first in range(D) for last in range(first, D)]
                 if D <= 4 else [(0, 0), (D-1, D-1), (0, D-1), (D//4, 3*D//4-1)])
    for first, last in intervals:
        layers = layer_records[first:last+1]
        product = F(1)
        for layer in layers:
            product *= 1-layer['conditional_hazard']
        require(product >= layers[-1]['exit_survival'], 'interior survival product')
        require(sum((x['hazard'] for x in layers), F(0))
                <= layers[0]['enter_survival']-layers[-1]['exit_survival'],
                'interior paired-hazard budget')
        interval_T = layers[-1]['T_before']-layers[0]['T_before']
        interval_B = sum((x['enter_rank'] for x in layers), F(0))
        require(interval_T >= gamma*interval_B+layers[-1]['exit_rank']-layers[0]['enter_rank'],
                'interior retained-exit rank ledger')

    jensen = 'not_applicable'
    if A and H:
        if len(positive_ranks) == 1:
            rr = next(iter(positive_ranks))
            require(B == A*rr and H == A*gamma**rr, 'constant-rank Jensen identity')
            jensen = 'exact_equality'
        else:
            log_gamma = ln_interval(1/gamma)
            log_ratio = ln_interval(A/H)
            require(B*log_gamma[0] >= A*log_ratio[1], 'weighted Jensen rational certificate')
            jensen = 'strict_rational_certificate'
    if D and 0 < delta < 1:
        lower = reliability_lower(D, M, gamma, delta)
        require(T >= lower, 'finite reliability bound fails on exact envelope schedule')
        b = 1-delta
        root_complement = one_minus_root_interval(b, D)
        product_lower = D*b*ln_interval(1/root_complement[1])[0]/ln_interval(1/gamma)[1]
        require(U >= product_lower, 'exact-product cumulative r-s lower bound failed')
        require(T >= gamma*(Q+product_lower)+final_rank-M,
                'symplectic product and retained-exit ledger failed')
    else:
        lower = None
        product_lower = None
    return {
        'n':n, 'gamma':str(gamma), 'D':D, 'prescribed_T_slots':T,
        'initial_rank':M, 'final_rank_expectation':str(final_rank),
        'extinction_probability':str(delta), 'alive_cut_time':str(A),
        'rank_cut_time':str(B), 'rank_hazard_sum':str(H),
        'symplectic_half_rank_cut_time':str(Q),
        'rank_minus_symplectic_cut_time':str(U),
        'paired_pivot_hazard_sum':str(paired_H),
        'survival_product_checked':True,
        'conditional_paired_Jensen_cases':conditional_jensen_cases,
        'exact_product_rank_bound_checked':product_lower is not None,
        'interval_product_and_exit_cases':len(intervals),
        'noise_first_extinction':str(extinction_in_noise),
        'free_first_extinction':str(extinction_in_free),
        'weighted_Jensen':jensen,
        'finite_bound_checked':lower is not None,
        'finite_bound_strictly_positive':lower is not None and lower > 0,
        'interpretation':'Exact sampled-envelope process, not a certified successful physical storage channel.'
    }


def check_schedules():
    specs = [
        (1, tuple(range(4)), [('N',), ('H',0), ('T',2), ('N',), ('S',0), ('N',)]),
        (1, span([1]), [('T',2), ('N',), ('H',0), ('T',2), ('N',), ('R',1), ('N',)]),
        (2, tuple(range(16)), [('N',), ('CX',0,1), ('T',4), ('N',), ('R',1), ('H',1), ('T',8), ('N',), ('M',3), ('N',)]),
        (2, span([1,8]), [('H',0), ('N',), ('T',4), ('CX',1,0), ('N',), ('M',5), ('T',8), ('N',)]),
        (3, tuple(range(64)), [('N',), ('CX',0,2), ('T',8), ('N',), ('R',2), ('H',1), ('T',16), ('N',), ('CX',2,1), ('N',)]),
        (3, span([1,10,48]), [('T',8), ('N',), ('CX',1,2), ('H',0), ('N',), ('T',16), ('R',4), ('N',), ('M',5), ('N',)]),
    ]
    # A Pauli-axis pi/8 rotation is a single T conjugated by free Clifford
    # operations. A fixed sweep of all axes replenishes a nonzero envelope.
    # This long fixture exercises a strictly positive reliability lower bound;
    # it is still only an envelope process, not an encoded-memory construction.
    sweep = [('N',), *(('T',axis) for axis in range(1,16))]
    specs.append((2, tuple(range(16)), sweep*100))
    result = []
    for n, start, schedule in specs:
        for gamma in (F(1,100), F(1,20), F(1,3)):
            result.append(check_schedule(n, gamma, start, schedule))
    return result


def check_positive_part():
    count = 0
    for D in range(1,31):
        for delta in (F(1,100), F(1,10), F(1,2), F(3,4), F(99,100)):
            A0 = D*(1-delta)
            for j in range(21):
                A = A0 + (D-A0)*j/20
                if A0 <= delta:
                    require(positive_log_lower(A0/delta) == 0, 'positive-part zero regime')
                else:
                    lo_A, hi_A = ln_interval(A/delta)
                    lo_0, hi_0 = ln_interval(A0/delta)
                    if A == A0:
                        require(lo_A == lo_0 and hi_A == hi_0, 'same logarithm mismatch')
                    else:
                        require(A*lo_A >= A0*hi_0, 'positive-regime monotonicity')
                count += 1
    return count


def check_schmidt_secants():
    count = 0
    for k in range(1, 17):
        baseline = F(1, 2**k)
        for s in range(2*k+1):
            branch_cap = min(F(1), F(2)**(s-k))
            secant = baseline+(1-baseline)*F(s, k)
            require(branch_cap <= secant, 'Bell Schmidt-number secant failed')
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Optional copy of the full JSON result; stdout is authoritative.')
    args = parser.parse_args()
    paired_pivots = check_paired_pivots()
    subspaces = check_subspaces()
    schedules = check_schedules()
    monotonicity = check_positive_part()
    positive_cases = sum(x['finite_bound_strictly_positive'] for x in schedules)
    require(positive_cases > 0, 'all finite schedule inequalities were vacuous')
    result = {
        'schema_version':2, 'status':'pass',
        'arithmetic':'Exact binary spaces and Fraction probabilities; exact rational logarithm enclosures.',
        'subspaces':subspaces, 'paired_pivot_subspaces':paired_pivots,
        'Bell_Schmidt_secant_cases':check_schmidt_secants(),
        'exact_schedules':len(schedules),
        'strictly_positive_schedule_bounds':positive_cases,
        'schedule_results':schedules, 'positive_part_cases':monotonicity,
        'limits':'Finite checks supplement the written proof. The schedules verify envelope dynamics, first-extinction accounting, survival products and the exit ledger, not Bell fidelity, physical feasibility, capacity or an optimal T coefficient. Paired pivots are exhaustive through three qubits and sampled on four; full noise/gate dynamics are exhaustive through two qubits and sampled on three.'
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
