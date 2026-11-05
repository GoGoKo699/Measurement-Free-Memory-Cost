#!/usr/bin/env python3
"""Exact-rational outward certificates for the reliability ledger.

This evaluates necessary bounds only. The displayed physical parameters are not
claimed feasible, and numerical certificates do not establish the theorem.
"""
from fractions import Fraction as F
import json
from pathlib import Path


def ln_interval(x, terms=64):
    """Enclose ln(x) by rational atanh sums after binary range reduction."""
    x = F(x)
    if x <= 0:
        raise ValueError("logarithm requires a positive rational")
    exponent = 0
    while x >= 2:
        x /= 2
        exponent += 1
    while x < 1:
        x *= 2
        exponent -= 1

    def unit(y):
        t = (y-1)/(y+1)
        power = t
        total = F(0)
        for j in range(terms):
            total += 2*power/(2*j+1)
            power *= t*t
        remainder = 2*power/((2*terms+1)*(1-t*t))
        return total, total+remainder

    low, high = unit(x)
    two_low, two_high = unit(F(2))
    if exponent >= 0:
        return low+exponent*two_low, high+exponent*two_high
    return low+exponent*two_high, high+exponent*two_low


def quotient(numerator, denominator):
    nl, nh = numerator
    dl, dh = denominator
    assert nl >= 0 and dl > 0
    return nl/dh, nh/dl


def round_interval(interval, digits=40):
    """Shorten rational endpoints outward before subsequent series work."""
    scale = 10**digits
    low, high = interval
    return (F((low.numerator*scale)//low.denominator, scale),
            F(-((-high.numerator*scale)//high.denominator), scale))


def one_minus_root_interval(b, D, terms=24):
    """Enclose 1-b**(1/D) without cancellation or floating-point powers.

    For the certified examples x=-ln(b)/D is in (0,1]. The alternating
    expansion of 1-exp(-x) has decreasing terms there. An even partial sum
    is a lower bound; adding the next positive term is an upper bound.
    Monotonicity in x permits separate outward log endpoints.
    """
    assert 0 < b < 1 and D >= 1 and terms > 0 and terms % 2 == 0
    log_b = ln_interval(1/b)
    low_x, high_x = round_interval(tuple(x/D for x in log_b))
    assert 0 < low_x <= high_x <= 1, 'alternating-series certificate domain'

    def endpoint(x):
        term, total = x, F(0)
        for j in range(1, terms+1):
            total += term if j % 2 else -term
            term *= x/(j+1)
        return total, total+term

    low = endpoint(low_x)[0]
    high = endpoint(high_x)[1]
    result = round_interval((low, high))
    assert 0 < result[0] <= result[1] < 1
    return result


def outward_decimal(interval, digits=12):
    scale = 10**digits
    low, high = interval
    a = (low.numerator*scale)//low.denominator
    b = -((-high.numerator*scale)//high.denominator)
    def fixed(integer):
        sign = '-' if integer < 0 else ''
        integer = abs(integer)
        return f"{sign}{integer//scale}.{integer%scale:0{digits}d}"
    assert F(a, scale) <= low <= high <= F(b, scale)
    return {"lower":fixed(a), "upper":fixed(b), "decimal_places":digits,
            "outward_denominator":str(scale)}


def ceiling(x):
    return -((-x.numerator)//x.denominator)


def certify(D, n, k=1, gamma=F(1,100), epsilon=F(1,100)):
    q = 1-gamma
    delta = epsilon/(1-F(1,4**k))
    A = D*(1-delta)
    assert A/delta > 1
    log_ratio = ln_interval(A/delta)
    log_noise = ln_interval(1/gamma)
    reliability = quotient(tuple(gamma*A*x for x in log_ratio),log_noise)
    reliability = tuple(x-2*n for x in reliability)
    a = 2*k*(1-epsilon-F(1,4**k))/(1-F(1,4**k))
    b = 1-delta
    quantum_floor = k*max(F(0), (1-epsilon-F(1,2**k))/(1-F(1,2**k)))
    root_complement = one_minus_root_interval(b, D)
    product_log = (ln_interval(1/root_complement[1])[0],
                   ln_interval(1/root_complement[0])[1])
    product_floor = quotient(tuple(b*x for x in product_log), log_noise)
    rank_floor = tuple(max(a, quantum_floor+x) for x in product_floor)
    sharpened = tuple(gamma*D*x+a-2*n for x in rank_floor)
    assert ceiling(sharpened[0]) == ceiling(sharpened[1])
    assert sharpened[0] > reliability[1]
    kappa = tuple(-x for x in reversed(ln_interval(q)))
    reserve = ln_interval(2*n/a)
    numerator=(D*kappa[0]-reserve[1],D*kappa[1]-reserve[0])
    assert numerator[0] > 0
    logarithmic = tuple(a*x for x in numerator)
    integer_relaxed = quotient(numerator,ln_interval(1+1/a))
    assert ceiling(integer_relaxed[0]) == ceiling(integer_relaxed[1])
    assert ceiling(reliability[0]) == ceiling(reliability[1])

    # Exact old scalar ledger. Its greedy optimal schedule coasts until the last
    # integer t with M*q**t >= a. Binary search avoids summing every candidate.
    M=F(2*n)
    lo,hi=0,1
    while hi <= D and M*q**hi >= a:
        lo,hi=hi,min(D+1,2*hi)
    while lo+1 < hi:
        mid=(lo+hi)//2
        if M*q**mid >= a:
            lo=mid
        else:
            hi=mid
    t0=lo
    if t0 >= D:
        scalar=F(0)
    else:
        scalar=a/q-M*q**t0+(D-t0-1)*(1-q)*a/q

    return {
        "parameters":{"k":k,"n":n,"D":D,"gamma":str(gamma),
                      "epsilon":str(epsilon),"delta":str(delta),"A":str(A)},
        "reliability_lower":outward_decimal(reliability),
        "integer_gate_count":ceiling(reliability[0]),
        "sharpened_joint_lower":outward_decimal(sharpened),
        "sharpened_integer_gate_count":ceiling(sharpened[0]),
        "rank_floor":str(a),
        "symplectic_half_rank_floor":str(quantum_floor),
        "strictly_beats_rank_only_integer_count":(
            ceiling(sharpened[0]) > ceiling(reliability[1])),
        "old_payload_logarithmic_lower":outward_decimal(logarithmic),
        "old_individual_gate_relaxation":outward_decimal(integer_relaxed),
        "old_individual_integer_gate_count":ceiling(integer_relaxed[0]),
        "old_scalar_ledger_lower":outward_decimal((scalar,scalar)),
        "old_scalar_integer_gate_count":ceiling(scalar),
        "old_scalar_coasting_rounds":t0,
        "strictly_beats_both_old_integer_counts":(
            ceiling(reliability[0]) > max(ceiling(integer_relaxed[1]),ceiling(scalar))),
        "proof_of_decimal_enclosures":{
            "method":"64-term exact Fraction atanh series after binary range reduction; 24-term alternating series for 1-exp(-x), with 40-place outward rational rounding",
            "tail":"2*t^(2*N+1)/((2*N+1)*(1-t^2))",
            "one_minus_b_root":outward_decimal(root_complement, digits=30),
            "log_inverse_root_complement":outward_decimal(product_log),
            "exact_product_rank_floor":outward_decimal(product_floor),
            "log_A_over_delta":outward_decimal(log_ratio),
            "log_inverse_gamma":outward_decimal(log_noise),
            "log_inverse_q":outward_decimal(kappa),
            "log_width_reserve":outward_decimal(reserve),
            "log_one_plus_inverse_a":outward_decimal(ln_interval(1+1/a))},
        "scope":"Necessary-bound comparison only. No claim of feasible storage parameters, optimal T coefficient, or theorem proof."
    }


def main():
    examples=[certify(100000,200),certify(1000000,200)]
    assert all(x['strictly_beats_both_old_integer_counts'] for x in examples)
    assert all(x['strictly_beats_rank_only_integer_count'] for x in examples)
    result={"status":"pass","arithmetic":"Exact Fraction throughout; decimal output rounded outward",
            "examples":examples,
            "limits":"Finite arithmetic certificates supplement, rather than replace, the written extinction theorem and physical resource accounting."}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
