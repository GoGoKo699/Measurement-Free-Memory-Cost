#!/usr/bin/env python3
"""Exact scalar inversions and rational-interval guessing-decay certificates.

All certificate decisions use Fraction arithmetic. Finite checks supplement the
all-parameter proofs in science/guessing_decay.md; they do not certify capacity.
No third-party dependencies. Prints deterministic JSON; writes no files.
"""
from fractions import Fraction as Q
from functools import lru_cache
import json
import hashlib
from math import comb, isqrt, factorial

BITS = 80
ZERO, ONE = Q(0), Q(1)


def outward(pair, bits=BITS):
    lo, hi = pair
    scale = 1 << bits
    low = (lo*scale).numerator // (lo*scale).denominator
    high = -((-hi*scale).numerator // (-hi*scale).denominator)
    return Q(low, scale), Q(high, scale)


def add(a, b): return a[0]+b[0], a[1]+b[1]
def neg(a): return -a[1], -a[0]
def sub(a, b): return add(a, neg(b))
def mul(a, b):
    vals = [x*y for x in a for y in b]
    return outward((min(vals), max(vals)))
def div(a, b):
    assert b[0] > 0
    return mul(a, (1/b[1], 1/b[0]))
def point(x): return Q(x), Q(x)


@lru_cache(None)
def log_interval(x):
    assert x > 0
    if x < 1:
        return neg(log_interval(1/x))
    count = 0
    while x >= 2:
        x /= 2
        count += 1
    def unit(y):
        t = (y-1)/(y+1)
        terms = 40
        low = 2*sum((t**(2*j+1)/Q(2*j+1) for j in range(terms)), ZERO)
        tail = 2*t**(2*terms+1)/(Q(2*terms+1)*(1-t*t))
        return low, low+tail
    low, high = unit(x)
    if count:
        l2, h2 = unit(Q(2))
        low += count*l2
        high += count*h2
    return outward((low, high))


def log_range(x):
    return log_interval(x[0])[0], log_interval(x[1])[1]


def alpha_interval(q, z):
    t = mul(point(2), log_range(z))
    a = log_range(add(point(1-q), mul(point(q), mul(z, z))))
    return div(a, t)


def derivative_interval(q, rho, z):
    # d/d(lambda) [lambda*(2-rho*alpha/(1-alpha))].
    z2 = mul(z, z)
    base = add(point(1-q), mul(point(q), z2))
    t = mul(point(2), log_range(z))
    a = log_range(base)
    h = div(mul(point(q), z2), base)
    numerator = sub(mul(mul(t, t), h), mul(a, a))
    gap = sub(t, a)
    assert gap[0] > 0
    hessian_factor = div(numerator, mul(gap, gap))
    return sub(point(2), mul(point(rho), hessian_factor))


def objective_interval(q, rho, z):
    alpha = alpha_interval(q, z)
    ratio = div(alpha, sub(point(1), alpha))
    return mul(log_range(z), sub(point(2), mul(point(rho), ratio)))


def seventh_root_two(bits=70):
    scale = 1 << bits
    lo, hi = scale, 2*scale
    while hi-lo > 1:
        middle = (lo+hi)//2
        if middle**7 <= 2*scale**7:
            lo = middle
        else:
            hi = middle
    assert Q(lo, scale)**7 <= 2 < Q(hi, scale)**7
    return Q(lo, scale), Q(hi, scale)


def small(pair, digits=12):
    lo, hi = pair
    scale = 10**digits
    low = (lo*scale).numerator//(lo*scale).denominator
    high = -((-hi*scale).numerator//(-hi*scale).denominator)
    return [str(Q(low, scale)), str(Q(high, scale))]


def certified_optimum(q, rho):
    critical = 2*(1-q)/q
    if rho >= critical:
        return {'rho': str(rho), 'case': 'zero', 'E_interval': ['0', '0']}
    if not rho:
        value = mul(point(Q(2,7)), log_interval(Q(2)))
        return {'rho': '0', 'case': 'right_endpoint', 'E_interval': small(value)}
    endpoint = seventh_root_two()
    endpoint_derivative = derivative_interval(q, rho, endpoint)
    if endpoint_derivative[0] >= 0:
        # Strict concavity proves the endpoint is the global maximizer.
        value = objective_interval(q, rho, endpoint)
        return {'rho': str(rho), 'case': 'right_endpoint',
                'E_interval': small(value), 'endpoint_derivative': small(endpoint_derivative)}
    assert endpoint_derivative[1] < 0
    lo, hi = Q(1), endpoint[1]
    for _ in range(65):
        middle = (lo+hi)/2
        derivative = derivative_interval(q, rho, point(middle))
        if derivative[0] > 0:
            lo = middle
        elif derivative[1] < 0:
            hi = middle
        else:
            break  # Previous endpoint signs still certify the root bracket.
    assert lo > 1
    assert derivative_interval(q, rho, point(lo))[0] > 0
    assert derivative_interval(q, rho, point(hi))[1] < 0
    value = objective_interval(q, rho, (lo, hi))
    assert value[0] > 0
    return {'rho': str(rho), 'case': 'unique_interior_maximum',
            'z_interval': small((lo, hi)), 'E_interval': small(value)}


def injection(D, width, floor, retention):
    return max([ZERO]+[floor/retention+(D-ell)*(1-retention)*floor/retention
                       -width*retention**(ell-1) for ell in range(1,D+1)])


def inverse(D, width, budget, retention):
    if D == 0:
        return width
    return min(retention*(budget+width*retention**(ell-1)) /
               (1+(D-ell)*(1-retention)) for ell in range(1,D+1))


def exact_inversion_checks():
    checks = 0
    for retention in (Q(1,3), Q(2,3), Q(99,100)):
        for D in range(1,9):
            for width in (Q(2), Q(7), Q(24)):
                for budget in (ZERO, Q(1,2), Q(2), Q(17)):
                    u = inverse(D,width,budget,retention)
                    assert injection(D,width,u,retention) == budget
                    for floor in (u/2, u, u+Q(1,100)):
                        assert (injection(D,width,floor,retention) <= budget) == (floor <= u)
                        checks += 1
    return checks


def finite_u(D, width, budget, alpha):
    # Monotone in retention. Round each power outward to prevent denominator
    # growth while preserving rigorous upper/lower bounds for all 2,000 cuts.
    lower, upper = width*alpha[0], width*alpha[1]
    power_low = power_high = ONE
    selected = None
    for ell in range(1,D+1):
        low = alpha[0]*(budget+width*power_low)/(1+(D-ell)*(1-alpha[0]))
        high = alpha[1]*(budget+width*power_high)/(1+(D-ell)*(1-alpha[1]))
        lower = min(lower, low)
        if high < upper:
            upper, selected = high, ell
        power_low = outward((power_low*alpha[0],power_low*alpha[0]))[0]
        power_high = outward((power_high*alpha[1],power_high*alpha[1]))[1]
    assert 0 <= lower <= upper
    return outward((lower,upper)), selected


def exp_upper(x):
    assert 0 <= x < 20
    total = term = ONE
    for j in range(1,101):
        term = term*x/j
        total += term
    next_term = term*x/101
    return outward((total,total+next_term/(1-x/102)))[1]


def root_upper(x, degree, bits=80):
    assert x >= 0
    if not x:
        return ZERO
    scale = 1 << bits
    if degree == 2:
        low = isqrt((x.numerator*scale**2)//x.denominator)
    else:
        low, high = 0, scale*max(1,-((-x.numerator)//x.denominator))
        while high-low > 1:
            middle = (low+high)//2
            if middle**degree*x.denominator <= x.numerator*scale**degree:
                low = middle
            else:
                high = middle
    answer = Q(low+1,scale)
    assert Q(low,scale)**degree <= x <= answer**degree
    return answer


def upper_hull(points):
    result = []
    for item in points:
        while len(result) >= 2:
            a,b = result[-2:]
            old = (b[1]-a[1])/(b[0]-a[0])
            new = (item[1]-b[1])/(item[0]-b[0])
            if old > new:
                break
            result.pop()
        result.append(item)
    for x,y in points:
        assert y <= polygon_value(result,x)
    for a,b in zip(result,result[1:]):
        assert b[1] >= a[1]
    return result


def polygon_value(vertices, y):
    if y >= vertices[-1][0]:
        return ONE
    for a,b in zip(vertices,vertices[1:]):
        if y <= b[0]:
            return min(ONE,a[1]+(y-a[0])*(b[1]-a[1])/(b[0]-a[0]))
    raise AssertionError('outside polygon')


def integer_profiles(k,n,z):
    d = 2**k
    beta = (2+root_upper(Q(2),2))/4
    pair_rows = []
    points = {name:[] for name in ('integer_direct','specified_recovery','combined')}
    for r in range(2*k+1):
        rank_caps = {name:ZERO for name in points}
        for s in range(r//2+1):
            c = r-2*s
            if s+c > n:
                continue
            dimension = min(ONE,Q(2)**(s+c-k))
            lifting = min(ONE,Q(2)**s*beta**k)
            direct = min(ONE,dimension,lifting,root_upper(Q(2)**(s-k),6))
            fidelity = min(ONE,Q(2)**(r-2*k),Q(2)**(s-k))
            values = []
            for ell in range(k+1):
                tau = Q(1,2**(ell+1)) if ell < k else ZERO
                coeff = sum((Q(comb(k,j))*(1-2**j*tau) for j in range(ell+1)),ZERO)
                values.append(tau+coeff*fidelity)
            u = min(values)
            assert u >= Q(1,d)
            renes = (1+root_upper((d-1)*(d*u-1),2))/d
            recovery = min(ONE,dimension,lifting,renes)
            if r == 0:
                # Rank zero is an input-independent channel, stronger than
                # any recovery relaxation at this endpoint.
                direct = recovery = Q(1,d)
            caps = {'integer_direct':direct,'specified_recovery':recovery,
                    'combined':min(direct,recovery)}
            pair_rows.append({'r':r,'s':s,'c':c,**{a:str(b) for a,b in caps.items()}})
            for name,cap in caps.items():
                rank_caps[name] = max(rank_caps[name],cap)
        for name in points:
            points[name].append((z**r,rank_caps[name]))
    assert len(pair_rows) == 169
    hulls = {name:upper_hull(items) for name,items in points.items()}
    digest = hashlib.sha256(json.dumps(pair_rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return hulls, {'pairs_checked':len(pair_rows),'pair_caps_sha256':digest,
                   'hull_vertex_ranks':{name:[next(r for r in range(2*k+1) if z**r == x)
                                             for x,_ in vertices] for name,vertices in hulls.items()},
                   'scope':'Outer caps for the specified published-recovery reductions and direct dimension theorem, with existing dimension and dimension-lifting caps. No exhaustive prior theorem or optimal operational score claim.'}


def finite_examples():
    k, n, D, q, z = 12, 120, 2000, Q(99,100), Q(27,25)
    assert z**7 <= 2
    alpha = alpha_interval(q,point(z))
    logz = log_interval(z)
    rows = []
    profiles, profile_receipt = integer_profiles(k,n,z)
    # The existing exact moment certificate proves this stronger polygon.
    vertices = [(ONE,Q(1,4096)),(z**11,Q(1,4)),(z**24,ONE)]
    slopes = [(b[1]-a[1])/(b[0]-a[0]) for a,b in zip(vertices,vertices[1:])]
    assert slopes[0] > slopes[1] > 0
    assert profiles['integer_direct'] == profiles['combined'] == vertices
    for budget in (0,171,200,243,300):
        u, ell = finite_u(D,Q(2*n),Q(budget),alpha)
        decay = max(ZERO,logz[0]*(2*k-u[1]))
        y_upper = exp_upper(u[1]*logz[1])
        polygon_upper = polygon_value(vertices,y_upper)
        profile_upper = {name:polygon_value(hull,y_upper) for name,hull in profiles.items()}
        simple_upper = ONE if decay == 0 else 1/(
            sum((decay**j/Q(factorial(j)) for j in range(101)),ZERO))
        # exp(-decay) <= reciprocal of a positive truncated exp series.
        assert polygon_upper <= simple_upper
        rows.append({'m_T':budget,'coasting_index':ell,'U_interval':small(u),
                     'simple_negative_log_score_lower':small((decay,decay)),
                     'simple_score_upper':small((simple_upper,simple_upper))[1],
                     'certified_polygon_score_upper':small((polygon_upper,polygon_upper))[1],
                     **{name+'_score_upper':small((value,value))[1] for name,value in profile_upper.items()}})
    return {'parameters':{'k':k,'n':n,'D':D,'q':str(q),'z':str(z)},'integer_profile_certificate':profile_receipt,'rows':rows}


def main():
    checks = exact_inversion_checks()
    q = Q(99,100)
    critical = 2*(1-q)/q
    optima = [certified_optimum(q,rho) for rho in
              (ZERO,Q(1,200),Q(1,100),Q(3,200),Q(1,50),critical)]
    print(json.dumps({'status':'pass','exact_inversion_checks':checks,
        'critical_rho':str(critical),'certified_asymptotic_optima':optima,
        'finite_examples':finite_examples(),
        'scope':'Exact scalar identities and rational-interval objective/finite certificates. All-parameter statements use the written proof. No achievable capacity or physical optimality is certified.'},indent=2))


if __name__ == '__main__':
    main()
