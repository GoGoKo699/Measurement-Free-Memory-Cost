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


def power_interval(x, exponent):
    assert exponent >= 0
    result, base = point(1), point(x)
    while exponent:
        if exponent & 1:
            result = mul(result,base)
        base = mul(base,base)
        exponent //= 2
    return result


@lru_cache(None)
def adaptive_log_interval(x):
    """Atanh-series enclosure with an explicit adaptive rational tail."""
    assert x >= 1
    original = x
    count = 0
    while x >= 2:
        x /= 2
        count += 1
    t = (x-1)/(x+1)
    power, total, terms = t, ZERO, 0
    while True:
        total += 2*power/Q(2*terms+1)
        terms += 1
        power *= t*t
        tail = 2*power/(Q(2*terms+1)*(1-t*t))
        if tail <= Q(1,1 << 105):
            break
    result = total, total+tail
    if count:
        ln2 = log_interval(Q(2))
        result = result[0]+count*ln2[0], result[1]+count*ln2[1]
    assert original >= 1
    return outward(result)


def adaptive_log_range(x):
    return adaptive_log_interval(x[0])[0], adaptive_log_interval(x[1])[1]


def pgf_series(q,z,start=1,derivative=False):
    """Enclose sum_{j=start}^infty phi_j or its lambda derivative.

    Direct head terms end once b=(z^2-1)*q^J <=1/4. The remaining
    logarithmic/geometric series alternates with strictly decreasing terms.
    Even truncation is a lower sum and the next odd truncation an upper sum.
    Every operation, including q^J, is enclosed outward with Fractions.
    """
    assert 0 < q < 1 and z[0] >= 1 and start >= 1
    z2 = mul(z,z)
    a = sub(z2,point(1))
    qpower = power_interval(q,start)
    result = point(0)
    head = 0
    while a[1]*qpower[1] > Q(1,4):
        product = mul(a,qpower)
        if derivative:
            term = div(mul(z2,qpower),add(point(1),product))
        else:
            term = mul(point(Q(1,2)),adaptive_log_range(add(point(1),product)))
        result = add(result,term)
        qpower = mul(qpower,point(q))
        head += 1
    b = mul(a,qpower)
    assert b[1] <= Q(1,4)
    tail = point(0)
    qsmall = ONE
    bpower = point(1)
    prefactor = mul(z2,qpower)
    for index in range(1,62):
        qsmall *= q
        if derivative:
            term = div(mul(prefactor,bpower),point(1-qsmall))
        else:
            bpower = mul(bpower,b)
            term = div(bpower,point(2*index*(1-qsmall)))
        if index == 61:
            # After 60 terms the exact remainder is nonnegative and no larger
            # than the next positive term.
            tail = tail[0],tail[1]+term[1]
            break
        tail = add(tail,term if index & 1 else neg(term))
        if derivative:
            bpower = mul(bpower,b)
    result = add(result,tail)
    return (max(ZERO,result[0]),result[1]), {'head_terms':head,'alternating_tail_terms':60}


def full_intercept(z):
    lnz = log_range(z)
    ln2 = log_interval(Q(2))
    candidates = [mul(point(2),lnz),
                  mul(point(Q(1,6)),add(ln2,mul(point(5),lnz))),ln2]
    return min(item[0] for item in candidates),min(item[1] for item in candidates)


def full_objective(q,rho,z):
    series,receipt = pgf_series(q,z)
    return sub(full_intercept(z),mul(point(rho),series)),receipt


def certified_full_optimum(q,rho):
    critical = 2*(1-q)/q
    if rho >= critical:
        return {'rho':str(rho),'case':'zero','E_interval':['0','0']}
    if not rho:
        return {'rho':'0','case':'right_endpoint_z2','E_interval':small(log_interval(Q(2)))}
    kink = seventh_root_two()
    derivative_kink,_ = pgf_series(q,kink,derivative=True)
    derivative_end,_ = pgf_series(q,point(2),derivative=True)
    kink_left = sub(point(2),mul(point(rho),derivative_kink))
    kink_right = sub(point(Q(5,6)),mul(point(rho),derivative_kink))
    endpoint_left = sub(point(Q(5,6)),mul(point(rho),derivative_end))
    if endpoint_left[0] >= 0:
        z,case = point(2),'right_endpoint_z2'
    elif kink_left[0] >= 0 and kink_right[1] <= 0:
        z,case = kink,'profile_kink_z_2power1over7'
    else:
        if kink_left[1] < 0:
            lo,hi,slope = ONE,kink[1],Q(2)
            case = 'unique_root_below_profile_kink'
        else:
            assert kink_right[0] > 0 and endpoint_left[1] < 0
            lo,hi,slope = kink[0],Q(2),Q(5,6)
            case = 'unique_root_above_profile_kink'
        for _ in range(50):
            middle = (lo+hi)/2
            derivative,_ = pgf_series(q,point(middle),derivative=True)
            sign = sub(point(slope),mul(point(rho),derivative))
            if sign[0] > 0:
                lo = middle
            elif sign[1] < 0:
                hi = middle
            else:
                break
        left,_ = pgf_series(q,point(lo),derivative=True)
        right,_ = pgf_series(q,point(hi),derivative=True)
        assert sub(point(slope),mul(point(rho),left))[0] > 0
        assert sub(point(slope),mul(point(rho),right))[1] < 0
        z = lo,hi
    value,receipt = full_objective(q,rho,z)
    assert value[0] > 0
    return {'rho':str(rho),'case':case,'z_interval':small(z),
            'E_interval':small(value),'series_certificate':receipt}


def pgf_schedule_checks():
    cases = 0
    # A gate in backward gap h has age d-h at cut d, when d>h. Verify every
    # column and every selected terminal window, not one particular schedule.
    for D in range(1,25):
        weights = [Q(3,4)**j for j in range(D+1)]
        for L in range(1,D+1):
            T = D-L+1
            bound = sum(weights[1:T+1],ZERO)
            for h in range(D+1):
                column = sum((weights[d-h] for d in range(L,D+1) if d>h),ZERO)
                assert column <= bound
                cases += 1
    identities = 0
    for q in (Q(1,3),Q(2,3),Q(99,100)):
        for z in (Q(27,25),Q(3,2),Q(2)):
            value = z*z
            for j in range(1,25):
                value = 1-q+q*value
                assert value == 1+q**j*(z*z-1)
                identities += 1
    profiles = 0
    for k in range(1,25):
        for r in range(2*k+3):
            # Raising to42 removes the sixth/seventh roots in the inherited
            # profile and each affine support line.
            if 6*r <= 5*k:
                cap42 = Q(2)**(42*(r-k))
            elif r <= 2*k:
                cap42 = Q(2)**(6*(r-2*k))
            else:
                cap42 = ONE
            for z in (ONE,Q(27,25),Q(3,2),Q(2),Q(3)):
                support42 = max(z**(42*(r-2*k)),
                                Q(2)**(-7*k)*z**(42*r-35*k),
                                Q(2)**(-42*k)*z**(42*r))
                assert cap42 <= support42
                profiles += 1
    return {'all_gate_column_windows_checked':cases,
            'semigroup_squared_identities':identities,'full_profile_support_checks':profiles}


def full_finite_examples():
    k,n,D,q = 12,120,2000,Q(99,100)
    rows = []
    for z in (Q(27,25),Q(2)):
        prefix = [point(0)]
        qpower = point(1)
        a = z*z-1
        for _ in range(D):
            qpower = mul(qpower,point(q))
            term = mul(point(Q(1,2)),adaptive_log_range(add(point(1),mul(point(a),qpower))))
            prefix.append(add(prefix[-1],term))
        intercept = full_intercept(point(z))
        polygon = [(ONE,Q(1,4096)),(z**11,Q(1,4)),(z**24,ONE)] if z==Q(27,25) else None
        for budget in (0,171,200,243,300):
            best = None
            lower = None
            selected = None
            for L in range(1,D+1):
                T = D-L+1
                width = sub(prefix[D],prefix[L-1])
                moment_bound = mul(point(Q(1,T)),add(mul(point(2*n),width),mul(point(budget),prefix[T])))
                lower = moment_bound[0] if lower is None else min(lower,moment_bound[0])
                if best is None or moment_bound[1] < best:
                    best,selected = moment_bound[1],L
            decay = max(ZERO,k*intercept[0]-best)
            score_upper = ONE if not decay else 1/sum((decay**j/Q(factorial(j)) for j in range(101)),ZERO)
            item = {'z':str(z),'m_T':budget,'selected_L':selected,
                    'log_moment_bound_interval':small((lower,best)),
                    'negative_log_score_lower':small((decay,decay)),
                    'score_upper':small((score_upper,score_upper))[1]}
            if polygon is not None:
                polygon_upper = polygon_value(polygon,exp_upper(best))
                item['existing_polygon_score_upper'] = small((polygon_upper,polygon_upper))[1]
            rows.append(item)
    return {'parameters':{'k':k,'n':n,'D':D,'q':str(q)},'rows':rows}


def full_generating_function_certificate():
    q = Q(99,100)
    rates = (ZERO,Q(1,1000),Q(1,200),Q(1,100),Q(3,200),Q(1,50),2*(1-q)/q)
    optima = [certified_full_optimum(q,rho) for rho in rates]
    comparisons = []
    for item,rho in zip(optima,rates):
        old = certified_optimum(q,rho)
        new_lower = Q(item['E_interval'][0])
        old_upper = Q(old['E_interval'][1])
        if rho < 2*(1-q)/q:
            assert new_lower > old_upper
        comparisons.append({'rho':str(rho),'old_E_interval':old['E_interval'],
                            'full_E_interval':item['E_interval']})
    return {'status':'pass','schedule_and_profile_checks':pgf_schedule_checks(),
            'certified_optima':optima,'strict_improvement_over_previous_bound':comparisons,
            'finite_examples':full_finite_examples(),
            'series_scope':'Rational outward head terms plus 60-term alternating tail after (z^2-1)q^J<=1/4; derivative tails certified likewise. No floating-point certificate decisions.'}

def main():
    checks = exact_inversion_checks()
    q = Q(99,100)
    critical = 2*(1-q)/q
    optima = [certified_optimum(q,rho) for rho in
              (ZERO,Q(1,200),Q(1,100),Q(3,200),Q(1,50),critical)]
    print(json.dumps({'status':'pass','exact_inversion_checks':checks,
        'critical_rho':str(critical),'certified_asymptotic_optima':optima,
        'finite_examples':finite_examples(),
        'full_generating_function':full_generating_function_certificate(),
        'scope':'Exact scalar identities and rational-interval objective/finite certificates. All-parameter statements use the written proof. No achievable capacity or physical optimality is certified.'},indent=2))


if __name__ == '__main__':
    main()
