#!/usr/bin/env python3
"""Finite checks and exact certificates for moment processing bounds.

No third-party dependencies. All certificate decisions use integers/Fractions.
The complete continuum noise theorem is proved in the scientific account;
finite subspace checks supplement it. No physical achievability is inferred.
Run this file directly; it prints one deterministic JSON object and writes none.
"""
from fractions import Fraction as Q
from math import comb, isqrt
from collections import deque
from functools import lru_cache
import hashlib
import json

F = Q

BITS=80

def log_unit(x,N=40):
    assert 1<=x<=2
    t=(x-1)/(x+1)
    s=2*sum((t**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0))
    e=2*t**(2*N+1)/(Q(2*N+1)*(1-t*t))
    return s,s+e

def log_interval(x):
    assert x>0
    if x<1:
        a,b=log_interval(1/x);return -b,-a
    m=0
    while x>=2:x/=2;m+=1
    a,b=log_unit(x);c,d=log_unit(Q(2))
    return a+m*c,b+m*d

def outward(interval,bits=BITS):
    lo,hi=interval;scale=1<<bits
    l=(lo*scale).numerator//(lo*scale).denominator
    h=-((-hi*scale).numerator//(-hi*scale).denominator)
    result=Q(l,scale),Q(h,scale)
    assert result[0]<=lo<=hi<=result[1]
    return result

def ratio(x,y):
    assert x[0]>=0 and y[0]>0
    return outward((x[0]/y[1],x[1]/y[0]))

def fidelity_lower(k,p):
    d=2**k;g=p*p+(1-p)**2/Q(d-1)
    candidates=[('measure-and-prepare',p/d)]
    for ell in range(k+1):
        tau=Q(1,2**(ell+1)) if ell<k else Q(0)
        A=sum((Q(comb(k,j))*(1-2**j*tau) for j in range(ell+1)),Q(0))
        candidates.append((f'spectral ell={ell}',max(Q(0),g-tau)/A))
    return max(candidates,key=lambda x:x[1])

def gates(a,M,D,q):
    # Increasing in a and decreasing in q whenever numerator is positive.
    kap=log_interval(1/q);width=log_interval(M/a);step=log_interval(1+1/a)
    return max(Q(0),D*kap[0]-width[1])/step[1], max(Q(0),D*kap[1]-width[0])/step[0]

def injection(a,M,D,q):
    power=Q(1);best=Q(0)
    for ell in range(1,D+1):
        best=max(best,a/q+(D-ell)*(1-q)*a/q-M*power)
        power*=q
    return best

def small(interval,scale=10**12):
    lo,hi=interval
    l=Q((lo*scale).numerator//(lo*scale).denominator,scale)
    h=Q(-((-hi*scale).numerator//(-hi*scale).denominator),scale)
    assert l<=lo<=hi<=h
    return [str(l),str(h)]

def spaces(n):
    dim=2*n
    found={frozenset([0])}
    queue=deque(found)
    while queue:
        v=queue.popleft()
        for w in range(1,1<<dim):
            if w not in v:
                u=v|frozenset(x^w for x in v)
                if u not in found:
                    found.add(u);queue.append(u)
    return sorted(found,key=lambda v:(len(v),sorted(v)))


def rank(v):
    a=len(v);assert a and not a&(a-1)
    return a.bit_length()-1


def pivot_positions(v,n):
    rows=[]
    for x in sorted(v):
        for y in rows:
            x=min(x,x^y)
        if x:
            rows.append(x);rows.sort(reverse=True)
    pivots=[x.bit_length()-1 for x in rows]
    assert len(pivots)==rank(v)
    # Independent exhaustive injection check.
    images={sum(((x>>p)&1)<<i for i,p in enumerate(pivots)) for x in v}
    assert len(images)==len(v)
    return pivots


def symp(v,w,n):
    return sum(((v>>(2*j))&1)*((w>>(2*j+1))&1)+((v>>(2*j+1))&1)*((w>>(2*j))&1) for j in range(n))%2


def finite_envelope_checks():
    summaries=[];checks=0;pivot_checks=0;gate_checks=0;reset_checks=0
    for n,expected in [(1,5),(2,67),(3,2825)]:
        vs=spaces(n);assert len(vs)==expected
        for v in vs:
            r=rank(v);pivots=pivot_positions(v,n)
            subset_ranks=[]
            for s in range(1<<n):
                mask=sum(3<<(2*j) for j in range(n) if (s>>j)&1)
                surviving=frozenset(x for x in v if x&~mask==0)
                rr=rank(surviving);pp=sum((s>>(p//2))&1 for p in pivots)
                assert rr<=pp;pivot_checks+=1
                subset_ranks.append(rr)
                # Backward input-independent reset: project labels onto outside.
                projected=frozenset(x&mask for x in v)
                assert rank(projected)<=r;reset_checks+=1
            for q in [F(1,10),F(1,2),F(9,10)]:
                for t in [F(3,2),F(2),F(3),F(4)]:
                    lhs=sum(q**s.bit_count()*(1-q)**(n-s.bit_count())*t**subset_ranks[s] for s in range(1<<n))
                    rhs_base=1-q+q*t*t
                    assert lhs*lhs<=rhs_base**r,(n,sorted(v),q,t,lhs)
                    checks+=1
            if n<=2:
                for w in range(1,1<<(2*n)):
                    updated=v|frozenset(x^w for x in v) if any(symp(x,w,n) for x in v) else v
                    assert rank(updated)<=r+1;gate_checks+=1
        summaries.append({'n':n,'binary_subspaces':len(vs)})
    q=F(1,10);t=F(2)
    # Stochastic domination by fully paired coordinate survival is false.
    classical_survives=2*q-q*q
    paired_survives=q
    assert classical_survives>paired_survives
    # A parity-aware packed kernel need not preserve convexity of test functions.
    a=1-q+q*t*t;b=1-q+q*t
    second_difference=a*b-2*a+b
    assert second_difference<0
    result={'status':'pass','subspace_sets':summaries,'rational_moment_cases':checks,
            'pivot_coupling_checks':pivot_checks,'reset_rank_checks':reset_checks,
            'T_growth_checks':gate_checks,
            'moment_statement':'E[t^r_after | V]^2 <= (1-q+q*t^2)^r_before, t>=1',
            'stochastic_domination_counterexample':{'q':str(q),'span_Z1_Z2_survival':str(classical_survives),'span_X1_Z1_survival':str(paired_survives)},
            'packed_kernel_convexity_counterexample':{'q':str(q),'t':str(t),'second_difference_at_rank2':str(second_difference)},
            'limits':'Finite exact checks support the pivot/pairing proof; no physical-achievability or universal branch-adaptive-controller claim.'}
    return result

def recovery_u(k,s,c):
    r=2*s+c;F=min(Q(1),Q(2)**(r-2*k),Q(2)**(s-k))
    vals=[]
    for ell in range(k+1):
        tau=Q(1,2**(ell+1)) if ell<k else Q(0)
        A=sum((Q(comb(k,j))*(1-2**j*tau) for j in range(ell+1)),Q(0))
        vals.append((tau+A*F,ell))
    u,ell=min(vals)
    return u,ell,F


def combined_certificate():
    k,n,p,q,z=12,120,Q(2,5),Q(99,100),Q(27,25)
    d=2**k;b0=Q(1,d);Y0,Y1,Y2=Q(1),z**11,z**24
    points=[(Y0,b0),(Y1,Q(1,4)),(Y2,Q(1))]
    slopes=[(b[1]-a[1])/(b[0]-a[0]) for a,b in zip(points,points[1:])]
    assert slopes[0]>slopes[1]>0
    # Upper dyadic sqrt2 enclosure proves the inherited classical beta^k cap.
    scale=1<<100;sqrt2_hi=Q(isqrt(2*scale*scale)+1,scale)
    assert sqrt2_hi**2>2
    beta_hi=(2+sqrt2_hi)/4;classical_cap=beta_hi**k
    witnesses=[]
    for r in range(2*k+1):
        Y=z**r
        a,b=points[:2] if Y<=Y1 else points[1:]
        target=a[1]+(b[1]-a[1])*(Y-a[0])/(b[0]-a[0])
        assert b0<=target<=1
        for s in range(r//2+1):
            c=r-2*s
            if s+c>n:continue
            method=None;detail={}
            dim=min(Q(1),Q(2)**(s+c-k))
            if dim<=target:method='output_dimension'
            elif target**6>=Q(2)**(s-k):method='direct_quantum_dimension_sixth_power'
            else:
                u,ell,F=recovery_u(k,s,c)
                gd=target*target+(1-target)**2/(d-1)
                if u<=gd:
                    method='normalized_DFW_Renes';detail={'ell':ell,'Fcap':str(F),'u':str(u),'g_of_target':str(gd)}
                elif s==0 and classical_cap<=target:method='known_classical_BB84'
            assert method is not None,(r,s,c,target)
            witnesses.append({'r':r,'s':s,'c':c,'target':str(target),'method':method,**detail})
    # Above rank2k, score<=1 and the polygon's constant continuation is enoug
    Yfloor=Y1+(p-Q(1,4))*(Y2-Y1)/Q(3,4)
    assert Yfloor==(4*z**11+z**24)/5
    alpha=ratio(log_interval(1-q+q*z*z),
                  (2*log_interval(z)[0],2*log_interval(z)[1]))
    floor=ratio(log_interval(Yfloor),log_interval(z))
    old=Q(2*k)*(p-b0)/(1-b0)
    cases=[]
    for D,count in [(2000,244),(700,58)]:
        new=(gates(floor[0],Q(2*n),D,alpha[1])[0],
             gates(floor[1],Q(2*n),D,alpha[0])[1])
        assert count-1<new[0]<=new[1]<count
        cases.append({'D':D,'combined_moment_individual_T_interval':small(new),'rounded_necessary_count':count})
    result={'status':'pass','parameters':{'k':k,'n':n,'p':str(p),'q':str(q),'z':str(z)},
            'polygon_ranks':[0,11,24],'polygon_scores':['1/4096','1/4','1'],
            'concavity_slopes':[str(v) for v in slopes],
            'complete_pair_witnesses':witnesses,'pairs_checked':len(witnesses),
            'moment_lower_exact':str(Yfloor),'moment_lower_formula':'(4*z^11+z^24)/5',
            'log_moment_floor_interval':small(floor),'effective_retention_interval':small(alpha),
            'cases':cases,
            'scope':'Valid combined branch-profile upper polygon and conditional necessary count. No globally optimal z, exact physical score curve or feasible physical circuit is asserted.'}
    methods={}
    for witness in witnesses:methods[witness['method']]=methods.get(witness['method'],0)+1
    result.pop('complete_pair_witnesses')
    result['witness_methods']=methods
    result['pair_witness_sha256']=hashlib.sha256(json.dumps(witnesses,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return result

k,n,p,q=12,120,Q(2,5),Q(99,100)
b0=Q(1,2**k)
cap_lower=Q(525654108194,10**12)
w=(p-b0)/(cap_lower-b0)


@lru_cache(None)
def alpha(z):
    if z==1:return q,q
    return ratio(log_interval(1-q+q*z*z),
                   (2*log_interval(z)[0],2*log_interval(z)[1]))


@lru_cache(None)
def floor_upper(z):
    if z==1:return 16*w
    return ratio(log_interval(1+w*(z**16-1)),log_interval(z))[1]


def injection_upper(a,retention,D):
    # A feasible continuous scalar injection schedule: coast for t rounds,
    # then inject delta each round. It upper-bounds the optimized scalar bound.
    width=log_interval(Q(2*n)/a)
    decay=log_interval(1/retention)
    t=min(D,max(0,(width[0]/decay[1]).numerator//(width[0]/decay[1]).denominator))
    return (D-t)*(1-retention)*a/retention


def recovery_family_ceiling():
    # The profile at r16 permits (s,c)=(4,8). Its cap is the Renes inverse
    # of the normalized spectral minimum u=283/1024. The lower cap below is
    # rational and strictly below that exact inverse.
    from math import comb
    # The standard dimension-lift cap 2^s*beta^k is inactive at this
    # witness, so including it cannot remove the feasible profile mixture.
    assert Q(7,5)**2 < 2 and 16*Q(17,20)**12 > 1
    F=Q(1,256)
    us=[]
    for ell in range(k+1):
        tau=Q(1,2**(ell+1)) if ell<k else Q(0)
        A=sum((Q(comb(k,j))*(1-2**j*tau) for j in range(ell+1)),Q(0))
        us.append((tau+A*F,ell))
    u,ell=min(us)
    assert ell==1 and u==Q(283,1024)
    g=cap_lower**2+(1-cap_lower)**2/(2**k-1)
    assert b0<p<cap_lower<1 and g<u
    assert 0<w<1
    targets={2000:223,700:52}
    leaves=[]
    def visit(lo,hi,depth=0):
        assert depth<30
        a_hi=floor_upper(hi);alpha_lo=alpha(lo)[0]
        bounds={D:max(gates(a_hi,Q(2*n),D,alpha_lo)[1],
                      injection_upper(a_hi,alpha_lo,D)) for D in targets}
        if all(bounds[D]<targets[D] for D in targets):
            leaves.append({'z_lo':str(lo),'z_hi':str(hi),
                           'ledger_upper':{str(D):small((Q(0),v))[1] for D,v in bounds.items()}})
        else:
            mid=(lo+hi)/2
            visit(lo,mid,depth+1);visit(mid,hi,depth+1)
    visit(Q(1),Q(2))
    # For z>=2 the feasible mixture is supported on ranks0and16, so its
    # log moment is <=16; alpha(z)>=alpha(2). Both are monotonicity bounds.
    tail={D:max(gates(Q(16),Q(2*n),D,alpha(Q(2))[0])[1],
                injection_upper(Q(16),alpha(Q(2))[0],D)) for D in targets}
    assert all(tail[D]<targets[D] for D in targets)
    result={'status':'pass','parameters':{'k':k,'n':n,'p':str(p),'q':str(q)},
            'specified_family':'Normalized DFW/Renes branch cap + measure-and-prepare + finite output dimension + known classical-memory bound, then rank maximization and exact concavification in z^r.',
            'recovery_rank16_pair':{'s':4,'c':8,'Fcap':str(F),'ell':ell,'u':str(u),'score_cap_lower':str(cap_lower),'g_of_cap_lower':str(g)},
            'feasible_profile_mixture':{'ranks':[0,16],'weight_at_rank16':str(w)},
            'coverage':'Adaptive exact rational subintervals cover1<z<=2; monotone analytic tail coversz>=2.',
            'intervals':leaves,'tail_z_at_least_2':{str(D):small((Q(0),v))[1] for D,v in tail.items()},
            'largest_rounded_ledger_count_at_most':{str(D):v for D,v in targets.items()},
            'scope':'An upper ceiling on the strength of the specified recovery-only moment formulas, including individual-T and continuous-injection ledgers. Not a physical gate upper bound or a comparison with all possible recovery theorems.'}
    result.pop('intervals')
    result['interval_count']=len(leaves)
    result['interval_certificate_sha256']=hashlib.sha256(json.dumps(leaves,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    result['finite_interval_maximum_upper']={str(D):str(max(Q(leaf['ledger_upper'][str(D)]) for leaf in leaves)) for D in targets}
    return result


def main():
    result={
        'status':'PASS',
        'moment_noise':finite_envelope_checks(),
        'combined_branch_certificate':combined_certificate(),
        'specified_recovery_family_ceiling':recovery_family_ceiling(),
        'scope':'Finite algebra and exact scalar certificates supplement the written proofs; no physical-feasibility or model-verification claim.'
    }
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
