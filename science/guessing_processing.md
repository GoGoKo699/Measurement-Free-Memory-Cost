# Delayed-basis guessing and continuing processing

[Model](model.md) · [Processing budgets](control_cost.md) · [Lattice Memory Certification](https://github.com/GoGoKo699/Lattice-Memory-Certification)

Exact recovery of a hidden BB84 string supplies an operational premise for the processing and interval bounds. In the device model below, any unconditional average score above random guessing forces a positive expected Pauli-envelope rank at every circuit cut. For $`k\ge8`$ challenge qubits, the convenient floor is

```math
R_t\ge a_k(p):=2k\left[\frac{p-2^{-k}}{1-2^{-k}}\right]_+.
\qquad\text{(G1)}
```

The full-output envelope and dynamical ledger belong to this project. The finite quantum-memory guessing theorem used below belongs to [Lattice Memory Certification](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/dcff3be0246b27b07fb325c31293199e6fddbcf7/science/guessing-bound.md); its original input-length and dimension notation is translated here to $`k`$ challenge qubits and $`2^s`$ retained quantum dimension. The channel factorization and score-based composition are shared results. Established recovery reductions already yield both (G1) and the small-block mean-rank floor below. The transfer supplies a direct operational proof of their consequences. The [moment-processing chapter](moment_processing.md) retains the actual noise distribution and gives a stronger dynamical ledger from the direct profile and complementary recovery bounds. Section 5 and Appendices A–B give the comparison, including scores at or below one half.

## 1. Challenge, chronology, and counted device

Choose $`X,\Theta\in\{0,1\}^k`$ independently and uniformly. Give the device only $`|X\rangle_\Theta=\bigotimes_{j=1}^kH^{\Theta_j}|X_j\rangle`$. Keep $`X`$ secret throughout. Disclose only $`\Theta`$ after every charged storage operation has ended. The score is the unconditional exact whole-string probability

```math
p=\Pr[\widehat X=X].
\qquad\text{(G2)}
```

An abort is a failure. An accepted run or a score conditioned on a favorable event is not (G2).

Use the [memory model](model.md): $`n`$ is the total noisy persistent width, $`D`$ the number of independent product-depolarizing rounds, and $`q=1-\gamma\in(0,1)`$. Every persistent data, syndrome, controller, auxiliary, and resource record occupies counted slots and undergoes that noise. Between layers the externally prescribed circuit is fixed. Clifford operations and input-independent replacement resets are free; interior T/T-dagger gates are charged. Other primitives require their actual implementation or a proved rank-growth charge. The encoder and basis-dependent final decoding may be arbitrary CPTP endpoints, with no protected side channel between them.

Counting controller records does not make arbitrary feedback free. A controlled Clifford can require a non-Clifford coherent implementation. Nonlinear decisions and outcome-dependent reset states must be represented by the permitted charged primitives. The conclusion applies to this fixed-circuit model, not automatically to an arbitrary adaptive measurement-based controller.

## 2. One envelope and a quantum-plus-classical CPTP bottleneck

Let $`\Phi`$ be the channel to the full physical output **before** basis-dependent readout. Initialize one simultaneous backward envelope at $`W_n=\mathbb F_2^{2n}`$, with $`M=2n`$. It covers every terminal observable for every basis and endpoint decoder. The number of challenge or guessed bits does not justify replacing $`M`$ by $`2k`$.

At a primitive cut $`t`$, average the earlier noise into a CPTP prefix $`\Pi_t`$ from the challenge input to the memory. A realization $`b`$ of future identity/replacement choices selects a CPTP tail $`S_{b,t}`$ to the full physical output. Independent Markov noise makes $`\Pi_t`$ independent of $`b`$. The [envelope invariant](pauli_envelopes.md) and [charged-gate updates](control_cost.md#3-rank-growth-charges-for-gates) give

```math
\mathrm{range}S_{b,t}^\dagger\subseteq\mathcal A_{V_{b,t}}.
\qquad\text{(G3)}
```

Write $`r=\dim V`$, let $`2s`$ be the rank of its restricted symplectic form, and put $`c=r-2s`$. In Clifford coordinates $`U`$, [Lemma C1](control_cost.md#lemma-c1--canonical-algebra-and-conditional-expectation) identifies

```math
\mathcal A_V\simeq M_{2^s}\otimes\mathrm{Diag}_{2^c}\otimes I_G,
\qquad \dim G=2^{n-s-c}.
\qquad\text{(G4)}
```

The trace-preserving conditional expectation $`\mathcal E_V`$ onto this algebra is CPTP, self-adjoint, and satisfies $`S_{b,t}=S_{b,t}\mathcal E_V`$. Its explicit factorization is

```math
\begin{aligned}
F_V(\rho)&=\sum_{z\in\{0,1\}^c}
\mathrm{Tr}_G\!\left[\langle z|_C U\rho U^\dagger|z\rangle_C\right]
\otimes|z\rangle\langle z|_C,\\
J_V(\sigma_{QC})&=U^\dagger
\left(\sigma_{QC}\otimes\frac{I_G}{\dim G}\right)U,
\qquad \mathcal E_V=J_VF_V.
\end{aligned}
\qquad\text{(G5)}
```

The classical blocks in $`F_V(\rho)`$ are positive and their traces sum to one. Both maps are CPTP, and $`F_V`$ has a genuine output of $`s`$ quantum qubits and a classical register with $`2^c`$ values. For an actual endpoint decoder $`\mathcal D_\theta`$, define

```math
K_{b,t}=F_V\Pi_t,
\qquad L_{b,t,\theta}=\mathcal D_\theta S_{b,t}J_V.
\qquad\text{(G6)}
```

These are normalized CPTP encoding and decoding maps whose composition reproduces every actual outcome distribution. This is a channel factorization, not a conditioning on a subnormalized state or a purification with uncounted quantum memory.

The flags $`b`$ are only independent random **channel choices** in a convex decomposition. They are neither supplied information nor outcomes to postselect. One may bound each sampled channel by the optimum decoder that knows $`b`$; doing so only enlarges its score. Here $`p_{b,t}`$ is the actual fixed endpoint decoder's score on the sampled channel; optimized branch scores only upper-bound it. Linearity gives

```math
p=\mathbb E_b p_{b,t}
\qquad\text{(G7)}
```

at every primitive cut, on either side of every gate and noise layer. Coupling the future flags across cuts supplies the same adjacent-cut envelope recurrence as in the original processing proof.

## 3. From branch scores to an all-cut rank floor

Apply two separate bounds to the original CPTP encoder $`K_{b,t}`$:

```math
p_{b,t}\le\min\left\{1,\ 2^{(s-k)/6},\ 2^{s+c-k}\right\}.
\qquad\text{(G8)}
```

The quantum-dimension term is Lattice's finite BB84 theorem, which permits unlimited classical information. For the last term, fix a basis. Its $`d=2^k`$ uniform inputs give normalized output states $`\rho_x\le I_{QC}`$. Every guessing POVM therefore obeys

```math
\frac1d\sum_x\mathrm{Tr}(M_x\rho_x)
\le\frac1d\sum_x\mathrm{Tr}M_x
=\frac{2^{s+c}}d.
\qquad\text{(G9)}
```

Use (G9) before the optional classical pure-state refinement in Lattice's proof: those extra mathematical labels need not fit in the physical $`c`$-bit register. Uniform $`X`$ and the absence of an extra supplied label are essential here.

Since $`s+c=r-s`$, maximize the right side of (G8) over the relaxed range $`0\le s\le r/2`$. Below $`r=5k/6`$ the maximum is at $`s=0`$; above it, the two exponential bounds meet at $`s=(6r-5k)/7`$. This yields

```math
p_{b,t}\le h_k(r):=
\begin{cases}
2^{r-k},&0\le r\le5k/6,\\
2^{(r-2k)/7},&5k/6\le r\le2k,\\
1,&r\ge2k.
\end{cases}
\qquad\text{(G10)}
```

Integer $`s,c`$ and the constraint $`s+c\le n`$ can only lower this upper curve. It is not claimed to be an attainable optimum.

Put $`b_0=2^{-k}`$, $`H_k=2^{-k/6}`$, and $`t_0=5k/6`$. Each exponential piece is convex. For $`k\ge8`$, its shared endpoint obeys

```math
H_k\le b_0+(1-b_0)\frac5{12}.
\qquad\text{(G11)}
```

It holds directly at $`k=8`$ and persists: with $`x=2^{-k/6}`$, the ratio $`(H_k-b_0)/(1-b_0)=(x-x^6)/(1-x^6)`$ increases with $`x\in(0,1)`$. Its derivative has numerator $`1-6x^5+5x^6\ge0`$. Convexity on the two pieces now puts $`h_k`$ below the chord from $`(0,b_0)`$ to $`(2k,1)`$; for $`r\ge2k`$ that chord is at least one. Average this affine bound using (G7) to obtain (G1) at **every** cut.

For $`1\le k\le7`$, the concave polygon through $`(0,b_0),(t_0,H_k),(2k,1)`$, extended constantly thereafter, majorizes $`h_k`$. Jensen's inequality and inversion give the all-cut floor $`R_t\ge a_k(p)`$ with

```math
a_k(p)=
\begin{cases}
0,&p\le b_0,\\
t_0\dfrac{p-b_0}{H_k-b_0},&b_0\lt p\le H_k,\\
t_0+(2k-t_0)\dfrac{p-H_k}{1-H_k},&H_k\le p\le1.
\end{cases}
\qquad\text{(G12)}
```

Its decreasing slopes are equivalent to the reversed inequality in (G11), which holds for these seven integers. A simple affine alternative, valid for all $`k`$, is

```math
R_t\ge
\frac{[p-b_0]_+}
{\max\{(H_k-b_0)/t_0,\ (1-b_0)/(2k)\}}.
\qquad\text{(G13)}
```

Thus every above-chance score supplies a positive rank floor, including scores at or below one half. Averaging an inverted exponential, or demanding the floor on each noise branch individually, would not establish this conclusion.

## 4. Gate and contiguous-interval consequences

Take any proved floor $`a\gt0`$ from Section 3. A backward noise layer gives $`R'\le qR`$, a free primitive gives $`R'\le R`$, and a charged gate gives $`R'\le R+b_g`$. With $`B=\sum_gb_g`$ and $`\kappa=-\ln q`$, the existing [processing laws](control_cost.md#4-logarithmic-and-additive-resource-laws) therefore give

```math
B\ge a[D\kappa-\ln(2n/a)]_+,
\qquad
m_T\ge\left\lceil
\frac{[D\kappa-\ln(2n/a)]_+}{\ln(1+1/a)}
\right\rceil.
\qquad\text{(G14)}
```

The second inequality applies to an actual Clifford+T implementation and counts both T and T-dagger. The stronger [discrete injection bound](control_cost.md#5-exact-finite-depth-bound-from-discrete-injections) is

```math
\begin{aligned}
B&\ge\mathcal L(D,2n,a,q),\\
\mathcal L(D,M,a,q)
&=\max\left\{0,\max_{1\le\ell\le D}
\left[\frac aq+(D-\ell)\frac{(1-q)a}{q}-Mq^{\ell-1}\right]\right\}.
\end{aligned}
\qquad\text{(G15)}
```

Set $`\mathcal L(0,M,a,q)=0`$. For integer T count, round this lower bound upward and retain the stronger applicable inequality. The function is the exact optimum of necessary continuous scalar injection constraints, not a sufficient processing budget or a matching physical construction.

For any contiguous interval $`I`$ between primitive cuts, the **same global floor** holds throughout. Replace $`D,B,2n`$ in (G14)–(G15) by $`D_I,B_I,M_I`$, where $`M_I\le2n`$ is an upper bound on expected envelope rank at its later cut. The later suffix can determine that entering random envelope, but noise within $`I`$ is independent of the suffix. Hence the contraction and charge recurrence hold unchanged. Assign boundary gates consistently. If $`D_I\ge1`$, feasibility also requires $`a\le2nq`$.

For $`k=12,p=0.4,n=120,\gamma=0.01,D=2000`$, (G1) gives $`a=9.5964835\ldots`$. The individual-gate bound requires at least 171 T/T-dagger gates; (G15) alone requires 163. For a 700-round contiguous interval, the individual-gate bound requires at least 39; (G15) alone requires 37. These are conditional necessary statements, not an assertion that these parameters are achievable.

## 5. Comparison with existing fidelity and recovery reductions

Let $`F_{\rm rec}`$ be the optimum squared Bell-overlap fidelity after CPTP recovery of the input-to-output channel. Let $`F_{\rm pg}\le F_{\rm rec}`$ use pretty-good recovery, and let $`u`$ be basis-averaged pretty-good guessing success. The full sampled branch channel remains CPTP, and every recovered logical channel factors through $`\mathcal E_V`$. The existing superoperator and coherent-block bounds therefore give

```math
F_{\rm rec}\le\min\{1,2^{r-2k},2^{s-k}\}.
\qquad\text{(G16)}
```

For clarity, the following comparison uses an optimal average guessing score $`p_*\ge\max\{p,1/d\}`$, with $`d=2^k`$. Replacing it by an actual lower score $`p\ge1/d`$ preserves the inequalities below. [Renes (2017), Eq. (12)](https://arxiv.org/abs/1707.01114), applied basis by basis and then averaged by convexity, gives

```math
g_d(p_*):=p_*^2+\frac{(1-p_*)^2}{d-1}\le u.
\qquad\text{(G17)}
```

Barnum–Knill's weaker $`p_*^2\le u`$ suffices for the proposition below. With $`S_\ell=\sum_{j=0}^{\ell}\binom kj`$, [Dupuis–Fawzi–Wehner (2015), Eq. (27) and Appendix C](https://arxiv.org/abs/1305.1316v3) imply

```math
u\le S_\ell F_{\rm pg}
+2^{-\ell-1}\mathbf1_{\ell\lt k},
\qquad 0\le\ell\le k.
\qquad\text{(G18)}
```

Here $`2^{-H_2(R|B)}=dF_{\rm pg}`$ for the normalized channel Bell state. The fidelity in the cited operational interpretation is a root fidelity; $`F_{\rm pg}`$ here is its square.

Keeping normalization in the same proof gives a sharper finite bound. Set

```math
\begin{aligned}
\tau_\ell&=2^{-\ell-1}\mathbf1_{\ell\lt k},\\
A_\ell&=\sum_{j=0}^{\ell}\binom kj(1-2^j\tau_\ell).
\end{aligned}
\qquad\text{(G18a)}
```

The spectral weights in [Dupuis–Fawzi–Wehner, Eqs. (13)–(15)](https://arxiv.org/abs/1305.1316v3), divided by $`d`$, give numbers $`w_v\ge0`$ with $`\sum_vw_v=1`$ and $`w_v\le F_{\rm pg}`$. For product BB84 measurements, their spectrum gives

```math
u=\sum_{v\in\{0,1,3\}^k}2^{-|v|}w_v
\le\tau_\ell+A_\ell F_{\rm pg}.
\qquad\text{(G18b)}
```

Here $`|v|`$ counts nonidentity Pauli labels; labels containing $`2`$ have spectral coefficient zero. To prove the inequality, subtract $`\tau_\ell\sum_vw_v=\tau_\ell`$ from the spectral sum. Discard nonpositive terms, and apply $`w_v\le F_{\rm pg}`$ to the remaining $`2^j\binom kj`$ coefficients at each weight $`j\le\ell`$. This retains the probability mass already assigned to those terms. Since $`A_\ell\le S_\ell`$, it strengthens (G18), with equality at $`\ell=k`$. This is a consequence of the cited spectral constraints, not a separately assumed uncertainty relation. Equations (G17) and (G18b) yield

```math
F_{\rm rec}\ge
\max_{0\le\ell\le k}
\frac{[g_d(p_*)-\tau_\ell]_+}{A_\ell}.
\qquad\text{(G19)}
```

This comparison is useful below one half. At $`p_*=1/2`$, for example, $`\ell=2`$ gives a positive lower bound of order $`k^{-2}`$. For fixed positive score, a fixed $`\ell`$ with $`2^{-\ell-1}\lt p_*^2`$ gives a polynomial-in-$`k`$ fidelity lower bound.

Two simple consequences complete the comparison. For each basis, applying its guessing POVM and preparing the guessed basis vector gives Bell fidelity $`p_\theta/d`$. A recovery channel that samples the basis internally therefore gives $`F_{\rm rec}\ge p_*/d`$. Combining this with (G16) yields $`p_*\le\min\{1,2^{r-k}\}`$. In addition, (G18b) with $`\ell=0`$ gives

```math
u\le\frac{1+F_{\rm pg}}2,
\qquad
F_{\rm rec}\ge[2g_d(p_*)-1]_+.
\qquad\text{(G20)}
```

The first inequality also agrees with the complementary-pair bound in [Berta–Coles–Wehner (2014), Eq. (B23)](https://arxiv.org/abs/1302.5902v2). The spectral proof above applies directly to product BB84 measurements and needs no completion to a mutually unbiased set. The entire product-BB84 collection is not a complete mutually unbiased set or a 2-design.

**Prior-recovery comparison.** For every integer $`k\ge8`$, (G16), (G18b), and $`p_*^2\le u`$, together with the measure-and-prepare reduction, already imply, branch by branch,

```math
p_*\le2^{-k}+(1-2^{-k})\min\left\{\frac r{2k},1\right\}.
\qquad\text{(G21)}
```

Appendix A proves the statement analytically for $`k\ge11`$ and supplies exact rational certificates for $`k=8,9,10`$. Averaging it gives precisely (G1). Appendix B gives 49 exact rational certificates that the same prior-recovery family lies below the concave polygon used for (G12), at every interior integer rank for $`1\le k\le7`$. Thus prior recovery already implies **both displayed all-cut floors**, and therefore all gate and interval consequences obtained from them.

The comparison must be performed before averaging future noise: using only a recovery bound for the already averaged score can lose strength. The direct quantum-plus-classical proof supplies a short operational route. This conclusion concerns the convenient floors (G1) and (G12), not optimal guessing or every stronger bound obtainable by optimizing integer quantum/classical block sizes.

### Distribution-aware strengthening

The comparison above concerns the convenient mean-rank floors (G1) and (G12). It does not exhaust the information in the noisy circuit. The [exponential-moment theorem](moment_processing.md) keeps $`\mathbb E[z^r]`$ for $`1\lt z\le2`$ rather than replacing the sampled ranks by an arbitrary mixture with the same mean. A pivot-coordinate argument contracts $`X_z=\log_z\mathbb E[z^r]`$ by $`\alpha_z=\ln(1-q+qz^2)/(2\ln z)`$ at each noise layer; a T/T-dagger gate adds at most one.

That theorem fixes one recovery of the averaged full channel, then uses its branchwise Bell fidelities and their linear average at every cut. This differs from applying an averaged-score recovery estimate only to the mean-rank ledger. The recovery lower bound is inherited; the stronger consequence comes from the dynamical moment estimate. Applying the direct quantum-plus-classical branch curve in the moment variable gives a complementary floor and keeps every above-chance score useful. The full-output initial potential is $`2n`$, and an interval reserve must bound $`\log_z\mathbb E[z^r]`$, rather than merely $`\mathbb E[r]`$.

## 6. Exact boundaries and interpretation

At $`p=2^{-k}`$, discard the challenge and guess randomly at the free endpoint. This uses zero interior processing for any duration, so no positive budget follows at or below chance. An above-chance score can also reflect classical information; the processing conclusion alone does not establish quantum storage.

At fixed width, a mathematical noise branch replaces all $`n`$ registers in one round with probability $`\gamma^n`$. Thereafter the output is independent of $`X`$, even given $`\Theta`$, and its score is at most $`2^{-k}`$. Other branches have score at most one. Consequently every gate budget obeys

```math
p\le2^{-k}+(1-2^{-k})(1-\gamma^n)^D.
\qquad\text{(G22)}
```

The continuing budget is a necessary law along feasible families, not proof of indefinitely reliable storage at fixed width. The flags in this argument remain mathematical proof aids.

A guessing score does not verify the noise strength, independence, controller restrictions, or absence of side channels. Those remain model premises. Estimating $`p`$ from experimental trials is a separate statistical task. No LWE or remote-preparation wrapper is imported: such a compilation would need its own joint-state, chronology, controller, and resource-preservation reduction. These conclusions establish neither optimal T count nor a hardware certification procedure.

## A. Proof that prior recovery already gives the large-block floor

It suffices to show the stronger $`p_*\le r/(2k)`$ at every integer $`1\le r\le2k`$ for $`k\ge11`$. Rank zero is a replacement channel with score at most $`2^{-k}`$; ranks at least $`2k`$ are trivial.

We use an endpoint criterion. For $`G(x)=A2^x+c-x^2/(4k^2)`$ with $`A\gt0`$, $`G''`$ is strictly increasing. If $`G'(L)\le0`$, then $`G'`$ either stays negative or changes from negative to positive after its possible minimum. Thus $`G`$ has no interior maximum on $`[L,U]`$, and $`G(x)\le\max\{G(L),G(U)\}`$.

**First range, $`1\le r\le4k/5`$.** The measure-and-prepare reduction gives $`p_*\le2^{r-k}`$. The convex function $`2^{x-k}-x/(2k)`$ is nonpositive at $`x=1`$ for $`k\ge4`$ and at $`x=4k/5`$ because $`2^{-k/5}\le2/5`$ for $`k\ge11`$. It is therefore nonpositive throughout this range.

**Second range, $`4k/5\le r\le k`$.** Set $`S_2=(k^2+k+2)/2`$. Equations (G16) and (G18) with $`\ell=2`$ give $`p_*^2\le S_2 2^{r-2k}+1/8`$. The endpoint criterion with $`L=4k/5,U=k`$ applies if

```math
\begin{aligned}
S_2 2^{-6k/5}&\le7/200,\\
S_2 2^{-k}&\le1/8,\\
(\ln2)kS_2 2^{-6k/5}&\le2/5.
\end{aligned}
\qquad\text{(A1)}
```

At $`k=11`$, $`S_2=67`$ and $`2^{-66/5}\lt2^{-13}`$ verify all three inequalities, using $`\ln2\lt1`$. They persist for larger $`k`$ because $`(\ln S_2)'\le2/k`$, so each left side decreases on $`[11,\infty)`$. Thus $`p_*^2\le r^2/(4k^2)`$.

**Third range, $`k+1\le r\le3k/2`$.** With $`\ell=1`$, $`p_*^2\le(k+1)2^{r-2k}+1/4`$. Apply the endpoint criterion with $`L=k+1,U=3k/2`$. Sufficient conditions are

```math
\begin{aligned}
(k+1)2^{1-k}&\le\frac1{2k}\lt\frac{2k+1}{4k^2},\\
(k+1)2^{-k/2}&\le5/16,\\
(\ln2)k^2 2^{1-k}&\le1/2.
\end{aligned}
\qquad\text{(A2)}
```

At $`k=11`$, the first and third follow from $`132/1024\lt1/2`$ and $`121/1024\lt1/2`$. The middle one is $`12/(32\sqrt2)\le5/16`$, equivalent after squaring to $`36\le50`$. Their left sides, after multiplying the first by $`k`$, decrease thereafter by their logarithmic derivatives. These conditions bound both endpoint values by zero and ensure the left derivative is nonpositive. No integer rank is omitted between the second and third ranges.

**Fourth range, $`3k/2\le r\le2k`$.** Equations (G16) and (G20) give $`p_*\le H(r):=\sqrt{(1+2^{r-2k})/2}`$. With $`z=2^{r-2k}`$,

```math
H''(r)=\frac{(\ln2)^2z(z+2)}{4\sqrt2(1+z)^{3/2}}\gt0.
\qquad\text{(A3)}
```

The right endpoint satisfies $`H(2k)=1`$. At the left, $`H(3k/2)\le3/4`$ follows from $`2^{-k/2}\le1/8`$, valid for $`k\ge6`$. Convexity puts the entire curve below $`r/(2k)`$. The four ranges finish the analytic proof.

For $`k=8,9,10`$, the following exact rational calculation verifies all 54 remaining integer ranks. For each, either measure-and-prepare directly bounds the score by $`r/(2k)`$, or the normalized spectral bound (G18b) bounds its square by $`r^2/(4k^2)`$:

```python
from fractions import Fraction as F
from math import comb

for k in (8, 9, 10):
    for r in range(1, 2*k + 1):
        f = F(2**r, 2**(2*k))
        candidates = []
        for ell in range(k+1):
            tau = F(1, 2**(ell+1)) if ell < k else F(0)
            A = sum(comb(k, j)*(1-2**j*tau)
                    for j in range(ell+1))
            candidates.append(tau + A*f)
        z = min(candidates)
        line = F(r, 2*k)
        assert F(2**r, 2**k) <= line or z <= line*line
```

These finite certificates complete (G21); they do not replace the cited general recovery theorems.

## B. Small-block recovery comparison

For $`1\le k\le7`$, let $`P_k(r;H)`$ be the polygon from Section 3 with vertices $`(0,b_0),(t_0,H),(2k,1)`$. We show that the same recovery inequalities already give $`p_*\le P_k(r;H_k)`$ at every integer rank. Rank zero is a replacement channel; ranks at least $`2k`$ have the trivial upper bound one. Only $`\sum_{k=1}^7(2k-1)=49`$ interior cases remain.

These cases admit exact rational certificates. Put $`N=2^{40}`$ and choose the integer $`m`$ by bisection so that

```math
m^6 2^k\le N^6\lt(m+1)^6 2^k.
\qquad\text{(B1)}
```

Then $`H_-=m/N\le H_k\lt(m+1)/N`$. At fixed interior rank, both polygon pieces increase with $`H`$: their derivatives are $`r/t_0`$ and $`(2k-r)/(2k-t_0)`$. Hence the rational value $`L=P_k(r;H_-)`$ is no greater than $`P_k(r;H_k)`$ and lies in $`[1/d,1]`$.

For each of the 49 cases, exact fraction arithmetic verifies at least one of

```math
2^{r-k}\le L,
\qquad\text{or}\qquad
\min_{0\le\ell\le k}
\left(\tau_\ell+A_\ell2^{r-2k}\right)\le g_d(L).
\qquad\text{(B2)}
```

The first is the measure-and-prepare bound. The second, together with (G16)–(G18b), gives $`g_d(p_*)\le g_d(L)`$. Since $`p_*\ge1/d`$ and $`g_d`$ increases on $`[1/d,1]`$, it implies $`p_*\le L`$. Averaging the resulting polygon bound and applying Jensen's inequality reproduces (G12).

The [canonical verifier](../verification/verify_guessing_processing.py) implements the integer bisection and all 49 checks in `exact_small_block_recovery_certificates`; the [reference output](../verification/expected/guessing_processing.json) records every rational endpoint and witness. No floating-point tolerance enters these certificates. They establish this finite comparison of the displayed bounds, not an optimal operational score curve.

## Finite guessing and quantified decay refinements

For three challenges and one retained quantum qubit, the [Lattice finite refinement](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/900b54da514d462c7969da80abff0bf0b7c84e49/science/finite-guessing.md) gives $`p\le(2+7\sqrt2)/16`$, including unrestricted classical records and arbitrary encoders. It uses a rank-sensitive two-copy purity constraint, followed by the inherited finite-alphabet discrimination inequality. The stated normalized recovery-only specialization gives the weaker $`(1+2\sqrt7)/8`$. This finite improvement supplies a sharper branch cap at $`k=3,s=1`$; form the concave majorant in $`z^r`$ before averaging. It does not tensorize against arbitrary joint encoders.

The [decay corollary](guessing_decay.md) proves a finite inverse of the existing injection ledger and an explicit asymptotic exponential score guarantee. Its critical processing rate is a necessary-converse threshold; no achievable capacity or optimal numerical coefficient is inferred.


## 9. Reliability-sensitive processing at fixed challenge size

The [rank-extinction reliability law](reliability_processing.md#6-delayed-basis-bb84-scores) uses the same common full-output envelope before basis-dependent readout. Its zero event gives input-independent output and therefore score at most chance, including aborts. For every challenge size and every strictly above-chance unconditional score, it yields a finite T/interval bound without a mean-rank chord approximation. At fixed challenge size, positive noise and fixed above-chance score, logarithmic physical width requires order D ln D continuing processing. This includes retained classical input dependence and does not supply a quantum-dimension certificate or an LWE preparation wrapper. The full proof and parameter domains are owned by the linked chapter.
