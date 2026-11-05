# Explicit decay from the full noisy generating function

[Moment dynamics](moment_processing.md) · [Discrete ledger](control_cost.md) · [BB84 interface](guessing_processing.md)

Keeping the entire generating function through each noise layer strengthens the earlier fixed-parameter logarithmic-moment converse. The resulting finite theorem covers every prescribed placement of interior T/T-dagger slots. Its asymptotic consequence gives a strictly larger whole-string decay exponent at every subcritical processing rate, including the random-guessing exponent $`\ln2`$ when the normalized rate tends to zero. These are necessary resource statements; they do not establish an achievable processing capacity.

## 1. Premises, the complete profile and its transport

Use the unchanged [memory model](model.md). The device receives only $`k`$ uniformly random BB84 challenge qubits; $`X`$ stays secret and only the basis string is disclosed after storage. The score $`p`$ is unconditional exact-string success, with aborts failing. All $`n`$ persistent physical qubits, including controller, syndrome and workspace registers, undergo $`D`$ independent product-depolarizing rounds. Put $`q=1-\gamma\in(0,1)`$. Encoding and final decoding remain arbitrary free CPTP endpoints. The budget $`m_T`$ counts prescribed interior T/T-dagger slots of the fixed circuit. Mathematical noise flags remain unavailable to the controller.

At each primitive cut use the same full-output envelope from [moment dynamics](moment_processing.md), before basis-dependent readout. For its sampled rank $`r`$ define

```math
\mathcal M_t(z)=\mathbb E[z^{r_t}],\qquad
H_t(\lambda)=\ln\mathcal M_t(e^\lambda),\qquad \lambda\ge0.
```

The envelope starts at deterministic rank $`2n`$. The inherited branch score bound is

```math
p_{b,t}\le h_k(r_{b,t}),\qquad
h_k(r)=\begin{cases}
2^{r-k},&0\le r\le5k/6,\\
2^{(r-2k)/7},&5k/6\le r\le2k,\\
1,&r\ge2k.
\end{cases}
\qquad\text{(D1)}
```

For every $`\lambda\ge0`$, its complete exponential support is

```math
h_k(r)\le e^{\lambda r-kA(\lambda)},\qquad
A(\lambda)=\min\left\{2\lambda,
\frac{\ln2+5\lambda}{6},\ln2\right\}.
\qquad\text{(D2)}
```

Indeed, $`\lambda r-\ln h_k(r)`$ is affine on each of the three regions and attains its minimum at $`r=0,5k/6,2k`$. This proves (D2) for all real $`r\ge0`$, and in particular for integer envelope ranks. The three slopes of $`A`$ are $`2,5/6,0`$, with breakpoints $`\ln2/7`$ and $`\ln2`$. Actual common-decoder score averaging, $`p=\mathbb E_b p_{b,t}`$, now gives the all-cut floor

```math
H_t(\lambda)\ge\ln p+kA(\lambda).
\qquad\text{(D3)}
```

It is not a score premise on individual noise branches. For $`p=0`$ all subsequent score upper bounds hold trivially; use $`-\ln0=+\infty`$ in limits.

### Lemma D1 — Full generating-function transition

Define

```math
f_q(z)=\sqrt{1-q+qz^2},\qquad
\phi_j(\lambda)=\frac12\ln\!\left[1+q^j(e^{2\lambda}-1)\right],
\quad j\ge0.
\qquad\text{(D4)}
```

One backward noise layer obeys

```math
\mathcal M'(z)\le\mathcal M(f_q(z)),\qquad
H'(\lambda)\le H(\phi_1(\lambda)).
\qquad\text{(D5)}
```

Moreover $`f_q^{\circ j}(z)^2=1+q^j(z^2-1)`$, or equivalently $`\phi_1^{\circ j}=\phi_j`$. A T/T-dagger slot obeys $`H'(\lambda)\le H(\lambda)+\lambda`$; free Clifford operations and input-independent resets do not increase $`H`$.

**Proof.** The conditional moment lemma gives $`\mathbb E[z^{r'}\mid V]\le(1-q+qz^2)^{r/2}=f_q(z)^r`$. Average over the incoming suffix envelope. Independence of the next noise layer from that envelope is exactly the existing model premise. No Jensen step is needed. The iterate identity follows from the affine recurrence for $`f_q(z)^2-1`$. A T gate increases envelope rank by at most one on every branch; the remaining primitive statements are inherited. $`\square`$

This retains more information than replacing $`\mathcal M(f_q(z))`$ by $`\mathcal M(z)^{\alpha_z}`$, where $`\alpha_z=\phi_1(\lambda)/\lambda`$. Both inequalities remain valid, but their finite relaxations differ.

## 2. Finite width and time for every fixed gate schedule

### Theorem D2 — A generating-function cut ledger

For $`D\ge1`$, choose any integer $`1\le L\le D`$ and put $`T_L=D-L+1`$. Then every permitted circuit satisfies, for every $`\lambda\ge0`$,

```math
-\ln p\ge
\left[kA(\lambda)-
\frac{2n\sum_{d=L}^{D}\phi_d(\lambda)
+m_T\sum_{j=1}^{T_L}\phi_j(\lambda)}{T_L}\right]_+.
\qquad\text{(D6)}
```

Take the strongest bound over $`L`$ and $`\lambda`$. It suffices to use $`0\le\lambda\le\ln2`$: above $`\ln2`$, $`A`$ is constant and both cost terms increase. For $`D=0`$ this theorem gives no nontrivial decay statement.

**Proof.** Traverse the circuit backward. Cut $`d`$ is immediately after crossing the $`d`$th noise layer and before following backward primitives, for $`1\le d\le D`$. Let $`a_h`$ count T/T-dagger slots in backward gap $`h`$, after noise layer $`h`$ and before layer $`h+1`$, with $`h=0`$ denoting the gap before the first noise layer. Slots after the last backward noise layer affect none of these cuts and can be omitted. Thus $`\sum_{h=0}^{D-1}a_h\le m_T`$.

Applying Lemma D1 successively, and evaluating every preceding inequality at the transported argument, gives

```math
H_d(\lambda)\le2n\phi_d(\lambda)
+\sum_{h=0}^{d-1}a_h\phi_{d-h}(\lambda).
\qquad\text{(D7)}
```

This induction allows arbitrary fixed Clifford/reset operations and any placement of charged slots within each gap. It counts each prescribed slot regardless of input or coherent control.

Sum (D7) over the selected cuts $`d=L,\ldots,D`$. For a slot in gap $`h`$, its contributing ages $`d-h`$ are distinct positive integers, with at most $`T_L`$ terms. Since $`\phi_j(\lambda)`$ is nonnegative and decreases with $`j`$,

```math
\sum_{\substack{L\le d\le D\\d\gt h}}
\phi_{d-h}(\lambda)
\le\sum_{j=1}^{T_L}\phi_j(\lambda).
```

Consequently the average $`H_d`$ is at most the fraction in (D6). Its average is also at least $`\ln p+kA(\lambda)`$ by (D3). Rearrange and use $`p\le1`$ to obtain the positive part. $`\square`$

The term containing $`2n`$ is the explicit initial-width allowance. The cut averaging makes the theorem independent of the gate schedule. There is no random selection of storage circuits, expected executed charge, free nonlinear controller, or supplied noise flag.

A stronger finite branch majorant can be composed with this ledger. If $`P_z`$ is any established nondecreasing concave majorant of the branch score in $`Y=z^r`$, set $`B_L`$ equal to the fraction in (D6). At least one selected cut has $`H_d\le B_L`$, so

```math
p\le P_z(e^{B_L}),\qquad z=e^\lambda.
\qquad\text{(D8)}
```

This preserves the existing full polygon and integer-algebra/recovery refinements. They need not be reduced to the simpler support (D2).

## 3. A stronger asymptotic exponential guarantee

### Corollary D3 — Full-profile necessary-rate exponent

Fix $`q\in(0,1)`$ and finite $`\rho\ge0`$. For any sequence of permitted circuits with

```math
k_i\to\infty,\quad D_i\to\infty,\quad n_i\ge k_i,\quad
\ln(n_i/k_i)=o(D_i),\quad
\limsup_i\frac{m_{T,i}}{k_iD_i}\le\rho,
\qquad\text{(D9)}
```

the unconditional scores satisfy

```math
\liminf_i\frac{-\ln p_i}{k_i}\ge\mathcal E_q(\rho),
\qquad
\mathcal E_q(\rho)=
\max_{0\le\lambda\le\ln2}
\left[A(\lambda)-\rho S_q(\lambda)\right]_+,
```

```math
S_q(\lambda)=\sum_{j=1}^{\infty}\phi_j(\lambda)
=\frac12\sum_{j=1}^{\infty}
\ln\!\left[1+q^j(e^{2\lambda}-1)\right].
\qquad\text{(D10)}
```

**Proof.** First fix $`\lambda\in[0,\ln2]`$ independently of $`i`$. The inequality $`\ln(1+x)\le x`$ gives $`\phi_j(\lambda)\le(e^{2\lambda}-1)q^j/2`$, so $`S_q(\lambda)`$ converges. Choose

```math
L_i=1+\left\lceil\frac{\ln(2n_i/k_i)}{-\ln q}\right\rceil.
```

Then $`L_i=o(D_i)`$, eventually $`L_i\le D_i`$, and $`2n_iq^{L_i}\le k_i`$. The normalized initial-width term in (D6) obeys

```math
\frac{2n_i}{k_iT_{L_i}}
\sum_{d=L_i}^{D_i}\phi_d(\lambda)
\le\frac{e^{2\lambda}-1}{2(1-q)T_{L_i}}\to0.
```

Also $`T_{L_i}/D_i\to1`$, and the normalized gate term has limsup at most $`\rho S_q(\lambda)`$. Apply (D6) to this fixed $`\lambda`$, then take the supremum of the resulting lower bounds on the same liminf. This uses no interchange of an optimized finite parameter with a limit. Continuity and compactness make the supremum a maximum. Scores equal to zero are covered by the extended logarithm convention. $`\square`$

For every $`\eta\gt0`$, a positive exponent implies the eventual bound $`p_i\le\exp[-k_i(\mathcal E_q(\rho)-\eta)]`$. The noise parameter is fixed. Simultaneous weak-noise limits need uniform estimates in addition to this corollary.

## 4. Positivity, endpoints and certified optimization

Termwise differentiation is justified uniformly on the compact parameter interval by geometric bounds. It gives

```math
S_q'(\lambda)=\sum_{j\ge1}
\frac{q^je^{2\lambda}}{1+q^j(e^{2\lambda}-1)},
\qquad
S_q'(0)=\frac q{1-q},\qquad S_q''(\lambda)\gt0.
```

For $`\rho\gt0`$, $`A-\rho S_q`$ is strictly concave, including across the two downward slope jumps of $`A`$. Its derivative at zero is $`2-\rho q/(1-q)`$. Therefore

```math
\mathcal E_q(\rho)\gt0
\quad\Longleftrightarrow\quad
0\le\rho\lt\frac{2(1-q)}q,
\qquad
\mathcal E_q(0)=\ln2.
\qquad\text{(D11)}
```

At and above the critical rate, strict convexity gives $`S_q(\lambda)\gt\lambda q/(1-q)`$ for $`\lambda\gt0`$, while $`A(\lambda)\le2\lambda`$, so the exponent is zero. Below it, the positive derivative at zero gives a positive value. At zero normalized rate, the maximum is the right endpoint $`\lambda=\ln2`$. This recovers the random-guessing exponent, which is compatible with discarding the challenge and guessing uniformly.

For a positive subcritical rate the untruncated objective has one positive global maximum. It is either a derivative zero in one of the smooth regions, the kink $`\lambda=\ln2/7`$, or the right endpoint. Its smooth derivative is $`2-\rho S_q'`$ below the kink and $`5/6-\rho S_q'`$ above it. Certified signs of these monotone derivatives therefore determine and bracket the global maximizer.

The [verifier](../verification/verify_guessing_decay.py) uses rational outward enclosures throughout. After a finite head, choose $`J`$ so that $`b=(e^{2\lambda}-1)q^J\le1/4`$. The remaining tails are

```math
\begin{aligned}
\sum_{j=J}^{\infty}\phi_j(\lambda)
&=\frac12\sum_{\ell=1}^{\infty}
\frac{(-1)^{\ell+1}b^\ell}{\ell(1-q^\ell)},\\
\sum_{j=J}^{\infty}\phi_j'(\lambda)
&=e^{2\lambda}q^J\sum_{\ell=1}^{\infty}
\frac{(-1)^{\ell+1}b^{\ell-1}}{1-q^\ell}.
\end{aligned}
```

In each series successive absolute terms have ratio at most $`b`$. The 60-term even partial sum and the next odd partial sum enclose the exact tail. Logarithmic head terms have an explicit atanh-series remainder; the algebraic kink is enclosed by exact seventh-power comparisons. Thus the reported optimizer and sum enclosures do not rely on floating-point convergence or a numerical global-search assertion.

At $`q=99/100`$, the certified exponent intervals are:

| $`\rho`$ | Maximizer | $`\mathcal E_q(\rho)`$ interval |
|---:|---|---:|
| $`0`$ | Right endpoint | $`[0.693147180559,0.693147180560]`$ |
| $`0.001`$ | Right endpoint | $`[0.597010325033,0.597010325034]`$ |
| $`0.005`$ | Interior above the kink | $`[0.217453493626,0.217453493627]`$ |
| $`0.01`$ | Kink | $`[0.095025987871,0.095025987872]`$ |
| $`0.015`$ | Kink | $`[0.043517956012,0.043517956013]`$ |
| $`0.02`$ | Interior below the kink | $`[0.000100281684,0.000100281685]`$ |
| $`2/99`$ | Zero objective | $`0`$ |

These certify the displayed converse objective, rather than the best possible operational exponent.

## 5. Comparison with the earlier fixed-parameter ledger

The earlier scalar result remains valid. For $`1\lt z\le2^{1/7}`$, put $`\lambda=\ln z`$, $`\alpha_z=\phi_1(\lambda)/\lambda`$ and $`M=2n`$. Its finite inverse is

```math
U_z=\min\left\{M\alpha_z,\;
\min_{1\le\ell\le D}
\frac{\alpha_z(m_T+M\alpha_z^{\ell-1})}
{1+(D-\ell)(1-\alpha_z)}\right\},
\qquad p\le\min\{1,z^{U_z-2k}\}.
\qquad\text{(D12)}
```

This is obtained by rearranging every discrete injection inequality at the all-cut floor $`[2k+\ln p/\lambda]_+`$, together with the immediate post-noise width cap $`M\alpha_z`$. If the floor is zero, the displayed score bound follows directly. The initial allowance is retained. The same inversion with a proved concave profile gives $`p\le P_z(z^{U_z})`$. It is exact only as an inverse of the necessary scalar constraints; it is not a physical construction.

Its asymptotic exponent was

```math
E_q^{\rm old}(\rho)=
\max_{0\le\lambda\le\ln2/7}
\left[2\lambda-\rho\lambda\frac{\alpha_{e^\lambda}}
{1-\alpha_{e^\lambda}}\right]_+,
```

with the continuous zero value at $`\lambda=0`$. The improvement in (D10) is rigorous. Strict convexity of $`\phi_1`$, $`\phi_1(0)=0`$ and $`\phi_1(\lambda)=\alpha\lambda`$ imply

```math
\phi_j(\lambda)\le\lambda\alpha^j,\qquad
S_q(\lambda)\lt\lambda\frac\alpha{1-\alpha}
\quad(\lambda\gt0),
```

where strictness already holds at $`j=2`$. On the old parameter interval $`A(\lambda)=2\lambda`$. For every positive subcritical $`\rho`$, evaluating the new objective at a positive maximizer of the old one gives a strictly larger value. At $`\rho=0`$, the larger admissible parameter range instead raises $`2\ln2/7`$ to $`\ln2`$. At and above the common critical rate both guarantees are zero. Thus

```math
\mathcal E_q(\rho)\gt E_q^{\rm old}(\rho)
\quad\text{for }0\le\rho\lt2(1-q)/q.
\qquad\text{(D13)}
```

Finite bounds remain complementary because (D6) uses a schedule-independent cut average. For example, at $`k=12,n=120,q=99/100,D=2000`$:

| Prescribed $`m_T`$ | Full generating function, specified $`z`$ | Earlier scalar inverse with existing polygon |
|---:|---|---:|
| $`0`$ | $`z=2:\ p\le0.000244140789`$ | $`z=27/25:\ p\le0.000244168541`$ |
| $`171`$ | $`z=27/25:\ p\le0.252821100596`$ using (D8) | $`p\le0.248198005521`$ |
| $`200`$ | $`z=27/25:\ p\le0.316524376419`$ using (D8) | $`p\le0.310300082242`$ |
| $`243`$ | $`z=27/25:\ p\le0.427189188587`$ using (D8) | $`p\le0.417900524354`$ |

All numbers are outward upper bounds. Retain the strongest applicable finite inequality; an asymptotic improvement does not imply that every finite relaxation improves at each fixed parameter. The existing individual-gate ledger still excludes $`p=2/5`$ at $`m_T=243`$.

The finite verifier also retains the complete comparison of all 169 integer algebra pairs $`(s,c)`$ with $`r=2s+c\le24`$ and $`s+c\le n`$. It includes output dimension, dimension lifting $`p\le2^s\beta^k`$ with $`\beta=(1+1/\sqrt2)/2`$, the direct cap $`2^{(s-k)/6}`$, and the specified normalized DFW/Renes recovery caps. At $`z=27/25`$ the integer direct and combined hulls both equal the existing three-vertex polygon; the specified recovery hull has rank vertices $`0,13,16,19,24`$. For $`m_T=200`$ that recovery hull and the scalar inverse give $`p\le0.367329468646`$, compared with $`0.310300082242`$ from the direct polygon. This is a complete fixed-parameter comparison of those specified bound families, with no additional integer improvement there.

The historical $`223/52`$ ceilings remain scoped to their specified recovery-formula comparison at $`p=2/5`$. They do not determine a full decay exponent, every possible recovery refinement, or a physical upper processing cost. Recovery bounds may still be combined with either dynamical ledger. A fixed spectral cutoff justified for a nonvanishing score cannot simply be reused when the score decays exponentially.

## 6. Proof ownership, verification and limits

The branch profile, the conditional moment lemma, fixed primitive charges and CPTP factorization are inherited from the existing science. The addition is propagation of the full generating function without Jensen compression, its arbitrary-schedule cut ledger, and the stronger explicit decay exponent. The old finite inverse and its certified profile comparisons remain useful finite consequences.

The verifier retains 864 exact scalar inversion cases and all 169 finite algebra pairs, and adds 5,200 gate-column/window checks, 216 generating-function semigroup identities, 3,360 full-profile support checks, certified infinite sums and global optimizer brackets, and reproducible finite examples. These checks supplement the universal written proofs. Independent checks of the same transport and schedule steps are recorded in the private research handoff.

The common critical value is a necessary-rate threshold for exponential failure under the stated growing-family conditions. Neither exponent establishes processing capacity, the optimal numerical T coefficient, physical optimality, or practical constants. A score verifies none of the noise, controller or chronology premises. No protected workspace, hidden adaptive lookup, LWE preparation wrapper, hardware gate-time noise, locality condition or altered challenge distribution is introduced.
