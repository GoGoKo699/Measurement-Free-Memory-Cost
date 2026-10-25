# A finite-depth processing budget for quantum memory

This chapter proves necessary gate budgets for the memory model in [model.md](model.md). The binary envelope and its noise contraction are established in [pauli_envelopes.md](pauli_envelopes.md). All gates charged here lie inside the chosen storage interval. An arbitrary fixed encoder and final decoder may remain free endpoints.

Write $`q=1-\gamma\in(0,1)`$, $`\kappa=-\ln q`$, and let $`D`$ be the number of full noise layers. The width bound $`n`$ includes every persistent physical register. A gate count is a necessary resource in this precise model, not an assertion that every implementation of quantum error correction requires a particular non-Clifford resource.

## 1. From a binary subspace to a quantum bottleneck

For $`V\le\mathbb F_2^{2n}`$, let $`r=\dim V`$. Its Pauli span $`\mathcal A_V`$ is a unital $`*`$-algebra: multiplication adds binary labels, up to a phase, and each Pauli is Hermitian. Let the symplectic form restricted to $`V`$ have rank $`2s`$, with radical dimension $`c=r-2s`$.

### Lemma C1 — Canonical algebra and conditional expectation

There is a Clifford change of coordinates under which


```math
\mathcal A_V
\simeq M_{2^s}\otimes\operatorname{Diag}_{2^c}
\otimes I_{2^{n-s-c}}.
\tag{C1}
```


The trace-preserving conditional expectation onto this algebra is


```math
\mathcal E_V(X)
=\frac{1}{|V^\perp|}\sum_{w\in V^\perp}P_wXP_w,
\qquad
V^\perp=\{w:[v,w]=0\text{ for all }v\in V\}.
\tag{C2}
```


It is CPTP, self-adjoint, and idempotent. It retains exactly the Pauli terms whose labels belong to $`V`$.

**Proof.** Symplectic Gram–Schmidt supplies $`s`$ anticommuting pairs and $`c`$ mutually commuting radical generators. Extend these to a symplectic basis of the ambient space. A Clifford maps the pairs to $`X_j,Z_j`$ on the first $`s`$ qubits and the radical generators to $`Z`$ on the next $`c`$ qubits. Their generated algebra has the form (C1). This uses the usual binary symplectic description of Clifford operations [Gottesman, 1998](https://arxiv.org/abs/quant-ph/9807006).

Equation (C2) is an average of unitary conjugations, so it is CPTP. For a Pauli $`P_v`$, its coefficient after twirling is the average of the character $`(-1)^{[v,w]}`$ over $`w\in V^\perp`$. This average is one if $`v\in(V^\perp)^\perp=V`$, and zero otherwise. Thus $`\mathcal E_V`$ is the orthogonal projection onto $`\mathcal A_V`$ in the Hilbert–Schmidt inner product. It is consequently self-adjoint and idempotent. $`\square`$

### Lemma C2 — Factorization and entanglement bound

Suppose a CPTP tail channel $`\Theta`$, including its fixed terminal decoder, obeys


```math
\operatorname{range}\Theta^\dagger\subseteq\mathcal A_V.
\tag{C3}
```


Then


```math
\Theta=\Theta\circ\mathcal E_V.
\tag{C4}
```


For every state shared between the memory and an arbitrary reference, the output of $`\Theta`$ has Schmidt number at most $`2^s`$ across that bipartition. In particular, its overlap with a maximally entangled state of Schmidt rank $`K=2^k`$ satisfies


```math
F\le\min\{1,2^{s-k}\}
\le\min\{1,2^{r/2-k}\}.
\tag{C5}
```


Here $`F=\langle\Omega_K|\rho|\Omega_K\rangle`$ is the squared fidelity with the pure target.

**Proof.** By (C3), $`\mathcal E_V\circ\Theta^\dagger=\Theta^\dagger`$. Taking adjoints and using self-adjointness of $`\mathcal E_V`$ gives (C4).

In the coordinates of (C1), write the memory as $`QCG`$, with $`Q`$ containing $`s`$ qubits, $`C`$ containing $`c`$ qubits, and $`G`$ the remaining factor. For every joint input $`\rho_{QCGR}`$,


```math
(\mathcal E_V\otimes\operatorname{id}_R)(\rho)
=\sum_{z\in\{0,1\}^c}
\rho^{(z)}_{QR}\otimes|z\rangle\langle z|_C
\otimes\frac{I_G}{\dim G},
\tag{C6}
```


up to tensor ordering, where each $`\rho^{(z)}_{QR}`$ is positive and their traces sum to one. Every $`QR`$ block admits a pure-state decomposition of Schmidt rank at most $`\dim Q=2^s`$. Tensoring the memory-side factors and mixing the blocks preserves that bound. A local CPTP map also preserves an upper bound on Schmidt number: apply each Kraus operator to the memory side of each pure-state decomposition, which cannot increase the corresponding Schmidt rank. These are the standard Schmidt-number properties introduced by [Terhal and Horodecki, 2000](https://arxiv.org/abs/quant-ph/9911117).

For completeness, any normalized pure state $`|\psi\rangle`$ of Schmidt rank at most $`S`$ has $`|\langle\Omega_K|\psi\rangle|^2\le S/K`$. In matrix form its coefficient matrix $`A`$ has Frobenius norm one and rank at most $`S`$, so $`|\operatorname{Tr}A|\le\|A\|_*\le\sqrt S`$. Dividing the squared overlap by $`K`$ gives the claim. The same inequality holds for a convex mixture of such states. Set $`S=2^s`$ and use $`2s\le r`$. $`\square`$

### Lemma C2a — Superoperator rank and Bell entanglement fidelity

Let $`\Lambda`$ be a CPTP channel from a $`K`$-dimensional logical system to a $`K`$-dimensional logical system. Regard it as a complex-linear superoperator on the $`K^2`$-dimensional operator space. Then


```math
F_{\mathrm e}(\Lambda)
=\frac{\operatorname{Tr}_{\mathrm{sup}}\Lambda}{K^2}
\le\frac{\operatorname{rank}_{\mathbb C}\Lambda}{K^2}.
\tag{C5a}
```


If $`\Lambda=\Theta\circ\Pi`$, where $`\Pi`$ is any CPTP prefix and the tail obeys (C3), then


```math
F_{\mathrm e}(\Lambda)
\le\min\{1,2^{r}/K^2\}.
\tag{C5b}
```


**Proof.** For Kraus operators $`A_\alpha`$ of $`\Lambda`$, its matrix in the vectorized operator representation is $`\sum_\alpha\overline{A_\alpha}\otimes A_\alpha`$. Therefore


```math
\operatorname{Tr}_{\mathrm{sup}}\Lambda
=\sum_\alpha|\operatorname{Tr}A_\alpha|^2
=K^2F_{\mathrm e}(\Lambda)\ge0.
\tag{C5c}
```


Every eigenvalue of a CPTP superoperator has modulus at most one. One elementary justification is that all powers of the channel are trace-norm contractions on Hermitian operators; decomposing an arbitrary operator into Hermitian real and imaginary parts makes those powers uniformly bounded on the full finite-dimensional operator space. An eigenvalue of modulus greater than one would contradict that bound.

The number of nonzero eigenvalues, counted with algebraic multiplicity, is at most the matrix rank. This remains true for a nondiagonalizable superoperator: its generalized eigenspaces for nonzero eigenvalues are subspaces on which the operator is invertible. Hence


```math
0\le\operatorname{Tr}_{\mathrm{sup}}\Lambda
\le\sum_{\lambda\ne0}|\lambda|
\le\operatorname{rank}_{\mathbb C}\Lambda.
\tag{C5d}
```


This proves (C5a). By (C3), $`\operatorname{rank}_{\mathbb C}\Theta\le\dim_{\mathbb C}\mathcal A_V=2^r`$. Composition with $`\Pi`$ cannot increase rank, which gives (C5b). $`\square`$

The superoperator bound is stronger than (C5) at small binary rank and gives the main rank floor below. The Schmidt-number statement remains complementary: for example, a large commuting algebra has $`s=0`$, so (C5) still limits Bell fidelity to $`1/K`$ even when its binary dimension makes (C5b) trivial. Unlike the state-uniform Schmidt-number bound, (C5b) is used on a complete logical CPTP channel with the maximally entangled input specified in its entanglement fidelity.

## 2. High fidelity enforces a rank floor at every cut

Encode $`k\ge1`$ logical qubits by an arbitrary CPTP encoder. Pass the physical register through the storage circuit and a fixed arbitrary decoder. Assume that the resulting channel has Bell entanglement fidelity


```math
F\ge1-\epsilon,
\qquad 0\le\epsilon<1-4^{-k}.
\tag{C7}
```


At each primitive circuit cut $`t`$, construct the sampled backward envelope $`V_{b,t}`$ of the entire decoded tail. Let


```math
R_t=\mathbb E_b\dim V_{b,t},
\qquad
a=\frac{2k(1-\epsilon-4^{-k})}{1-4^{-k}}>0.
\tag{C8}
```


### Theorem C3 — Uniform expected-rank floor

At every primitive cut, including cuts immediately before and after each gate or noise layer,


```math
R_t\ge a.
\tag{C9}
```


**Proof.** Fix a cut and any realization $`b`$ of the future noise. Average the earlier noise to obtain a CPTP prefix $`\Pi_t`$ from the logical input to the physical register at the cut. The Markov noise model makes this prefix independent of $`b`$. Let $`\Theta_{b,t}`$ be the sampled tail including its fixed decoder. Lemma E1 of [pauli_envelopes.md](pauli_envelopes.md) puts $`\operatorname{range}\Theta_{b,t}^\dagger`$ inside $`\mathcal A_{V_{b,t}}`$. Thus the full logical channel $`\Lambda_b=\Theta_{b,t}\circ\Pi_t`$ is CPTP and satisfies (C5b).

For every $`r\ge0`$,


```math
\min\{1,2^{r-2k}\}
\le 4^{-k}+(1-4^{-k})\frac{r}{2k}.
\tag{C10}
```


On $`[0,2k]`$, this is the secant-line bound for a convex exponential; for $`r\ge2k`$, its right side is at least one. Average (C5b) and (C10) over the sampled tail. Bell fidelity is linear in the channel, so


```math
1-\epsilon\le F
\le 4^{-k}+(1-4^{-k})\frac{R_t}{2k}.
\tag{C11}
```


Rearrangement proves (C9). $`\square`$

The lower bound applies to the expectation of the envelope rank. The prefix may be any CPTP encoder-and-evolution channel, including a noisy, nonisometric one. No special form of its actual state at the cut is required. The theorem does not require every sampled branch to retain rank $`a`$, and it does not identify the envelope rank with the actual amount of transmitted quantum information. A larger envelope only weakens the ensuing necessary bound.

## 3. Rank-growth charges for gates

A rank-growth charge $`b_g\ge0`$ for a gate is valid if an envelope can be updated across that gate with


```math
\dim V'\le\dim V+b_g
\tag{C12}
```


for every incoming envelope. This is a structural upper bound on support growth, not an assumed physical cost metric. The total charge is $`B=\sum_g b_g`$.

### Lemma C4 — T gates and arbitrary local gates

A physical $`T`$ or $`T^\dagger`$ gate has charge one. More generally, a fixed rotation about any Pauli axis has charge at most one. An arbitrary unitary acting on $`s_g`$ specified qubits has charge at most $`2s_g`$. Clifford unitaries and replacement resets have charge zero.

**Proof.** Let $`U=\exp(i\theta P_w)`$. A Pauli commuting with $`P_w`$ is unchanged under conjugation. An anticommuting Pauli maps into the linear span of $`P_v`$ and $`P_wP_v`$, whose binary labels are $`v`$ and $`v+w`$. Thus


```math
V'\subseteq V+\langle w\rangle,
\qquad \dim V'\le\dim V+1.
\tag{C13}
```


For $`T`$ and $`T^\dagger`$, $`w`$ is the single-qubit $`Z`$ label, up to an irrelevant global phase. A rotation about a large Pauli has the same envelope charge, but is not thereby assigned the physical implementation cost of one T gate.

For a unitary on a set $`J`$ of $`s_g`$ qubits, conjugation changes only the Pauli factor on $`J`$. The outside label is $`\pi_{J^c}v`$, and an arbitrary label can appear inside $`J`$. An envelope is therefore


```math
V'=\pi_{J^c}V+W_J,
\qquad
\dim V'\le\dim V+2s_g,
\tag{C14}
```


where $`W_J`$ is the $`2s_g`$-dimensional coordinate subspace supported on $`J`$. The zero charges follow from Lemma E1. $`\square`$

The local bound (C14) also holds for any fixed CPTP primitive supported on $`J`$, since the outside Pauli factor is unchanged. Such a primitive must still be represented honestly as part of the storage dynamics; a noiseless external controller acting on an uncounted memory is not a local primitive on $`J`$.

For Clifford+T circuits, $`B=m_T`$, the number of actual $`T`$ and $`T^\dagger`$ gates. For arbitrary local unitaries one may take $`B=2\sum_g s_g`$, summing only the charged gates, or use a smaller gate-specific charge when (C12) has been proved.

## 4. Logarithmic and additive resource laws

Initialize the decoded-tail envelope with dimension $`M\le2n`$, taking $`M=2n`$ when no decoder-specific information is used. The rank floor implies $`M\ge a`$; otherwise (C7) is impossible even before a noise layer is considered.

### Theorem C5 — General logarithmic budget

Every memory satisfying (C7) obeys


```math
B\ge a\left[D\kappa-\ln\!\left(\frac Ma\right)\right]_+,
\qquad [x]_+=\max\{0,x\}.
\tag{C15}
```


For a Clifford+T implementation, the stronger individual-gate version is


```math
m_T\ge
\frac{\left[D\kappa-\ln(M/a)\right]_+}
{\ln(1+1/a)}.
\tag{C16}
```


The right side may be rounded up because $`m_T`$ is an integer.

**Proof.** Traverse the circuit backward. At every cut $`R\ge a`$. A noise layer gives $`R'\le qR`$ by Lemma E2. A free operation gives $`R'\le R`$. A gate of charge $`b_g`$ gives $`R'\le R+b_g`$, hence


```math
\ln R'-\ln R
\le\ln(1+b_g/R)
\le b_g/a.
\tag{C17}
```


Summing the logarithmic increments over the entire backward circuit gives


```math
\ln(a/M)
\le\ln(R_{\mathrm{end}}/M)
\le D\ln q+B/a.
\tag{C18}
```


This proves (C15). When each charged gate is T or $`T^\dagger`$, replace its increment by the sharper $`\ln(1+1/a)`$ before summing, which gives (C16). $`\square`$

The dimension reserve $`M`$ contributes only the finite term $`\ln(M/a)`$. For $`M=2n`$, long storage needs a total charge of order $`k\gamma D`$ at weak noise and fixed sufficiently small error. The next theorem gives a stronger exact optimization of the scalar rank inequalities, including the finite initial reserve.

## 5. Exact finite-depth bound from discrete injections

Define, for an integer $`D\ge1`$, $`M\ge a>0`$, and $`q\in(0,1)`$,


```math
\mathcal L(D,M,a,q)
=\max\!\left\{0,
\max_{1\le\ell\le D}
\left[
\frac aq+(D-\ell)\frac{(1-q)a}{q}
-Mq^{\ell-1}
\right]\right\}.
\tag{C19}
```


Set $`\mathcal L(0,M,a,q)=0`$.

### Theorem C6 — Discrete processing budget

Every circuit satisfying (C7) and the rank-growth assumptions satisfies


```math
B\ge\mathcal L(D,M,a,q).
\tag{C20}
```


For Clifford+T circuits one may replace $`B`$ by $`m_T`$ and round the right side up. A simpler consequence is


```math
B\ge\left[\frac{(1-q)a}{q}D-(M-a)\right]_+.
\tag{C21}
```


The bound (C19) is exactly the optimum of the continuous scalar injection problem induced by the rank inequalities. It is not a claim of a matching physical T-count construction.

**Proof.** Number the noise layers in backward order, $`j=1,\ldots,D`$. Let $`b_j`$ be the total charge between the previous backward noise layer and layer $`j`$, with $`b_1`$ including all charged gates between the terminal cut and the first backward noise. Let $`R_j`$ be the expected envelope rank immediately after processing backward through layer $`j`$, and set $`R_0=M`$. Combining gate growth, free rank loss, and one-layer contraction gives


```math
R_j\le q(R_{j-1}+b_j),
\qquad R_j\ge a.
\tag{C22}
```


Charged gates after the last backward noise layer are included in $`B`$ but not in $`\sum_{j=1}^D b_j`$. This can only weaken the lower bound. Such gates cannot eliminate the fidelity bottleneck at the cut immediately after the last backward noise layer. Thus no assumption about gate placement at the forward input boundary is needed.

Define the dominating recurrence $`x_0=M`$, $`x_j=q(x_{j-1}+b_j)`$. Induction gives $`x_j\ge R_j`$, so feasibility requires


```math
\sum_{i=1}^j q^{j-i+1}b_i\ge a-Mq^j,
\qquad j=1,\ldots,D.
\tag{C23}
```


Minimize $`\sum_i b_i`$ over nonnegative real $`b_i`$ subject to (C23). Every physical circuit supplies a feasible point in this relaxation. For a chosen $`\ell\in\{1,\ldots,D\}`$, multiply the constraint at $`j<D`$ by $`\lambda_j=(1-q)/q`$ if $`j\ge\ell`$, and zero otherwise; multiply the constraint at $`j=D`$ by $`\lambda_D=1/q`$. The total coefficient of $`b_i`$ is


```math
\sum_{j=i}^D\lambda_jq^{j-i+1}
=\begin{cases}
q^{\ell-i},&i<\ell,\\
1,&i\ge\ell.
\end{cases}
\tag{C24}
```


Every coefficient is at most one. The weighted right sides telescope to


```math
\sum_{j=1}^D\lambda_j(a-Mq^j)
=\frac aq+(D-\ell)\frac{(1-q)a}{q}-Mq^{\ell-1}.
\tag{C25}
```


Therefore $`B\ge\sum_i b_i`$ is at least every expression in (C25) and at least zero. This proves (C20). Taking $`\ell=1`$ gives (C21).

To prove exactness for the scalar relaxation, choose the greedy injections


```math
b_j=\left[\frac aq-x_{j-1}\right]_+,
\qquad x_j=q(x_{j-1}+b_j).
\tag{C26}
```


Let $`t_0=\lfloor\ln(M/a)/\kappa\rfloor`$. If $`D\le t_0`$, the zero-injection trajectory $`x_j=Mq^j`$ remains at least $`a`$, so the optimum is zero. If $`D>t_0`$, the trajectory coasts through layer $`t_0`$, then uses


```math
b_{t_0+1}=\frac aq-Mq^{t_0},
\qquad
b_j=\frac{(1-q)a}{q}\quad(j>t_0+1).
\tag{C27}
```


Its cost equals the dual value (C25) at $`\ell=t_0+1`$, namely


```math
\mathcal L(D,M,a,q)
=\frac aq-Mq^{t_0}
+(D-t_0-1)\frac{(1-q)a}{q}
\quad(D>t_0).
\tag{C28}
```


The primal feasible point and matching dual certificate establish optimality. Exact equality $`Mq^{t_0}=a`$ causes no ambiguity: the first positive injection is still at layer $`t_0+1`$. $`\square`$

### Dimension cap and feasibility

The physical envelope rank never exceeds $`C=2n`$. Thus if $`D\ge1`$, the rank floor and the contraction at any noise layer require


```math
a\le qC=2nq.
\tag{C29}
```


If $`a>qC`$, the fidelity requirement is impossible within this model, regardless of the gate budget. When $`a\le M\le C`$ and $`a\le qC`$, the optimal scalar schedule (C26) respects the cap even immediately before each noise layer: during coasting its rank is at most $`M\le C`$, and during replenishment its pre-noise rank is $`a/q\le C`$. Consequently adding the cap to the scalar relaxation does not change (C19) in its feasible parameter regime.

This observation does not make the scalar schedule a physical construction. Actual envelope ranks are random integers, each injection arises from a realizable gate, and preservation of a large envelope is only a necessary condition for transmission. The theorem proves the cost required by these necessary conditions.

For fixed $`M,a,q`$, the necessary-bound function has the asymptotic slope


```math
\lim_{D\to\infty}\frac{\mathcal L(D,M,a,q)}{D}
=\frac{(1-q)a}{q}.
\tag{C30}
```


This describes the scalar lower bound, without asserting that arbitrarily long high-fidelity storage is physically possible at fixed width. Rare complete depolarization of the entire finite memory imposes an additional ultimate limit; see [interpretation.md](interpretation.md). Along a family with fixed $`a,q`$ and $`\ln M=o(D)`$, the same asymptotic lower-bound slope follows from (C28).

The coefficient $`(1-q)/q`$ is stronger than $`-\ln q`$. They agree to first order in $`\gamma`$ as $`\gamma\to0`$. Bounds (C16) and (C20) may both be applied; the integer-gate logarithmic refinement and the continuous injection optimization use different relaxations.

## 6. Every contiguous storage interval has its own budget

### Corollary C7 — Interval resource law

Choose any two primitive cuts and the contiguous interval between them. Let $`D_I`$ be the number of full noise layers in that interval and $`B_I`$ the sum of charges of gates inside it. Assign boundary gates according to the chosen cuts so that no gate is counted twice. Let $`M_I`$ be any deterministic upper bound on the expected decoded-tail envelope rank at the later cut. Then


```math
B_I\ge\mathcal L(D_I,M_I,a,q),
\qquad
B_I\ge a\left[D_I\kappa-\ln(M_I/a)\right]_+.
\tag{C31}
```


One may always take $`M_I=2n`$. If the physical width at the later cut is smaller and all its persistent registers are represented there, the corresponding twice-width bound is available.

**Proof.** The global fidelity guarantee gives the same rank floor $`a`$ at every cut of the selected interval by Theorem C3. Start backward propagation at its later cut. The entering envelope can depend on noise and gates in the later suffix, but its expected rank is at most $`M_I`$, and noise inside the interval is independent of that suffix. The recurrence (C22) and logarithmic argument (C17) therefore hold unchanged. The earlier prefix can be any CPTP channel; the proof of Theorem C3 imposes no further condition on it or on the actual state it produces. $`\square`$

If only an upper bound $`M_I`$ is used, the dominating recurrence starts at that upper bound rather than the exact expectation; this preserves every inequality. The interval theorem shows that an arbitrarily long interior segment cannot be supported solely by gates outside that segment. The allowed initial reserve for each interval is still finite, at most $`2n`$.

## 7. The final decoder and the meaning of the terminal rank

For a fixed decoder $`\mathcal D`$, its Pauli-support envelope $`V_{\mathrm{dec}}`$ is defined by (E5). Starting at $`M=\dim V_{\mathrm{dec}}`$ is valid because every pulled-back terminal observable already lies in $`\mathcal A_{V_{\mathrm{dec}}}`$. The decoder is not subsequently charged as an internal operation. Its complete action is included in the initial support condition.

If $`\mathcal D`$ is a Clifford unitary followed by discarding all but $`k`$ physical output qubits, then the $`2k`$ independent output Pauli generators pull back to $`2k`$ independent Pauli labels. Hence $`M=2k`$. A decoder built solely from the free primitives after selecting $`k`$ output qubits has $`M\le2k`$, by the same backward rank argument. In that case the logarithmic reserve is $`\ln(2k/a)`$, independent of the physical storage width.

An arbitrary free decoder does not admit this simplification. For a concrete example, measure three physical qubits in the computational basis, compute their majority, and prepare that classical output bit. With $`Z`$ eigenvalues $`z_j\in\{\pm1\}`$, the majority observable is


```math
\mathcal D^\dagger(Z)
=\frac12(Z_1+Z_2+Z_3-Z_1Z_2Z_3).
\tag{C32}
```


Its Pauli-support span has dimension three, although the output is only one qubit. Even a unitary decoder followed by discarding one qubit can have the full two-qubit support rank four. Take


```math
U=\operatorname{CNOT}_{1\to2}(I\otimes THT),
\qquad \mathcal D(\rho)=\operatorname{Tr}_2(U\rho U^\dagger).
\tag{C33}
```


Then $`\mathcal D^\dagger(Z)=Z_1`$, while


```math
\mathcal D^\dagger(X)
=X_1\left(\frac{Z_2}{\sqrt2}+\frac{X_2+Y_2}{2}\right).
\tag{C34}
```


The four displayed nonidentity Pauli labels span $`\mathbb F_2^4`$. Thus one cannot set $`M=2k`$ merely from a $`k`$-qubit decoder output. The dimension of its output operator space is not a bound on the binary span of its Pauli supports.

The conservative theorem uses $`M=2n`$ and permits an unrestricted fixed final decoder. This is the relevant formulation for a decoder that performs nonlinear processing of accumulated syndrome records only after storage ends. If a construction instead uses a free final Clifford decoder, then all preceding nonlinear coherent recoveries remain inside the charged interval and must be included in $`B`$. Moving a nonlinear terminal recovery into the free decoder is allowed only if its effect is included when computing $`V_{\mathrm{dec}}`$. Both conventions are consistent; they define different storage boundaries.

Arbitrary final local recovery cannot invalidate the Schmidt-number bottleneck. A protected side channel carrying part of the logical input would invalidate the assumed model, because the state reaching the output would then have a component that did not pass through the cut algebra.

## 8. Interpretation and verification boundary

These theorems establish a necessary processing budget beyond the free Clifford-and-replacement primitives. Their strength is that resets may inject arbitrary fixed magic states, the encoder and final decoder may be arbitrary, and logical entanglement with a reference is retained in the analysis. The obstruction is the loss of the algebra capable of carrying a sufficiently high-fidelity logical channel through every temporal cut. Its complex dimension supplies the main fidelity bound; its coherent block size gives the complementary Schmidt-number restriction.

The theorems do not establish a physically optimal T-count, a universal gate-depth threshold, or a lower bound on measurement-based QEC with a noiseless external classical decoder. A measured controller must be included in a separate operational model. A coherent implementation of its nonlinear processing can be covered by charging its actual gate decomposition or a proved local rank-growth budget.

The program [verification/verify_pauli_algebra.py](../verification/verify_pauli_algebra.py) checks the finite Pauli-algebra identities by dense matrices, verifies the superoperator fidelity identity and rank bound for explicit logical CPTP channels, and verifies primal/dual equality of (C19)–(C28) over an exact-rational parameter grid. These tests address support propagation, the CPTP conditional expectation, the terminal-decoder distinction, the logical-channel bottleneck, and the discrete placement optimization. They do not replace the analytic proof or claim asymptotic evidence from small systems.
