# Uniform width and gate bounds at weak noise

The [storage model](model.md) counts every persistent qubit in the physical width and every prescribed interior T or T-dagger slot. Fix a target error ceiling $`0<\epsilon_0<1/2`$. The width and gate converses use recovered Bell fidelity $`F_{\mathrm e}\ge1-\epsilon`$, and the construction gives half-diamond error at most $`\epsilon`$. Throughout, $`k,D\ge1`$ are integers, $`0<\epsilon\le\epsilon_0`$ and $`0<\gamma\le1/20`$. Write $`q=1-\gamma`$.

The minimum width, allowing unrestricted gate count, is

```math
n_{\min}=\Theta_{\epsilon_0}\!\left(k+\frac{\ln(D/\epsilon)}{\ln(1/\gamma)}\right).
```

The achievable width constant also depends on the fixed Tanner family. At long-batch width proportional to $`k+\ln((1+\gamma D)/\epsilon)`$, the optimal gate order is $`\gamma D[k+\ln((1+\gamma D)/\epsilon)]`$ for sufficiently large $`\gamma D`$. The constants in these comparisons are uniform in $`\gamma`$. The two width scales can differ when $`\gamma`$ tends to zero.

## 1. Minimum width and achievable batch sizes

Factoring the recovered channel through $`n`$ qubits bounds its Bell fidelity by $`2^{n-k}`$ when $`n\le k`$, so $`F_{\mathrm e}>1/2`$ requires $`n\ge k`$. At each round all $`n`$ slots are independently replaced with probability $`\gamma^n`$. A trajectory containing a complete replacement has logical Bell fidelity $`4^{-k}`$. The probability of at least one complete replacement must therefore be at most $`\delta=\epsilon/(1-4^{-k})`$. This gives

```math
1-(1-\gamma^n)^D\le\delta,
\qquad D\gamma^n\le-\ln(1-\delta)\le4\epsilon.
```

The last inequality uses $`\delta\le4\epsilon/3<2/3`$ and $`-\ln(1-\delta)\le\delta/(1-\delta)`$. The resulting lower bound $`n\ge\ln(D/(4\epsilon))/\ln(1/\gamma)`$ combines with $`n\ge k`$ to give the uniform order above. The argument applies to any prescribed gate count and to any endpoint decoder. For $`D=0`$ the minimum is exactly $`k`$ and zero interior gates suffice.

For the upper bound, choose any integer batch length $`1\le h\le\min\{D,\lfloor1/(20\gamma)\rfloor\}`$, and write $`J=\lceil D/h\rceil`$ and $`\tau=\gamma h`$. Use the [coherent Tanner circuit with fixed inner protection](tanner_memory.md#8-constant-inner-concatenation-extends-the-physical-noise-range). Every epoch uses the same fixed circuit, with independent input-free resets before syndrome and work registers are reused. All of these registers occupy counted noisy slots between epochs.

The total fresh nonidentity probability on an inner physical site per batch is at most $`3\tau/4`$. After inner correction, conditional on the earlier history, the fresh residual errors on outer positions are independent. Each outer position has a fresh nonidentity error probability at most

```math
\pi(\tau)=\frac{((15/2)\tau)^{2^{L_*}}}{10}\le u/2,
```

where $`u`$ is the fixed outer decoder threshold and $`L_*`$ the fixed inner concatenation level. A CSS component exceeds weight $`uN`$ with probability at most $`\exp[-N\mathrm{KL}(u\Vert\pi(\tau))]`$. Here $`N`$ is the Tanner block length, and $`\mathrm{KL}`$ is the binary relative entropy using natural logarithms. The two components can be correlated at the same site, so the proof uses their two marginal tails and a union bound. The ratio $`\mathrm{KL}(u\Vert\pi(\tau))/\ln(1/\tau)`$ is positive and continuous over $`0<\tau\le1/20`$, with limit $`u2^{L_*}>0`$ at zero. It has a positive lower bound depending only on the fixed family.

Take the first available block length above the payload and tail thresholds. The [bounded ratios of available lengths](tanner_memory.md#7-available-lengths-padding-and-diamond-error) let us choose

```math
N=O\!\left(k+\frac{\ln(2J/\epsilon)}{\ln(1/\tau)}\right)
```

with total bad-trajectory probability at most $`\epsilon`$. The fixed residual invariant and final recovery give the identity on every good trajectory, including an attached reference. Any remaining trajectory contributes at most its probability to half-diamond error. The coherent circuit requires $`O(N)`$ total noisy slots and $`O(N)`$ prescribed gates per recovery. The last recovery can form part of the unrestricted endpoint decoder. Thus this family achieves

```math
n=O\!\left(k+\frac{\ln(2J/\epsilon)}{\ln(1/\tau)}\right),
\qquad
m_T=O\!\left((J-1)\left[k+\frac{\ln(2J/\epsilon)}{\ln(1/\tau)}\right]\right).
```

Taking $`h=1`$ proves the uniform width upper bound with $`O((D-1)n)`$ prescribed gates using the allocated width. For fixed $`k,\epsilon`$ and $`\gamma=D^{-1/2}`$, this choice uses bounded width as $`D`$ grows. Longer batches reduce the recovery frequency and can increase the width bound for the same target error.

## 2. Extinction over a block of noise rounds

Let $`r`$ be the rank of the backward Pauli envelope at an entry cut. The zero envelope is absorbing. A T or T-dagger gate may add at most one binary generator. Clifford gates preserve rank and replacement resets cannot increase it. The [pivot-column noise estimate](pauli_envelopes.md#lemma-e2--one-noise-layer-contracts-expected-rank) bounds the rank after a noise layer by the number of surviving pivot coordinates. For $`0\le z\le1`$ it gives

```math
\mathbb E[z^{r'}\mid V]\ge(1-q+qz)^r.
```

Indeed, a site that holds $`d`$ of the $`r`$ pivot coordinates, with $`d\in\{1,2\}`$, contributes $`1-q+qz^d\ge(1-q+qz)^d`$ to the generating function for the number of surviving pivots. Site flags are independent conditional on the incoming envelope. The inequality extends to $`z=0`$ with $`0^0=1`$. Across a charged gate the corresponding lower bound is $`z\,z^r`$.

Group $`H=\lceil1/\gamma\rceil`$ successive noise layers into one backward block. Assign each intervening gate to its following noise layer in backward time. Let $`M_j`$ be the charged gate count in block $`j`$, and let $`\ell_g`$ be the number of noise layers in that block preceding gate $`g`$. Let $`G_t`$ be the gates in this block with $`\ell_g<t`$. Iterating the one-step generating function bounds at $`z=0`$ gives, just after layer $`t`$,

```math
\Pr(r_t=0\mid V_{\rm entry})
\ge(1-q^t)^r\prod_{g\in G_t}
       (1-q^{t-\ell_g}).
```

Choose a target among $`t=\lceil H/2\rceil,\ldots,H`$. Its entry-rank factor is at least $`e^{-r}`$ because $`1-q^t\ge1-e^{-1/2}>e^{-1}`$. A gate contributes zero negative-log penalty on earlier targets. Its mean penalty is at most

```math
\frac{2}{H}\sum_{\ell=1}^H-\ln(1-q^\ell)
\le 2[1-\ln(1-e^{-1})]<3.
```

Here $`1-q^\ell\ge1-e^{-\ell/H}\ge(1-e^{-1})\ell/H`$ and $`\sum_{\ell=1}^H\ln(H/\ell)\le H`$. Some target therefore incurs total gate penalty at most $`3M_j`$. Absorption then gives a conditional probability at least $`e^{-r-3M_j}`$ that the envelope is zero by the end of the block. The argument needs no uniform spacing of the gates.

## 3. Gate lower bound and matching scale

Put $`B=\lfloor D/H\rfloor`$ and $`b=1-\delta`$. Complete extinction makes the logical channel a replacement on that trajectory, so its total probability is at most $`\delta`$. For any integer $`0\le S<B`$, discard the first $`S`$ complete backward blocks from the survival calculation. If $`R_j`$ is expected rank at block $`j`$ entry, then $`R_{j+1}\le q^H R_j+M_j`$ and $`R_1\le 2n`$. Since $`q^H\le e^{-1}`$,

```math
\sum_{j=S+1}^B R_j
\le \frac{2ne^{-S}+m_T}{1-q^H}
\le 4ne^{-S}+2m_T .
```

Skipped gates contribute to the new rank reserve, and the upper bound counts each gate once in the numerator of the sum. Write $`A=\sum_{j=S+1}^B\Pr(r_j>0)\ge b(B-S)`$. First extinction from a surviving block entry gives disjoint events. Applying Jensen to the conditional hazards, over remaining blocks and surviving entry histories, gives

```math
\delta\ge A \exp[-(4ne^{-S}+5m_T)/A]
\ge b(B-S)\exp[-(4ne^{-S}+5m_T)/(b(B-S))].
```

The second inequality uses the monotonicity of $`A e^{-v/A}`$ for fixed $`v\ge0`$. Consequently

```math
m_T\ge\frac{1}{5}\left[b(B-S)\ln\frac{b(B-S)}{\delta}-4ne^{-S}\right]_+ .
```

The finite bound applies when $`B\ge1`$. It also applies to unconditional delayed-basis BB84 score $`2^{-k}<p<1`$ in the same independent product-noise model with $`\delta=(1-p)/(1-2^{-k})`$, since complete extinction limits that branch's score to chance. Perfect score is impossible at finite width for $`D\ge1`$ because complete physical replacement has positive probability. The [one-layer reliability bounds](reliability_processing.md) remain available for finite comparisons and use more of the commutation structure.

Write $`x=\gamma D`$ and $`L=\ln((1+x)/\epsilon)`$. At width $`n\le C(k+L)`$ for fixed $`C>0`$, the new bound and the [payload law](control_cost.md#4-logarithmic-and-additive-resource-laws) together imply

```math
m_T=\Omega_{C,\epsilon_0}\bigl(x(k+L)\bigr)
```

for sufficiently large $`x`$, with a threshold depending only on $`C,\epsilon_0`$. If $`k\le L`$, use $`S=0`$. One has $`B\ge20x/21-1`$, $`b`$ bounded away from zero, and $`\ln(bB/\delta)\ge L-O_{\epsilon_0}(1)`$. The reliability term is of order $`xL`$ and the width reserve is $`O_C(L)`$. If $`k>L`$, the payload floor $`a=2kb`$ demands $`m_T\ge a[D(-\ln q)-\ln(2n/a)]_+`$. Here $`D(-\ln q)\ge x`$ and the width logarithm is bounded, giving order $`xk`$. The longest allowed batch in Section 1 has $`\tau`$ bounded below by a positive constant for large $`x`$. It achieves the same gate order within a sufficiently large fixed width factor. At fixed positive noise this recovers the fixed-noise matching result.

The same gate lower order holds even for $`n\le e^{x/8}(k+L)`$ once $`x`$ is sufficiently large. In the case $`k\le L`$, use $`S=\lfloor B/2\rfloor`$ to make $`ne^{-S}=O(L)`$. If $`k>L`$, the payload width logarithm is at most $`x/8+O_{\epsilon_0}(1)`$. The complete width–gate frontier at a width near the uniform minimum, where the reliability term can be smaller than $`L`$, remains open.

## 4. Relation to prior work

The code, residual-error invariant and coherent recovery circuits are the existing [Tanner construction and its references](tanner_memory.md#9-attribution-interpretation-and-finite-checks). The variable-batch upper bound keeps the noise-dependent Chernoff exponent. The gate lower bound uses the stochastic Pauli representation of Nelson, Rajakumar and Gullans and the complete-erasure survival mechanism adapted in the [one-layer reliability converse](reliability_processing.md). Grouping rounds and averaging a target cut give a lower bound uniform in weak noise under the stated fixed circuit model.
