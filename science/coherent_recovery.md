# Coherent recovery with an explicit gate budget

This construction complements the lower bounds in [Control cost](control_cost.md). It stores unknown quantum states using only the operations and noise schedule of the [memory model](model.md), with every nonlinear correction performed coherently and charged. It gives a constructive upper bound within a power of a logarithm of the long-time lower bound. It does not establish optimal $`T`$ count.

The ingredients are the standard five-qubit perfect code, concatenation, and an exact seven-$`T`$ implementation of a Toffoli gate; their attribution is in [Related work](../references/RELATED_WORK.md). The syndrome polynomials and resource accounting below make their use explicit.

## 1. A five-qubit recovery using 42 T gates

Use four commuting independent stabilizer generators

```math
g_1=XZZXI,\quad g_2=IXZZX,\quad
g_3=XIXZZ,\quad g_4=ZXIXZ,
```

and logical representatives

```math
\overline X=XXXXX,\qquad \overline Z=ZZZZZ.
```

A Clifford decoding unitary can simultaneously map the four generators to the $`Z`$ operators on four syndrome wires and map the logical representatives to $`X`$ and $`Z`$ on the remaining wire. This follows by completing the commuting generators and the anticommuting logical pair to a symplectic basis. Signs can be fixed by Pauli conjugations. This decoding is free in the stated gate budget.

For a physical Pauli error $`E`$, let $`s_i`$ be its commutation bit with $`g_i`$: zero for commuting and one for anticommuting. Under the decoder, the syndrome wires are the computational basis state $`s_1s_2s_3s_4`$. The residual Pauli on the logical wire has bits

```math
x(E)=[E,\overline Z]_{\mathrm{sp}},
\qquad
z(E)=[E,\overline X]_{\mathrm{sp}}.
```

For each single-site $`X`$, $`Y`$, or $`Z`$ error, these bits are respectively $`(1,0)`$, $`(1,1)`$, and $`(0,1)`$. The fifteen nontrivial single-site errors have distinct nonzero syndromes:

| Syndrome | Error | Logical X bit | Logical Z bit |
|---|---|---:|---:|
| 0000 | Identity | 0 | 0 |
| 1000 | X on qubit 2 | 1 | 0 |
| 0100 | Z on qubit 5 | 0 | 1 |
| 1100 | X on qubit 3 | 1 | 0 |
| 0010 | Z on qubit 3 | 0 | 1 |
| 1010 | Z on qubit 1 | 0 | 1 |
| 0110 | X on qubit 4 | 1 | 0 |
| 1110 | Y on qubit 3 | 1 | 1 |
| 0001 | X on qubit 1 | 1 | 0 |
| 1001 | Z on qubit 4 | 0 | 1 |
| 0101 | Z on qubit 2 | 0 | 1 |
| 1101 | Y on qubit 2 | 1 | 1 |
| 0011 | X on qubit 5 | 1 | 0 |
| 1011 | Y on qubit 1 | 1 | 1 |
| 0111 | Y on qubit 5 | 1 | 1 |
| 1111 | Y on qubit 4 | 1 | 1 |

These tables have the following algebraic normal forms over $`\mathbb F_2`$:

```math
f_X(s)=s_1\oplus s_4\oplus s_1s_3
       \oplus s_2s_3\oplus s_2s_4,
```

```math
f_Z(s)=s_2\oplus s_3\oplus s_1s_2
       \oplus s_1s_4\oplus s_3s_4.
```

Apply the conditional logical correction $`X^{f_X(s)}`$ and $`Z^{f_Z(s)}`$, then reset the four syndrome wires. The order changes only a syndrome-dependent phase, which does not change the recovered logical channel after those wires are discarded.

The two linear terms in each polynomial require CNOT or controlled-$`Z`$ gates. The three quadratic terms in $`f_X`$ require three Toffoli gates; the three in $`f_Z`$ require three controlled-controlled-$`Z`$ gates. A controlled-controlled-$`Z`$ is Clifford-equivalent to Toffoli. Thus an exact implementation uses at most

```math
6\times 7=42
```

$`T`$ or $`T^\dagger`$ gates. There are no additional work qubits: the four controls are the decoded syndrome wires, and the target is the decoded logical wire. All syndrome bits remain quantum registers during the operation. No external measurement record, classical lookup, or data-dependent address operation is supplied for free.

On every error of weight at most one this operation restores the logical input exactly, including entanglement with an external reference. For an arbitrary Pauli error it leaves a definite residual logical Pauli. That last property permits recursive analysis.

## 2. Concatenation and error suppression

Concatenate this code for $`L\ge0`$ levels. A block holding one logical qubit uses $`5^L`$ physical qubits. Decode and correct the level-one blocks, reset their syndrome wires, then decode and correct the next level, continuing to the root. Finally reencode the corrected root through the inverse Clifford encoding tree.

The number of recovery nodes in one top-level block is

```math
1+5+\cdots+5^{L-1}=\frac{5^L-1}{4}.
```

The operation uses only the original $`5^L`$ wires. The released syndrome wires become the clean inputs for reencoding. A global Clifford layer may route the surviving wires as needed. Gate duration and geometric locality are not charged in this model.

Suppose each physical qubit has an independent Pauli error with probability at most $`p_0`$. Different child blocks depend on disjoint physical-error samples. Their residual logical Paulis are therefore independent, even when the encoded logical state is entangled across blocks or with a reference. A parent block certainly succeeds if at most one child has a nonidentity residual Pauli. Consequently,

```math
p_{\ell+1}\le {5\choose2}p_\ell^2=10p_\ell^2,
\qquad
p_L\le \frac{(10p_0)^{2^L}}{10}.
```

The pair union bound is deliberately conservative. It neither assumes that every two-child error fails nor uses a fitted threshold.

For $`k`$ separately encoded logical qubits and $`J`$ recovery intervals, a union bound gives entanglement fidelity at least

```math
F_{\mathrm e}\ge1-kJp_L.
```

An earlier logical Pauli error does not change the probability estimate for a later interval: after each operation the system is back in the encoded space, and the next physical noise samples are independent. The good event of no logical error in any block or interval acts as the identity on the complete $`k`$-qubit logical state. Its probability is input independent. The recovered channel is therefore a convex combination of the identity channel and a CPTP bad-event channel, with bad-event probability at most $`kJp_L`$. Consequently the same bound also gives half-diamond error at most $`kJp_L`$.

## 3. Batching weak noise

Let $`0<\gamma\le1/20`$ and let storage last $`D\ge1`$ fully noisy rounds. Set

```math
\tau=\left\lfloor\frac1{20\gamma}\right\rfloor,
\qquad
J=\left\lceil\frac D\tau\right\rceil.
```

Perform coherent recovery after every $`\tau`$ noise rounds and after the shorter last interval. Since $`\tau\ge1/(40\gamma)`$,

```math
J\le40\gamma D+1.
```

An interval of length $`t\le\tau`$ has effective depolarizing strength

```math
\gamma_t=1-(1-\gamma)^t\le t\gamma\le\frac1{20}.
```

In the Pauli-mixture convention, its nonidentity error probability is $`p_0=3\gamma_t/4\le3/80`$. Therefore

```math
p_L\le\frac1{10}\left(\frac38\right)^{2^L}.
```

The recovery circuit is ideal within its processing interval, exactly as in the converse model. If a physical implementation spreads it over additional noisy layers, those layers and their error propagation require a different analysis.

## 4. Explicit resource guarantee

For $`k,D\ge1`$, $`0<\epsilon<1/2`$, and $`0<\gamma\le1/20`$, define

```math
A=\max\!\left\{1,\,
\frac{\ln\!\bigl(kJ/(10\epsilon)\bigr)}{\ln(8/3)}
\right\},
\qquad
L=\left\lceil\log_2 A\right\rceil.
```

Then $`2^L\ge A`$ and $`kJp_L\le\epsilon`$. This defines a memory with

```math
n=k5^L,
\qquad
m_T\le \frac{42kJ}{4}(5^L-1),
\qquad
\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon,\qquad F_{\mathrm e}\ge1-\epsilon.
```

For $`L=0`$ the recovery tree is empty and the gate count is zero; the same probability bound covers that case. Counting the recovery after the final noise interval makes the stated upper bound conservative. The free endpoint decoder can then be chosen Clifford.

With $`\alpha=\log_2 5`$, the integer choice obeys $`5^L\le5A^\alpha`$. Uniformly in the weak-noise range above, this gives

```math
n=O\!\left(
k\left[1+\ln\frac{k(\gamma D+1)}{\epsilon}\right]^\alpha
\right),
```

```math
m_T=O\!\left(
k(\gamma D+1)
\left[1+\ln\frac{k(\gamma D+1)}{\epsilon}\right]^\alpha
\right).
```

Constants do not depend on $`k,D,\gamma,\epsilon`$ within the stated ranges. The upper bound pays for error correction coherently and uses no protected workspace. Its logarithmic overhead is an achieved bound, not a lower bound on the optimum.

## 5. Relation to the converse

The [gate-budget theorem](control_cost.md) gives an extensive necessary cost proportional to $`k\gamma D`$ after its initial width allowance. This construction realizes that dependence with a power of a logarithm overhead in both width and gate count. For its width, the logarithmic width allowance in the simple converse is only of order $`\gamma^{-1}\ln\ln(k(\gamma D+1)/\epsilon)`$ at fixed fidelity.

A complete optimum would require closing the remaining logarithmic factors and controlling sharper dependence on fidelity and code rate. Neither ordinary classical decoder running time nor the use of a code with good distance automatically supplies such a coherent circuit bound.

The [zero-non-Clifford construction](clifford_lifetime.md) realizes the alternative resource choice: an exponentially large noisy record of syndrome history and nonlinear recovery deferred to the endpoint. These constructions serve different resource regimes.

## Verification

The executable [five-qubit checker](../verification/verify_five_qubit.py) checks the stabilizer commutation and rank, all sixteen correctable-error syndromes, the two Boolean polynomials, every Pauli error on a five-qubit block, and the exact logical-error enumerator of this recovery rule. It also checks the chosen integer batching parameters and the conservative recursion. These finite calculations support the explicit construction; the unbounded concatenation statement follows from the independence and induction arguments above.
