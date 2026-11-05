# Coherent Tanner memory with linear space and processing per epoch

## 1. Result and scope

Existing good quantum Tanner codes and their constant-iteration parallel decoder give a stronger sufficient resource bound in the [same memory model](model.md). This chapter proves the fixed-circuit compilation, accounts for its total width and prescribed T count, and combines it with the published residual-error invariant. The addition is a resource consequence of inherited codes, decoders and measurement-free correction mechanisms.

Fix an admissible Tanner family as specified below. There is a family-dependent constant $`\gamma_*\gt0`$ such that, for $`0\lt\gamma\le\gamma_*`$, $`k\ge1`$, $`D\ge0`$ and $`0\lt\epsilon\lt1/2`$, a memory exists with

```math
n=O\!\left(k+\ln\frac{J_+}{\epsilon}\right),
\qquad m_T=O(nJ),
\qquad J=O(\gamma D+1),
\qquad J_+=\max\{1,J\},
```

and

```math
\frac12\|\Lambda-\mathrm{id}\|_\diamond\le\epsilon.
\qquad\text{(T1)}
```

Here $`n`$ includes data, syndrome, controller and all temporary work qubits. The constants are independent of $`k,D,\epsilon,\gamma`$ in the stated range. They can be very large, and $`\gamma_*`$ can be very small. This noise range is not the full range of the [five-qubit construction](coherent_recovery.md).

When $`k=\Omega(\ln(J_+/\epsilon))`$ and $`\gamma D`$ is sufficiently long for the converse's initial-width allowance, (T1) gives $`n=O(k)`$ and $`m_T=O(k\gamma D)`$. It closes the earlier polylogarithmic overhead in asymptotic order in that regime. It does not determine the optimal coefficient, practical constants, locality or a threshold with faulty processing gates.

## 2. Inherited code and decoder statements

[Gu, Tang, Caha, Choe, He and Kubica](https://doi.org/10.1007/s00220-024-04951-6), [arXiv:2306.12470v2](https://arxiv.org/abs/2306.12470v2), supply the code and decoding ingredients. Theorem numbers below refer to that inspected preprint; the corresponding statements must be identified by content when another version is used.

Write $`N`$ for Tanner data-block length, distinct from total noisy width $`n`$. Fix local-code parameters $`\Delta,\rho,d_r,\kappa_c`$ satisfying their Theorem 2.3, with $`\rho\ne1/2`$, and put

```math
b=\Delta^2,
\qquad R_0=(1-2\rho)^2\gt0,
\qquad d_0=\frac{d_r^2\kappa_c^2}{256\Delta}.
\qquad\text{(T2)}
```

The code has at least $`R_0N`$ logical qubits and distance $`d\ge d_0N`$. Each local view has $`b`$ physical positions. Each physical position lies in exactly one view of each of the four vertex types $`V_{ij}`$. Local check generators have bounded weight and bounded incidence degree. Admissible fixed local codes are an inherited existence result; no concrete practical generator matrices or useful numerical threshold are asserted here. They can be chosen once and compiled into the family of fixed circuits.

The notation $`\kappa_c`$ denotes local product expansion and differs from $`\kappa=-\ln q`$ elsewhere in this repository. Choose $`0\lt\delta_c\lt1/18`$. Theorem 4.20 states, for each CSS component, that

```math
A|e|_R+B|d_s|_V\le\min\{C_cN,d/\Delta\}
\quad\Longrightarrow\quad
|e+\widehat f^{(t)}|_R\le a_t|e|_R+\beta|d_s|_V,
\qquad\text{(T3)}
```

where $`|e|_R`$ is minimum weight modulo the corresponding stabilizer subgroup, $`d_s`$ is syndrome error, and

```math
A=\frac{48}{\kappa_c\Delta},\qquad
B=\frac{6\Delta}{\kappa_c},\qquad
C_c=\frac{d_r^2\delta_c^3\kappa_c}{4096\Delta^2},
```

```math
a_t=\frac{24}{5\kappa_c}(1-g)^t,
\qquad g=\frac{1-18\delta_c}{16},
\qquad \beta=\frac{6\Delta^2}{\kappa_c\delta_c}.
```

Choose once and for all an integer $`t\ge1`$ with $`a:=a_t\le1/2`$. Algorithms 3–4 run for exactly $`t`$ iterations and use only the current syndrome. The repeated-round invariant is also established in their Theorem 3.5; Section 5 below specializes and checks it with $`d_s=0`$, including endpoint correctness.

Perfect current syndrome extraction is allowed here because all processing between noisy layers is ideal. Scratch is reset immediately before extraction. Thus no inherited threshold for noisy extraction gates or measurement errors is needed, and none is silently claimed.

## 3. From the parallel decoder to a fixed Boolean circuit

This is a circuit-size argument, not an identification of RAM runtime with size. Fix an ordering of local physical positions, generators, codewords and decompositions. Resolve every tie lexicographically in these fixed orders.

The initialization at each vertex in $`V_1`$ maps at most $`b`$ local syndrome bits to a minimum-weight compatible correction, or to zero if none exists. This is a finite truth table.

In each vertex-type substep of Algorithm 4, the local input is the current $`b`$-bit mismatch on $`Q(v)`$. Eligibility of a local word $`x`$ is exactly

```math
4|Z\cap x|\ge3|x|,
```

because $`|Z|-|Z+x|=2|Z\cap x|-|x|`$. Choose an eligible nonzero word of maximum weight, or zero if none exists. Then choose $`x=c+r`$ minimizing the specified column/row count, with fixed ties. This is a map from $`b`$ bits to $`2b`$ bits. Its truth table includes local minimizations, comparisons, tie-breaking and the no-update case.

A map of $`s\le b`$ bits is implemented by computing all $`2^s`$ minterms and XORing the appropriate minterms to each output. Negations are NOT gates. With at most $`2b`$ outputs, safe bounds for each table are

```math
A_b=(b-1)2^b
\quad\text{binary AND gates},
\qquad S_b=4b2^b+4b
\quad\text{total Boolean gates}.
```

Thus lookup is physically compiled and charged; it is not an available free primitive. The constants are exponential in the fixed local length.

There are $`2N/b`$ initialization tables and $`4tN/b`$ update tables per CSS component. For both components together,

```math
L_N=\frac{(4+8t)N}{b},
\qquad A_N\le L_NA_b,
\qquad S_N\le L_NS_b+32(1+t)N.
\qquad\text{(T4)}
```

The last term safely covers the XORs and copies for mismatch initialization, accumulated corrections, successive array updates and final correction masks. Fresh records for successive arrays give a Boolean directed acyclic circuit of this size. Within one vertex type the local views are disjoint, so updates have no write conflicts. The four types are processed in a fixed order and the loop is unrolled a fixed number of times.

All inter-view connections are hardwired incidence wires. They have no data-dependent addresses. Wire permutations and fanout are Clifford operations in the global control model. The circuit description and local tables are compiled choices, rather than persistent external controller data. No geometric-routing or hardware-depth conclusion follows.

## 4. Coherent compilation, reset and exact channel

Reset the scratch pool to a fixed zero state immediately before each recovery. Extract the commuting stabilizer eigenvalues coherently using Clifford gates. A local basis has at most $`b`$ generators, so both CSS components together require at most $`4N`$ syndrome qubits, including global redundancy.

Evaluate the Boolean circuit reversibly into fresh work bits while retaining its syndrome input. An AND is one Toffoli; XOR and NOT are Clifford. Apply the two output masks to the data by controlled Pauli gates, then reverse the Boolean evaluation. A standard exact Toffoli uses seven T/T-dagger slots, giving

```math
T_{\mathrm{epoch}}\le14A_N
\le14\frac{(4+8t)(b-1)2^b}{b}\,N.
\qquad\text{(T5)}
```

This counts both computation and uncomputation. Every slot is prescribed, including for zero syndrome or no successful update. It is not an expected executed cost.

All work returns to zero. Reset the syndrome pool on the same fixed schedule after every epoch. A conservative width bound, including correction masks and reversible gate work, is

```math
n\le C_WN,\qquad C_W=64(1+t)b2^b.
```

One may also use $`C_T=112(1+t)b2^b`$ in $`T_{\mathrm{epoch}}\le C_TN`$. These loose constants emphasize the asymptotic scope of the result. All idle scratch remains in the counted noisy pool between epochs; resetting it before use discards its noise and retains no earlier information.

To verify the channel, let $`\Pi_s`$ denote the joint syndrome projector and $`P_{f(s)}`$ the deterministic Pauli correction. Extraction, reversible evaluation, correction and inverse evaluation send an arbitrary pure data input to

```math
\sum_s P_{f(s)}\Pi_s|\psi\rangle
\otimes|s\rangle\otimes|0_{\mathrm{work}}\rangle.
```

The original syndrome record is unchanged by correcting the data, so inverse Boolean evaluation remains valid. Resetting that record implements exactly

```math
\mathcal R(\omega)=
\sum_sP_{f(s)}\Pi_s\omega\Pi_sP_{f(s)}^\dagger.
\qquad\text{(T6)}
```

This equality holds for arbitrary operators and therefore after tensoring with any untouched reference. The reset implements the intended syndrome-measurement-and-correction channel; it does not assume coherent superpositions of different syndromes survive. On an encoded input following a Pauli trajectory, the syndrome is definite and carries no unknown logical amplitude. No syndrome, decision or mathematical noise flag is retained or supplied to the controller between epochs.

## 5. A correctable repeated-round invariant

Set

```math
\theta=\min\{C_c,d_0/\Delta\},
\qquad u=\min\left\{\frac{\theta}{4A},\frac{d_0}{8},\frac18\right\}\gt0,
\qquad v=\frac{au}{1-a}\le u.
\qquad\text{(T7)}
```

Suppose each fresh batch error has each CSS component of weight at most $`uN`$, and the previous residual components have reduced weight at most $`vN`$. The next incoming reduced weight is at most $`(v+u)N\le2uN`$, which satisfies (T3) with $`d_s=0`$. Its new residual is at most

```math
a(v+u)N=vN.
```

The invariant starts at zero and holds at every epoch. Minimum-weight representatives are only proof devices: stabilizer-equivalent errors have the same syndrome, and the deterministic decoder returns the same correction, preserving equivalence of the residuals. No such representative is physically computed.

The two component representatives have combined Pauli support at most

```math
2vN\le2uN\le d_0N/4\lt d/2.
```

One fixed endpoint decoder therefore corrects every trajectory satisfying the invariant, including logical states entangled with a reference. This uses the code's ordinary distance guarantee and ample strict slack. The endpoint may use arbitrary computation after charged storage has ended.

## 6. Independent depolarizing batches and reliability

A batch of $`h`$ noise layers has product-depolarizing retention $`q^h`$. Each data position's X component and Z component separately have probability

```math
p_h=\frac{1-q^h}{2}.
```

They are independent across positions. They are correlated with each other at a position. Only the two marginal binomial tails and a union bound will be used.

For the proved noise range, choose

```math
0\lt\gamma\le\gamma_*:=u/2,
\qquad h=\left\lfloor\frac{u}{\gamma}\right\rfloor,
\qquad J=\left\lceil\frac D h\right\rceil.
\qquad\text{(T8)}
```

Then $`h\ge u/(2\gamma)`$, and Bernoulli's inequality gives $`1-q^h\le h\gamma\le u`$. Thus $`p_h\le u/2`$, also for a shorter final batch. Let

```math
c=\mathrm{KL}(u\|u/2)
=u\ln2+(1-u)\ln\frac{1-u}{1-u/2}\gt0.
```

Chernoff's bound gives probability at most $`e^{-cN}`$ that a specified CSS component exceeds $`uN`$ in an epoch. For both components and all epochs,

```math
P_{\mathrm{bad}}\le2Je^{-cN}.
\qquad\text{(T9)}
```

This event concerns fresh noise, so dependence of residual errors on prior rounds causes no difficulty. No independence between X and Z components or between final residuals is assumed. A less sharp explicit sufficient choice may use the Hoeffding exponent $`u^2/2`$ in place of $`c`$.

## 7. Available lengths, padding and diamond error

The needed length family has bounded multiplicative gaps. Section 2.3 of Gu et al. uses fixed-degree Ramanujan Cayley graphs on $`\mathrm{PSL}_2(Q^i)`$ with fixed odd prime power $`Q`$ and $`\Delta=Q+1`$. [Dinur, Evra, Livne, Lubotzky and Mozes, arXiv:2111.04808](https://arxiv.org/abs/2111.04808), Claim 6.7, gives these for every integer $`i\gt2`$. With the fourfold-cover face convention,

```math
N_i=\frac{\Delta^2Q^i(Q^{2i}-1)}2,
\qquad
\frac{N_{i+1}}{N_i}
=Q^3\frac{1-Q^{-2i-2}}{1-Q^{-2i}}
\lt2Q^3.
```

Choose the first available length satisfying

```math
N\ge\max\left\{N_0,\frac{k}{R_0},
\frac1c\ln\frac{2J_+}{\epsilon}\right\}.
\qquad\text{(T10)}
```

Here $`N_0`$ is a fixed family starting length. Bounded gaps imply $`N=O(k+\ln(J_+/\epsilon))`$. The free encoder pads unused logical axes with fixed states; the free final decoder returns only the desired $`k`$ logical qubits. All physical padding remains counted.

On every good Pauli trajectory the same endpoint decoder is exactly the identity logical channel. Trajectory probabilities are independent of the logical input. Thus, for a CPTP channel $`\Xi`$,

```math
\Lambda=(1-P_{\mathrm{bad}})\mathrm{id}
+P_{\mathrm{bad}}\Xi,
\qquad
\frac12\|\Lambda-\mathrm{id}\|_\diamond
\le P_{\mathrm{bad}}\le\epsilon.
```

The decomposition proves reference-stable worst-case correctness, rather than only averaged codeword fidelity. Finally,

```math
J\le\frac{2\gamma D}{u}+1,
\qquad m_T\le C_TNJ,
\qquad n\le C_WN,
```

proving (T1). For $`D=0`$, take $`J=0`$ and $`J_+=1`$; no interior recovery is needed.

## 8. Attribution, interpretation and finite checks

Good Tanner codes, the parallel decoder, its contraction and its repeated-round use are inherited. Coherent syndrome extraction followed by measurement-free nonlinear correction is also established prior work, including [Hong, Guo and Lucas](https://doi.org/10.1038/s41467-024-55570-7), who study measurement-free error correction and passive LDPC memory constructions. This result adds the explicit fixed-circuit translation and total noisy-width/prescribed-T consequence in this model. It is not a claim of the first coherent LDPC memory or of new code/decoder performance.

The [related-work account](../references/RELATED_WORK.md) separates these ingredients. The [finite verifier](../verification/verify_tanner_memory.py) checks all two-input/two-output lookup tables on clean inputs, a nonlinear correction channel's full Choi state with a reference, exact scalar invariant and batching margins, length ratios, resource inequalities and a sufficient-size scalar certificate. The scalar fixture does not instantiate an admissible Tanner code or prove a practical threshold. The universal theorem relies on the argument above and the inherited code/decoder results; it is not certified by a small numerical simulation.
