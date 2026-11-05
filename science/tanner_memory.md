# Coherent Tanner memory with linear space and processing per epoch

## 1. Result and scope

Existing good quantum Tanner codes and their constant-iteration parallel decoder give a stronger sufficient resource bound in the [same memory model](model.md). This chapter proves the fixed-circuit compilation, accounts for its total width and prescribed T count, and combines it with the published residual-error invariant. A fixed number of inner five-qubit concatenation levels extends the physical noise range. These are resource consequences of inherited codes, decoders and measurement-free correction mechanisms.

Fix an admissible Tanner family as specified below. For $`0\lt\gamma\le1/20`$, $`k\ge1`$, $`D\ge0`$ and $`0\lt\epsilon\lt1/2`$, a memory exists with

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

Here $`n`$ includes data, syndrome, controller and all temporary work qubits. The constants depend on the fixed Tanner family and are independent of $`k,D,\epsilon,\gamma`$ in the stated range. They can be very large. Sections 2–7 retain the original direct Tanner construction, which requires a possibly tiny family-dependent $`\gamma_*\gt0`$. Section 8 proves the extended range by explicit composition with the [five-qubit construction](coherent_recovery.md); it does not assume that the Tanner decoder itself has that physical noise threshold.

At fixed $`\gamma\in(0,1/20]`$, write $`L=\ln(D/\epsilon)`$. The [reliability converse](reliability_processing.md), together with the payload floor, proves that every family with $`D\to\infty`$, $`0\lt\epsilon\le\epsilon_0\lt1/2`$ and $`n\le C(k+L)`$ needs $`m_T=\Omega_{\gamma,C,\epsilon_0}(D(k+L))`$. Equation (T1) attains this order at width $`O(k+L)`$, including fixed $`k`$. For $`k=\Omega(L)`$ it reduces to the payload-dominated order $`k\gamma D`$. This sharpens the resource comparison without altering the inherited decoder or this construction. It does not determine a numerical coefficient, uniform weak-noise law, practical constants, locality or a threshold with faulty processing gates.

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

For the direct Tanner construction, choose

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

proving (T1) in the original range $`0\lt\gamma\le\gamma_*`$. For $`D=0`$, take $`J=0`$ and $`J_+=1`$; no interior recovery is needed.

## 8. Constant inner concatenation extends the physical noise range

The [five-qubit recovery](coherent_recovery.md) gives an exact coherent circuit which decodes one concatenated block, corrects it and resets its syndrome wires. For independent physical Pauli errors with nonidentity probability at most $`p_0`$, its fresh residual logical Pauli has probability at most

```math
p_L\le\frac{(10p_0)^{2^L}}{10}.
```

This circuit uses $`W=5^L`$ physical wires and $`42(W-1)/4`$ prescribed T/T-dagger slots, with no additional work qubits. The syndrome wires are reset and reused by the Clifford reencoder. These are the existing five-qubit circuit and concatenation bounds, rather than new decoding statements.

Choose once, using the outer decoder's constant $`u`$ from (T7),

```math
L_* =\min\left\{\ell\in\mathbb Z_{\ge0}:
\eta_\ell:=\frac1{10}\left(\frac38\right)^{2^\ell}
\le\frac u2\right\},
\qquad W=5^{L_*}.
\qquad\text{(T11)}
```

The minimum exists and depends only on the fixed outer family. In particular, the inner concatenation depth does not grow with payload, duration or target error. For $`0\lt\gamma\le1/20`$, replace (T8) by

```math
h=\left\lfloor\frac1{20\gamma}\right\rfloor,
\qquad J=\left\lceil\frac D h\right\rceil,
\qquad J\le40\gamma D+1.
\qquad\text{(T12)}
```

Encode each of the $`N`$ outer Tanner positions into a separate $`W`$-qubit inner block. Following every full batch, and the shorter final batch, execute the following fixed schedule within one ideal processing interval:

1. Decode and correct all inner blocks through their fixed $`L_*`$ levels; reset the released inner syndrome wires.
2. Reset the outer scratch pool, extract the current outer syndrome, execute the fixed Tanner circuit from Sections 3–4, and reset its syndrome records.
3. Reencode every outer data wire into its inner block by Clifford operations using the released zero wires.

During the intervening noisy layers, all $`WN`$ inner data wires and the entire outer scratch pool remain noisy and counted. No record survives in protected storage. The temporarily decoded outer data are exposed only during ideal processing, as permitted by the existing model. Each reset has a prescribed, input-independent output and timing.

### Fresh inner errors and pre-existing outer residuals

A batch of at most $`h`$ physical noise layers has nonidentity Pauli probability

```math
p_0=\frac34(1-q^h)\le\frac{3h\gamma}{4}\le\frac3{80}
```

per physical position. For a fixed inner block, let $`E`$ be its complete fresh physical Pauli pattern. The inner correction circuit sends an arbitrary encoded logical operator to a definite logical Pauli conjugation $`\omega\mapsto P_{\ell(E)}\omega P_{\ell(E)}^\dagger`$. This holds on all logical inputs, including those entangled with other blocks or a reference; syndrome resets implement the associated correction channel exactly.

A pre-existing outer residual Pauli acts as a logical Pauli $`P`$ of the inner block. Its encoded representative commutes with the inner stabilizers, so it changes no extracted inner syndrome. Commuting it through $`E`$ changes only an irrelevant Pauli phase. Consequently the corrected logical operator is

```math
P_{\ell(E)}P\,\omega\,P^\dagger P_{\ell(E)}^\dagger.
\qquad\text{(T13)}
```

Thus the old outer residual is followed by a fresh error whose label depends only on the new physical noise. It does not change the distribution of that fresh label or the validity of the inner error estimate. The same conclusion follows level by level in the concatenation tree.

Different inner blocks use disjoint fresh physical-error samples. Conditional on any prior noise history, their fresh logical Pauli labels are therefore independent, with nonidentity probability at most $`\eta_{L_*}\le u/2`$. The resulting logical Pauli channel need not be depolarizing. Its X and Z components at the same outer position can be correlated. Each component separately is an independent collection across positions with marginal probabilities at most $`u/2`$. Its weight is stochastically bounded by $`\mathrm{Bin}(N,u/2)`$, even if the marginal probabilities differ between positions or epochs.

Apply the two marginal Chernoff bounds and a union bound, using the same $`c=\mathrm{KL}(u\|u/2)`$ as before. With probability at least $`1-2Je^{-cN}`$, every fresh outer error has each CSS component of weight at most $`uN`$. Section 5's residual invariant then applies without change. In particular, inner decoding failures are corrected by the outer code when they lie in this region; no unnecessary union bound requiring every inner block to succeed is imposed.

The exact inner correction channel, the outer channel (T6), and the Clifford reencoder compose as CPTP maps on the entire system and reference. Every trajectory in the good event has the same correctable final outer residual. The free endpoint first decodes the inner encoding and then applies the fixed outer endpoint decoder. It is the identity logical channel on every good trajectory. The input-independent trajectory mixture proves the same half-diamond bound as Section 7.

### Width, prescribed processing and reliability

Use (T10) with the epoch count from (T12). The outer code length still obeys $`N=O(k+\ln(J_+/\epsilon))`$. A conservative simultaneous allocation reserves separate inner blocks and the complete outer circuit pool; it gives

```math
n\le(W+C_W)N,
\qquad
m_T\le\left[\frac{42(W-1)}4+C_T\right]NJ.
\qquad\text{(T14)}
```

The width bound deliberately permits overlap to be counted twice rather than assuming unproved scratch reuse. All circuits and every T/T-dagger slot are fixed in advance, including epochs whose syndromes are zero. The inner encoder, all syndrome extraction, routing and reencoding are Clifford operations already permitted by the model. The count includes the final recovery, although some final operations could instead be absorbed into the free endpoint.

Equations (T10)–(T14) prove (T1) uniformly for $`0\lt\gamma\le1/20`$. At $`D=0`$, set $`J=0`$; the endpoint encoding and decoding require no interior recovery. The extension can increase the already large constants. It supplies neither a practical numerical threshold for good Tanner codes nor tolerance of noise inside the processing circuit. In the smaller range $`\gamma\le\gamma_*`$, the earlier direct construction remains available and can have lower overhead.

## 9. Attribution, interpretation and finite checks

Good Tanner codes, the parallel decoder, its contraction and its repeated-round use are inherited. Coherent syndrome extraction followed by measurement-free nonlinear correction is also established prior work, including [Hong, Guo and Lucas](https://doi.org/10.1038/s41467-024-55570-7), who study measurement-free error correction and passive LDPC memory constructions. This result adds the explicit fixed-circuit translation and total noisy-width/prescribed-T consequence in this model. It is not a claim of the first coherent LDPC memory or of new code/decoder performance.

The fixed inner concatenation in Section 8 combines that theorem with the existing five-qubit circuit. Its addition is the proved physical-noise range and explicit combined resource accounting, not a new concatenated code or threshold mechanism.

The [related-work account](../references/RELATED_WORK.md) separates these ingredients. The [finite verifier](../verification/verify_tanner_memory.py) checks all two-input/two-output lookup tables on clean inputs, a nonlinear correction channel's full Choi state with a reference, exact scalar invariant and batching margins, length ratios, resource inequalities and sufficient-size scalar certificates. It also exhaustively checks that every five-qubit fresh Pauli pattern preserves the inner syndrome and acquires precisely the expected residual when any previous logical Pauli is present, checks biased Pauli-error distributions with exact rational arithmetic, and checks the combined schedule and constant-depth choice. The scalar fixtures do not instantiate an admissible Tanner code or prove a practical threshold. The universal theorem relies on the argument above and the inherited code/decoder results; it is not certified by a small numerical simulation.
