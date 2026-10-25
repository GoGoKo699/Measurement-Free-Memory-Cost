# Quantum storage with Clifford control and a noisy syndrome archive

## Result and scope

An exponentially growing archive of syndrome information gives a quantum memory with zero non-Clifford gates during storage. For sufficiently small independent depolarizing noise, the construction stores any $`k`$ logical qubits, including their entanglement with a reference, for $`D`$ noise rounds with recovered half-diamond error at most $`\epsilon`$. Its total noisy width is

```math
N=\exp\!\bigl(O(\gamma D)\bigr)\,
 O\!\left((k+\log(J/\epsilon))
 \log\frac{2J(k+\log(J/\epsilon))}{\epsilon}\right),
\qquad J=O(1+\gamma D).
```

All data, live records, retired record registers, and temporary fanout targets are included in $`N`$. The final decoder is the only place where majority voting or a nonlinear classical syndrome decoder is evaluated. During storage, the controller follows a fixed schedule of Clifford unitaries and resets to $`|0\rangle`$.

The construction uses the **global ideal control model**: an arbitrary Clifford unitary can be completed between two consecutive noise rounds, with no noise inside that unitary. This is the unrestricted model of the converse. It does not imply the same lifetime for a bounded-depth, geometrically local, or noisy-gate implementation. Initial encoding and final recovery are noiseless endpoints outside the storage interval; their computational complexity is not bounded here. No measurement result, external classical register, or noiseless ancillary memory is available during storage.

The code-existence argument and the classical broadcast process are established ingredients. The contribution of this chapter is their integration into a precisely scheduled, entirely Clifford, measurement-free quantum memory, with finite error and workspace accounting. Together with the logarithmic Clifford-memory converse, this proves the order of the zero-non-Clifford lifetime in the stated model. It does not establish the optimal width constant or a matching minimum non-Clifford gate rate at polynomial width.

## 1. Noise, endpoints, and parameters

A physical noise round applies

```math
\mathcal N_\gamma^{\otimes N},\qquad
\mathcal N_\gamma(\rho)
=(1-\gamma)\rho+\frac{\gamma}{2}I\mathrm{Tr}\rho,
```

to every allocated physical qubit. Noise is independent between qubits and between rounds. Equivalently, the one-qubit Pauli probabilities are

```math
\Pr(I)=1-\frac{3\gamma}{4},\qquad
\Pr(X)=\Pr(Y)=\Pr(Z)=\frac{\gamma}{4}.
```

On a computational-basis record the bit-flip probability is $`\gamma/2`$. A phase flip on a record is immaterial once that record is disconnected from the data.

Let $`\mathcal V`$ encode $`k`$ logical qubits into the memory, let $`\Phi_D`$ be the storage channel, and let $`\mathcal R`$ be the final recovery. We prove

```math
\frac12\|\mathcal R\Phi_D\mathcal V-\mathrm{id}_{2^k}\|_\diamond
\le\epsilon. \tag{1}
```

In particular, with the squared-overlap convention,

```math
F_e(\Lambda)=\langle\Omega_{2^k}|(\mathrm{id}\otimes\Lambda)
(|\Omega_{2^k}\rangle\langle\Omega_{2^k}|)|\Omega_{2^k}\rangle,
```

the recovered entanglement fidelity is at least $`1-\epsilon`$.

Assume initially that

```math
k,D\in\mathbb N_{\ge1},\qquad
0<\epsilon<1,\qquad 0<\gamma\le\frac1{16}.
```

Batch the physical rounds into intervals of length at most

```math
\tau=\left\lfloor\frac1{16\gamma}\right\rfloor,\qquad
J=\left\lceil\frac D\tau\right\rceil. \tag{2}
```

Put $`t_j=\min\{j\tau,D\}`$, $`t_0=0`$, and $`\ell_j=t_j-t_{j-1}`$. The effective depolarization in block $`j`$ is

```math
\Gamma_j=1-(1-\gamma)^{\ell_j}\le\Gamma,
\qquad
\Gamma=1-(1-\gamma)^\tau\le\tau\gamma\le\frac1{16}.
\tag{3}
```

Write $`Q=1-\Gamma`$, so $`Q\ge15/16`$, and define

```math
\alpha=\frac1{16},\qquad
R_0=1-h_2(1/8)-\frac18\log_2 3>0,
```

```math
p_* = \frac{3\Gamma}{4}<\alpha,
\qquad
c_\Gamma=D_{\mathrm{bin}}(\alpha\|p_*)>0, \tag{4}
```

where $`h_2`$ is the binary entropy in bits and

```math
D_{\mathrm{bin}}(a\|p)
=a\ln\frac ap+(1-a)\ln\frac{1-a}{1-p}
```

uses natural logarithms. Set

```math
n_0=
\left\lceil\max\left\{
\frac{k+1}{R_0},\,
\frac{\ln(2J/\epsilon)}{c_\Gamma}
\right\}\right\rceil,
\qquad
s=n_0-k,\qquad
w=\lfloor n_0/16\rfloor. \tag{5}
```

Here $`n_0`$ is the data-block width, $`s`$ the number of stabilizer generators, and $`w`$ the guaranteed number of correctable errors per block. The choice implies $`s>0`$.

When $`J\ge2`$, put $`M=s(J-1)`$ and choose $`B`$ as the smallest positive odd integer satisfying

```math
B\ge\frac8{Q^2}\ln\frac{2M}{\epsilon}. \tag{6}
```

When $`J=1`$, set $`B=0`$ and omit the syndrome archive.

### Theorem 1 — finite-parameter Clifford quantum memory

For the above parameters, an encoding, a fixed measurement-free Clifford/reset storage schedule, and a final recovery satisfying (1) exist with

```math
N=
\begin{cases}
n_0,&J=1,\\[2mm]
n_0+Bs\bigl(2\cdot3^{J-2}-1\bigr),&J\ge2.
\end{cases} \tag{7}
```

Exactly $`D`$ physical noise rounds occur. The storage schedule uses no T gates, no other non-Clifford gates, and no magic-state ancillas. Each control interval requires only resets to $`|0\rangle`$ followed by one ideal global Clifford unitary. Final recovery is allowed to discard the archive and its work registers.

The rest of this chapter proves the theorem.

## 2. A finite stabilizer-code existence lemma

The symplectic stabilizer construction and its Gilbert–Varshamov counting argument are standard; see CRSS, Theorem 2 [1], and the more developed stabilizer formalism in [2]. We give the finite bound needed here, including the integer slack, so the theorem does not rely on unspecified asymptotic code parameters.

**Lemma 2.** If $`n\ge1`$, $`k\ge1`$, and $`R_0n\ge k+1`$, there is a stabilizer code on $`n`$ qubits encoding exactly $`k`$ qubits, with distance at least $`2\lfloor n/16\rfloor+1`$.

**Proof.** Choose uniformly an isotropic subspace

```math
S\subset\mathbb F_2^{2n},\qquad \dim S=n-k,
```

for the usual nondegenerate symplectic form. Such a subspace defines a stabilizer code encoding $`k`$ qubits. Since the symplectic group is transitive on nonzero vectors and preserves the uniform distribution of $`S`$, every nonzero Pauli label $`v`$ has the same probability of belonging to $`S^\perp`$. Counting all nonzero vectors in that orthogonal space gives exactly

```math
\Pr\{v\in S^\perp\}
=\frac{2^{n+k}-1}{2^{2n}-1}
\le\frac43\,2^{k-n}. \tag{8}
```

Let $`r=1/8`$, $`y=r/[3(1-r)]`$, and $`w=\lfloor n/16\rfloor`$. Since $`0<y<1`$ and $`2w\le rn`$,

```math
\sum_{j=0}^{2w}3^j\binom nj
\le y^{-rn}(1+3y)^n
=2^{n[h_2(r)+r\log_2 3]}
=2^{n(1-R_0)}. \tag{9}
```

Consequently the probability that $`S^\perp`$ contains any nonzero Pauli label of weight at most $`2w`$ is at most

```math
\frac43\,2^{k-n}\sum_{j=1}^{2w}3^j\binom nj
\le\frac43\,2^{k-R_0n}
\le\frac23<1. \tag{10}
```

Some $`S`$ therefore has no such label. In particular, every nontrivial logical Pauli, whose label lies in $`S^\perp\setminus S`$, has weight at least $`2w+1`$. The resulting code corrects every Pauli error of weight at most $`w`$. This proves the lemma. The argument actually gives a pure code, which is stronger than required. $`\square`$

Apply Lemma 2 with $`n=n_0`$. Let $`E_j`$ be the net Pauli error on the data during noise block $`j`$. Its weight has distribution $`\mathrm{Bin}(n_0,3\Gamma_j/4)`$, and the increments are independent between blocks. A binomial Chernoff bound and (3)–(5) give

```math
\Pr\{\mathrm{wt}(E_j)>w\}
\le \exp(-n_0 c_\Gamma)
\le\frac{\epsilon}{2J}. \tag{11}
```

Indeed, an integer exceeding $`\lfloor\alpha n_0\rfloor`$ exceeds $`\alpha n_0`$, and the tail is increasing in its Bernoulli parameter. The probability that any data increment is uncorrectable by this weight guarantee is at most $`\epsilon/2`$.

## 3. Exact coherent syndrome extraction

Choose independent stabilizer generators $`g_1,\ldots,g_s`$. For a Pauli $`A`$, define $`\mathrm{syn}(A)\in\mathbb F_2^s`$ through its commutation signs with these generators. The code has a Clifford encoder $`U`$ with

```math
V|\psi\rangle=U\bigl(|\psi\rangle\otimes|0\rangle^{\otimes s}\bigr),
\qquad U Z_{k+a}U^\dagger=g_a.
```

One global Clifford operation extracts $`B`$ replicas of each syndrome bit using only the data and $`Bs`$ fresh record targets:

1. Apply $`U^\dagger`$ to the data.
2. Apply CNOTs from the $`s`$ decoded syndrome coordinates to their respective $`B`$ fresh targets.
3. Apply $`U`$ to the data.

For every Pauli $`A`$ and every logical state, including a state entangled with an untouched reference, the transformation is

```math
A V|\psi\rangle\otimes|0\rangle^{\otimes Bs}
\longmapsto
A V|\psi\rangle\otimes|\mathrm{syn}(A)\rangle^{\otimes B}.
\tag{12}
```

To see this, $`U^\dagger A U`$ is a Pauli. Its action on the final $`s`$ input coordinates produces a definite computational-basis syndrome; these coordinates are unentangled with the logical input in each Pauli branch. Copying them therefore leaves the data unchanged after re-encoding. There is no temporary syndrome register besides the $`Bs`$ targets themselves.

Throughout storage, record qubits are never used as controls or targets of a gate acting on data after their extraction. A live record is used only as a CNOT control on fresh record children. Conditional on the complete data Pauli history, all such records are computational-basis bits subjected to independent classical bit flips. This is an exact description of the conditional quantum channel. It does not assume the unknown logical state is classical and does not copy quantum logical information.

## 4. Schedule and archive dynamics

Let

```math
A_j=E_jE_{j-1}\cdots E_1,
\qquad \sigma_j=\mathrm{syn}(A_j),
\qquad \sigma_0=0.
```

At the beginning of noise block $`j`$, the current data syndrome is $`\sigma_{j-1}`$. The fixed schedule is as follows.

| Stage | Operation |
|---|---|
| Start of block 1 | Keep the encoded data; no archive is needed for the known zero syndrome. |
| Start of block $`j\ge2`$ | Reset all designated new targets to $`|0\rangle`$. Extract $`B`$ copies of each bit of $`\sigma_{j-1}`$. For every old live archive leaf, copy its bit into three fresh children by CNOTs. |
| Remaining control in the block | Retire the old parent slots from the live archive; they are ignored and may be reset for reuse at a later block. All other physical control is the identity. |
| Noise in block $`j`$ | Apply $`\ell_j`$ ordinary all-qubit noise rounds. All allocated slots, including inactive and retired ones, are noisy. |
| After physical round $`D`$ | Use the unrestricted final recovery to read archive leaves and measure the data syndrome $`\sigma_J`$. |

The syndrome extraction and every fanout in a given block compose to one global Clifford unitary. The retire operation is only a change in the predetermined designation of registers; no physical removal, noiseless storage, or measurement is needed. Slots from earlier generations can remain noisy until they are reset as fresh targets.

The archived syndromes are exactly $`\sigma_1,\ldots,\sigma_{J-1}`$. The final syndrome $`\sigma_J`$ is measured only at the allowed recovery endpoint. Thus every noise block has a syndrome before and after it, without a hidden extra extraction or uncounted noise round.

### Lemma 3 — uniform reconstruction of one archived bit

For every age of an archive, recursive ternary majority on its leaves reconstructs each noisy root with error at most $`1/4`$. Including the root's initial noise, each replica reconstructs its true extracted syndrome bit with error at most

```math
\frac12-\frac Q4. \tag{13}
```

**Proof.** Put $`e_*=\Gamma/2`$. Initially regard the root value immediately after the noise in its extraction block as the bit to be reconstructed. Let $`\eta_h`$ upper-bound the error after $`h`$ further generations. With no descendants, the final leaf is the root itself and $`\eta_0=0`$. Each child first receives its parent's bit, then suffers a bit flip of probability at most $`e_*`$. Its recursively reconstructed value differs from the parent with probability at most

```math
e_*+(1-2e_*)\eta_h=e_*+Q\eta_h.
```

The three child subtrees are independent conditional on the parent. If the final shorter block has smaller noise, replacing that noise by $`e_*`$ only increases this bound while $`\eta_h\le1/2`$. Therefore

```math
\eta_{h+1}\le f_3(e_*+Q\eta_h),
\qquad f_3(u)=3u^2-2u^3. \tag{14}
```

If $`\eta_h\le1/4`$, then

```math
e_*+Q\eta_h\le\frac{1+\Gamma}{4}\le\frac{17}{64},
\qquad
f_3(17/64)=\frac{22831}{131072}<\frac14.
```

Induction proves the first claim. The initially extracted root also undergoes noise before its first descendants are made. Restoring this initial noise gives error at most

```math
e_*+Q/4=\frac12-Q/4,
```

which proves (13). $`\square`$

The initial root noise must be included: $`1/4`$ is a bound for reconstructing an already noisy root, not for reconstructing the syndrome before that root's first noisy block.

For each syndrome bit, the $`B`$ initial replicas and all their descendant trees use disjoint physical registers at any common time. Reused registers are reset, and their future noise is independent of past noise. Conditional on the data history, the $`B`$ reconstruction errors are consequently independent. Hoeffding's inequality and (13) imply

```math
\Pr\{\text{majority of the }B\text{ replicas is wrong}\}
\le\exp(-BQ^2/8). \tag{15}
```

There are $`M=s(J-1)`$ archived syndrome bits. Equation (6) and a union bound show that all are correctly reconstructed with probability at least $`1-\epsilon/2`$, uniformly in the data history.

The broadcast process is a classical, established model [3]. The elementary recursion (14) is sufficient here; neither an optimal tree-reconstruction threshold nor optimality of recursive majority is used or claimed.

## 5. Final recovery and the channel guarantee

The final recovery measures the surviving archive leaves in the computational basis, performs the tree and replica majorities, and measures the final data syndrome $`\sigma_J`$ noiselessly. It then forms

```math
\widehat\delta_j=\widehat\sigma_j\oplus\widehat\sigma_{j-1},
\qquad j=1,\ldots,J,
```

with the exact endpoint values $`\widehat\sigma_0=0`$ and $`\widehat\sigma_J=\sigma_J`$. For each $`\widehat\delta_j`$, choose a deterministic minimum-weight Pauli representative $`C_j`$ with that syndrome. Syndrome representatives always exist because the stabilizer generators are independent. Apply the inverse of the estimated cumulative Pauli, then Clifford-decode the data and retain the $`k`$ logical output qubits.

Consider the good event consisting of:

1. $`\mathrm{wt}(E_j)\le w`$ for every data increment;
2. every archived syndrome bit is reconstructed correctly.

On this event, syndrome linearity gives

```math
\widehat\delta_j=\mathrm{syn}(E_j).
```

The minimum-weight representative satisfies $`\mathrm{wt}(C_j)\le w`$. Thus $`C_j^{-1}E_j`$ has zero syndrome and weight at most $`2w`$, strictly smaller than the code distance. It is a stabilizer up to a scalar phase. Since Pauli products commute up to scalar phase, the estimated cumulative Pauli differs from $`A_J`$ by a stabilizer and a global phase. Its inverse exactly restores the encoded logical state, including any entanglement with a reference.

Equations (11) and (15) give a total bad-event probability $`\beta\le\epsilon`$. This probability does not depend on the logical input. More explicitly, all physical depolarization can be sampled as classical Pauli labels independently of the input. Syndrome extraction, record outcomes, and the good-event test are then deterministic functions of those labels. Every good branch induces precisely the identity channel on the logical input. Therefore the complete recovered channel can be written

```math
\mathcal R\Phi_D\mathcal V
=(1-\beta)\mathrm{id}+\beta\Lambda_{\rm bad},
\qquad 0\le\beta\le\epsilon, \tag{16}
```

for a CPTP channel $`\Lambda_{\rm bad}`$ when $`\beta>0`$. The diamond distance between two channels is at most two, so (16) proves (1). Its maximally entangled overlap is at least $`1-\beta`$, giving the asserted entanglement-fidelity guarantee.

Deferred Pauli correction is valid here because the intervening data operation is nondemolition extraction of its current syndrome. Applying physical feedback after each increment is unnecessary. All nonlinear processing of syndrome information has been postponed until the recovery endpoint.

## 6. Width, work, and lifetime

Let $`L_j`$ denote the number of live archive leaves immediately after block $`j`$. For $`J\ge2`$,

```math
L_1=0,\qquad
L_j=3L_{j-1}+Bs
=\frac{Bs}{2}(3^{j-1}-1),\quad j\ge2. \tag{17}
```

At the start of block $`j`$, the $`L_{j-1}`$ old leaves coexist with all $`L_j`$ fresh targets: three children per parent, plus the $`Bs`$ new roots. This parent/child overlap is essential to the workspace count. A pool of

```math
L_{J-1}+L_J=Bs(2\cdot3^{J-2}-1) \tag{18}
```

archive slots suffices, because $`L_{j-1}+L_j`$ increases with $`j`$. At each block, reset $`L_j`$ slots outside the current live set, use them as targets, and then designate only those targets as live. Retired slots need no further access until a later reset. Equation (18), together with the $`n_0`$ data qubits, is exactly (7).

For completeness, the number of fresh record-target preparations over the entire storage interval is

```math
\sum_{j=2}^J L_j
=\frac{Bs}{4}\bigl(3^J-2J-1\bigr), \tag{19}
```

and the total number of record-to-record CNOTs is $`3\sum_{j=2}^J L_{j-1}`$. There are $`Bs(J-1)`$ data-to-record CNOTs in the extraction circuits, in addition to the Clifford encoder and its inverse for each extraction. These are finite predetermined circuits. Their elementary gate depth is not charged in the global-control model. They introduce no additional quantum work registers beyond (7).

The two endpoints may use arbitrary computational work as permitted by their channel definition. In particular, minimum-weight decoding need not be computationally efficient. Endpoint work is not hidden noiseless memory used during the $`D`$ physical noise rounds. The archive is allowed to be discarded after recovery; a requirement that every register be returned clean would define a different endpoint condition.

The constants can be chosen uniformly for $`0<\gamma\le1/16`$: since $`\Gamma\le1/16`$,

```math
c_\Gamma\ge D_{\rm bin}(1/16\|3/64)>0,
\qquad Q\ge15/16.
```

Thus

```math
n_0=O(k+\ln(J/\epsilon)),\qquad
B=O\!\left(\ln\frac{2J(k+\ln(J/\epsilon))}{\epsilon}\right),
\tag{20}
```

and

```math
N=O\!\left[
(k+\ln(J/\epsilon))
\ln\frac{2J(k+\ln(J/\epsilon))}{\epsilon}\;3^J
\right]. \tag{21}
```

Since $`\lfloor x\rfloor\ge x/2`$ for $`x\ge1`$, equation (2) implies

```math
J\le32\gamma D+1. \tag{22}
```

For fixed $`k`$ and $`0<\epsilon<1`$, the family therefore has $`\ln N=O(1+\gamma D)`$, with logarithmic corrections already absorbed for large $`D`$. Conversely, for every sufficiently large available width one can select $`J`$ of order $`\log N`$, pad unused slots if necessary, and obtain a lifetime of order $`\gamma^{-1}\log N`$. For example, selecting $`J\le(\ln N)/(2\ln3)`$ leaves room for the factors in (21) once $`N`$ is large enough. Fix a fidelity target $`4^{-k}<F_0<1`$ and take $`\epsilon\le1-F_0`$. Together with the Clifford-memory converse, the optimal lifetime at this fixed logical size and target fidelity has order

```math
D_{\max}(N)=\Theta(\gamma^{-1}\log N)
```

throughout the stated small-noise regime, with constants depending on the fixed logical size and error tolerance. This is a comparison of lifetime order in the same endpoint and global-control model. The exact best exponent relating width and lifetime is not determined.

## 7. Boundary cases and limitations

- If $`D=0`$, use $`N=k`$ and the identity channel. No code or archive is needed.
- If $`\gamma=0`$, $`N=k`$ suffices for every finite $`D`$; the batching formula is unnecessary.
- If $`J=1`$, the ordinary code and final recovery handle the single noise block, with no archived syndrome and no $`\log0`$ parameter.
- If $`k=0`$, the storage task has a one-dimensional logical system and is vacuous; zero data qubits suffice under the convention allowing an empty memory.
- The theorem concerns $`0<\epsilon<1`$. Exact error zero under nonzero depolarization is not asserted. For $`\epsilon\ge1`$, the stated approximation target is vacuous.
- The explicit constant $`1/16`$ is a convenient sufficient noise range, not a threshold claim. Better codes, larger branching factors, or better tree decoders can change constants.
- The result is an existence construction. It uses a finite stabilizer search and a possibly exponential-time final syndrome decoder. It does not supply efficient code search, a local gate schedule, or a decoder available during storage.
- This construction shows that exponential noisy width is a possible zero-T escape from the converse. Optimal non-Clifford resource rates for polynomial width remain a separate problem.

## 8. Verification and provenance

The accompanying [archive checker](../verification/verify_archive.py) is a finite validation of delicate steps, not a replacement for the proof. It checks the complete Pauli-history correction identity for the five-qubit code over three increments, exact small ternary broadcast distributions against the recursion, the initial-noise correction, integer resource scheduling, and finite coding/error inequalities for representative parameters. Its reproducible output is [saved archive output](../verification/expected/archive.json).

The stabilizer/Gilbert–Varshamov ingredients and tree-broadcast model are inherited. Equations (8)–(10) give a self-contained finite specialization of the standard random-subspace argument. Equations (12)–(19) establish the particular memory schedule and its quantum-channel guarantee. No claim of novelty is made for standard code existence, majority concentration, deferred Pauli frames, or classical broadcasting individually.

## References

1. A. R. Calderbank, E. M. Rains, P. W. Shor, and N. J. A. Sloane, “Quantum Error Correction and Orthogonal Geometry,” *Physical Review Letters* **78**, 405–408 (1997), Theorems 1–2. [DOI](https://doi.org/10.1103/PhysRevLett.78.405); [corrected arXiv v3 and full text](https://arxiv.org/html/quant-ph/9605005v3). The corrected theorem gives the asymptotic stabilizer rate $`1-h_2(2\delta)-2\delta\log_2 3`$.
2. A. R. Calderbank, E. M. Rains, P. W. Shor, and N. J. A. Sloane, “Quantum Error Correction via Codes over GF(4),” *IEEE Transactions on Information Theory* **44**, 1369–1387 (1998). [DOI](https://doi.org/10.1109/18.681315); [arXiv v5 full text](https://arxiv.org/html/quant-ph/9608006v5). See Sections 2–3 for the additive/stabilizer construction and Clifford encoding. Earlier code-existence work is A. R. Calderbank and P. W. Shor, “Good Quantum Error-Correcting Codes Exist,” *Physical Review A* **54**, 1098–1105 (1996), [DOI](https://doi.org/10.1103/PhysRevA.54.1098), [arXiv](https://arxiv.org/abs/quant-ph/9512032).
3. W. Evans, C. Kenyon, Y. Peres, and L. J. Schulman, “Broadcasting on Trees and the Ising Model,” *Annals of Applied Probability* **10**, 410–433 (2000). [DOI](https://doi.org/10.1214/aoap/1019487349); [author-hosted paper](https://www.cs.ubc.ca/~will/papers/noisytrees.pdf); [author publication list](https://www.cs.ubc.ca/~will/papers/index.html). The present argument uses only the elementary independent-edge broadcast model and proves its sufficient recursive-majority bound directly.

Primary sources and corrected versions checked on 25 October 2026. The finite code lemma is proved above rather than relying on the incomplete Gilbert–Varshamov proof in the original version of Gottesman's thesis; see the author's [Section 7.1 erratum](https://personal.perimeterinstitute.ca/personal/dgottesman/thesis-errata.html).
