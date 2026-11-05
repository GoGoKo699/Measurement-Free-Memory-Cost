# Literature, contribution, and operational comparison

**Project:** Measurement-Free-Memory-Cost  
**Supporting-literature audit:** 25 October 2026  
**Manuscript citation check:** 3 November 2026; repository alignment 4 November 2026  
**Reopened construction audit:** 5 November 2026, covering good Tanner-code decoding, code-length availability and passive measurement-free LDPC memories.
**Scope:** the construction audit supplements the earlier manuscript-citation check. It is a targeted primary-source comparison, not an exhaustive novelty search. Bibliography keys refer to `references.bib`.

## 1. Contribution in relation to prior work

The main result is a quantitative requirement for **sustained coherent non-Clifford processing** in a noisy quantum memory. A finite amount of physical width supplies a logarithmic time allowance; maintaining a recoverable quantum payload beyond that allowance requires non-Clifford operations at an extensive rate. The statement remains meaningful when input-independent resets, even to magic states, are free.

For $`n`$ physical qubits, $`D`$ independent noise rounds, $`q=1-\gamma\in(0,1)`$, and a recovered $`k`$-qubit channel of entanglement fidelity at least $`1-\epsilon`$, the gate-budget theorem is

```math
m_T\ge a\left[D(-\ln q)-\ln\frac{2n}{a}\right]_+,
\qquad
a = \frac{2k(1-\epsilon-4^{-k})}{1-4^{-k}}.
```

The nontrivial regime is $`\epsilon\lt 1-4^{-k}`$. At fixed acceptable infidelity and beyond the width allowance, this gives $`m_T=\Omega(k\gamma D)`$. The companion zero-$`T`$ statement gives replacement-channel convergence on the scale $`\gamma^{-1}\log(n/\epsilon)`$, including a diamond-norm formulation. These statements are proved in the project's Pauli-envelope and control-cost chapters.

The new step is contraction of the simultaneous binary Pauli-envelope rank and its conversion into reference-stable forgetting, fidelity-dependent gate counts, and temporal interval laws. Nelson, Rajakumar, and Gullans already proved a quantitative forgetting theorem for noisy Clifford processing with resets. The present result sharpens its width dependence and relates the rank bottleneck to continuing internal work. Pauli symplectic coordinates, operator-algebra decomposition, Schmidt-number bounds, concatenated codes, and seven-$`T`$ Toffoli circuits are standard ingredients.

## 2. Direct prior: Nelson, Rajakumar, and Gullans

The essential predecessor is **Nelson, Rajakumar, and Gullans, “Non-Clifford Gates Are Required for Long-Term Memory,” CCC 2026**, DOI [10.4230/LIPIcs.CCC.2026.38](https://doi.org/10.4230/LIPIcs.CCC.2026.38), first posted as [arXiv:2510.08451](https://arxiv.org/abs/2510.08451) in October 2025. The comparison here uses the final proceedings paper (`nelson2026nonclifford`).

Its Theorem 1 gives loss of pairwise trace distinguishability after $`O(\gamma^{-1}\log n\log(n/\delta))`$ layers at tolerance $`\delta`$. Lemma 6 controls simultaneous survival of nontrivial Pauli labels. Section 2.3, Eq. (4), already permits replacement by **arbitrary single-qubit states**; Eqs. (8)–(11) specify the reset adjoint. Fresh magic states are therefore already allowed. Backward Pauli evolution and stochastic complete-depolarization trajectories also belong to this predecessor. The project's expected-rank contraction sharpens the width dependence to a single logarithm, gives a reference-stable channel bound, and yields a fidelity-dependent gate budget at every temporal cut. Allowing joint replacement states is a model extension, not the headline novelty.

The public arXiv record still listed only v1 at the 3 November manuscript check. A related four-author QCTiP 2026 abstract adds Dominik Hangleiter and uses the title “Limitations on Measurement-free Fault-tolerant Protocols using Clifford Circuits”; it states broader qualitative results but no quantitative replacement for the published bound. See the [conference abstract, p. 38](https://qctipconf.github.io/qctip_2026/QCTiP2026_Book_of_Abstracts_002.pdf). This audit cannot exclude an unpublished refinement.

## 3. Model boundaries that matter for comparison

The project's depolarizing layer is

```math
\mathcal D_\gamma=(1-\gamma)\mathrm{id}
   +\gamma\,R_{I/2},
```

independently across persistent qubits and rounds. Its processing intervals can contain arbitrary ideal Clifford gates, replacement resets, and counted $`T/T^\dagger`$ gates. Consequently, $`D`$ is a number of noisy intervals, **not automatically the depth of an elementary noisy gate circuit**.

These conventions determine which comparisons with other memory models are valid.

| Issue | Meaning in this project | Consequence |
|---|---|---|
| Reset | A replacement state independent of the unknown input and the replaced state | An external source retaining input-dependent memory, or temporally correlated program registers that carry the payload, is outside the model. |
| Classical records | Every record kept between noise rounds belongs to the noisy width $`n`$ | A protected adaptive classical controller would change the allowed operations. |
| Magic | Fresh input-independent magic states may be free | The lower bound measures coherent non-Clifford gate work; it is not a lower bound on the number of supplied magic states. |
| Endpoints | Initial encoding and final decoding may be arbitrary | The cost is for sustaining the memory between endpoints, not for preparing its initial codeword. |
| Recovery work | All coherent nonlinear feedback inside the memory is charged | A RAM decoder's classical running time cannot be silently equated with the $`T`$-count of a fixed coherent implementation. |
| Noise convention | $`\gamma`$ is the probability of complete depolarization in the channel decomposition | It differs by a factor from conventions that use total nonidentity Pauli-error probability. |

A necessary gate budget is not a sufficiency statement at fixed width. For $`n`$ qubits, each complete noise layer has probability $`\gamma^n`$ of replacing every register simultaneously, after which input-independent processing cannot restore the unknown state. The project therefore also records the elementary ceiling $`F_{\mathrm e}\le4^{-k}+(1-4^{-k})(1-\gamma^n)^D`$, regardless of the counted gate supply.

The rank argument uses an envelope of Pauli labels. An arbitrary reset can leave a nonlinear set of surviving labels, so a proof that presumes the actual surviving set is always a vector space is insufficient. The project instead propagates a linear span that safely contains it. Likewise, extending the theorem to arbitrary channels that map individual Paulis to scalar multiples of Paulis requires additional structure; the verified proof relies on the linear action on binary Pauli labels.

## 4. Existing positive results for measurement-free correction

The existence of scalable quantum computation without intermediate measurement is longstanding. Aharonov and Ben-Or's original threshold theorem explicitly describes measurement-free computation with polylogarithmic overhead ([arXiv:quant-ph/9611025](https://arxiv.org/abs/quant-ph/9611025), `aharonov1997faulttolerant`). It does not give the resource law sought here.

Several recent papers supply concrete architectures:

| Primary reference | Established contribution relevant here | What it does not establish for this comparison |
|---|---|---|
| Paz-Silva, Brennen, and Twamley, PRL **105**, 100501 (2010), [DOI](https://doi.org/10.1103/PhysRevLett.105.100501), `pazsilva2010fault` | Fault-tolerant measurement-free Bacon–Shor circuits and a threshold analysis allowing slow/noisy measurements | A lower bound, or an asymptotically optimal coherent $`T`$ budget for an arbitrary memory |
| Heußen, Locher, and Müller, PRX Quantum **5**, 010333 (2024), [DOI](https://doi.org/10.1103/PRXQuantum.5.010333), `heussen2024measurementfree` | Fault-tolerant coherent feedback for low-distance CSS codes, including auxiliary preparation and hardware-oriented benchmarks | A uniform $`O(k\gamma D)`$ memory construction |
| Butt et al., Science Advances **11**, eadv2590 (2025), [full primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155506/), `butt2025measurementfree` | A scalable universal measurement-free architecture with code switching and concatenation | Removal of the reliability overhead in this project's work bound |
| Veroni, Paler, and Giudice, PRX Quantum **6**, 040337 (2025), [DOI](https://doi.org/10.1103/lkk1-v6wp), `veroni2025universal` | Universal Bacon–Shor processing and a disposable Toffoli feedback gadget that enables measurement-free concatenation | A theorem equating its physical resource cost to the present lower bound |

In the last paper, Sec. IV.1 explicitly gives the usual recursively increasing correctable-error order under concatenation; its circuit-level threshold study is a separate physical result. These constructions support a polylogarithmic reliability overhead at sufficiently low fixed noise. Their numerical threshold values should not be substituted for the present ideal-processing model's $`\gamma`$.

### Good-code decoding and existing passive measurement-free memories

**Leverrier and Zémor** introduced the good quantum Tanner codes and sequential/parallel mismatch decoders (`leverrier2022tanner`, `leverrier2023decoding`). **Gu, Tang, Caha, Choe, He and Kubica**, CMP **405**, 85 (2024), [published paper](https://doi.org/10.1007/s00220-024-04951-6), [arXiv:2306.12470v2](https://arxiv.org/abs/2306.12470v2), `gu2024singleshot`, prove the single-shot and repeated-round residual mechanism. The inspected v2 and Springer versions both state it as Theorems 3.5 and 4.20, using Algorithms 3–4. A fixed number of parallel iterations preserves a correctable residual region using only the current syndrome. Their good-code theorem is inherited from the code/decoder sources. **Dinur, Evra, Livne, Lubotzky and Mozes**, [arXiv:2111.04808v2, Claim 6.7](https://arxiv.org/abs/2111.04808v2), `dinur2022locallytestable`, supply the Ramanujan family used to verify bounded-ratio code-length availability. These are coding and decoding contributions, not results originating in this repository.

**Hong, Guo and Lucas**, Nature Communications **16**, 316 (2025), [DOI](https://doi.org/10.1038/s41467-024-55570-7), `hong2025memory`, already give passive LDPC memories and coherent correction with resettable ancillas. The inspected [arXiv v3](https://arxiv.org/abs/2403.10599v3), Appendix E, explicitly stores stabilizer information coherently, applies feedback, and resets the ancillas. Its concrete HGP theorem (B.8) gives a lifetime exponential in the square root of data width, and the local implementation uses a constant number of ancillas per data qubit. Thus coherent syndrome feedback, fixed local processing and measurement-free LDPC memory are prior mechanisms. The present Tanner compilation is not claimed as their first realization. This project's ideal interior operations and charged T slots require a separate resource translation; the prior mechanism is credited even where its physical noise and lifetime statement differ.

### A promising but nonmatching low-overhead comparison

**Guo, Hong, Kaufman, and Lucas**, “Toward Self-Correcting Quantum Codes for Neutral Atom Arrays,” PRX Quantum **7**, 010301 (2026), [final paper](https://journals.aps.org/prxquantum/pdf/10.1103/mfmt-fwkg), `guo2026toward`, is particularly useful. Appendix D implements local greedy decisions using fixed-size reversible majority gadgets and a bounded-color schedule. **Our inference:** after compiling those fixed gadgets, a sweep uses $`O(n)`$ non-Clifford gates. This is a concrete possible route to low work overhead.

It is not a matching theorem here. The sustained-decoder results are numerical and phenomenological, and the paper leaves full noisy coherent-decoder analysis for future work. It does not prove a suitable finite-rate code family with arbitrary-duration storage in the precise model above. Even a proved $`O(kD)`$ cost at a fixed $`\gamma\gt 0`$ would not establish a bound uniform in weak noise of the form $`O(k\gamma D)`$; the dependence of constants on $`\gamma`$ matters.

## 5. Same-model constructive comparison supplied by this project

The project includes an explicit upper construction rather than relying on an informal threshold-theorem conversion. Its ingredients are the standard perfect $`[[5,1,3]]`$ code, concatenation, coherent correction, and periodic batching. The perfect code was introduced by **Laflamme, Miquel, Paz, and Zurek**, PRL **77**, 198–201 (1996), [DOI](https://doi.org/10.1103/PhysRevLett.77.198), and **Bennett, DiVincenzo, Smolin, and Wootters**, PRA **54**, 3824–3851 (1996), [DOI](https://doi.org/10.1103/PhysRevA.54.3824); keys `laflamme1996perfect` and `bennett1996mixed`.

The explicitly derived syndrome correction polynomials use three Toffoli and three CCZ gates per internal decoding node. Their 42-$`T`$ total uses exact seven-$`T`$ decompositions. **Selinger**, PRA **87**, 042302 (2013), Eq. (1), p. 1 of the [author manuscript](https://www.mathstat.dal.ca/~selinger/papers/toffoli.pdf), explicitly displays the three-wire seven-$`T`$ Toffoli; Eq. (3) credits Amy et al.'s depth improvement. Selinger's separate $`T`$-depth-one construction uses four extra ancillas and is not the construction invoked here. See `selinger2013tdepth` and `amy2013meet`.

For $`0\lt \gamma\le1/20`$, take

```math
\tau=\left\lfloor\frac{1}{20\gamma}\right\rfloor,
\qquad J=\left\lceil\frac{D}{\tau}\right\rceil.
```

The effective physical error probability over a batch is bounded by a constant below the elementary concatenation threshold. Disjoint child blocks give the recursion

```math
p_{\ell+1}\le {5\choose2}p_\ell^2=10p_\ell^2.
```

With $`L`$ concatenation levels, the number of internal nodes per code block is $`(5^L-1)/4`$. The explicit witness therefore has

```math
n=k5^L,
\qquad
m_T\le42kJ\frac{5^L-1}{4},
```

up to the convention about including an unnecessary final recovery inside the counted memory. Choosing $`L`$ to make $`kJp_L\le\epsilon`$, and using $`J\le40\gamma D+1`$, gives

```math
\begin{aligned}
n&=O\!\left(k\left[1+\log\frac{k(\gamma D+1)}{\epsilon}\right]^{\log_2 5}\right),\\
m_T&=O\!\left(k(\gamma D+1)\left[1+\log\frac{k(\gamma D+1)}{\epsilon}\right]^{\log_2 5}\right).
\end{aligned}
```

This calculation is an inference from the project's explicit circuit and elementary recursion, not a theorem attributed to the above papers. The detailed syndrome table and the 42-$`T`$ circuit are verified separately in the project. The perfect code, concatenation, and seven-$`T`$ synthesis themselves are prior ingredients. The demonstrated rate matches the lower bound **up to a polylogarithmic factor**, not up to constants.

Ideal processing between noisy intervals is essential to this upper comparison: a full recursive decode/correct/re-encode operation is not exposed to extra noise proportional to its elementary gate depth. This is the same processing model as the lower bound, so the comparison is valid internally; it does not by itself establish a hardware circuit-level threshold or a physical gate-time cost.

### Fixed-circuit Tanner resource consequence

The [Tanner construction](../science/tanner_memory.md) now supplies the previously missing resource translation in this repository. It hardwires local syndrome minimization, eligible-codeword selection, deterministic ties and row/column decompositions as bounded-size Boolean circuits. Graph incidence fixes routing; constant iterations give a linear circuit and counted resettable workspace. Exact coherent compilation charges every prescribed T/T-dagger slot, with no protected controller. Resetting scratch immediately before ideal extraction supplies the current syndrome within the existing model.

Combining this compilation with the inherited residual invariant, correlated X/Z tail accounting, weak-noise batching and bounded-ratio code lengths gives, after a fixed number of inner five-qubit levels, throughout $`0\lt\gamma\le1/20`$,

```math
n=O\!\left(k+\ln\frac{J}{\epsilon}\right),
\qquad m_T=O(nJ),
\qquad J=O(\gamma D+1),
\qquad\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

The mechanism is largely suggested by combining established good-code decoding and standard coherent simulation. The proved addition is the fixed-circuit compilation/accounting and its matching-order consequence for this processing converse. It closes the multiplicative logarithmic gap when $`k=\Omega(\ln(J/\epsilon))`$ and duration exceeds the width allowance. It does not determine the numerical coefficient, practical resource size or hardware threshold. For smaller payloads an additive reliability logarithm remains. The extension to $`\gamma\le1/20`$ is a separate composition proof: fixed-level inner correction bounds each fresh outer Pauli error without assuming isotropy, then ideal outer correction restores the residual invariant. Ordinary constant concatenation and standard coherent correction remain inherited tools; the addition is the exact same-model resource composition.

### Quantum storage with zero interior non-Clifford gates

The separate Clifford-only construction in this project stores syndrome histories in exponentially growing noisy branching records and performs the nonlinear decoding at the free final endpoint. Its guarantee is for an unknown quantum state, including an external reference. For fixed payload and fidelity, it realizes lifetime proportional to $`\gamma^{-1}\log n`$ at sufficiently small noise, establishing the order of the zero-$`T`$ width/lifetime tradeoff in the same ideal-processing model. It does not supply a polynomial-width memory with the optimal non-Clifford rate.

The construction combines established ingredients. The stabilizer existence and binary-geometric code framework come from Calderbank, Rains, Shor, and Sloane, [PRL **78**, 405–408 (1997)](https://doi.org/10.1103/PhysRevLett.78.405) and [IEEE Transactions on Information Theory **44**, 1369–1387 (1998)](https://doi.org/10.1109/18.681315), keys *calderbank1997orthogonal* and *calderbank1998gf4*. The 1997 paper's [corrected arXiv v3](https://arxiv.org/abs/quant-ph/9605005v3) is the relevant version for the rate statement. Classical independent-edge broadcasting is an established model; see Evans, Kenyon, Peres, and Schulman, [Annals of Applied Probability **10**, 410–433 (2000)](https://doi.org/10.1214/aoap/1019487349), *evans2000broadcasting*. The project's particular finite memory schedule, quantum-channel guarantee, and resource accounting are what require its own proof. Classical broadcast reconstruction and random stabilizer-code existence individually are not new.

## 6. Standard methods and other nearby work

**Operator algebras and channel fidelity.** The finite-dimensional matrix/direct-sum structure behind subsystem correction is standard operator quantum error correction; see Kribs, Laflamme, Poulin, and Lesosky, *Quantum Information & Computation* **6**, 382–399 (2006), [publisher record](https://www.rintonpress.com/journals/doi/QIC6.4-5-6.html), `kribs2006operator`. Schmidt number and its behavior under local operations are established by Terhal and Horodecki, PRA **61**, 040301(R) (2000), [DOI](https://doi.org/10.1103/PhysRevA.61.040301), `terhal2000schmidt`.

The project uses the dimension $`2^r`$ of a binary-rank-$`r`$ Pauli algebra as a superoperator-rank bottleneck. For a logical channel on dimension $`K=2^k`$, its proof derives $`F_{\mathrm e}=\mathrm{Tr}_{\mathrm{super}}(\Lambda)/K^2\le\mathrm{rank}(\Lambda)/K^2`$.

The fidelity identity is standard: see Nielsen, [Physics Letters A **303**, 249–252 (2002), Eq. (16)](https://arxiv.org/pdf/quant-ph/0205035), *nielsen2002simple*, and the original entanglement-fidelity formulation of Schumacher, [PRA **54**, 2614–2628 (1996)](https://doi.org/10.1103/PhysRevA.54.2614), *schumacher1996sending*. The subsequent rank inequality uses that a channel has spectral radius at most one. This yields the sharper branch bound $`F_{\mathrm e}\le\min\{1,2^{r-2k}\}`$ and the rank floor stated above.

The Pauli algebra's coherent matrix blocks also have size at most $`2^{r/2}`$, giving a complementary Schmidt-number bound. The channel-rank argument supplies the $`4^{-k}`$ replacement baseline and the theorem's scalar rank floor. Retaining the coherent block size also detects purely classical algebras, even when their total rank is large. General algebra and Schmidt-number theory are established ingredients; the rank contraction and its dynamical resource consequences are the new application.

**Magic-resource lower bounds.** Beverland, Campbell, Howard, and Kliuchnikov, *Quantum Science and Technology* **5**, 035009 (2020), [arXiv:1904.01124](https://arxiv.org/abs/1904.01124), `beverland2020lower`, study resource-state conversion, synthesis, and computational lower bounds with stabilizer operations free. Their setting allows adaptive measurement and charges nonstabilizer resource states. Because the present model permits arbitrary replacement states but restricts coherent dynamics, its memory bound is not a direct application of those state-magic monotones.

**Recent neighboring questions.** Tang et al.'s [arXiv:2610.12098](https://arxiv.org/abs/2610.12098), `tang2026purification`, bounds the $`T`$-count of a specified symmetry-based purification task. It does not state a universal noisy-memory lifetime/work law. Baumann, Stein, and Linnhoff-Popien's [arXiv:2610.26228](https://arxiv.org/abs/2610.26228), `baumann2026ceiling`, studies memory ceilings with a limited supply of fresh qubits and more general dynamics; that reset-budget problem differs from the unlimited-reset, restricted-processing problem here. They address related resource questions with different operational constraints.

**Autonomous is not automatically measurement-free in this sense.** Midha and Abanin's [arXiv:2609.38139](https://arxiv.org/abs/2609.38139), `midha2026autonomous`, studies local inference for quantum LDPC codes with idealized syndrome information and an offline error-decoding setting. Its classical processing and noise assumptions do not supply a fully noisy coherent memory construction for this comparison.

The phrase “memory-magic exchange” also occurs in [arXiv:2609.37368](https://arxiv.org/abs/2609.37368), but the memory there is classical storage in a streaming compiler for rotation angles. It is not a direct prior on preserving an unknown noisy quantum state. Its resource tradeoff therefore has a different operational meaning.

## 7. Delayed-basis guessing and existing recovery reductions

[Lattice Memory Certification](https://github.com/GoGoKo699/Lattice-Memory-Certification) bounds exact BB84 whole-string guessing with bounded retained quantum dimension and unlimited classical records. Its finite theorem gives $`p\le(2^s/2^k)^{1/6}`$ for $`k`$ inputs and $`s`$ retained quantum qubits. This project owns the noisy full-output envelope and dynamical charge ledger. The [shared transfer](../science/guessing_processing.md) combines them through a genuine CPTP factorization into $`s`$ quantum qubits and $`c`$ classical bits at each sampled cut, with binary rank $`r=2s+c`$. The finite classical alphabet also gives $`p\le2^{s+c-k}`$. Their concave majorant yields an above-chance expected-rank floor at every cut and hence the existing gate and contiguous-interval laws.

The relevant prior comparison includes **Dupuis, Fawzi, and Wehner**, “Entanglement sampling and applications,” IEEE Transactions on Information Theory **61**, 1093–1112 (2015), [arXiv:1305.1316v3, Eqs. (13)–(15), (27) and Appendix C](https://arxiv.org/abs/1305.1316v3), `dupuis2015sampling`; **Renes**, “Better bounds on optimal measurement and entanglement recovery, with applications to uncertainty and monogamy relations,” PRA **96**, 042328 (2017), [arXiv:1707.01114, Eq. (12)](https://arxiv.org/abs/1707.01114), `renes2017better`; and **Berta, Coles, and Wehner**, “Entanglement-assisted guessing of complementary measurement outcomes,” PRA **90**, 062127 (2014), [arXiv:1302.5902v2, Eq. (B23)](https://arxiv.org/abs/1302.5902v2), `berta2014guessing`.

These reductions are useful at scores at or below one half. Write $`u`$ for pretty-good guessing success and $`F_{\rm pg}`$ for squared pretty-good recovery fidelity. Retaining the normalization of DFW's spectral weights gives $`u\le\tau_\ell+A_\ell F_{\rm pg}`$, where $`\tau_\ell=2^{-\ell-1}`$ for $`\ell\lt k`$, $`\tau_k=0`$, and $`A_\ell=\sum_{j=0}^\ell\binom kj(1-2^j\tau_\ell)`$. This is a direct consequence of their spectral argument and improves the simpler threshold estimate in their Eq. (27).

Combine it with Renes's $`p^2+(1-p)^2/(2^k-1)\le u`$ and the channel-rank cap **before** averaging noise. [The complete comparison](../science/guessing_processing.md#5-comparison-with-existing-fidelity-and-recovery-reductions) proves that the prior route already gives the same floor $`R_t\ge2k[p-2^{-k}]_+/(1-2^{-k})`$ for every $`k\ge8`$: an analytic argument covers $`k\ge11`$ and 54 exact rational certificates cover $`k=8,9,10`$. For $`1\le k\le7`$, another 49 exact interior-rank certificates recover the displayed polygon floor; its endpoints follow directly. Thus the convenient mean-rank score floors and their gate and interval consequences already follow from this recovery route. The direct factorization gives a shorter operational proof. This comparison does not determine optimal guessing or compare every possible sharper integer-algebra bound.

The [exponential-moment extension](../science/moment_processing.md) combines the direct quantum-plus-classical profile and complementary inherited recovery bounds with a stronger dynamical estimate. Independent product replacement and a pivot-coordinate injection give $`\mathbb E[z^{r'}\mid V]\le(1-q+qz^2)^{r/2}`$; Jensen then contracts $`\log_z\mathbb E[z^r]`$ by a fixed factor. This preserves actual noise probabilities that a mean-rank secant discards. A fixed recovery of the averaged channel supplies the all-cut potential floor. The resulting finite and long-duration processing consequences are not additional pointwise guessing inequalities, and the convenient floors compared above remain recovery-implied. The fixed-score leading coefficient improves while the weak-noise order remains $`k\gamma D`$.

All these consequences retain the device-model assumptions. Neither a score nor the proof certifies independent noise, controller restrictions, or the absence of protected records. The LWE-based remote-preparation compilation in the sister project is not imported into this physical processing statement.


The finite moment comparison keeps the recovery inequalities branchwise before concavification in $`z^r`$. Its exact certificate bounds the strength of the explicitly specified recovery-only moment family over all $`z\gt 1`$: at the twelve-qubit example its two ledgers give rounded counts at most 223 and 52. The direct quantum-dimension/output-dimension/classical certificate yields 244 and 58. These are comparisons of converse formulas, not physical gate upper bounds or universal recovery optimality statements. The known classical cap used in the finite witness has a self-contained Gram-row-sum proof in the moment chapter.


The [quantitative guessing-decay corollary](../science/guessing_decay.md) is an inversion and controlled asymptotic consequence of this existing moment ledger. It retains the finite width reserve, gives an explicit exponent positive below $`2(1-q)/q`$, and compares the simple floor with the full polygon, integer-algebra profiles and the specified recovery reductions. Propagating the full generating-function argument supplies a stronger explicit exponent with the same threshold. This uses the existing branchwise guessing polygon and noise pivot bound; the additional step is the profile semigroup and fixed-slot convolution ledger. It is not a new pointwise uncertainty theorem or a processing-capacity achievability claim.
