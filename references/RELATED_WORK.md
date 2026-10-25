# Literature, contribution, and operational comparison

**Project:** Measurement-Free-Memory-Cost  
**Literature checked:** 25 October 2026  
**Scope:** primary literature checked through this date. Bibliography keys refer to `references.bib`.

## 1. Contribution in relation to prior work

The main result is a quantitative requirement for **sustained coherent non-Clifford processing** in a noisy quantum memory. A finite amount of physical width supplies a logarithmic time allowance; maintaining a recoverable quantum payload beyond that allowance requires non-Clifford operations at an extensive rate. The statement remains meaningful when input-independent resets, even to magic states, are free.

For $`n`$ physical qubits, $`D`$ independent noise rounds, $`q=1-\gamma\in(0,1)`$, and a recovered $`k`$-qubit channel of entanglement fidelity at least $`1-\epsilon`$, the gate-budget theorem is

```math
m_T\ge a\left[D(-\ln q)-\ln\frac{2n}{a}\right]_+,
\qquad
a = \frac{2k(1-\epsilon-4^{-k})}{1-4^{-k}}.
```

The nontrivial regime is $`\epsilon<1-4^{-k}`$. At fixed acceptable infidelity and beyond the width allowance, this gives $`m_T=\Omega(k\gamma D)`$. The companion zero-$`T`$ statement gives replacement-channel convergence on the scale $`\gamma^{-1}\log(n/\epsilon)`$, including a diamond-norm formulation. These statements are proved in the project's Pauli-envelope and control-cost chapters.

The new element is the rank contraction and the resulting quantitative lifetime/work laws. The idea that Clifford processing with resets cannot sustain memory is prior work. Pauli symplectic coordinates, operator-algebra decomposition, Schmidt-number bounds, concatenated codes, and seven-$`T`$ Toffoli circuits are standard ingredients.

## 2. Direct prior: Nelson, Rajakumar, and Gullans

The essential predecessor is **Nelson, Rajakumar, and Gullans, “Non-Clifford Gates Are Required for Long-Term Memory,” CCC 2026**, DOI [10.4230/LIPIcs.CCC.2026.38](https://doi.org/10.4230/LIPIcs.CCC.2026.38), first posted as [arXiv:2510.08451](https://arxiv.org/abs/2510.08451) in October 2025. The comparison here uses the final proceedings paper (`nelson2026nonclifford`).

Its Theorem 1 and Lemma 6 give loss of pairwise trace distinguishability after $`O(\gamma^{-1}\log n\log(n/\epsilon))`$ layers. Section 2.3, Eq. (4), already permits replacement by **arbitrary single-qubit states**; Eqs. (8)–(11) specify the reset adjoint. Fresh magic states are therefore already allowed. Backward Pauli evolution and stochastic complete-depolarization trajectories also belong to this predecessor. The project sharpens the dependence to a single logarithm and obtains an explicit extensive gate-budget lower bound. Allowing joint replacement states is a model extension, not the headline novelty.

The public arXiv record still listed only v1 when checked. A related four-author QCTiP 2026 abstract adds Dominik Hangleiter and uses the title “Limitations on Measurement-free Fault-tolerant Protocols using Clifford Circuits”; it states broader qualitative results but no quantitative replacement for the published bound. See the [conference abstract, p. 38](https://qctipconf.github.io/qctip_2026/QCTiP2026_Book_of_Abstracts_002.pdf). This audit cannot exclude an unpublished refinement.

## 3. Model boundaries that matter for comparison

The project's depolarizing layer is

```math
\mathcal D_\gamma=(1-\gamma)\mathrm{id}
   +\gamma\,R_{I/2},
```

independently on every persistent qubit. Its processing intervals can contain arbitrary ideal Clifford gates, replacement resets, and counted $`T/T^\dagger`$ gates. Consequently, $`D`$ is a number of noisy intervals, **not automatically the depth of an elementary noisy gate circuit**.

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

### A promising but nonmatching low-overhead comparison

**Guo, Hong, Kaufman, and Lucas**, “Toward Self-Correcting Quantum Codes for Neutral Atom Arrays,” PRX Quantum **7**, 010301 (2026), [final paper](https://journals.aps.org/prxquantum/pdf/10.1103/mfmt-fwkg), `guo2026toward`, is particularly useful. Appendix D implements local greedy decisions using fixed-size reversible majority gadgets and a bounded-color schedule. **Our inference:** after compiling those fixed gadgets, a sweep uses $`O(n)`$ non-Clifford gates. This is a concrete possible route to low work overhead.

It is not a matching theorem here. The sustained-decoder results are numerical and phenomenological, and the paper leaves full noisy coherent-decoder analysis for future work. It does not prove a suitable finite-rate code family with arbitrary-duration storage in the precise model above. Even a proved $`O(kD)`$ cost at a fixed $`\gamma>0`$ would not establish a bound uniform in weak noise of the form $`O(k\gamma D)`$; the dependence of constants on $`\gamma`$ matters.

## 5. Same-model constructive comparison supplied by this project

The project includes an explicit upper construction rather than relying on an informal threshold-theorem conversion. Its ingredients are the standard perfect $`[[5,1,3]]`$ code, concatenation, coherent correction, and periodic batching. The perfect code was introduced by **Laflamme, Miquel, Paz, and Zurek**, PRL **77**, 198–201 (1996), [DOI](https://doi.org/10.1103/PhysRevLett.77.198), and **Bennett, DiVincenzo, Smolin, and Wootters**, PRA **54**, 3824–3851 (1996), [DOI](https://doi.org/10.1103/PhysRevA.54.3824); keys `laflamme1996perfect` and `bennett1996mixed`.

The explicitly derived syndrome correction polynomials use three Toffoli and three CCZ gates per internal decoding node. Their 42-$`T`$ total uses exact seven-$`T`$ decompositions. **Selinger**, PRA **87**, 042302 (2013), Eq. (1), p. 1 of the [author manuscript](https://www.mathstat.dal.ca/~selinger/papers/toffoli.pdf), explicitly displays the three-wire seven-$`T`$ Toffoli; Eq. (3) credits Amy et al.'s depth improvement. Selinger's separate $`T`$-depth-one construction uses four extra ancillas and is not the construction invoked here. See `selinger2013tdepth` and `amy2013meet`.

For $`0<\gamma\le1/20`$, take

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

### What is still open

No inspected primary source supplied a proof of $`O(k\gamma D)`$ total $`T`$-count, uniformly in small $`\gamma`$, at long horizons with every persistent data and syndrome register included in the noisy width. This is a bounded literature finding, not a proof that no such construction exists. Eliminating the explicit polylogarithmic factor while preserving the exact resource accounting remains an operational tightness question. Constant-rate codes and efficient coherent decoding are plausible routes, but neither classical RAM decoding time nor ideal classical side processing may be treated as free coherent circuitry.

### Quantum storage with zero interior non-Clifford gates

The separate Clifford-only construction in this project stores syndrome histories in exponentially growing noisy branching records and performs the nonlinear decoding at the free final endpoint. Its guarantee is for an unknown quantum state, including an external reference. For fixed payload and fidelity, it realizes lifetime proportional to $`\gamma^{-1}\log n`$ at sufficiently small noise, establishing the order of the zero-$`T`$ width/lifetime tradeoff in the same ideal-processing model. It does not supply a polynomial-width memory with the optimal non-Clifford rate.

The construction combines established ingredients. The stabilizer existence and binary-geometric code framework come from Calderbank, Rains, Shor, and Sloane, [PRL **78**, 405–408 (1997)](https://doi.org/10.1103/PhysRevLett.78.405) and [IEEE Transactions on Information Theory **44**, 1369–1387 (1998)](https://doi.org/10.1109/18.681315), keys *calderbank1997orthogonal* and *calderbank1998gf4*. The 1997 paper's [corrected arXiv v3](https://arxiv.org/abs/quant-ph/9605005v3) is the relevant version for the rate statement. Classical independent-edge broadcasting is an established model; see Evans, Kenyon, Peres, and Schulman, [Annals of Applied Probability **10**, 410–433 (2000)](https://doi.org/10.1214/aoap/1019487349), *evans2000broadcasting*. The project's particular finite memory schedule, quantum-channel guarantee, and resource accounting are what require its own proof. Classical broadcast reconstruction and random stabilizer-code existence individually are not new.

## 6. Standard methods and other nearby work

**Operator algebras and channel fidelity.** The finite-dimensional matrix/direct-sum structure behind subsystem correction is standard operator quantum error correction; see Kribs, Laflamme, Poulin, and Lesosky, *Quantum Information & Computation* **6**, 382–399 (2006), [publisher record](https://www.rintonpress.com/journals/doi/QIC6.4-5-6.html), `kribs2006operator`. Schmidt number and its behavior under local operations are established by Terhal and Horodecki, PRA **61**, 040301(R) (2000), [DOI](https://doi.org/10.1103/PhysRevA.61.040301), `terhal2000schmidt`.

The project uses the dimension $`2^r`$ of a binary-rank-$`r`$ Pauli algebra as a superoperator-rank bottleneck. For a logical channel on dimension $`K=2^k`$, its proof derives $`F_{\mathrm e}=\mathrm{Tr}_{\mathrm{super}}(\Lambda)/K^2\le\mathrm{rank}(\Lambda)/K^2`$.

The fidelity identity is standard: see Nielsen, [Physics Letters A **303**, 249–252 (2002), Eq. (16)](https://arxiv.org/pdf/quant-ph/0205035), *nielsen2002simple*, and the original entanglement-fidelity formulation of Schumacher, [PRA **54**, 2614–2628 (1996)](https://doi.org/10.1103/PhysRevA.54.2614), *schumacher1996sending*. The subsequent rank inequality uses that a channel has spectral radius at most one. This yields the sharper branch bound $`F_{\mathrm e}\le\min\{1,2^{r-2k}\}`$ and the rank floor stated above.

The Pauli algebra's coherent matrix blocks also have size at most $`2^{r/2}`$, giving a complementary Schmidt-number bound. Keeping only the total rank yields a weaker scalar floor, while retaining the coherent block size also detects purely classical algebras. General algebra and Schmidt-number theory are established ingredients; the rank contraction and its dynamical resource consequences are the new application.

**Magic-resource lower bounds.** Beverland, Campbell, Howard, and Kliuchnikov, *Quantum Science and Technology* **5**, 035009 (2020), [arXiv:1904.01124](https://arxiv.org/abs/1904.01124), `beverland2020lower`, study resource-state conversion, synthesis, and computational lower bounds with stabilizer operations free. Their setting allows adaptive measurement and charges nonstabilizer resource states. Because the present model permits arbitrary replacement states but restricts coherent dynamics, its memory bound is not a direct application of those state-magic monotones.

**Recent neighboring questions.** Tang et al.'s [arXiv:2610.12098](https://arxiv.org/abs/2610.12098), `tang2026purification`, bounds the $`T`$-count of a specified symmetry-based purification task. It does not state a universal noisy-memory lifetime/work law. Baumann, Stein, and Linnhoff-Popien's [arXiv:2610.26228](https://arxiv.org/abs/2610.26228), `baumann2026ceiling`, studies memory ceilings with a limited supply of fresh qubits and more general dynamics; that reset-budget problem differs from the unlimited-reset, restricted-processing problem here. They address related resource questions with different operational constraints.

**Autonomous is not automatically measurement-free in this sense.** Midha and Abanin's [arXiv:2609.38139](https://arxiv.org/abs/2609.38139), `midha2026autonomous`, studies local inference for quantum LDPC codes with idealized syndrome information and an offline error-decoding setting. Its classical processing and noise assumptions do not supply a fully noisy coherent memory construction for this comparison.

The phrase “memory-magic exchange” also occurs in [arXiv:2609.37368](https://arxiv.org/abs/2609.37368), but the memory there is classical storage in a streaming compiler for rotation angles. It is not a direct prior on preserving an unknown noisy quantum state. Its resource tradeoff therefore has a different operational meaning.
