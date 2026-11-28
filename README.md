# Measurement-Free Memory Cost

**How much non-Clifford processing is needed to keep unknown quantum information alive?**

This project studies a memory whose data and syndrome records all remain noisy, with no protected external classical controller. It gives payload- and reliability-dependent lower bounds on continuing non-Clifford processing, a sharp order for Clifford-only quantum lifetime, and explicit storage constructions. The necessary and achieved gate orders match at the width used by long batches of coherent recovery, uniformly at weak noise and including a one-qubit payload.

The model permits arbitrary ideal global Clifford operations and input-independent replacement resets between independent depolarizing-noise layers. Fresh reset states may include magic states. Initial encoding and final recovery may be arbitrary. Every internal T/T-dagger gate is counted, and every persistent register belongs to the noisy physical width.

## Run the coherent recovery

1. **Encode.** Choose a code that fits the payload and target reliability. For [Tanner memory](science/tanner_memory.md#7-available-lengths-padding-and-diamond-error), initialize spare logical inputs and apply the outer encoder. The [extended noise-range construction](science/tanner_memory.md#8-constant-inner-concatenation-extends-the-physical-noise-range) also encodes every outer position in a fixed number of five-qubit levels. Data, syndrome slots and idle work all undergo independent depolarizing noise.
2. **Recover between batches.** At each scheduled ideal gap, including after a shortened final batch, decode and correct the inner blocks when present. The [five-qubit node](science/coherent_recovery.md#1-a-five-qubit-recovery-using-42-t-gates) Clifford-decodes a block, applies its two syndrome polynomials directly as logical X and Z corrections using four controlled Clifford gates, three Toffolis and three controlled-controlled Z gates, then resets the four syndrome wires before reencoding. For Tanner recovery, use the decoded outer wires. Reset outer workspace, extract stabilizer signs coherently, run the fixed [reversible decoder](science/tanner_memory.md#4-coherent-compilation-reset-and-exact-channel), and apply its controlled X and Z corrections. Uncompute decoder work from the retained syndrome, reset the syndrome and restore the inner encoding when used.
3. **Read the payload.** Apply the final decoder. The returned register preserves the unknown input and its reference correlations within the stated error bound. Internal recovery uses coherent gates and input-independent resets, without measurement or postselection.

Count every noisy slot and every prescribed internal T or T-dagger gate, including uncomputation and rounds with zero syndrome. The Tanner circuit counts fixed lookup tables and assumes ideal global processing in the gaps. The separate [Clifford-only archive](science/clifford_lifetime.md#4-schedule-and-archive-dynamics) stores a noisy syndrome history for final processing.

## Delayed-basis guessing and processing

[Lattice Memory Certification](https://github.com/GoGoKo699/Lattice-Memory-Certification) bounds exact BB84 whole-string guessing through a quantum-plus-classical bottleneck. Here that bottleneck is supplied by the Pauli algebra of a sampled storage tail. Its unconditional guessing score enforces an expected-rank floor at every cut, which can replace the fidelity premise in the gate and contiguous-interval budgets.

The [complete transfer and recovery comparison](science/guessing_processing.md) live in this repository; Lattice supplies the [finite guessing theorem](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/main/science/guessing-bound.md) and its [operational application](https://github.com/GoGoKo699/Lattice-Memory-Certification/blob/main/science/noisy-processing.md). The device receives only the quantum challenge; the hidden string stays secret and only the bases arrive after storage. The common terminal envelope is the full physical output, with $`M=2n`$. All persistent controller records are noisy and counted, and their processing must use the charged circuit model. Noise flags remain proof aids.

The transfer is useful below one-half success. Existing recovery inequalities reproduce the convenient mean-rank score floors for all challenge sizes. The [exponential-moment bound](science/moment_processing.md) then strengthens the dynamical ledger by retaining the actual product-noise probabilities. At fixed positive guessing score and fixed noise strength, its long-storage coefficient approaches that of a full logical payload as the challenge grows. The [correlated replacement extension](science/guessing_processing.md#correlated-replacement) covers conditional erasure probabilities in the mean-rank ledger. The [classical-verifier composition](science/guessing_processing.md#classical-verifier-composition) transfers a score after accounting for its preparation simulation error. Both require independent evidence for their device premises.

The [full generating-function ledger](science/guessing_decay.md) strengthens the explicit decay exponent by propagating its argument through every noise layer, rather than contracting one fixed moment. Its optimized exponent reaches $`\ln2`$ at zero T rate and is positive exactly below the same necessary rate $`2(1-q)/q`$. This remains a converse, rather than an achievable processing capacity.

## Principal results

Let $`n\ge k\ge1`$ be integer physical and logical widths, $`D\ge0`$ the integer number of noisy rounds, and $`0\lt\gamma\lt1`$ the depolarizing strength in the convention $`\mathcal N_\gamma=(1-\gamma)\mathrm{id}+\gamma\mathcal R_{I/2}`$. Write $`q=1-\gamma`$ and $`\kappa=-\ln q`$.

For recovered Bell entanglement fidelity at least $`1-\epsilon`$, set

```math
a=2k\,\frac{1-\epsilon-4^{-k}}{1-4^{-k}},
\qquad 0\le\epsilon\lt 1-4^{-k}.
```

The necessary processing budget includes

```math
m_T\ge
a\left[
D\kappa
-\ln\!\left(\frac{2n}{a}\right)
\right]_+.
```

Here $`[x]_+=\max\{x,0\}`$. At fixed high fidelity and sufficiently long storage, the required work is at least of order $`k\gamma D`$ after a finite logarithmic width allowance. [Control cost](science/control_cost.md) proves sharper finite-depth and integer-gate bounds, together with a constraint on every contiguous storage interval.

The [reliability converse](science/reliability_processing.md) also retains the probability that the envelope loses every nonconstant input direction. For $`D\ge1`$ and $`0\lt\epsilon\lt1-4^{-k}`$, put $`\delta=\epsilon/(1-4^{-k})`$ and $`b=1-\delta`$. It gives

```math
m_T\ge\left[
\frac{\gamma Db}{\ln(1/\gamma)}
\left[\ln\frac{Db}{\delta}\right]_+-2n
\right]_+.
```

A rank-$`r`$ envelope with symplectic rank $`2s`$ has conditional extinction probability at least $`\gamma^{r-s}`$. The displayed rank-only consequence is a simple comparison bound. The full finite theorem retains the coherent-block floor, the exact survival product and the exiting-rank requirement, giving stronger necessary counts and interval bounds. Keeping total extinction below $`\delta`$ requires cumulative rank that must be replenished by charged gates.

| Result | Guarantee | Authoritative account |
|---|---|---|
| Clifford/reset forgetting | Half-diamond distance to a replacement channel is at most $`1-[1-(1-\gamma)^D]^{2n}`$ | [Binary Pauli envelopes](science/pauli_envelopes.md) |
| Continuing processing | Extensive non-Clifford budget with fidelity, width, and timing dependence | [Control cost](science/control_cost.md) |
| Reliability-dependent processing | Rank extinction forces an additional reliability logarithm at modest width; fixed-noise order matching includes small payloads | [Reliability processing](science/reliability_processing.md) |
| Uniform weak-noise bounds | Minimum width, a family of batch sizes and a matching gate order at long-batch width | [Uniform width and gate bounds](science/weak_noise_tradeoff.md) |
| Guessing-score premise | Delayed-basis whole-string success supplies an all-cut rank floor | [Guessing and processing](science/guessing_processing.md) |
| Exponential-moment strengthening | Actual product-noise probabilities give a stronger family of gate and interval bounds | [Moment processing](science/moment_processing.md) |
| Quantitative guessing decay | Full generating-function propagation and an explicit stronger exponential score bound, with finite width retained | [Guessing decay](science/guessing_decay.md) |
| Sharp zero-T quantum lifetime | $`\Theta(\gamma^{-1}\ln n)`$ at fixed logical size and target fidelity $`4^{-k}\lt F_0\lt 1`$, for sufficiently weak noise | [Clifford quantum lifetime](science/clifford_lifetime.md) |
| Five-qubit coherent recovery | Explicit 42-T node and concatenated reliability overhead, for $`\gamma\le1/20`$ | [Coherent recovery](science/coherent_recovery.md) |
| Constant-iteration Tanner recovery | Counted width $`O(k+\ln(J_+/\epsilon))`$ and prescribed T count $`O(nJ)`$ for $`\gamma\le1/20`$, using fixed-level inner protection | [Tanner memory](science/tanner_memory.md) |

For the constructive recovery bound with $`k,D\ge1`$, let $`\alpha=\log_2 5`$ and $`\ell=1+\ln[k(\gamma D+1)/\epsilon]`$. For $`0\lt \gamma\le1/20`$ and $`0\lt \epsilon\lt 1/2`$, an explicit concatenated five-qubit construction achieves

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr),
\qquad
\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

Each decoding node uses at most **42 T/T-dagger gates and no additional work qubits**. The logarithmic overhead belongs to this five-qubit construction. It retains the stated noise range.

For a fixed good quantum Tanner-code family, a fixed number of inner five-qubit concatenation levels permits a second coherent construction throughout $`0\lt\gamma\le1/20`$. Write $`J_+=\max\{1,J\}`$, including $`D=J=0`$. Its guarantee is

```math
n=O\!\left(k+\ln\frac{J_+}{\epsilon}\right),
\qquad m_T=O(nJ),
\qquad J=O(\gamma D+1),
\qquad\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

The [proof and compilation ledger](science/tanner_memory.md) counts syndrome, local tables, tie-breaking, updates, routing and resettable workspace. Constant iterations keep a residual error inside a correctable region, and the free endpoint decoder completes recovery.

At weak noise $`0<\gamma\le1/20`$, the [uniform width and gate bounds](science/weak_noise_tradeoff.md) give

```math
n_{\min}=\Theta\!\left(k+\frac{\ln(D/\epsilon)}{\ln(1/\gamma)}\right)
```

for $`D\ge1`$ and target error bounded above by a fixed constant below one half. For a batch length $`h`$ with $`1\le h\le\min\{D,\lfloor 1/(20\gamma)\rfloor\}`$, put $`J=\lceil D/h\rceil`$. The coherent Tanner circuit achieves

```math
n=O\!\left(k+\frac{\ln(2J/\epsilon)}{\ln(1/(\gamma h))}\right),
\qquad m_T=O(n(J-1))
```

with its allocated width $`n`$ and half-diamond error at most $`\epsilon`$. For long batches, write $`x=\gamma D`$ and $`L=\ln((1+x)/\epsilon)`$. At width proportional to $`k+L`$ with a sufficient construction allowance, the necessary and sufficient gate order is $`\Theta(x(k+L))`$ for sufficiently large $`x`$. The gate lower bound extends to $`n\le e^{x/8}(k+L)`$. As noise tends to zero, the minimum width can lie below the long-batch width, leaving the optimal gate count there open.

Gu, Tang, Caha, Choe, He and Kubica supply the decoder and repeated-round invariant. Existing passive LDPC and measurement-free correction supply the correction mechanism. The addition here is its explicit fixed-circuit resource consequence in the present model. Nelson, Rajakumar, and Gullans proved a quantitative Clifford-memory forgetting bound and introduced the stochastic Pauli framework used here. Bharti, Haug and Tanggara established the complete-erasure/Jensen mechanism for a cumulative spacetime law. The reliability converse adapts that mechanism to backward Pauli rank and converts cumulative rank into specifically charged T/T-dagger processing. The [literature comparison](references/RELATED_WORK.md) explains the sharper reference-stable contraction, the processing and timing laws, and the inherited coding and circuit ingredients.

## Read the science

Start with the [model, assumptions, and exact claims](science/model.md). Follow the envelope proof to the fidelity bottleneck and processing budget, then compare the archive and the two coherent constructions. [Interpretation](science/interpretation.md) explains the matching-order regime and the remaining physical and finite-constant questions.

| Reader question | Reading route |
|---|---|
| Why do arbitrary fresh magic states not remove this bound? | [Envelope updates](science/pauli_envelopes.md) → [processing budget](science/control_cost.md) |
| How is the quantum payload related to binary Pauli rank? | [Channel-rank and algebra bottlenecks](science/control_cost.md) |
| Why does a small payload still pay for long reliable storage? | [Absorbing rank extinction and reliability-dependent processing](science/reliability_processing.md) |
| Can a delayed-basis guessing score replace the fidelity premise? | [CPTP factorization, score floors and recovery comparison](science/guessing_processing.md) |
| What is gained by retaining the noise distribution? | [Exponential moments, fixed recovery and strengthened budgets](science/moment_processing.md) |
| How fast must success decay below the necessary processing rate? | [Finite inverse ledger and explicit exponent](science/guessing_decay.md) |
| Can a quantum state survive logarithmic time with zero internal T count? | [Noisy syndrome archive](science/clifford_lifetime.md) |
| What does a fully charged recovery circuit look like? | [Five-qubit correction](science/coherent_recovery.md) and [constant-iteration Tanner compilation](science/tanner_memory.md) |
| Which resources remain necessary even with more gates? | [Width, endpoints, and interpretation](science/interpretation.md) |
| Which parts come from existing work? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

## Learn the project

The single educational anchor shared with Lattice Memory Certification is **John Watrous, *The Theory of Quantum Information* (2018)**. Follow the [learning guide](learning/README.md), [source-to-project bridge](learning/WATROUS_BRIDGE.md), and [30 project-intrinsic milestones](learning/MILESTONES.md): 12 background/preliminaries and 18 narrative/innovation targets. The final milestones connect the constructions and the guessing-score corollary to the resource frontier. Primary stabilizer, uncertainty and recovery sources remain separately attributed research references.

[The LLM reading guide](llms.txt) maps scientific questions to authoritative files.

## Reproduce the finite checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
python3 verification/check_learning.py
```

The scientific checks cover Pauli-algebra and channel-rank identities, the exact finite-depth scalar optimization, the syndrome archive and workspace schedule, the explicit five-qubit recovery, the guessing-score transfer with its prior-recovery comparison, and the exponential-moment bounds with exact finite certificates. The learning check verifies the 30-milestone route and local links. [Verification](verification/README.md) records the environment, expected outputs, and scope.

The [finite resource comparison](verification/Graph/README.md) provides the
Python-generated width and processing curves, pinned dependencies, and a check
against the committed plot and numerical certificates.

The [finite processing lower-bound figure](verification/Graph/README.md#finite-processing-lower-bounds)
compares the refined reliability, rank-only reliability and payload bounds at
the same physical width, with independent checks of its labelled values.

Finite computations validate the stated finite identities. General converses, concatenation guarantees, and asymptotic resource laws rely on the written proofs.

## Scope

All persistent qubits undergo independent depolarizing noise. Processing between noise layers is ideal and may be global. The constructions do not include extra noisy rounds proportional to elementary gate depth, geometric locality, or a protected external controller. Final recovery lies outside the charged storage interval.

All constructions protect arbitrary unknown inputs and their reference correlations. The zero-T construction pays exponential noisy width and uses unrestricted final nonlinear decoding. Five-qubit concatenation retains a polylogarithmic reliability overhead. Fixed-level inner protection extends the Tanner recovery bound to $`\gamma\le1/20`$. The width–gate frontier at short batches, optimal numerical T count, practical constants and hardware costs remain open.

## Purpose and contact

The project is currently under collaboration, and the manuscript is being written. For questions, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).


