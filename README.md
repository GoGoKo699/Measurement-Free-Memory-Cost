# Measurement-Free Memory Cost

**How much non-Clifford processing is needed to keep unknown quantum information alive?**

This project studies a memory whose data and syndrome records all remain noisy, with no protected external classical controller. It gives a quantitative lower bound on continuing non-Clifford processing, a sharp order for Clifford-only quantum lifetime, and explicit storage constructions.

The model permits arbitrary ideal global Clifford operations and input-independent replacement resets between independent depolarizing-noise layers. Fresh reset states may include magic states. Initial encoding and final recovery may be arbitrary. Every internal T/T-dagger gate is counted, and every persistent register belongs to the noisy physical width.

## Explanatory report

Read the [complete explanatory report (PDF)](learning/report/Measurement_Free_Memory_Cost_Report.pdf) for a self-contained explanation of the memory model, the processing lower bound, the two storage constructions, and the remaining resource gap.

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

| Result | Guarantee | Authoritative account |
|---|---|---|
| Clifford/reset forgetting | Half-diamond distance to a replacement channel is at most $`\min\{1,2n(1-\gamma)^D\}`$ | [Binary Pauli envelopes](science/pauli_envelopes.md) |
| Continuing processing | Extensive non-Clifford budget with fidelity, width, and timing dependence | [Control cost](science/control_cost.md) |
| Sharp zero-T quantum lifetime | $`\Theta(\gamma^{-1}\ln n)`$ at fixed logical size and target fidelity $`4^{-k}\lt F_0\lt 1`$, for sufficiently weak noise | [Clifford quantum lifetime](science/clifford_lifetime.md) |
| Explicit coherent recovery | Width and T count with a polylogarithmic overhead relative to the extensive lower bound | [Coherent recovery](science/coherent_recovery.md) |

For the constructive recovery bound with $`k,D\ge1`$, let $`\alpha=\log_2 5`$ and $`\ell=1+\ln[k(\gamma D+1)/\epsilon]`$. For $`0\lt \gamma\le1/20`$ and $`0\lt \epsilon\lt 1/2`$, an explicit concatenated five-qubit construction achieves

```math
n=O(k\ell^\alpha),
\qquad
m_T=O\bigl(k(\gamma D+1)\ell^\alpha\bigr),
\qquad
\frac12\|\Lambda-\mathrm{id}_{2^k}\|_\diamond\le\epsilon.
```

Each decoding node uses at most **42 T/T-dagger gates and no additional work qubits**. The logarithmic overhead is an upper bound for this construction; optimal T count is not established.

Nelson, Rajakumar, and Gullans proved a quantitative Clifford-memory forgetting bound and introduced the stochastic Pauli framework used here. The [literature comparison](references/RELATED_WORK.md) explains the sharper reference-stable contraction, the processing and timing laws, and the inherited coding and circuit ingredients.

## Read the science

Start with the [model, assumptions, and exact claims](science/model.md). Follow the envelope proof to the fidelity bottleneck and processing budget, then compare the two constructions. [Interpretation](science/interpretation.md) explains the physical resource choices and remaining gap.

| Reader question | Reading route |
|---|---|
| Why do arbitrary fresh magic states not remove this bound? | [Envelope updates](science/pauli_envelopes.md) → [processing budget](science/control_cost.md) |
| How is the quantum payload related to binary Pauli rank? | [Channel-rank and algebra bottlenecks](science/control_cost.md) |
| Can a quantum state survive logarithmic time with zero internal T count? | [Noisy syndrome archive](science/clifford_lifetime.md) |
| What does a fully charged recovery circuit look like? | [Five-qubit correction and concatenation](science/coherent_recovery.md) |
| Which resources remain necessary even with more gates? | [Width, endpoints, and interpretation](science/interpretation.md) |
| Which parts come from existing work? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

## Learn the project

The single educational anchor is **Daniel Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation*** ([arXiv:0904.2557](https://arxiv.org/abs/0904.2557)). Follow the [source-to-project bridge](learning/GOTTESMAN_BRIDGE.md) and the [30 project-intrinsic milestones](learning/MILESTONES.md): 12 background/preliminaries and 18 narrative/innovation targets. Personal coaching and learner progress are kept outside this scientific repository.

[The LLM reading guide](llms.txt) maps scientific questions to authoritative files.

## Reproduce the finite checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The suite checks Pauli-algebra and channel-rank identities, the exact finite-depth scalar optimization, the syndrome archive and workspace schedule, and the explicit five-qubit recovery. [Verification](verification/README.md) records the environment, expected outputs, and scope.

The [finite resource comparison](verification/Graph/README.md) provides the
Python-generated width and processing curves, pinned dependencies, and a check
against the committed plot and numerical certificates.

Finite computations validate the stated finite identities. General converses, concatenation guarantees, and asymptotic resource laws rely on the written proofs.

## Scope

All persistent qubits undergo independent depolarizing noise. Processing between noise layers is ideal and may be global. The constructions do not include extra noisy rounds proportional to elementary gate depth, geometric locality, or a protected external controller. Final recovery lies outside the charged storage interval.

Both constructions protect arbitrary unknown inputs and their reference correlations. The zero-T construction pays exponential noisy width and uses unrestricted final nonlinear decoding. Continuing recovery pays internal T count and a polylogarithmic reliability overhead. The complete optimum between them remains open.

## Purpose and contact

The project is currently under collaboration, and the manuscript is being written. For questions, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

