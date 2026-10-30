# Independent audit of OpenAI's October 2026 AI-generated mathematics

On 6 October 2026, OpenAI released a collection of AI-generated mathematics manuscripts ([openai/math](https://github.com/openai/math)); its catalogue now lists about 720. OpenAI reports that about 42% of the top-line results (300 of 719, per its `history.md`) have accompanying Lean formalizations. Its history also records the withdrawal of three manuscripts on 7 October 2026.

This repository records an independent audit of selected unformalized results, plus Lean formalizations of some of the combinatorial parts.

> **How this was produced.** The audits and the Lean code were written by Anthropic's Claude (Opus 5.5) acting as research agents, under the direction of the repository owner. The owner is not a professional mathematician. Every finding listed under "Issues" below was re-checked against the original LaTeX source by the coordinating agent, and then reviewed independently by OpenAI's Codex (GPT-6.1), a model from a different developer, against the same sources. That review overturned one of our original findings (see "Withdrawn finding"). None of the findings has been reviewed by a domain expert. Please treat them as leads for experts, not as verdicts.

## Method

- For the computational checks reported here we used only independently written scripts. We did not run the code or verification scripts shipped with openai/math. Reproducing the supplied computations would also have value, but it is a different check from ours.
- Depending on the manuscript, we:
  1. stated the main theorem precisely and compared it with the literature;
  2. recomputed selected constants, exponents and finite identities, using exact (rational or symbolic) arithmetic where possible and high-precision numerics or sampling otherwise;
  3. searched for counterexamples numerically where possible;
  4. read selected new arguments for gaps.
- Much of the reading used text extracted from the PDFs, and line numbers in the reports refer to that text unless a `.tex` file is named. PDF extraction can drop overlines and other marks (this caused the withdrawn finding below), so every finding under "Issues" was re-verified against the LaTeX source.
- Findings are classified as *definite error*, *gap / unjustified step*, *unclear exposition*, or *checked OK*.
- "No error found" means only that nothing was found in what we checked. Long probabilistic, analytic or arithmetic-geometric arguments were often checked only at the level of structure and bookkeeping.
- Reports are in Japanese, the working language of the project, in `audits/<result-number>/REPORT.md`. Result numbers follow openai/math's `CONTENTS.md`.

## Findings

### Issues

| Result | Manuscript | Issue |
|---|---|---|
| 074 | Every four-dimensional Kakeya set has full Hausdorff dimension | **Definite (minor) error, citation only.** In the version dated 24 September 2026, the paper cites "Lemma 2.3" of the 3D Kakeya maximal paper three times (`01-introduction.tex:135`, `02-framework.tex:930`, `05-scalar-tools.tex:79`, and the bibliography note). In the 3D paper, whose theorem environments share one counter per section, 2.3 is *Definition 2.3 (Critical exponent)*; the content used matches its *Lemma 3.3 (Weighted set estimate in a line chart)*. We identified no mathematical consequence. Most of the 4D argument was not audited. |
| 171 | The hypercube Ramsey number has linear order | **Compressed argument (not an error).** §9 ("adapted height construction") cannot apply the height-selection lemma of §3 directly, because the lemma's hypothesis b > b₀ + (D+1)ρ fails in §9's regime; the paper says so. It lists the modifications to Steps 2–4 of the lemma's proof (three new crowd tests, an overlap bound across slices, new parameter constraints) and then states that "the same scale induction" follows. The adapted induction is not written out. We found no point where it fails, but §9's conclusions depend on it, so writing it out would make the step checkable. The report also notes sketch-level or incompletely checked constructions in §§5, 16 and 18; much of the probabilistic machinery was checked only structurally. |

### Withdrawn finding

| Result | What we claimed | Why it was wrong |
|---|---|---|
| 201 | A Counterexample to Kurosh's Division-Ring Problem: §3.1 starts from k₀ = Q, which does not contain all roots of unity as Theorem 2.1 requires. | The LaTeX source has k₀ = Q̄ (`\overline{\mathbb Q}`), the algebraic closure of Q, which does contain them. Our text extraction of the PDF had dropped the overline. The paper is correct here. We thank the Codex review for catching this. |

### Additional audits: checks and open questions

No definite error was found in these. "No error found" is a statement about what we checked, not about the whole paper.

| Result | Topic | Main dependency or caveat |
|---|---|---|
| 004 | Hilbert's tenth problem over Q (Turing reduction) | Unconditional in form, but relies essentially on OpenAI's unreviewed companion paper on a pointwise 2-converse (audited as 004b). |
| 004b | Pointwise 2-converse, E(Q)[2] ≠ 0 | 004 needs the rank-one branch (ring-class / ultrafilter / graph machinery), which is the hardest and least checked part. The p = 2 use of Kato's Euler system needs expert confirmation. |
| 011 | Poisson–Dirichlet law for the prime factors of p − 1 (Granville's conjecture) | All 37 exponent and parameter checks pass, and the reduction, lifts and major-arc estimate check out. In Theorem 3.1 (a marked Type II estimate), F is any function bounded by 1 and invariant under multiplication by primes in the specified bands. We infer from this a marked Elliott–Halberstam-type consequence near level 1, with absolute values over prime moduli, well beyond known results; this inference needs expert confirmation. It is not a counterexample, but it calls for close scrutiny. The key new step, the coefficient-independent moment bound (Prop 4.1), was only partly checked, and our confidence in the full proof is low. |
| 029 | Artin's primitive-root conjecture, lower bound x/(log x)² | Rests on (i) a uniform zero-free region for Hecke L-functions over cyclotomic fields containing μ₁₂, which implies a quasi-Riemann hypothesis, and (ii) the Type II estimate of OpenAI's result 011 (Thm 3.1), imported verbatim. The imported statement matches 011's, so 029 inherits 011's main risk. Separately, the central zero-free-region argument was not verified. |
| 038 | Fujita's freeness conjecture, all dimensions | New variational method. A possible rigor gap in §5.2 concerns the admissibility of the test valuation v_ε; the report considers it likely repairable. Novelty is concentrated in Prop 5.1, which still needs expert review. |
| 048 | Lipman–Zariski conjecture counterexample | All 68 explicit formulas reproduced. The surface is not explicit, and a single sign-convention chain is a single point of failure. |
| 059 | Zariski multiplicity conjecture counterexample | Finite arithmetic reproduced with our own code. The application of the cited real-Seifert-form theorem needs a check against its source. For the companion paper in C⁴ (multiplicities four and five), the central integral Seifert-form "tube" computations were not checked. |
| 074 | 3D Kakeya maximal conjecture | Closing case analysis and exponents verified. Core §6–9 unaudited. |
| 077 | Fourier restriction, sphere, p > 3 | Exponents are within the necessary conditions, and the bookkeeping is consistent. By our reading, the core geometric theorem alone implies the 3D Kakeya set conjecture in Minkowski form. Its ~4,700 lines were largely unchecked, apart from selected estimates. |
| 101 | Simplex conjecture for isotropic constants, and asymmetric Mahler | About 220,000 sampled quadratic-form tests, the scalar tables and the counterexample searches found no violation. Some margins are extremely thin (down to 5·10⁻⁷). Transport regularity, the degree argument and the approximation and limiting steps still need expert review. |
| 164 | Hindman's finite sums and products conjecture | We found no issue in §2's reduction, which assumes the paper's two analytic principles. The analytic core (§3–8) remains largely unchecked. See the Lean section. |
| 201 | Kurosh's division-ring problem (counterexample) | The derivation of the main theorem (§3, §12) checks out, and an 11-block solvability test agreed numerically up to the obstruction dimension. The core §5–7 (response-space splitting, transfer, linkage, replication) was not checked in full. |
| 284 | Randomized vs quantum query complexity | For each fixed k ≥ 2 the paper claims R ≥ c_k Q^{4−2/k}/polylog(Q), so the exponent approaches 4. It does not give a single family with R ≥ Q⁴/polylog. The paper's title ("Nearly Quartic") states this accurately; the shorter description "optimal quartic separation" in the openai/math catalogue is stronger than what is claimed. |

We also inspected 12 results that OpenAI marks as formalized (including 003, 005, 013, 145, 161, 189). The question was whether each Comparator statement and its definitions in Lean say what the paper claims, allowing for disclosed restrictions (for example, 084 is formalized only for ratio 1/2). We found no undisclosed mismatch. We did not build or verify the proofs.

## Lean formalizations (`lean/`)

The project uses Lean v4.34.1 and Mathlib `d13f23b7`, the same versions as openai/math. All of it is our own code; nothing is copied from openai/math.

| File | Content | Status |
|---|---|---|
| `Koeru/Hindman164.lean` | Statements of Theorem 1.1, with the separation condition (1.1), and Corollary 1.2 of result 164. Proofs that Thm 1.1 ⇒ the standard finite-sums-and-products conjecture, and that **Thm 1.1 ⇒ Cor 1.2**. A non-vacuity check for m = 1. | Proved, except Thm 1.1 itself (`sorry`, OpenAI's claim). |
| `Koeru/FiniteSums.lean` | Finite form of the finite sums theorem (Folkman–Rado–Sanders), derived from Mathlib's Hindman theorem by compactness. | Proved. |
| `Koeru/Ramsey.lean` | Finite hypergraph Ramsey theorem with per-colour targets, and a simultaneous version for all sizes ≤ q. Mathlib has no Ramsey theorem. | Proved. |
| `Koeru/Lemma25.lean` | Lemma 2.5 of result 164 (selection of a product chain). Our block layout differs from the paper's; only the existence of n is claimed. The paper takes x to be a sequence of positive integers; we prove it for every x : ℕ → ℕ, which is stronger. | Proved. |

Except for `theorem_1_1`, every theorem depends only on `propext`, `Classical.choice` and `Quot.sound` (see the `#print axioms` lines).

To build:

```
cd lean
lake exe cache get
lake build
```

This code is not intended for Mathlib. Mathlib's policy on AI use closes AI-generated pull requests that lack expert supervision, and this code was not written under such supervision.

## Licence

Apache-2.0 (see `LICENSE`). The audited manuscripts belong to their authors and are not redistributed here. The reports quote them only briefly.
