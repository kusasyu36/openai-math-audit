# Independent audit of OpenAI's October 2026 AI-generated mathematics

On 6 October 2026, OpenAI released 720 AI-generated mathematics manuscripts ([openai/math](https://github.com/openai/math)). About 42% of the headline results come with Lean formalizations; the rest have not been machine-checked. Three manuscripts were withdrawn the next day, and all three were in the unformalized part.

This repository records an independent audit of selected unformalized results, plus Lean formalizations of some of the combinatorial parts.

> **How this was produced.** The audits and the Lean code were written by Anthropic's Claude (Opus 5.5) acting as research agents, under the direction of the repository owner. The owner is not a professional mathematician. Every finding labelled "confirmed" below was re-checked against the original source by the coordinating agent. None of the findings has been reviewed by a domain expert. Please treat them as leads for experts, not as verdicts.

## Method

- The audits use only our own code. **We never executed code from the openai/math repository**, including the verification scripts it ships, because checking a result with the author's own tools is not independent.
- For each manuscript we:
  1. stated the main theorem precisely and compared it with the literature;
  2. recomputed every explicit constant, exponent and finite identity with exact (rational or symbolic) arithmetic;
  3. searched for counterexamples numerically where possible;
  4. read the new arguments for gaps.
- Findings are classified as *definite error*, *gap / unjustified step*, *unclear exposition*, or *checked OK*.
- "No error found" means only that nothing was found in what we checked. Long probabilistic, analytic or arithmetic-geometric arguments were often checked only at the level of structure and bookkeeping.
- Reports are in Japanese, the working language of the project, in `audits/<result-number>/REPORT.md`. Result numbers follow openai/math's `CONTENTS.md`.

## Findings

### Confirmed issues

| Result | Manuscript | Issue |
|---|---|---|
| 201 | A Counterexample to Kurosh's Division-Ring Problem | **Definite (minor) error.** Theorem 2.1 assumes the base field contains all roots of unity. §3.1 starts from k₀ = Q and states "Thus it satisfies the hypotheses of Theorem 2.1". It does not. Starting from Q(μ∞) appears to repair this without affecting the main theorem. |
| 074 | Every four-dimensional Kakeya set has full Hausdorff dimension | **Definite (minor) error, citation only.** The paper cites "[24, Lemma 2.3]" of the 3D Kakeya maximal paper (l.188, l.931, bibliography entry [24]). In that paper 2.3 is *Definition 2.3 (Critical exponent)*; the content used matches its *Lemma 3.3*. The mathematics is not affected. |
| 171 | The hypercube Ramsey number has linear order | **Gap.** §9 ("adapted height construction") cannot apply the height-selection lemma of §3, because the lemma's radius hypothesis b > b₀ + (D+1)ρ fails in §9's regime. The paper instead adapts Steps 2–4 of the lemma's proof with new crowd tests, and asserts in one paragraph that "the same scale induction" holds. That induction is not written out, and §9's conclusions depend on it. |

### Other audited results (no error found)

| Result | Topic | Main dependency or caveat |
|---|---|---|
| 004 | Hilbert's tenth problem over Q (Turing reduction) | Unconditional in form, but relies essentially on OpenAI's unreviewed companion paper on a pointwise 2-converse (audited as 004b). |
| 004b | Pointwise 2-converse, E(Q)[2] ≠ 0 | 004 needs the rank-one branch (ring-class / ultrafilter / graph machinery), which is the hardest and least checked part. The p = 2 use of Kato's Euler system needs expert confirmation. |
| 029 | Artin's primitive-root conjecture, lower bound x/(log x)² | Rests on (i) a uniform zero-free region for Hecke L-functions over cyclotomic fields containing μ₁₂, which implies a quasi-Riemann hypothesis, and (ii) a sieve input imported from OpenAI's result 011. |
| 038 | Fujita's freeness conjecture, all dimensions | New variational method. A minor rigor gap in §5.2: admissibility of the test valuation v_ε. Novelty is concentrated in Prop 5.1. |
| 048 | Lipman–Zariski conjecture counterexample | All 68 explicit formulas reproduced. The surface is not explicit, and a single sign-convention chain is a single point of failure. |
| 059 | Zariski multiplicity conjecture counterexample | Finite arithmetic reproduced with our own code. The application of the cited real-Seifert-form theorem needs a check against its source. |
| 074 | 3D Kakeya maximal conjecture | Closing case analysis and exponents verified. Core §6–9 unaudited. |
| 077 | Fourier restriction, sphere, p > 3 | Exponents are within the necessary conditions, and the bookkeeping is consistent. The core geometric theorem alone would imply the 3D Kakeya set conjecture, and those ~4,700 lines are unaudited. |
| 101 | Simplex conjecture for isotropic constants, and asymmetric Mahler | About 220,000 quadratic-form checks, all scalar tables and counterexample searches passed. Some margins are extremely thin (down to 5·10⁻⁷). |
| 164 | Hindman's finite sums and products conjecture | §2 is sound. The analytic principles (§3–8) are unaudited. See the Lean section. |
| 284 | Randomized vs quantum query complexity | The paper proves R ≳ Q^{4−2/k}, i.e. it rules out every exponent below 4. The headline "optimal quartic separation" overstates this. |

We also checked 12 results that OpenAI marks as formalized (including 003, 005, 013, 145, 161, 189). The question was whether each Comparator statement in Lean says what the paper claims. We found no undisclosed mismatch.

## Lean formalizations (`lean/`)

The project uses Lean v4.34.1 and Mathlib `d13f23b7`, the same versions as openai/math. All of it is our own code; nothing is copied from openai/math.

| File | Content | Status |
|---|---|---|
| `Koeru/Hindman164.lean` | Statements of Theorem 1.1, with the separation condition (1.1), and Corollary 1.2 of result 164. Proofs that Thm 1.1 ⇒ the standard finite-sums-and-products conjecture, and that **Thm 1.1 ⇒ Cor 1.2**. A non-vacuity check for m = 1. | Proved, except Thm 1.1 itself (`sorry`, OpenAI's claim). |
| `Koeru/FiniteSums.lean` | Finite form of the finite sums theorem (Folkman–Rado–Sanders), derived from Mathlib's Hindman theorem by compactness. | Proved. |
| `Koeru/Ramsey.lean` | Finite hypergraph Ramsey theorem with per-colour targets, and a simultaneous version for all sizes ≤ q. Mathlib has no Ramsey theorem. | Proved. |
| `Koeru/Lemma25.lean` | Lemma 2.5 of result 164 (selection of a product chain). Our block layout differs from the paper's; only the existence of n is claimed. | Proved. |

Except for `theorem_1_1`, every theorem depends only on `propext`, `Classical.choice` and `Quot.sound` (see the `#print axioms` lines).

To build:

```
cd lean
lake exe cache get
lake build
```

This code is not intended for Mathlib. Mathlib's contribution rules require human-written, expert-supervised code.

## Licence

Apache-2.0 (see `LICENSE`). The audited manuscripts belong to their authors and are not redistributed here. The reports quote them only briefly.
