# Questions for experts

These are specific questions arising from our audit (revised 25 November 2026) of [openai/math](https://github.com/openai/math), not claims that the manuscripts contain errors. The audits were produced by AI agents under the direction of a non-mathematician (see the README). A short correction, reference, or assessment of any stated implication would be valuable; please open an issue on this repository.

In every case the supporting detail is in the report named at the end of the item. Reports are in Japanese; file and line references to the LaTeX sources of openai/math (commit `7ef5c0c`, 26 October 2026) are language-independent.

## Analytic number theory

**1. Does Theorem 3.1 of result 011 imply an Elliott–Halberstam-type estimate near level 1?**
Result 011 (Poisson–Dirichlet law for the prime factors of p − 1) proves a marked Type II estimate, Theorem 3.1, in which F is any function bounded by 1 and invariant under multiplication by primes in specified bands (`03-determinant.tex:37–58`). With α, β and the marked weight 𝒲 exactly as in Theorem 3.1, we infer, for fixed 0 < ε < 0.1,

Σ_{x^0.9 < ℓ < x^{1−ε}, ℓ prime} |Σ_{mn ≡ 1 (mod ℓ)} α_m β_n 𝒲(mn − 1)| ≪ X (log x)^{−D_*}.

Does this marked, fixed-residue-class consequence follow, and is a comparable estimate known? We did not identify one. The step that would carry it is the coefficient-independent moment bound, Proposition 4.1 (`04-memory.tex:11–35`), which we checked only in part. → `audits/011/REPORT.md`, §2.

**2. Can the zero-free input for results 029 and 142 be independently assessed?**
Theorem 1.2 of result 029 asserts that finite-order Hecke L-functions over cyclotomic fields containing μ₁₂ have no zeros in Re s > 1 − 10⁻⁶. Result 029 (Artin's primitive-root conjecture, lower bound x/(log x)²) and result 142 (deterministic polynomial factorization over prime fields "without GRH", which quotes it as Theorem 11.1 without proof, `10-analytic.tex:4–18`) both rest on it.
If correct, [Axiom Math's Lean formalization](https://github.com/AxiomMath/DedekindZetaZeroFree) (October 2026; Cureton–Thorner, with AxiomProver) of a 31/32 zero-free half-plane for all Dedekind zeta functions and finite-order Hecke L-functions supplies the input required by both papers. We have not built it and are not aware of a public review by a domain expert. It reuses analytic number theory code from openai/math, so it is not known to be independent of OpenAI's argument.
Separately, our audits have not verified result 029's own proof of Theorem 1.2. Neither our agent nor a second model from a different developer located an error in it. The places we could not check are the specialization of Kazhdan–Patterson to the paper's conventions, the unfolding and interchange steps, and the reflection step (Lemma 6.3). → `audits/rh_alpoge/REPORT.md`, `audits/029core/REPORT.md`.

**3. Correction concerning the q = 2 case of result 142.**
An earlier version of this page asked whether Proposition 11.3 of 142 (`10-analytic.tex:91–98`) with n = q = 2 gives a polylogarithmic bound for the least prime quadratic non-residue. Assuming result 003's formalized uniform 7/8 theorem for Dirichlet L-functions (`lean/ComparatorChallenges/DirichletSevenEighths.lean` in openai/math), the zero-free hypotheses needed for this case already follow, since Q(μ₂₄, √p) is abelian over Q. The stated polylogarithmic bound therefore requires no stronger zero-free input than 003. This does not verify the remaining argument of Proposition 11.3. We withdraw the question; we and our reviewing model had missed this. → `audits/rh_alpoge/REPORT.md`.

**4. Is the algebraic reduction of result 142 valid?**
Does Theorem 1.2 of result 142 correctly give deterministic factorization in time polynomial in B and the numerical sizes of the supplied auxiliary primes? This reduction is independent of result 029's zero-free proof and remains to be assessed even if Axiom Math's formalization supplies the analytic input. If valid, the resulting polynomial-time algorithm under GRH would already improve on Evdokimov's quasi-polynomial bound. → `audits/142/REPORT.md`.

## Arithmetic geometry

**5. Are the cyclotomic and ring-class inputs in result 004's companion justified?**
The companion paper to 004 (pointwise 2-converse for E(Q)[2] ≠ 0; 004 is Hilbert's tenth problem over Q by Turing reduction) has two branches, by the parity of the 2^∞-Selmer corank. For non-CM curves, the corank-zero argument uses Kato's rational Theorem 12.5(3), including its local H² correction, at p = 2 (`cyclotomic.tex:423–438`); it explicitly does not use the stronger integral statement 12.5(4). CM curves use a separate cited result. Is the passage from this rational divisibility to the determinant comparison valid?
Result 004 needs the corank-one argument, which also uses the even input. Can a specialist assess the simultaneous arithmetic interface identified in Remark `rd:interface-status` (`ring-determinants.tex:1760–1768`)? → `audits/004b/REPORT.md`.

## Harmonic analysis

**6. Does the geometric theorem of result 077 alone give 3D Kakeya?**
Theorem 2.2 of the restriction paper (`01-geometric.tex:107`, label `geo:main`), applied with s = D = 1 to one segment per δ-separated direction and independent uniform time, appears to give |⋃T| ≳_{κ,η} δ^{10κ+η} for every η > 0, after accounting for the retained sublaw. Letting κ → 0 gives full Minkowski dimension. Is this reduction correct? If so, §§2–5 (about 4,700 lines, largely unchecked by us) would provide an alternative proof of the known 3D Kakeya set theorem (Wang–Zahl) without taking Wang–Zahl or OpenAI's Kakeya theorem as inputs. → `audits/077/REPORT.md`.

## Complexity theory

**7. Is the Turing machine of result 109 oblivious, and what follows for network coding?**
Result 109 claims integer multiplication in O(n (log n)^{1−κ}) on a multitape Turing machine. We could not confirm that head movement depends only on n: the adapters to the cited numerical routines strip leading zeros (`07-resampling.tex`, around l.215–238). If an oblivious implementation with the same bound exists, then with Afshani–Freksen–Kamma–Larsen (2019) it would refute the undirected k-pairs network coding conjecture.
Separately, if the chunk-swap lemma (`04-swap.tex:285–301`) and its shape-determined implementation are correct, setting P = G = B = 1 transposes an m×m bit matrix in O(m²(log m)^τ) time, where τ < 1. Our reviewer derived an extension of Adler et al.'s (2006) two-tape simulation to any fixed number of tapes; together with their Theorem 6.2, this would refute the same conjecture. Is that extension valid? The paper mentions neither work. → `audits/109/REPORT.md`, `audits/codex_review/review_ncc109.md`.

## Group theory and topology

**8. Is there an independent reason why the base group of result 285 is not a-T-menable?**
Result 285 claims a torsion-free counterexample to the Kadison–Kaplansky projection conjecture. If the claimed matrix construction is correct, the torsion-free base group G_base has a projection in M_N(C_r*(G_base)) with nonintegral unnormalized trace (`setup.tex:6–17`). By Higson–Kasparov and trace integrality on the assembly image, G_base cannot be a-T-menable. Is there an independent geometric proof of this from its construction? → `audits/285/REPORT.md`.

**9. Corollary 2.3 of result 321.**
The counterexample to Wall's finite D(2) problem rests on Corollary 2.3 (`02-tower.tex:230–235`), which states that, for a finite connected 2-complex Y, H₂(Y; C_ρ) = 0 forces ρ to be trivial on every torsion element. Is the tower argument proving this statement valid? If it is, the constructed group has deficiency 0 and its presentation has a relation gap. Bartholdi–Bridson–Mikhailov–Nicholson (arXiv:2610.23637) derive an analogous deficiency conclusion for different groups, conditionally on this corollary. → `audits/321/REPORT.md`.

## What we can offer in return

- For results 266, 268, 138, 101, 048 and 059, the finite computations were reproduced with independently written code, which is in the corresponding `audits/` directory.
- For result 029, exact checks of the local lemma (Lemma 6.2) at one prime, and numerical tests of a classical analogue of the reflection step over Q(ω), are in `audits/029local/`, `audits/029reflect/` and `audits/003dual/`. They do not test Lemma 6.3 itself.
