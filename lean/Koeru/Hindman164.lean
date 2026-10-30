import Mathlib

/-!
# OpenAI 結果 164: Hindman の有限和・有限積予想（命題のみ）

出典: openai/math `preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026`
定理 1.1 と系 1.2 の命題を、論文の文言どおりに書く。定理 1.1 自体の証明は `sorry`（OpenAI の主張）。そこからの帰結は証明済み。

彩色は `χ : ℕ → Fin r`。論文の `χ : ℕ → [r]`（ℕ は正の整数）に対応する。
0 での値は使わない（a_j はすべて正）。
-/

namespace Koeru.Hindman164

open Finset

/-- 標準形の Hindman 予想（分離条件なし）:
有限彩色のもとで、任意の m について、m 個の相異なる正整数で、
空でない部分集合の和と積がすべて同じ色になるものがある。 -/
def FiniteSumsProductsConjecture : Prop :=
  ∀ (r : ℕ) (χ : ℕ → Fin r) (m : ℕ), ∃ (a : Fin m → ℕ) (c : Fin r),
    StrictMono a ∧ (∀ i, 0 < a i) ∧
    ∀ J : Finset (Fin m), J.Nonempty →
      χ (∑ j ∈ J, a j) = c ∧ χ (∏ j ∈ J, a j) = c

/-- 分離条件 (1.1): a₁ > R、かつ 2 ≤ d ≤ m で a_d > R (Σ_{k<d} a_k + Π_{k<d} a_k)^D。
添字は 0 始まり（論文の a_1 が `a 0`）。 -/
def Separated {m : ℕ} (R D : ℝ) (a : Fin m → ℕ) : Prop :=
  (∀ i : Fin m, i.val = 0 → R < a i) ∧
  ∀ d : Fin m, 1 ≤ d.val →
    R * (((∑ k ∈ univ.filter (· < d), a k : ℕ) : ℝ) +
         ((∏ k ∈ univ.filter (· < d), a k : ℕ) : ℝ)) ^ D < (a d : ℝ)

/-- 定理 1.1（OpenAI の主張）。r, m ≥ 1、R ≥ 2、D ≥ 1。 -/
def Theorem_1_1 : Prop :=
  ∀ (r m : ℕ), 1 ≤ r → 1 ≤ m → ∀ (R D : ℝ), 2 ≤ R → 1 ≤ D →
    ∀ χ : ℕ → Fin r, ∃ (a : Fin m → ℕ) (c : Fin r),
      StrictMono a ∧ (∀ i, 0 < a i) ∧
      (∀ J : Finset (Fin m), J.Nonempty →
        χ (∑ j ∈ J, a j) = c ∧ χ (∏ j ∈ J, a j) = c) ∧
      Separated R D a

/-- 有限和の集合 FS(A) と有限積の集合 FP(A)（A は有限集合）。 -/
def FS (A : Finset ℕ) : Finset ℕ :=
  (A.powerset.filter (·.Nonempty)).image (fun J => ∑ x ∈ J, x)

def FP (A : Finset ℕ) : Finset ℕ :=
  (A.powerset.filter (·.Nonempty)).image (fun J => ∏ x ∈ J, x)

/-- 系 1.2: 任意の有限彩色と m ≥ 1 に対し、m 元集合 A で FS(A) ∪ FP(A) が単色、
|FS(A)| = |FP(A)| = 2^m − 1、FS(A) ∩ FP(A) = A となるものがある。 -/
def Corollary_1_2 : Prop :=
  ∀ (r : ℕ) (χ : ℕ → Fin r) (m : ℕ), 1 ≤ m →
    ∃ (A : Finset ℕ) (c : Fin r), A.card = m ∧ (∀ x ∈ A, 0 < x) ∧
      (∀ x ∈ FS A ∪ FP A, χ x = c) ∧
      (FS A).card = 2 ^ m - 1 ∧ (FP A).card = 2 ^ m - 1 ∧
      FS A ∩ FP A = A

/-- 定理 1.1 は標準形の予想を含む（r = 0 と m = 0 の場合も込みで）。 -/
theorem conjecture_of_theorem (h : Theorem_1_1) : FiniteSumsProductsConjecture := by
  intro r χ m
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · exact (Fin.elim0 (χ 0))
  rcases Nat.eq_zero_or_pos m with rfl | hm
  · exact ⟨Fin.elim0, χ 0, fun i => i.elim0, fun i => i.elim0,
      fun J hJ => by
        obtain ⟨j, -⟩ := hJ
        exact j.elim0⟩
  obtain ⟨a, c, hmono, hpos, hcol, -⟩ := h r m hr hm 2 1 le_rfl le_rfl χ
  exact ⟨a, c, hmono, hpos, hcol⟩

/-- 命題が空でないことの確認: m = 1 の場合は自明に成り立つ
（a₁ を R より大きい正整数にとれば、和も積も a₁ 自身）。 -/
theorem theorem_1_1_m_one (r : ℕ) (hr : 1 ≤ r) (R D : ℝ) (χ : ℕ → Fin r) :
    ∃ (a : Fin 1 → ℕ) (c : Fin r),
      StrictMono a ∧ (∀ i, 0 < a i) ∧
      (∀ J : Finset (Fin 1), J.Nonempty →
        χ (∑ j ∈ J, a j) = c ∧ χ (∏ j ∈ J, a j) = c) ∧
      Separated R D a := by
  obtain ⟨n, hn⟩ := exists_nat_gt R
  refine ⟨fun _ => n + 1, χ (n + 1), fun i j hij => absurd hij (by omega),
    fun _ => Nat.succ_pos n, ?_, ?_⟩
  · intro J hJ
    have : J = univ := hJ.eq_univ
    subst this
    simp
  · refine ⟨fun _ _ => ?_, fun d hd => absurd hd (by omega)⟩
    push_cast
    linarith

/-! ## 第 2 段: 定理 1.1 ⇒ 系 1.2 -/

/-- 超増加条件: 各元 x は、x より小さい元の和と積の和の 2 倍より大きい。 -/
def SuperInc (A : Finset ℕ) : Prop :=
  ∀ x ∈ A, 2 * ((∑ y ∈ A.filter (· < x), y) + ∏ y ∈ A.filter (· < x), y) < x

lemma SuperInc.pos {A : Finset ℕ} (h : SuperInc A) {x : ℕ} (hx : x ∈ A) : 0 < x := by
  have := h x hx; omega

lemma SuperInc.three_le {A : Finset ℕ} (h : SuperInc A) {x : ℕ} (hx : x ∈ A) : 3 ≤ x := by
  have h1 := h x hx
  have h2 : 0 < ∏ y ∈ A.filter (· < x), y :=
    Finset.prod_pos (fun y hy => h.pos (Finset.mem_filter.1 hy).1)
  omega

lemma SuperInc.of_insert {a : ℕ} {s : Finset ℕ} (hlt : ∀ y ∈ s, y < a)
    (h : SuperInc (insert a s)) : SuperInc s ∧ 2 * ((∑ y ∈ s, y) + ∏ y ∈ s, y) < a := by
  constructor
  · intro x hx
    have hxa := hlt x hx
    have hfx : (insert a s).filter (· < x) = s.filter (· < x) := by
      rw [Finset.filter_insert, if_neg (by omega)]
    have := h x (Finset.mem_insert_of_mem hx)
    rwa [hfx] at this
  · have hfa : (insert a s).filter (· < a) = s := by
      rw [Finset.filter_insert, if_neg (lt_irrefl a), Finset.filter_true_of_mem hlt]
    have := h a (Finset.mem_insert_self a s)
    rwa [hfa] at this

lemma subset_of_insert_not_mem {a : ℕ} {s U : Finset ℕ} (hU : U ⊆ insert a s) (haU : a ∉ U) :
    U ⊆ s := fun y hy =>
  (Finset.mem_insert.1 (hU hy)).resolve_left (fun h => haU (h ▸ hy))

lemma erase_subset_of_insert {a : ℕ} {s U : Finset ℕ} (hU : U ⊆ insert a s) :
    U.erase a ⊆ s := fun y hy => by
  obtain ⟨hne, hyU⟩ := Finset.mem_erase.1 hy
  exact (Finset.mem_insert.1 (hU hyU)).resolve_left hne

lemma prod_le_of_subset {s U : Finset ℕ} (hs : SuperInc s) (hU : U ⊆ s) :
    ∏ y ∈ U, y ≤ ∏ y ∈ s, y :=
  Nat.le_of_dvd (Finset.prod_pos (fun y hy => hs.pos hy))
    (Finset.prod_dvd_prod_of_subset U s (fun y => y) hU)

/-- 超増加なら、部分集合の和は相異なる。 -/
lemma sum_inj (A : Finset ℕ) : SuperInc A →
    ∀ S ⊆ A, ∀ T ⊆ A, ∑ y ∈ S, y = ∑ y ∈ T, y → S = T := by
  refine Finset.induction_on_max A ?_ ?_
  · intro _ S hS T hT _
    rw [Finset.subset_empty.1 hS, Finset.subset_empty.1 hT]
  · intro a s hlt ih hA S hS T hT hST
    obtain ⟨hs, ha⟩ := SuperInc.of_insert hlt hA
    have hsum_lt : ∀ U ⊆ s, ∑ y ∈ U, y < a := fun U hU => by
      have := Finset.sum_le_sum_of_subset (f := fun y : ℕ => y) hU; omega
    by_cases haS : a ∈ S <;> by_cases haT : a ∈ T
    · have e1 : a + ∑ y ∈ S.erase a, y = ∑ y ∈ S, y := Finset.add_sum_erase S (fun y => y) haS
      have e2 : a + ∑ y ∈ T.erase a, y = ∑ y ∈ T, y := Finset.add_sum_erase T (fun y => y) haT
      have : S.erase a = T.erase a :=
        ih hs _ (erase_subset_of_insert hS) _ (erase_subset_of_insert hT) (by omega)
      calc S = insert a (S.erase a) := (Finset.insert_erase haS).symm
        _ = insert a (T.erase a) := by rw [this]
        _ = T := Finset.insert_erase haT
    · exfalso
      have e1 : a + ∑ y ∈ S.erase a, y = ∑ y ∈ S, y := Finset.add_sum_erase S (fun y => y) haS
      have := hsum_lt T (subset_of_insert_not_mem hT haT)
      omega
    · exfalso
      have e2 : a + ∑ y ∈ T.erase a, y = ∑ y ∈ T, y := Finset.add_sum_erase T (fun y => y) haT
      have := hsum_lt S (subset_of_insert_not_mem hS haS)
      omega
    · exact ih hs S (subset_of_insert_not_mem hS haS) T (subset_of_insert_not_mem hT haT) hST

/-- 超増加なら、部分集合の積は相異なる。 -/
lemma prod_inj (A : Finset ℕ) : SuperInc A →
    ∀ S ⊆ A, ∀ T ⊆ A, ∏ y ∈ S, y = ∏ y ∈ T, y → S = T := by
  refine Finset.induction_on_max A ?_ ?_
  · intro _ S hS T hT _
    rw [Finset.subset_empty.1 hS, Finset.subset_empty.1 hT]
  · intro a s hlt ih hA S hS T hT hST
    obtain ⟨hs, ha⟩ := SuperInc.of_insert hlt hA
    have hapos : 0 < a := hA.pos (Finset.mem_insert_self a s)
    have hlt_a : ∀ U ⊆ s, ∏ y ∈ U, y < a := fun U hU => by
      have := prod_le_of_subset hs hU; omega
    have hepos : ∀ U ⊆ insert a s, 0 < ∏ y ∈ U.erase a, y := fun U hU =>
      Finset.prod_pos (fun y hy => hs.pos (erase_subset_of_insert hU hy))
    by_cases haS : a ∈ S <;> by_cases haT : a ∈ T
    · have e1 : a * ∏ y ∈ S.erase a, y = ∏ y ∈ S, y := Finset.mul_prod_erase S (fun y => y) haS
      have e2 : a * ∏ y ∈ T.erase a, y = ∏ y ∈ T, y := Finset.mul_prod_erase T (fun y => y) haT
      have : S.erase a = T.erase a :=
        ih hs _ (erase_subset_of_insert hS) _ (erase_subset_of_insert hT)
          (Nat.eq_of_mul_eq_mul_left hapos (by rw [e1, e2, hST]))
      calc S = insert a (S.erase a) := (Finset.insert_erase haS).symm
        _ = insert a (T.erase a) := by rw [this]
        _ = T := Finset.insert_erase haT
    · exfalso
      have e1 : a * ∏ y ∈ S.erase a, y = ∏ y ∈ S, y := Finset.mul_prod_erase S (fun y => y) haS
      have h2 := hlt_a T (subset_of_insert_not_mem hT haT)
      have h3 : a * 1 ≤ a * ∏ y ∈ S.erase a, y := Nat.mul_le_mul_left a (hepos S hS)
      linarith
    · exfalso
      have e2 : a * ∏ y ∈ T.erase a, y = ∏ y ∈ T, y := Finset.mul_prod_erase T (fun y => y) haT
      have h2 := hlt_a S (subset_of_insert_not_mem hS haS)
      have h3 : a * 1 ≤ a * ∏ y ∈ T.erase a, y := Nat.mul_le_mul_left a (hepos T hT)
      linarith
    · exact ih hs S (subset_of_insert_not_mem hS haS) T (subset_of_insert_not_mem hT haT) hST

/-- 超増加なら、空でない部分集合の和と積が等しいのは、同じ 1 元集合のときだけ。 -/
lemma sum_eq_prod (A : Finset ℕ) : SuperInc A →
    ∀ S ⊆ A, S.Nonempty → ∀ T ⊆ A, T.Nonempty → ∑ y ∈ S, y = ∏ y ∈ T, y →
      ∃ x ∈ A, S = {x} ∧ T = {x} := by
  refine Finset.induction_on_max A ?_ ?_
  · intro _ S hS hSne
    rw [Finset.subset_empty.1 hS] at hSne
    exact absurd hSne Finset.not_nonempty_empty
  · intro a s hlt ih hA S hS hSne T hT hTne hST
    obtain ⟨hs, ha⟩ := SuperInc.of_insert hlt hA
    have hsum_lt : ∀ U ⊆ s, ∑ y ∈ U, y < a := fun U hU => by
      have := Finset.sum_le_sum_of_subset (f := fun y : ℕ => y) hU; omega
    have hprod_lt : ∀ U ⊆ s, ∏ y ∈ U, y < a := fun U hU => by
      have := prod_le_of_subset hs hU; omega
    have hepos : ∀ U ⊆ insert a s, 0 < ∏ y ∈ U.erase a, y := fun U hU =>
      Finset.prod_pos (fun y hy => hs.pos (erase_subset_of_insert hU hy))
    by_cases haS : a ∈ S <;> by_cases haT : a ∈ T
    · have e1 : a + ∑ y ∈ S.erase a, y = ∑ y ∈ S, y := Finset.add_sum_erase S (fun y => y) haS
      have e2 : a * ∏ y ∈ T.erase a, y = ∏ y ∈ T, y := Finset.mul_prod_erase T (fun y => y) haT
      have hS1 := hsum_lt _ (erase_subset_of_insert hS)
      rcases (T.erase a).eq_empty_or_nonempty with hTe | ⟨y, hy⟩
      · rw [hTe, Finset.prod_empty, mul_one] at e2
        have h0 : ∑ y ∈ S.erase a, y = 0 := by omega
        rw [Finset.sum_eq_zero_iff] at h0
        have hSe : S.erase a = ∅ := by
          rcases (S.erase a).eq_empty_or_nonempty with h | ⟨z, hz⟩
          · exact h
          · exfalso
            have := hs.pos (erase_subset_of_insert hS hz)
            have := h0 z hz
            omega
        refine ⟨a, Finset.mem_insert_self a s, ?_, ?_⟩
        · rw [← Finset.insert_erase haS, hSe]; simp
        · rw [← Finset.insert_erase haT, hTe]; simp
      · exfalso
        have hy3 : 3 ≤ y := hs.three_le (erase_subset_of_insert hT hy)
        have hP : y ≤ ∏ z ∈ T.erase a, z :=
          Nat.le_of_dvd (hepos T hT) (Finset.dvd_prod_of_mem (fun z => z) hy)
        have h3 : a * 3 ≤ a * ∏ z ∈ T.erase a, z := Nat.mul_le_mul_left a (le_trans hy3 hP)
        linarith
    · exfalso
      have e1 : a + ∑ y ∈ S.erase a, y = ∑ y ∈ S, y := Finset.add_sum_erase S (fun y => y) haS
      have := hprod_lt T (subset_of_insert_not_mem hT haT)
      omega
    · exfalso
      have e2 : a * ∏ y ∈ T.erase a, y = ∏ y ∈ T, y := Finset.mul_prod_erase T (fun y => y) haT
      have h2 := hsum_lt S (subset_of_insert_not_mem hS haS)
      have h3 : a * 1 ≤ a * ∏ y ∈ T.erase a, y := Nat.mul_le_mul_left a (hepos T hT)
      linarith
    · obtain ⟨x, hx, h1, h2⟩ := ih hs S (subset_of_insert_not_mem hS haS) hSne
        T (subset_of_insert_not_mem hT haT) hTne hST
      exact ⟨x, Finset.mem_insert_of_mem hx, h1, h2⟩

lemma card_nonempty_subsets (A : Finset ℕ) :
    (A.powerset.filter (·.Nonempty)).card = 2 ^ A.card - 1 := by
  have : A.powerset.filter (·.Nonempty) = A.powerset.erase ∅ := by
    ext S; simp [Finset.nonempty_iff_ne_empty, and_comm]
  rw [this, Finset.card_erase_of_mem (Finset.empty_mem_powerset A), Finset.card_powerset]

/-- 第 2 段: 定理 1.1 から系 1.2 が従う（R = 2, D = 1 の分離で十分）。 -/
theorem corollary_of_theorem (h : Theorem_1_1) : Corollary_1_2 := by
  intro r χ m hm
  unfold FS FP
  rcases Nat.eq_zero_or_pos r with rfl | hr
  · exact (χ 0).elim0
  obtain ⟨a, c, hmono, hpos, hcol, hsep⟩ := h r m hr hm 2 1 le_rfl le_rfl χ
  have hinj : Function.Injective a := hmono.injective
  set A := Finset.univ.image a with hAdef
  have hcard : A.card = m := by
    rw [hAdef, Finset.card_image_of_injective _ hinj]; simp
  have hfilt : ∀ d : Fin m, A.filter (· < a d) = (Finset.univ.filter (· < d)).image a := by
    intro d
    rw [hAdef, Finset.filter_image]
    congr 1
    ext k; simp [hmono.lt_iff_lt]
  have hSI : SuperInc A := by
    intro x hx
    obtain ⟨d, -, rfl⟩ := Finset.mem_image.1 hx
    rw [hfilt d, Finset.sum_image (fun i _ j _ e => hinj e),
      Finset.prod_image (fun i _ j _ e => hinj e)]
    rcases Nat.eq_zero_or_pos d.val with hd | hd
    · have hempty : Finset.univ.filter (· < d) = ∅ :=
        Finset.filter_false_of_mem (fun k _ hk => by have := Fin.lt_def.1 hk; omega)
      rw [hempty, Finset.sum_empty, Finset.prod_empty]
      have h2 : (2 : ℝ) < a d := hsep.1 d hd
      have : 2 < a d := by exact_mod_cast h2
      omega
    · have h2 := hsep.2 d hd
      rw [Real.rpow_one] at h2
      exact_mod_cast h2
  refine ⟨A, c, hcard, fun x hx => hSI.pos hx, ?_, ?_, ?_, ?_⟩
  · intro x hx
    rcases Finset.mem_union.1 hx with hx | hx
    · obtain ⟨S, hS, rfl⟩ := Finset.mem_image.1 hx
      obtain ⟨hSA, hSne⟩ := Finset.mem_filter.1 hS
      obtain ⟨J, -, rfl⟩ := Finset.subset_image_iff.1 (Finset.mem_powerset.1 hSA)
      show χ (∑ y ∈ J.image a, y) = c
      rw [Finset.sum_image (fun i _ j _ e => hinj e)]
      exact (hcol J (by simpa using hSne)).1
    · obtain ⟨S, hS, rfl⟩ := Finset.mem_image.1 hx
      obtain ⟨hSA, hSne⟩ := Finset.mem_filter.1 hS
      obtain ⟨J, -, rfl⟩ := Finset.subset_image_iff.1 (Finset.mem_powerset.1 hSA)
      show χ (∏ y ∈ J.image a, y) = c
      rw [Finset.prod_image (fun i _ j _ e => hinj e)]
      exact (hcol J (by simpa using hSne)).2
  · rw [Finset.card_image_of_injOn, card_nonempty_subsets, hcard]
    intro S hS T hT e
    exact sum_inj A hSI S (Finset.mem_powerset.1 (Finset.mem_filter.1 hS).1)
      T (Finset.mem_powerset.1 (Finset.mem_filter.1 hT).1) e
  · rw [Finset.card_image_of_injOn, card_nonempty_subsets, hcard]
    intro S hS T hT e
    exact prod_inj A hSI S (Finset.mem_powerset.1 (Finset.mem_filter.1 hS).1)
      T (Finset.mem_powerset.1 (Finset.mem_filter.1 hT).1) e
  · ext x
    simp only [Finset.mem_inter]
    constructor
    · rintro ⟨hx1, hx2⟩
      obtain ⟨S, hS, rfl⟩ := Finset.mem_image.1 hx1
      obtain ⟨T, hT, hTS⟩ := Finset.mem_image.1 hx2
      obtain ⟨hSA, hSne⟩ := Finset.mem_filter.1 hS
      obtain ⟨hTA, hTne⟩ := Finset.mem_filter.1 hT
      obtain ⟨y, hy, rfl, -⟩ := sum_eq_prod A hSI S (Finset.mem_powerset.1 hSA) hSne
        T (Finset.mem_powerset.1 hTA) hTne hTS.symm
      simpa using hy
    · intro hx
      have hmem : ({x} : Finset ℕ) ∈ A.powerset.filter (·.Nonempty) := by simp [hx]
      exact ⟨Finset.mem_image.2 ⟨{x}, hmem, by simp⟩, Finset.mem_image.2 ⟨{x}, hmem, by simp⟩⟩

/-- OpenAI の主張（第 1 段では未証明）。 -/
theorem theorem_1_1 : Theorem_1_1 := by
  sorry

end Koeru.Hindman164

#print axioms Koeru.Hindman164.corollary_of_theorem
#print axioms Koeru.Hindman164.conjecture_of_theorem
#print axioms Koeru.Hindman164.theorem_1_1_m_one
