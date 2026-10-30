import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Fintype.BigOperators

/-!
# 有限 Ramsey の定理（超グラフ版、全サイズ同時）

OpenAI 結果 164 の補題 2.5 が使う「繰り返し有限 Ramsey」:
各サイズ k ≤ q の部分集合の色がサイズだけで決まる、大きさ t の部分集合が取れる。

証明は独自に書いた（OpenAI の Lean コードは使っていない）。
色ごとの目標 T を 1 つずつ減らす古典的な帰納法（Erdős–Szekeres 型）で、順序は使わない。
-/

namespace Koeru.Ramsey

open Finset

/-- N 以上の大きさの集合 S と、k 元部分集合の r 彩色に対し、ある色 i と
大きさ T i の部分集合 H ⊆ S で、H の k 元部分集合がすべて色 i になるものがある。 -/
def RBound (k r : ℕ) (T : Fin r → ℕ) (N : ℕ) : Prop :=
  ∀ S : Finset ℕ, N ≤ S.card → ∀ χ : Finset ℕ → Fin r,
    ∃ i H, H ⊆ S ∧ H.card = T i ∧ ∀ A ⊆ H, A.card = k → χ A = i

/-- 有限 Ramsey の定理（色ごとに目標が異なる形）。 -/
theorem ramsey (k r : ℕ) : ∀ T : Fin r → ℕ, 0 < r → ∃ N, RBound k r T N := by
  induction k with
  | zero =>
    intro T hr
    refine ⟨univ.sup T, fun S hS χ => ?_⟩
    obtain ⟨H, hHS, hHc⟩ :=
      Finset.exists_subset_card_eq (le_trans (Finset.le_sup (f := T) (mem_univ (χ ∅))) hS)
    refine ⟨χ ∅, H, hHS, hHc, fun A _ hA => ?_⟩
    rw [Finset.card_eq_zero.1 hA]
  | succ k ih =>
    intro T0 hr
    suffices key : ∀ s, ∀ T : Fin r → ℕ, ∑ i, T i = s → ∃ N, RBound (k + 1) r T N from
      key _ T0 rfl
    intro s
    refine Nat.strong_induction_on s ?_
    intro s ihs T hT
    by_cases hz : ∃ i, T i = 0
    · obtain ⟨i, hi⟩ := hz
      refine ⟨0, fun S _ χ => ⟨i, ∅, empty_subset _, by simp [hi], fun A hA hAc => ?_⟩⟩
      rw [Finset.subset_empty.1 hA] at hAc
      simp at hAc
    · push Not at hz
      have hdec : ∀ j, ∃ N, RBound (k + 1) r (Function.update T j (T j - 1)) N := by
        intro j
        refine ihs _ ?_ _ rfl
        rw [← hT]
        apply Finset.sum_lt_sum
        · intro i _
          by_cases hij : i = j
          · rw [hij, Function.update_self]; omega
          · rw [Function.update_of_ne hij]
        · exact ⟨j, mem_univ _, by rw [Function.update_self]; have := hz j; omega⟩
      choose M hM using hdec
      obtain ⟨N', hN'⟩ := ih M hr
      refine ⟨N' + 1, fun S hS χ => ?_⟩
      obtain ⟨x, hx⟩ : S.Nonempty := Finset.card_pos.1 (by omega)
      have hS' : N' ≤ (S.erase x).card := by rw [Finset.card_erase_of_mem hx]; omega
      obtain ⟨j, H, hHS, hHc, hHcol⟩ := hN' (S.erase x) hS' (fun A => χ (insert x A))
      obtain ⟨i, G, hGH, hGc, hGcol⟩ := hM j H (le_of_eq hHc.symm) χ
      have hGS : G ⊆ S := fun y hy => Finset.mem_of_mem_erase (hHS (hGH hy))
      by_cases hij : i = j
      · subst hij
        have hxG : x ∉ G := fun h => by
          have := hHS (hGH h)
          simp at this
        refine ⟨i, insert x G, ?_, ?_, ?_⟩
        · intro y hy
          rcases Finset.mem_insert.1 hy with rfl | hy
          · exact hx
          · exact hGS hy
        · rw [Finset.card_insert_of_notMem hxG, hGc, Function.update_self]
          have := hz i
          omega
        · intro A hA hAc
          by_cases hxA : x ∈ A
          · have hAe : A.erase x ⊆ H := by
              intro y hy
              obtain ⟨hne, hyA⟩ := Finset.mem_erase.1 hy
              rcases Finset.mem_insert.1 (hA hyA) with h | h
              · exact absurd h hne
              · exact hGH h
            have := hHcol (A.erase x) hAe (by rw [Finset.card_erase_of_mem hxA, hAc]; omega)
            rwa [Finset.insert_erase hxA] at this
          · apply hGcol A _ hAc
            intro y hy
            rcases Finset.mem_insert.1 (hA hy) with h | h
            · exact absurd (h ▸ hy) hxA
            · exact h
      · refine ⟨i, G, hGS, ?_, hGcol⟩
        rw [hGc, Function.update_of_ne hij]

/-- H の中では、大きさ q 以下の部分集合の色が大きさだけで決まる。 -/
def Homog {β : Type*} (χ : Finset ℕ → β) (H : Finset ℕ) (q : ℕ) : Prop :=
  ∀ A ⊆ H, ∀ D ⊆ H, A.card = D.card → A.card ≤ q → χ A = χ D

/-- 全サイズ同時の有限 Ramsey の定理。 -/
theorem ramsey_simultaneous (r : ℕ) (hr : 0 < r) (q : ℕ) :
    ∀ t, ∃ N, t ≤ N ∧ ∀ S : Finset ℕ, N ≤ S.card → ∀ χ : Finset ℕ → Fin r,
      ∃ H ⊆ S, H.card = t ∧ Homog χ H q := by
  induction q with
  | zero =>
    intro t
    refine ⟨t, le_rfl, fun S hS χ => ?_⟩
    obtain ⟨H, hHS, hHc⟩ := Finset.exists_subset_card_eq hS
    refine ⟨H, hHS, hHc, fun A _ D _ hAD hA => ?_⟩
    have hA0 : A = ∅ := Finset.card_eq_zero.1 (by omega)
    have hD0 : D = ∅ := Finset.card_eq_zero.1 (by omega)
    rw [hA0, hD0]
  | succ q ih =>
    intro t
    obtain ⟨M, hM⟩ := ramsey (q + 1) r (fun _ => t) hr
    obtain ⟨N, hMN, hN⟩ := ih (max M t)
    refine ⟨N, le_trans (le_max_right _ _) hMN, fun S hS χ => ?_⟩
    obtain ⟨H, hHS, hHc, hHhom⟩ := hN S hS χ
    obtain ⟨i, G, hGH, hGc, hGcol⟩ := hM H (by rw [hHc]; exact le_max_left _ _) χ
    refine ⟨G, hGH.trans hHS, hGc, fun A hA D hD hAD hAq => ?_⟩
    rcases Nat.lt_or_ge A.card (q + 1) with h | h
    · exact hHhom A (hA.trans hGH) D (hD.trans hGH) hAD (by omega)
    · rw [hGcol A hA (by omega), hGcol D hD (by omega)]

end Koeru.Ramsey

#print axioms Koeru.Ramsey.ramsey_simultaneous
