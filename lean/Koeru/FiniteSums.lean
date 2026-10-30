import Mathlib.Combinatorics.Hindman
import Mathlib.Combinatorics.Compactness
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# 有限和の定理（有限形）

OpenAI 結果 164 の論文 §2.4 が使う「有限和の定理の有限形」（Folkman–Rado–Sanders）:
任意の m, r に対して F が存在し、ℕ の任意の r 彩色は、正整数 u₁, …, u_m で
空でない部分和がすべて F 以下・同色になるものを含む。

論文の説明どおり、Mathlib の Hindman の定理（無限形）とコンパクト性
（`Finset.rado_selection`）から導く。
-/

namespace Koeru

open Finset Hindman

/-- F 以下に収まる単色な有限和の配置が、彩色 χ の中にあること。 -/
def FSGood (m r F : ℕ) (χ : ℕ → Fin r) : Prop :=
  ∃ (u : Fin m → ℕ) (c : Fin r), (∀ d, 0 < u d) ∧
    ∀ J : Finset (Fin m), J.Nonempty → ∑ j ∈ J, u j ≤ F ∧ χ (∑ j ∈ J, u j) = c

/-- 正の項からなる列の有限和は正。 -/
lemma fs_pos : ∀ (a : Stream' ℕ) (m : ℕ), FS a m → (∀ n, 0 < a.get n) → 0 < m := by
  intro a m h
  induction h with
  | head a => intro hp; exact hp 0
  | of_tail a m _ ih => intro hp; exact ih (fun n => hp (n + 1))
  | cons a m _ ih =>
    intro hp
    have h0 : 0 < a.head := hp 0
    have h1 := ih (fun n => hp (n + 1))
    omega

/-- Hindman の定理（ℕ、正の項）: 有限和がすべて正で同色になる列がある。 -/
lemma infinite_fs (r : ℕ) (χ : ℕ → Fin r) :
    ∃ (b : Stream' ℕ) (c : Fin r), ∀ m, FS b m → 0 < m ∧ χ m = c := by
  let a : Stream' ℕ := fun n => n + 1
  obtain ⟨cl, hcl, b, hb⟩ :=
    FS_partition_regular a (Set.range fun c : Fin r => {n | 0 < n ∧ χ n = c})
      (Set.finite_range _) (by
        intro m hm
        have := fs_pos a m hm (fun n => Nat.succ_pos n)
        exact Set.mem_sUnion.2 ⟨_, ⟨χ m, rfl⟩, this, rfl⟩)
  obtain ⟨c, rfl⟩ := hcl
  exact ⟨b, c, fun m hm => hb m hm⟩

/-- 有限和の定理（有限形）。 -/
theorem finite_sums_finite (m r : ℕ) : ∃ F, ∀ χ : ℕ → Fin r, FSGood m r F χ := by
  by_contra hne
  push Not at hne
  choose g hg using hne
  obtain ⟨χ, hχ⟩ := Finset.rado_selection (β := fun _ : ℕ => Fin r) (fun s => g (s.sup id))
  obtain ⟨b, c, hb⟩ := infinite_fs r χ
  let u : Fin m → ℕ := fun d => b.get d
  have hsum : ∀ J : Finset (Fin m), J.Nonempty → FS b (∑ j ∈ J, u j) := by
    intro J hJ
    have := FS.finsetSum b (J.map Fin.valEmbedding) (by simpa using hJ)
    simpa [u, Finset.sum_map] using this
  let F := ∑ j : Fin m, u j
  have hle : ∀ J : Finset (Fin m), ∑ j ∈ J, u j ≤ F := fun J =>
    Finset.sum_le_sum_of_subset (Finset.subset_univ J)
  obtain ⟨t, hst, hagree⟩ := hχ (Finset.range (F + 1))
  have hFt : F ≤ t.sup id :=
    Finset.le_sup (f := id) (hst (Finset.mem_range.2 (Nat.lt_succ_self F)))
  apply hg (t.sup id)
  refine ⟨u, c, fun d => (hb _ (FS.singleton b d)).1, fun J hJ => ⟨le_trans (hle J) hFt, ?_⟩⟩
  have hmem : ∑ j ∈ J, u j ∈ Finset.range (F + 1) :=
    Finset.mem_range.2 (Nat.lt_succ_of_le (hle J))
  have h1 : χ (∑ j ∈ J, u j) = g (t.sup id) (∑ j ∈ J, u j) := hagree _ hmem
  rw [← h1]
  exact (hb _ (hsum J hJ)).2

end Koeru

#print axioms Koeru.finite_sums_finite
