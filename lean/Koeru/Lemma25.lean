import Koeru.FiniteSums
import Koeru.Ramsey
import Mathlib.Data.Finset.Sort
import Mathlib.Algebra.BigOperators.Ring.Finset

/-!
# OpenAI 結果 164 の補題 2.5（積の鎖の選択）

論文 §2.4 Lemma 2.5: 任意の m, r に対して n = n(m, r) が存在し、任意の彩色
χ : ℕ → [r] と任意の列 x₁, …, x_n に対して、長さ m の鎖 B₁, …, B_m で
x_{B₁}, …, x_{B_m} の空でない積がすべて同色になるものがある。

鎖 (2.2): 空でない尾 T_d とピボット i_d が T₁ < T₂ < ⋯ < T_m < i₁ < ⋯ < i_m を満たし、
B_d = T_d ∪ {i_d}。

証明は論文どおり、繰り返し有限 Ramsey（`Koeru.Ramsey.ramsey_simultaneous`）と
有限和の定理（`Koeru.finite_sums_finite`）による。位置の配置だけ論文と変え、
各ブロックを幅 W = 2F + 1 の区間に置く（n は大きくなるが、主張は n の存在だけ）。
-/

namespace Koeru.Lemma25

open Finset Koeru Koeru.Ramsey

/-- 鎖の条件 (2.2)。 -/
structure IsChain {m : ℕ} (T : Fin m → Finset ℕ) (piv : Fin m → ℕ) : Prop where
  ne : ∀ d, (T d).Nonempty
  tails : ∀ d e, d < e → ∀ a ∈ T d, ∀ b ∈ T e, a < b
  tail_piv : ∀ d e, ∀ a ∈ T d, a < piv e
  piv_mono : StrictMono piv

/-- ブロック B_d = T_d ∪ {i_d}。 -/
def block {m : ℕ} (T : Fin m → Finset ℕ) (piv : Fin m → ℕ) (d : Fin m) : Finset ℕ :=
  insert (piv d) (T d)

/-- 補題 2.5（積の鎖の選択）。 -/
theorem lemma_2_5 (m r : ℕ) (hr : 0 < r) : ∃ n, ∀ (χ : ℕ → Fin r) (x : ℕ → ℕ),
    ∃ (T : Fin m → Finset ℕ) (piv : Fin m → ℕ), IsChain T piv ∧
      (∀ d, block T piv d ⊆ range n) ∧
      ∃ c : Fin r, ∀ J : Finset (Fin m), J.Nonempty →
        χ (∏ d ∈ J, ∏ j ∈ block T piv d, x j) = c := by
  classical
  obtain ⟨F, hF⟩ := finite_sums_finite m r
  set W := 2 * F + 1 with hW
  set n0 := m * W + m with hn0
  obtain ⟨N, -, hN⟩ := ramsey_simultaneous r hr (2 * F) n0
  refine ⟨N, fun χ x => ?_⟩
  -- 部分集合の彩色 A ↦ χ(x_A) と、Ramsey による均質な集合 H
  let cχ : Finset ℕ → Fin r := fun A => χ (∏ j ∈ A, x j)
  obtain ⟨H, hHS, hHc, hHom⟩ := hN (range N) (by simp) cχ
  -- H の増加順の列挙 e'
  let en : Fin n0 ↪o ℕ := H.orderEmbOfFin hHc
  let e' : ℕ → ℕ := fun p => if h : p < n0 then en ⟨p, h⟩ else 0
  have he'mem : ∀ p, p < n0 → e' p ∈ H := fun p hp => by
    simp only [e', dite_eq_left hp]
    exact H.orderEmbOfFin_mem hHc _
  have he'lt : ∀ p p', p < p' → p' < n0 → e' p < e' p' := fun p p' h h' => by
    simp only [e', dite_eq_left h', dite_eq_left (lt_trans h h')]
    exact en.lt_iff_lt.2 (Fin.mk_lt_mk.2 h)
  have he'inj : ∀ p p', p < n0 → p' < n0 → e' p = e' p' → p = p' := by
    intro p p' hp hp' heq
    rcases lt_trichotomy p p' with h | h | h
    · exact absurd heq (ne_of_lt (he'lt p p' h hp'))
    · exact h
    · exact absurd heq.symm (ne_of_lt (he'lt p' p h hp))
  -- 各大きさの代表部分集合
  let pre : ℕ → Finset ℕ := fun k =>
    if hk : k ≤ H.card then Classical.choose (Finset.exists_subset_card_eq hk) else ∅
  have hpre : ∀ k, k ≤ H.card → pre k ⊆ H ∧ (pre k).card = k := fun k hk => by
    simp only [pre, dite_eq_left hk]
    exact Classical.choose_spec (Finset.exists_subset_card_eq hk)
  -- 有限和の定理を「大きさ 2a の色」に適用
  let col : ℕ → Fin r := fun a => cχ (pre (2 * a))
  obtain ⟨u, c, hupos, hu⟩ := hF col
  have hule : ∀ d, u d ≤ F := fun d => by simpa using (hu {d} (singleton_nonempty d)).1
  -- 位置の配置: 尾は [dW, dW + 2u_d − 1)、ピボットは mW + d
  let tpos : Fin m → Finset ℕ := fun d => Ico ((d : ℕ) * W) ((d : ℕ) * W + (2 * u d - 1))
  let pos : Fin m → Finset ℕ := fun d => insert (m * W + d) (tpos d)
  have hdW : ∀ d : Fin m, (d : ℕ) * W + W ≤ m * W := fun d => by
    have := Nat.mul_le_mul_right W (show (d : ℕ) + 1 ≤ m from d.isLt)
    rw [Nat.add_mul, Nat.one_mul] at this
    exact this
  have hdeW : ∀ d e : Fin m, d < e → (d : ℕ) * W + W ≤ (e : ℕ) * W := fun d e h => by
    have := Nat.mul_le_mul_right W (show (d : ℕ) + 1 ≤ e from h)
    rw [Nat.add_mul, Nat.one_mul] at this
    exact this
  have htpos : ∀ d p, p ∈ tpos d → (d : ℕ) * W ≤ p ∧ p < (d : ℕ) * W + W := by
    intro d p hp
    obtain ⟨h1, h2⟩ := Finset.mem_Ico.1 hp
    have := hule d
    omega
  have hpos_lt : ∀ d p, p ∈ pos d → p < n0 := by
    intro d p hp
    rcases Finset.mem_insert.1 hp with rfl | hp
    · have := d.isLt; omega
    · have := htpos d p hp; have := hdW d; omega
  -- 鎖
  let T : Fin m → Finset ℕ := fun d => (tpos d).image e'
  let piv : Fin m → ℕ := fun d => e' (m * W + d)
  have hblock : ∀ d, block T piv d = (pos d).image e' := fun d => by
    simp only [block, T, piv, pos, Finset.image_insert]
  have hchain : IsChain T piv := by
    refine ⟨fun d => ?_, fun d e h a ha b hb => ?_, fun d e a ha => ?_, fun d e h => ?_⟩
    · refine Finset.Nonempty.image (Finset.nonempty_Ico.2 ?_) e'
      have := hupos d
      omega
    · obtain ⟨p, hp, rfl⟩ := Finset.mem_image.1 ha
      obtain ⟨p', hp', rfl⟩ := Finset.mem_image.1 hb
      have h1 := htpos d p hp
      have h2 := htpos e p' hp'
      have h3 := hdeW d e h
      have h4 := hpos_lt e p' (Finset.mem_insert_of_mem hp')
      exact he'lt p p' (by omega) h4
    · obtain ⟨p, hp, rfl⟩ := Finset.mem_image.1 ha
      have h1 := htpos d p hp
      have h2 := hdW d
      have h3 := e.isLt
      exact he'lt p (m * W + e) (by omega) (by omega)
    · have h1 : (d : ℕ) < e := h
      have h2 := e.isLt
      exact he'lt (m * W + d) (m * W + e) (by omega) (by omega)
  -- ブロックの位置は互いに素
  have hdisj : ∀ d e : Fin m, d ≠ e → Disjoint (pos d) (pos e) := by
    intro d e hde
    rw [Finset.disjoint_left]
    intro p hpd hpe
    have hval : (d : ℕ) ≠ e := fun h => hde (Fin.ext h)
    rcases Finset.mem_insert.1 hpd with h1 | h1 <;> rcases Finset.mem_insert.1 hpe with h2 | h2
    · omega
    · have := htpos e p h2; have := hdW e; omega
    · have := htpos d p h1; have := hdW d; omega
    · have a1 := htpos d p h1
      have a2 := htpos e p h2
      rcases lt_or_gt_of_ne hde with h | h
      · have := hdeW d e h; omega
      · have := hdeW e d h; omega
  have hcard_pos : ∀ d, (pos d).card = 2 * u d := fun d => by
    have hnot : m * W + (d : ℕ) ∉ tpos d := fun h => by
      have := htpos d _ h; have := hdW d; omega
    rw [Finset.card_insert_of_notMem hnot, Nat.card_Ico]
    have := hupos d
    omega
  have hinjOn : ∀ S : Finset ℕ, (∀ p ∈ S, p < n0) → Set.InjOn e' S :=
    fun S hS p hp p' hp' heq => he'inj p p' (hS p hp) (hS p' hp') heq
  refine ⟨T, piv, hchain, fun d => ?_, c, fun J hJ => ?_⟩
  · rw [hblock d]
    intro y hy
    obtain ⟨p, hp, rfl⟩ := Finset.mem_image.1 hy
    exact hHS (he'mem p (hpos_lt d p hp))
  · -- 積を一つの集合 U 上の積に書き直す
    have hU_lt : ∀ p ∈ J.biUnion pos, p < n0 := by
      intro p hp
      obtain ⟨d, -, hp⟩ := Finset.mem_biUnion.1 hp
      exact hpos_lt d p hp
    have hprod : ∏ d ∈ J, ∏ j ∈ block T piv d, x j =
        ∏ j ∈ (J.biUnion pos).image e', x j := by
      rw [Finset.prod_image (hinjOn _ hU_lt), Finset.prod_biUnion (by
        intro d _ e _ h
        exact hdisj d e h)]
      refine Finset.prod_congr rfl (fun d _ => ?_)
      rw [hblock d, Finset.prod_image (hinjOn _ (hpos_lt d))]
    have hcardU : ((J.biUnion pos).image e').card = 2 * ∑ d ∈ J, u d := by
      rw [Finset.card_image_of_injOn (hinjOn _ hU_lt), Finset.card_biUnion (by
        intro d _ e _ h
        exact hdisj d e h)]
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl (fun d _ => hcard_pos d)
    have hUH : (J.biUnion pos).image e' ⊆ H := by
      intro y hy
      obtain ⟨p, hp, rfl⟩ := Finset.mem_image.1 hy
      exact he'mem p (hU_lt p hp)
    obtain ⟨hsumF, hcol⟩ := hu J hJ
    have hm : 1 ≤ m := by
      obtain ⟨d, -⟩ := hJ
      have := d.isLt
      omega
    have hWm : W ≤ m * W := Nat.le_mul_of_pos_left W hm
    have hsize : 2 * ∑ d ∈ J, u d ≤ H.card := by rw [hHc]; omega
    obtain ⟨hpreH, hprec⟩ := hpre _ hsize
    rw [hprod]
    show cχ ((J.biUnion pos).image e') = c
    rw [hHom _ hUH _ hpreH (by rw [hcardU, hprec]) (by rw [hcardU]; omega)]
    exact hcol

end Koeru.Lemma25

#print axioms Koeru.Lemma25.lemma_2_5
