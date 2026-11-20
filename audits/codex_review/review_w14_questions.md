公開前に直すべき点はあります。**確実な誤記は077の定理番号です。004bのKatoへの質問も、有理版と整係数版を混同した形になっています。** 原論文の主定理そのものを反証する誤りは、今回特定していません。

以下、`Q` は [public_repo/QUESTIONS.md](./public_repo/QUESTIONS.md)、`011/…` などは当該論文の `build/` 以下を指します。ファイル変更・収録コードの実行はしていません。

**(1) 全8項目の照合：判定＝一部正しい**

| ページの項目 | 判定 | 根拠と修正点 |
|---|---|---|
| 1：011 | 一部正しい | 定理3.1とProp 4.1の番号は正しい。ただし帰結にも **marked weight、特殊なα、roughなβ、固定剰余類1** が残る。Q:10の「distribution estimate」だけでは広く読める。`011/sections/03-determinant.tex:37–58`、`04-memory.tex:11–35`。差し替えは(2)。 |
| 2：029 | 正しい | 定理1.2、零点なし半平面、ζへの含意、142の引用は一致。011の輸入も本文で確認できる。`029/sections/01-introduction.tex:38–53`、`08-type-ii.tex:61–85`、`142/sections/10-analytic.tex:4–18`。Q:17には軽微な限定を補う。 |
| 3：142 | 正しい | 十分大きい全素数についての帰結は正しい。`142/sections/00-introduction.tex:46–55`、`10-analytic.tex:91–111`。詳細は(3)。 |
| 4：004b | 一部正しい | 枝の説明は基本的に正しい。ただしKatoへの問いが不適切で、CMの場合も別入力。`004b/sections/cyclotomic.tex:423–438`。詳細は(4)。 |
| 5：077 | 誤り | **Theorem 2.1ではなく2.2。** また厳密な帰結は部分測度と微小な追加損失を伴う。`077/preamble.tex:18–25`、`sections/01-geometric.tex:32–42, 107–119`。詳細は(5)。 |
| 6：109 | 正しい | 乗算全体のoblivious性を保留した修正版は妥当。局所補題からの経路には、補題の正しさと固定kテープへの拡張が必要。`109/sections/07-resampling.tex:215–238`、`04-swap.tex:285–301`、`audits/codex_review/review_ncc109.md:27–37`。 |
| 7：285 | 正しい | 非整数の**非正規化行列トレース**からbase groupについても従う。`285/sections/setup.tex:6–17`、`introduction.tex:40–52`。詳細は(6)。 |
| 8：321 | 一部正しい | 定理2.2・系2.3、条件付きのdeficiencyとrelation gapは正しい。ただし後続論文は**別の群**について同型の帰結を記す。Q:47はその区別を省略。`321/sections/02-tower.tex:103–109, 230–235`、`audits/321/REPORT.md:49–57,130`。 |

報告書・README以上に読める箇所は、主にQ:10の帰結の条件省略、Q:26のCMを含む「whole branch」、Q:34の損失の省略、Q:39の局所補題を無条件に成立済みと読む余地、Q:47の「same consequence」の対象です。ほかに具体的な過大主張は特定できませんでした。

直すべき文面：

Q:17末尾：

> The margins are 1/6000 on the reflection side and about 0.00086 before small power losses for the relevant non-principal rows. An additional loss of 10⁻³ at either of these estimates would invalidate the stated exponent bookkeeping.

これは `029/sections/06-zero-free.tex:650–683`、`audits/codex_review/review_w7_029core.md:29–48` に対応します。「証明を修復できない」という意味にはしません。

Q:39の局所補題以降：

> Separately, if the chunk-swap lemma (`04-swap.tex:285–301`) and its shape-determined implementation are correct, setting P = G = B = 1 transposes an m×m bit matrix in O(m²(log m)^τ) time, where τ < 1. Our reviewer derived an extension of Adler et al.’s two-tape simulation to any fixed number of tapes; together with their Theorem 6.2, this would refute the same conjecture. Is that extension valid?

既存論文の定理をそのまま引用しているのではない、という現在の区別は維持してください。[Adler et al., §6](https://www.cs.cornell.edu/~rdk/papers/AHJKRL05.pdf)、[AFKL, Theorem 1](https://drops.dagstuhl.de/storage/00lipics/lipics-vol132-icalp2019/LIPIcs.ICALP.2019.10/LIPIcs.ICALP.2019.10.pdf)との対応も妥当です。

Q:47：

> Corollary 2.3 (`02-tower.tex:230–235`) states that, for a finite connected 2-complex Y, H₂(Y; C_ρ) = 0 forces ρ to be trivial on every torsion element. Is the tower argument proving this statement valid? If it is, the constructed group has deficiency 0 and its presentation has a relation gap. Bartholdi–Bridson–Mikhailov–Nicholson derive an analogous deficiency conclusion for different groups, conditionally on this corollary.

後続論文の条件付けは明示されています。[同論文、Corollary B直前と証明](https://arxiv.org/html/2610.23637v1)。

**(2) 011のEH型帰結：判定＝一部正しい**

根拠：`011/sections/03-determinant.tex:37–58`、`audits/011/REPORT.md:38–51`。

含意自体は通ります。対象範囲では `mn−1=O(x)` なので、`x^0.9` を超える素因数は高々一つです。その素因数ℓに対応する複素位相をFで選べば、法についての絶対値を取り出せます。Fは、大きい対象素因数がちょうど一つならその位相、それ以外なら0と全域で定義すれば、帯の素数に対する不変性と有界性を満たします。

ただし通常のEHそのものでも、無重みの畳み込みの分布でもありません。「nothing comparable」は監査者の知識の申告として許容できますが、比較対象が曖昧です。同条件の既存評価を今回特定できなかった以上の断言は避けるべきです。

Q:10の帰結部分の差し替え：

> With α, β and the marked weight 𝒲 exactly as in Theorem 3.1, we infer, for fixed 0 < ε < 0.1,  
> Σ_{x^0.9 < ℓ < x^{1−ε}, ℓ prime} |Σ_{mn ≡ 1 (mod ℓ)} α_m β_n 𝒲(mn−1)| ≪ X(log x)^{−D_*}.  
> Does this marked, fixed-residue-class consequence follow, and is a comparable estimate known?

この形なら解析的整数論の専門家に十分明確です。「plausible」は主観的なので、具体式の成立と既知文献を尋ねる方が有用です。

**(3) 142のq=2：判定＝正しい**

根拠：`142/sections/00-introduction.tex:46–55`、`10-analytic.tex:91–98`。番号は `main.tex:18–26,60–70` のカウンタと入力順にも一致します。

n=2を固定すると `B=23+3⌈log₂p⌉`。十分大きいpでは `p>B^200000` が成立します。補助素数はℓ≡1 mod 24かつ `(p/ℓ)=−1`。二次相互法則により `(ℓ/p)=−1` です。したがって最小の素数平方非剰余に `≪(log p)^20000000` が出ます。

Q:20の差し替え：

> Taking n = q = 2 in Proposition 11.3 (`10-analytic.tex:91–98`) gives ℓ ≤ c₀(23 + 3⌈log₂p⌉)^{20000000}, with ℓ ≡ 1 (mod 24) and (p/ℓ) = −1. Quadratic reciprocity gives (ℓ/p) = −1. Thus, if the proposition is correct, the least prime quadratic non-residue is ≪ (log p)^{20000000} for every sufficiently large prime p. The proposition uses the zero-free theorem discussed in item 2.

現在の「unconditionally」は、**命題が正しければ**という条件下では正確です。

**(4) 004bとKato：判定＝一部正しい**

根拠：`004b/sections/cyclotomic.tex:423–438`、`introduction.tex:107–137`、`pointwise.tex:1029–1048`、`004/sections/03-parity.tex:457–494`。

- **Theorem 12.5(3)という引用番号は正しい。**
- 原典の12.5(3)はpを反転した高さ1の評価で、局所H²補正を伴います。
- **p≠2と「像がSL₂(Zp)を含む」を課すのは12.5(4)。** Remark 12.8(12.8.2)は非CMの場合に開部分群を含むことを述べています。したがってQ:26の問いは、主として原典に既に答えがあります。[Kato、印刷頁222–223、226、233](https://numdam.org/item/AST_2004__295__117_0.pdf)。
- CMの場合は004b自身がBurungale–Tianを使うと明記しています。「whole branch」は広すぎます。
- 004が必要とするcorank 1側という整理は妥当。ただしその証明もeven側の入力を使うため、二つの独立した証明経路ではありません。

Q:24–27の差し替え：

> **4. Are the cyclotomic and ring-class inputs in result 004’s companion justified?**  
> For non-CM curves, the corank-zero argument uses Kato’s rational Theorem 12.5(3), including its local H² correction, at p = 2 (`cyclotomic.tex:423–438`); it explicitly does not use the stronger integral statement 12.5(4). CM curves use a separate cited result. Is the passage from this rational divisibility to the determinant comparison valid?  
> Result 004 needs the corank-one argument, which also uses the even input. Can a specialist assess the simultaneous arithmetic interface identified in Remark `rd:interface-status` (`ring-determinants.tex:1760–1768`)?

**(5) 077とKakeya：判定＝一部正しい**

**番号は誤り、帰着は正しいです。** Q:34と報告書の双方を訂正してください。

根拠：[077/preamble.tex:18](./data/openai_math/preprints/Elliptic-capacity-propagation-and-Fourier-restriction-to-the-sphere-September-24-2026/build/preamble.tex:18) の共通カウンタ、`01-geometric.tex:32` のDefinition、`同:107` の定理から、**Definition 2.1、Theorem 2.2**です。107行という位置は正しい。

方向をδ間隔で一様に選び、各方向に一本の線分を置き、時間を独立一様にすると、s=1、D=1、K=O(1)で仮定を満たします。保持質量δ^{o(1)}と各δセルの上界から、管の和集合の体積下界はδ^{10κ+o(1)}、または任意のη>0について `≳_{κ,η}δ^{10κ+η}`。κを任意に小さくできるのでMinkowski次元3です。

Wang–ZahlとOpenAI Kakeyaを入力にしないという依存関係も支持できます。`00-introduction.tex:149–184`、`03-projection.tex:4–8,42–72`では、幾何部の外部入力はRen–Wangです。

Q:34の差し替え：

> Theorem 2.2 (`01-geometric.tex:107`, label `geo:main`), applied with s = D = 1 to one segment per δ-separated direction and independent uniform time, appears to give |⋃T| ≳_{κ,η} δ^{10κ+η} for every η > 0, after accounting for the retained sublaw. Letting κ → 0 gives full Minkowski dimension. Is this reduction correct? If so, §§2–5 would provide an alternative proof of the known 3D Kakeya set theorem without taking Wang–Zahl or OpenAI’s Kakeya theorem as inputs.

**(6) 285の非a-T-menable性：判定＝正しい**

根拠：`285/sections/setup.tex:6–17`、`introduction.tex:40–52`、`conclusion.tex:8–13`。

base groupはねじれなしで、行列射影の**非正規化**トレースが0と1/2の間です。a-T-menableならHigson–Kasparovによってassemblyが全射となり、その像上の整数トレース性と矛盾します。[Higson–Kasparov、Theorem 9.1](https://nigel.higson.ca/uploads/1/2/1/4/121496570/higson_kasparov_-_2001_-_e-theory_and_kk-theory_for_groups_which_act_properly_and_isometrically_on_hilbert_space.pdf)。

単なる「非自明な行列射影」からは従いませんが、この論文には必要な非整数トレースがあります。独立した幾何的理由を求める質問は適切です。その理由の不記載自体は証明の欠落ではありません。

Q:44の差し替え：

> If the claimed matrix construction is correct, the torsion-free base group G_base has a projection in M_N(C_r*(G_base)) with nonintegral unnormalized trace (`setup.tex:6–17`). By Higson–Kasparov and trace integrality on the assembly image, G_base cannot be a-T-menable. Is there an independent geometric proof of this from its construction?

**(7) 既知・自明な問い、追加候補：判定＝一部正しい**

- **項目4のKatoの仮定の問いは、原典の有理版／整係数版の区別でほぼ解消します。** 現状のまま専門家に投げるのは避けるべきです。
- 項目3の含意は二次相互法則による短い確認です。確認質問として残せますが、証明済みの帰着として項目2にまとめてもよいでしょう。
- 3D Kakeyaの集合版はWang–Zahlが既に証明しています。ただし077の**別証明としての独立性**を問うことは有用です。[Wang–Zahl](https://arxiv.org/abs/2502.17655)。
- 321の後続論文は系2.3を条件としており、塔の証明の検証を代行してはいません。

追加候補は一つです。README:66、`audits/142/REPORT.md:145`、`142/sections/00-introduction.tex:61–80`、`10-analytic.tex:172–178`にある、**029から独立した代数的帰着**は重要です。

追加できる英文：

> Does Theorem 1.2 of result 142 correctly reduce deterministic polynomial-time factorization to auxiliary primes of polynomial numerical size, independently of result 029? In particular, is its resulting polynomial-time algorithm under GRH valid?

**(8) 語調とREADME：判定＝一部正しい**

根拠：Q:3、10、17、27、34、47と [README:9](./public_repo/README.md:9)。

READMEの新文は適切で、必須修正はありません。Q全体も「誤りの告発ではない」という姿勢は保たれています。ただし、具体的な問いを置ける所で「plausible」「deserve reading」「Has anyone read it?」とするのは、専門家の時間を使う目的に合いません。

Q:3の差し替え：

> These are specific questions arising from our audit, not claims that the manuscripts contain errors. The audits were produced by AI agents under the direction of a non-mathematician. A short correction, reference, or assessment of any stated implication would be valuable; please open an issue.

READMEをさらに簡潔にするなら：

> **For specialists:** [QUESTIONS.md](QUESTIONS.md) collects eight specific questions from the audit, with source references.

全体として、公開の障害は多数の数学的誤りではなく、**077の番号、004bの既知事項への質問、帰結に残る条件の省略**です。この三点を直せば、専門家向けの入口として十分使えます。