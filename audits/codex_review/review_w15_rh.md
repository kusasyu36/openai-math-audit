**依存関係の更新は概ね妥当です。** Axiom の命題は、029・142 が必要とする零点なし領域を含みます。ただし、Alpöge の命題の説明などに具体的な誤記があります。原論文の結論を崩す誤りは、今回特定していません。

以下、`R`＝[対象 REPORT](audits/rh_alpoge/REPORT.md)、`S/`＝`audits/rh_alpoge/sources/`、`029/`・`142/`＝指定論文の `build/` とします。ファイル変更、ビルド、収録コードの実行、ネットワーク利用はしていません。

**(1) Axiom の命題忠実性：判定＝正しい。証明の検証は確認できず。**

`main` は、すべての数体 K について、ℂ∖{1} 上で正則で、Re s > 1 では Dedekind ゼータ関数と一致し、Re s > 31/32、s ≠ 1 で非消滅となる関数の存在を述べています。解析接続の一意性により、通常の解析接続についての主張です。s = 1 の極の形自体は、この命題には書かれていません。根拠：`S/axiom_DedekindZetaChallenge.lean:30–32`。

`RayClassCharacter` に、挙げられた抜け道は認めません。

- `m = ⊤` は有限部分が法 1 の狭義類群に対応し、問題ありません。`m = ⊥` は系の `hm` が除外します。
- `I = ⊥` では値が 0。その項は級数を変えません。
- 整イデアル上の非零値は、法と素な分数イデアルの群へ乗法的に延長できます。
- 整な生成元だけへの条件でも十分です。分数の射類主元 α に対し、分母を払う全正の β ≡ 1 mod m を中国剰余定理と近似で選べば、αβ と β の二つの整な主元への条件から、(α) への条件が得られます。
- 全実埋め込みでの正値条件は、すべての実素点を法に含めた狭義射類群を与えます。符号成分のある有限位数 Hecke 指標も含みます。複素素点の連結部分には非自明な連続有限像指標はありません。

したがって、固定した非零 m では「導手が m と全実素点からなる法を割る指標」の通常のイデアル表示です。射類群の有限性から値の有限位数性も従います。異なる m による表示の重複はありますが、余分な指標ではありません。根拠：同 `14–28,34–41`。

ただし、保存された challenge の証明は `sorry` です。これは課題文の書式であり、solution に欠陥がある証拠ではありません。README の証明済みという記述と comparator の設定だけから、実際の証明成功までは確認できません。根拠：`S/axiom_README.md:24–25,48–52`、`S/axiom_comparator.json:2–12`。

**(2) 029・142 への含意：判定＝正しい。ただし「解析全体が確認された」は不可。**

029 定理1.2は `heckeCorollary` の範囲内です。31/32 < 1 − 10⁻⁶ であり、原始指標と非原始指標の差である有限個の Euler 因子は Re s > 0 に零点を作りません。有理型接続との同定も、一致する右半平面からの一意性で足ります。根拠：`029/sections/01-introduction.tex:32–53`、Axiom challenge `34–41`。

142 の K・M はともに数体なので、`main` が直接、補題 `analytic:explicit` の零点なし仮定を供給します。Hecke 分解を経由する必要もありません。根拠：`142/sections/10-analytic.tex:39–47,101–111`。

正確な言い換えは、

> Axiom の形式化が正しければ、142 が必要とする一様な零点なし領域は029の証明を使わずに得られる。明示公式、一様零点計数、判別式評価、および142自身の帰着の正しさは別途必要である。

です。後者は同 `50–86,113–157`。R:8 の「危険が置き換わる」は、**零点なし入力について**と限定してください。

**(3) Alpöge の範囲：判定＝一部正しい。**

五命題が ζ と Dirichlet L 関数についてだけである点、および一般の Hecke L 関数を直接含まない点は正しいです。根拠：`S/alpoge_Statements.lean:3–16`、`S/mirror_tamnd_README.md:40–50`。

ただし、具体的な訂正が二つあります。

- **R:17「ディリクレの4つ」は誤記。** 独自定義四つの内訳は、ζ 版一つ、Dirichlet 版三つです。
- **R:20「ディリクレ版の命題は 0 < Re s < 1 の零点だけを対象」は誤り。** 非消滅命題は **Re s > 1/2 全体**を対象とします。帯内の零点を対象にするのは残る二命題です。

142 の M = Q(μ₁₂q, p^{1/q}) が q ≥ 3 で Q 上非アーベル、という理由は正しいです。例えば、含まれる実部分体 Q(p^{1/q}) は Eisenstein により次数 q で、非実な共役根を含まず非正規です。アーベル拡大の部分体なら正規になるため矛盾します。142 の適用条件も `10-analytic.tex:93,103–108` で確認できます。

ただし **029定理1.2の基礎体 F 自体は円分体でアーベル**です。R:9 は次のように直すのが正確です。

> 029 は円分体上の一般の有限位数 Hecke L 関数を扱う。その対応する類体や、142で必要となる Kummer 体は Q 上非アーベルになり、Dirichlet L 関数についての五命題の直接の対象ではない。

ERH を今後の期待としているという要約は正確です。根拠：`S/alpoge_remark.txt:6–10`。共同声明が認定対象として RH を述べ、Dirichlet 版には言及しないという R:18 も正確です：`S/lean_fro_joint_statement.txt:3`。

**(4) q = 2：判定＝体の包含と含意は正しい。Ankeny・Bach の精密な引用は確認できず。**

奇素数 p に対して、p* = (−1)^((p−1)/2)p と置けば √p* ∈ Q(μₚ)。i ∈ Q(μ₂₄) なので、

\[
\mathbb Q(\mu_{24},\sqrt p)\subseteq\mathbb Q(\mu_{24p})
\]

です。したがって K・M は Q 上アーベルで、そのゼータ関数の Dirichlet L 関数への分解から Dirichlet GRH が必要な領域を与えます。Axiom の `main` からも直接得られます。142 の残りの議論を用いれば、同じ q = 2 の補助素数評価が従います。根拠：`142/sections/10-analytic.tex:101–157`、補助素数の条件は `142/sections/00-introduction.tex:46–54`。

**Ankeny 1952、Bach 1990、係数 2、適用条件や定理番号を原典と照合できる資料は、今回の保存資料にはありません。** R:10 の「未確認」は維持し、精密な書誌や `2 log² p` を確認済みとして公開しないでください。

**(5) 後続論文：判定＝概ね正しい。ただし直接依存と手法利用を分けるべき。**

| 行 | 判定・根拠 |
|---|---|
| Liu 2610.12234 | **正しい**：003の議論のパラメータ変更であることは `S/arxiv_followups.txt:16–19` に明記。ただし、この論文自身の参考文献による「9月30日版」の確認まではできません。 |
| Xu–Sang | **正しい**：OpenAI の解析入力を仮定すると明記。同 `31–34`。 |
| Liu 2610.26835 | **正しい**：[Ope26a]＝003、[Ope26b]＝029 は `S/liu_2610.26835v2_abstract_ack_bib.tex:125–135` で確定。029の構成・局所計算の利用も同 `23–25` に明記。 |
| Chang–Wang | **一部正しい**：引用された一般体の入力への条件付き結果という説明は支持されます。ただし、各入力と003・029との対応は、この要旨だけでは確定できません。同 `arxiv_followups.txt:1–4`。 |
| Shunia | **正しい**：直接使用すると書かれているのは Cureton–Thorner。同 `26–29`。029を直接仮定する論文、とまでは言えません。 |
| Li–Zhong | **正しい**：OpenAI の Dirichlet 零点なし定理の使用が明記。同 `21–24`。 |

R:36 の一括した「いずれも、OpenAIの議論が正しいことを前提」は、次に置き換えるのが安全です。

> 保存した要旨には、OpenAIの結果の利用、議論の改変、または関連する入力への条件付き依存が認められる。これらは元の全証明の独立検証を示す資料ではない。

**検証と誤記している表の行はありません。** R:36 は明確に「検証ではない」としています。

ただし **R:42 の「Codex が下書き・検算」は謝辞の正確な要約ではありません。** 原文は下書き、参考文献確認、改良境界の探索。厳密な記号計算には SymPy を使用したとあります。根拠：Liu の保存 tex `31–35`。

Chang–Wang の行は、直前の数値からの単純な差と誤読されないよう、「要旨の B* から 1/160000 改良、終点 0.87446205769353…」とするのが明瞭です。

**(6) 029 の残る依存：判定＝概ね正しい。ただし範囲の列挙には訂正が必要。**

定理1.2の証明を構成するファイルは、

`02-arithmetic`、`02a-coefficient-family`、`03-theta`、`04-reflection`、`05-poisson`、`06-zero-free`、付録 `10-whittaker-calculation`

です。組版上は **§§3–8 と付録A**。根拠：`029/main.tex:76–89`、`029/sections/01-introduction.tex:144–145`、証明の終端は `06-zero-free.tex:754–764`。

したがって、**ご質問の列挙にある `02-arithmetic` は証明側**です。これは R:11 自体の誤記ではありません。同行にはファイル列挙はありません。

Axiom を受け入れた後に残る主な議論は、`01a-reduction`、`07-splitting` の零点なし入力以外の部分、`08-type-ii`、`11-sieve`、`12-construction` とその入力です。実際、`main` を K_q に直接適用すれば、`07-splitting.tex:70–98` の Hecke 分解による零点なし領域の移送も迂回できます。

また「011定理3.1だけを外部入力とする」と狭めてはいけません。係数移送の議論は011の**証明内部**にも依存しています。根拠：`08-type-ii.tex:92–97,125–127,187–202,210–240`。R:11 の「その他すべてとその入力」という広い整理なら、この点も含められます。

**(6b) 003だけで q = 2：判定＝正しい。ただし003の結果を受け入れる条件付き。**

決定的な根拠は、

`data/openai_math/lean/ComparatorChallenges/DirichletSevenEighths.lean:8–11`

です。法 q は任意の正の自然数、χ は任意の Dirichlet 指標で、原始性・導手・高さの制限はありません。除外するのは主指標の s = 1 だけです。`lean/docs/003.md:10–12` の説明とも一致します。

よって003の一様な 7/8 半平面を K・M の Dirichlet 分解に適用すれば、142が要求する 1 − 10⁻⁶ の領域は既に得られます。**新たな解析入力が必要な帰結として QUESTIONS 項目3を扱ったことは訂正すべきです。w14でもこの点を見落としていました。**

ただし R:10 の「もう専門家に聞く意味がない」は強すぎます。適切なのは、

> q = 2 の帰結に、003を超える零点なし定理は必要ない。ただし003の証明や、142の明示公式から補助素数への議論の検証とは別である。

です。

**(7) 公開英文：判定＝§4の修正方針は一部正しい。以下を提案します。**

QUESTIONS の項目2：

> **2. Can the zero-free input for results 029 and 142 be independently assessed?**  
> If correct, Axiom Math’s Lean formalization of a 31/32 zero-free half-plane for all Dedekind zeta functions and finite-order Hecke L-functions supplies the input required by both papers. We have not built it and are not aware of a public review by a domain expert. Separately, our audits have not verified result 029’s own proof of Theorem 1.2. → `audits/rh_alpoge/REPORT.md`, `audits/029core/REPORT.md`.

項目3は質問から訂正注記へ：

> **3. Correction concerning the q = 2 case of result 142.**  
> Assuming result 003’s formalized uniform 7/8 theorem for Dirichlet L-functions, the zero-free hypotheses needed for this case already follow, since Q(μ₂₄, √p) is abelian over Q. The stated polylogarithmic bound therefore requires no stronger zero-free input than 003. This does not verify the remaining argument of Proposition 11.3.

項目4：

> **4. Is the algebraic reduction of result 142 valid?**  
> Does Theorem 1.2 correctly give deterministic factorization in time polynomial in B and the numerical sizes of the supplied auxiliary primes? This reduction is independent of result 029’s zero-free proof and remains to be assessed even if Axiom Math’s formalization supplies the analytic input. → `audits/142/REPORT.md`.

README の029行：

> The zero-free input would follow from Axiom Math’s stated 31/32 Lean formalization, if correct. Result 029’s reduction, splitting estimates, marked Type II argument using result 011, sieve and prime construction still require expert review. We have not verified the main theorem.

029 (core) 行：

> Neither audit located a concrete error in the portions examined, but neither verified Theorem 1.2’s proof. Its conclusion is also implied by Axiom Math’s stated stronger formalization, if correct; this does not validate 029’s argument.

029 / 142 (applications) 行：

> Axiom Math’s stated formalization would supply the required zero-free hypotheses independently of 029’s proof. The applications also use degree and discriminant bounds, uniform zero counting and explicit formulas. Our reading found no concrete error in the portions examined; it was not a verification of the full arguments.

142 行：

> The proof as written cites result 029 for its zero-free input. Axiom Math’s stated Dedekind-zeta formalization would supply that input independently, if correct. The algebraic reduction and the deduction of polynomial-size auxiliary primes still require expert review. The q = 2 zero-free input already follows from result 003’s stated Dirichlet theorem.

これらの行に対する共通注記：

> We have not built Axiom Math’s formalization and are not aware of a public review by a domain expert.

003の三行は既存の検算記録を残し、共通して次を追記：

> Alpöge’s stated RH and Dirichlet GRH are stronger than 003’s zeta and Dirichlet conclusions. General Hecke L-functions over Q(√−3) are not among the five statements checked; those statements do not verify the reflection or coefficient calculations audited here.

**(8) 日付・版・帰属：判定＝一部正しい。**

- Alpöge の11月20日公開、共同声明の11月21日、Sutherland ページの11月24日更新は保存資料と一致します。`S/lean_fro_joint_statement.txt:3`、`S/sutherland_grh_page.txt:3`。Comparator v4.33.0 成功と v4.36.0 実行中の区別も正確です。
- **R:24 のリポジトリ作成日、最終コミット日、CI成功は、対応する日時・CI記録が保存資料にないため独立確認できません。** Liu の参考文献は同じコミットを10月15日参照としていますが、10月13日の公開日そのものの証拠ではありません。`S/liu_2610.26835v2_abstract_ack_bib.tex:52–56`。
- Axiom README は論文 **v1** にリンクし、Liu の [CT26] は **v2** を指定しています。この版の違いは明記すべきです。v1の主張が異なったという話は確認できません。`S/axiom_README.md:5–6`、Liu tex `45–50`。
- Cureton–Thorner が研究と数学的判断を担当し、AxiomProver が証明と形式化を生成したという帰属は正確です。`S/axiom_README.md:8–10`。
- **R:54 の「独立の第2の証明」は独立性未確認。** Axiom README は openai/math の解析数論コードの利用を明記しています。「別の形式化による含意」としてください。同 `39–44`。

総評として、Axiom を条件付きの代替入力として依存図に追加し、q = 2 の入力は003で既に足りたと訂正する方針を支持します。これは029・142の全証明の認定ではありません。修正すべき具体的な誤記は、Alpöge の定義の内訳・非消滅命題の範囲、およびLiu謝辞におけるCodexの役割です。