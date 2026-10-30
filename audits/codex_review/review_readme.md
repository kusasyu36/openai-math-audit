**公開前に、主に検証範囲と断定の強さを修正すべきです。** ファイルは変更していません。

1. **必須** — “released 720 … About 42% … the rest have not been machine-checked”
   → 初回コミットには原稿ディレクトリが722個あり、現在の目録は721本です。42%は上流の自己申告で、`history.md` は10月7日時点の300/719と説明しています。「Lean形式化なし」から「機械検証なし」は断定できません。
   → **“OpenAI released this collection on 6 October 2026. The catalogue updated on 26 October lists 721 manuscripts. OpenAI reports that approximately 42% of top-line results have accompanying Lean formalizations.”**

2. **必須** — “Three manuscripts were withdrawn the next day, and all three were in the unformalized part.”
   → 翌日の3本撤回は履歴で確認できます。ただし、撤回当時の3本すべての形式化状況は、取得済みローカル資料では確認できません。
   → **“OpenAI’s history records the withdrawal of three manuscripts on 7 October 2026.”**

3. **必須** — “For each manuscript … recomputed every explicit constant, exponent and finite identity with exact … arithmetic”
   → 全件・全項目・厳密演算という記述は報告書を超えます。048では高精度数値計算、101ではサンプリングも使っています。
   → **“Depending on the manuscript, we checked selected constants, exponents and finite identities using exact arithmetic or high-precision numerical calculations, searched for counterexamples, and examined selected arguments.”**

4. **任意** — “We never executed code … because checking a result with the author’s own tools is not independent.”
   → 各報告の実行方針とは整合しますが、“never”という履歴全体の否定はファイルだけでは証明できません。また、付属ツールの再実行にも再現性確認として価値があります。
   → **“For the computational checks reported here, we used independently written scripts rather than the supplied verification scripts. Reproducing supplied computations and independently checking the mathematics serve complementary purposes.”**

5. **必須** — 074：“version dated 23 September 2026”
   → 4D論文のREADME・本文・TeXの日付はいずれも**9月24日**です。引用番号の誤り自体は報告と整合します。また、4D本文の大部分は未監査です。
   → **“In the version dated 24 September 2026, … This is a citation-number error; we identified no mathematical consequence from it. Most of the 4D argument was not audited.”**

6. **必須** — “Other audited results (no error found)” 内の038：“A minor rigor gap”
   → 見出しと分類が紛らわしく、修復可能性も報告者の判断です。
   → 見出し：**“Additional audits: checks and unresolved questions”**
   → 038：**“A possible rigor gap in §5.2 concerns admissibility of the test valuation. The report considers it likely repairable; Proposition 5.1 still needs expert review.”**

7. **必須** — 171：“Compressed argument (not an error).”
   → 訂正後の§9評価には忠実ですが、報告に残る§5・§16・§18の概略的構成や未検証部分が落ちています。
   → 追記：**“The report also identifies sketch-level or incompletely checked constructions in §§5, 16 and 18. Much of the probabilistic machinery was checked only structurally.”**

8. **必須** — 059：“Finite arithmetic reproduced … real-Seifert-form theorem needs a check”
   → 主論文の要約としては妥当ですが、同じ報告が対象とするC⁴・重複度4と5の姉妹論文では、核心のtube計算が未検証です。
   → 追記：**“For the companion paper in C⁴ with multiplicities four and five, the central integral Seifert-form ‘tube’ computations were not checked.”**

9. **必須** — 101：“About 220,000 … checks … passed. Some margins are extremely thin”
   → 報告の最大の留保である輸送の正則性・写像度・極限操作が省かれ、数値検査で全体を裏付けたように読めます。
   → **“Approximately 220,000 sampled quadratic-form tests and the reported scalar checks found no violation. Transport regularity, the degree argument, approximation and limiting steps still need expert review.”**

10. **必須** — 011：“allows an arbitrary bounded function F”; 029：“inherits 011’s main risk”
    → Fには帯の素数に関する不変性条件があります。EH型帰結は監査側の推論で、専門家未確認です。011・029の全体への低い確信度と、029自身の零点排除定理の未検証性も明示すべきです。
    → 011：**“F is bounded and invariant under multiplication by primes in the specified bands. We infer a strong marked Elliott–Halberstam-type consequence, which needs expert confirmation. Proposition 4.1 was only partly checked; confidence in the full proof remains low.”**
    → 029追記：**“The central zero-free-region argument was not verified, independently of the unresolved Type II input from 011.”**

11. **必須** — 077：“3D Kakeya set conjecture … those ~4,700 lines are unaudited”
    → 報告の帰結は**Minkowski形**です。またhairbrush等の一部は確認済みなので、“unaudited”は確認実績を弱めています。
    → **“The report derives a Minkowski-dimension Kakeya consequence from the core geometric theorem. Its approximately 4,700 lines were largely unchecked, apart from selected estimates.”**

12. **必須** — 284：“The paper proves R ≳ Q^{4−2/k} … ‘optimal quartic separation’ overstates this.”
    → 対数因子が抜けています。“proves”や“overstates”も査読前の監査として断定的です。「最適指数4」と「単一族で四乗を対数因子まで達成」は区別できます。
    → **“The paper claims, for each fixed k ≥ 2, a separation R ≥ cₖQ^{4−2/k}/(log Q)^{Cₖ}. The exponent approaches four; the reported bounds do not establish a single quartic-up-to-polylogarithms family.”**

13. **必須** — “We checked 12 results … no undisclosed mismatch.”
    → 結論自体はfaithfulness報告と一致します。ただし、084のdyadic限定などの開示済み制限と、証明本体をビルド・検証していないことが欠落しています。
    → **“We inspected Comparator statements and definitions for 12 selected results, allowing for disclosed scope restrictions, such as the dyadic-only formalization of 084. We found no undisclosed mismatch in that inspection. We did not build or verify the proof implementations.”**

14. **必須** — “Mathlib’s contribution rules require human-written, expert-supervised code.”
    → 一律のAI生成禁止という要約は不正確です。公式方針はAI利用の開示、人間による理解・説明責任等を求めています。[公式方針](https://leanprover-community.github.io/contribute/index.html#use-of-ai)
    → **“This project is not intended for submission to Mathlib. Mathlib permits some AI-assisted contributions, subject to its contribution standards, disclosure requirements and human responsibility.”**

15. **任意** — 164：“§2 is sound. The analytic principles (§3–8) are unaudited.”
    → 報告に沿いますが、§2が二原理を仮定した帰着であることを明記すると公平です。§6.1等には部分検査もあります。
    → **“We found no issue in §2’s reduction assuming the Prediction and Alignment principles. The analytic core remains largely unchecked.”**

確認できた点：公開日、Lean v4.34.1／Mathlib `d13f23b7`の上流との一致、Issuesの074・171の原典再確認記録は整合しています。Lean表の4ファイルの内容も概ね正確で、唯一の実コード上の`共通sorry`は`theorem_1_1`です。他の主要定理の既存ビルド記録は記載の3公理のみを示し、二つの帰結は定理1.1を**仮定した含意**です。本レビューでは再ビルドしていません。

公開フォルダには原稿の`.tex`・`.pdf`・`paper.txt`はなく、長い英文転載も検出されませんでした。再配布説明に反する具体物は見当たりません。

**総合判定：修正後なら公開可能。現状は、検証範囲の過大表示と重要な留保の省略があるため、上記「必須」を直してからの公開を勧めます。**