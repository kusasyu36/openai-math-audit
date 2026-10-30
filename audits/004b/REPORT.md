# 監査レポート 004b: A pointwise 2-converse for elliptic curves with rational two-torsion（OpenAI、2026-10-07、91 ページ）

- 監査日: 2026-10-30（監査役: Claude サブエージェント。報告書は主担当が代筆保存）
- 位置づけ: 結果 004（ℚ 上の Hilbert 第 10 問題）の Lemma 3.4 が全面的に依存している土台の論文。
- 自作スクリプト: `audits/004b/descent2.py`（独自の完全 2-descent）、`checks.py`（明示的な計算の検算）。リポジトリのコードは実行していない。

## 結論

**明確な誤りは見つからなかった。** ただし、004 が実際に依存しているのは、この論文の中で検証できなかった部分（rank 1 の枝）である。

## 1. 主定理（Thm 1.1）

- **仮定**: E(Q)[2] ≠ 0、かつ s_2(E) := corank Sel_{2^∞}(E/Q) ∈ {0,1}（有限 Sel_2 の次元ではない）。
- **結論**: ord_{s=1} L(E,s) = rank E(Q) = s_2(E)、かつ Sha は有限。
- **証明の骨格**: Monsky の 2-parity → 解析的位数 = s_2（偶数側 pw:even-base、奇数側 pw:odd-base）→ Gross–Zagier–Kolyvagin。

## 2. 004 との接続

- 接続は正しい。004 は「C(Q)[2] ≠ 0 かつ corank ≤ 1 ⇒ Sha 有限」だけを使っている。
- **004 が依存しているのは奇数側（rank 1）の枝。** Cassels–Tate ペアリングにより、dim Sel_2 − 2 ≡ s_2 (mod 2)。よって Sel_2 の次元が 3 なら s_2 = 1 となり、004 は ring-class、odd 係数検出子、グラフの機構に全面的に依存する。この論文で最も重く、未検証の部分である。
- E_l は CM ではない（j = 35152/9）ので、CM 用の入力（Burungale–Tian 2026）は 004 には無関係。判別式が正なので、Δ < 0 の問題も起きない。

## 3. 構成と新規性

- **corank 0（cyclotomic）**: Kato のゼータ元を、二進立方体 G = (Z/2)^b 上の群環の元として扱う。Ferrero–Washington、群環の整数性、Fourier 交叉、中心値公式、Chevalley–Warning で、欠けた頂点での非消滅を示す。
- **係数テスト**: Waldspurger、重さ 3/2 の形式の偶奇（symbol）、trace と inertia の恒等式。
- **グラフ**: forest と prime network、address lemma。
- **rank 1（ring-class）**: 補助素数による変形、有限 Artin モデル、ultrafilter の極限、Kolyvagin の derivative、outer DVR による一様下界、missing-vertex 転送。
- **p = 2 での古典的な障害**（R = T、Ihara 型補題、Eisenstein イデアル、Heegner 点の mod 2 非消滅）は一切使わず、設計段階で迂回している。その結果、負担はすべて新しい構成に移っている。
- 証明は非構成的（定数は計算できない）。
- 006（Goldfeld）は「入力ではない」と明記されていて、循環依存はない。

## 4. 確認した項目（すべて整合）

- cy:binary（再構成と乱数 200 試行）、cy:fourier、符号単数、smoothing 因子。
- co:operators（1492 例）、co:symbol-interface（37128 例）、pw:good-simple（361 例）。
- rd:exact-frobenius、rd:residual-euler（記号計算）、根数の符号。
- **独自の 2-descent**:

| 曲線 | dim Sel_2 | s_2 の上界 | 既知の rank | 判定 |
|---|---|---|---|---|
| 合同数曲線 n = 1, 2, 3 | 2 | 0 | 0 | 整合 |
| n = 5, 6, 7, 14 | 3 | 1 | 1 | 整合 |
| n = 34, 41 | 4 | 2 | 2 | 定理の範囲外 |
| E_l（7 例） | 3 | 1 | — | 004 の監査と一致。5 例は有理点から rank ≥ 1 |

- 読んで整合を確認したもの: 局所指数の勘定、rd:active-local、rd:sha-comparison、2 での降下、Δ < 0 を回避する局所モデル。

## 5. 懸念点（明確な誤りはなし）

1. **【要専門家・中〜高】p = 2 での Kato 12.5(3) の有理整除性**（cyclotomic.tex）。偶数側の全体がこれに依存している。Kato の証明が p の奇数性を仮定していないか、Remark 12.8 などで「像 ⊇ SL₂(Z_p)」を「open image」に置き換えてよいか、を確認する必要がある。
2. **【要専門家・高、004 にとって決定的】ring-class の機構全体**（ring-limits.tex 121–981、ring-determinants.tex）。論文自身も Remark rd:interface-status で、既存の定理だけでは得られないと認めている。一様下界 C(E,A) は非構成的。
3. **【要専門家・中】odd 係数検出子**（coefficients.tex 254–564）。
4. **【要専門家・中】address lemma と packet・reserve の組合せ**（graphs.tex、pointwise.tex）。
5. **【低、004 には無関係】** CM の場合は Burungale–Tian 2026 に依存している。
6. **【説明不足・低】** bib のキー `Howard2012` が 2004 年の論文を指している。cy:missing-vertex の定数が、本文と証明で書き方が異なる（吸収されるので問題ない）。

## 6. 要約

接続は正しいが、004 は最も重い rank 1 の枝に依存している。p = 2 の古典的な障害を、新しい構成で迂回するという設計。偶数側は筋が通っていて、明示的な計算と 2-descent の標本はすべて整合した。p = 2 での Kato 12.5(3) と、ring-class・ultrafilter・グラフの機構は未検証。確信度は、偶数側が中、004 に必要な奇数側は低〜中。
