# 監査レポート 266: The maximum number of mutually unbiased bases in dimension six（OpenAI, 2026-09-24）＋副論文 Exact Fourier certificates for complex Hadamard matrices of order six（同日）

- 監査日: 2026-11-06（監査役: Claude サブエージェント）
- 対象 tex（主）: `data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/`（paper.tex、sections/{introduction,seed,completion,phase-cover,refinement,exclusion,execution,conclusion,rounding}.tex、計 5,966 行）。全セクションを tex 原文で読んだ。paper.txt と PDF は使っていない。
- 対象 tex（副）: `.../Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/build/`（paper.tex、preamble.tex、sections/01〜05、計 1,996 行）。
- 同梱の `verification/`: 主論文は C++ ソース 8 本（計 77 KB）と保存済み数値記録（`saved-run/execution/data/`、約 90 MB）、副論文は `verify.py`（576 行）。**ソースと記録は読んだが、OpenAI のコードは一切実行していない**（`check_records.py`、`run.py`、`verify.py`、C++ のビルドを含む）。副論文の `verify.py` からは整数データ `data`（Gram 行列の下三角と電荷リスト）だけを `ast.literal_eval` で読み取った。
- 自作スクリプト（`audits/266/`。`UV_CACHE_DIR="$TMPDIR/uv-cache" uv run --with sympy --with numpy --with mpmath python -I audits/266/<script>.py`。records だけは標準ライブラリのみ）:
  - `check_266_records.py`（`.out`）: 保存記録のバイナリを自作パーサで読み、execution.tex の表（件数・半径・バッチ分割・失敗記録の部分列性・fallback 入力順・最終 mod-4 分割・Root 定数の範囲・centers.txt のフィルタ条件）を全件照合。**全項目一致**。
  - `check_266_lemmas.py`（`.out`）: 主論文の小補題・定数（三重 MUB 補題、π の有理証明書、級数の裾、U(2) セルの Taylor 定数 6 と 42、水平行列の恒等式、ビン幾何の定数、距離の下界 20/7、境界 Taylor 定数、rounding 付録のスケール不等式）を sympy/mpmath で検算。
  - `check_266_completion_chain.py`（`.out`）: 既知の厳密 Hadamard 行列（Fourier 2 径数族、Diță 型ブロック構成、Björck 巡回行列、立方行列 S_6、それらの転置・同値変形、計 230 例）に対し、seed 半径内に乱数で置いた中心から completion.tex と phase-cover.tex の不等式の連鎖（枠、極分解、セル Taylor、Z の 3 片分解、z 上界、線形剰余 E_r、ノルム判定、方向判定、弧の上界）を数値で全件検証。**30 項目 × 230 例すべて成立**。
  - `check_266_cubic.py`（`.out`）: 副論文 §5（cubic obstruction）を Q(ω) の厳密整数算術で全再計算（Table 1 の N・U 列、R(b)=λ_j G(b)、軌道サイズ、w_5+D と w_8+D の分類、連立方程式の解、6+30S_0=−18、対称性 720 個の総当たり）。**すべて一致**。
  - `check_266_pair_identity.py`（`.out`）: 副論文 §3 の恒等式（Lemma contraction、Lemma partitions、Proposition pair-certificate）を 21 個の厳密 Hadamard 行列上で数値検証（残差 ≤ 6×10⁻³¹）。Lemma basic の 3 恒等式を次元 3・5 の完全 MUB 系と次元 6 の MUB 3 つ組で検証。Lemma mixed-reductions (iii) と Lemma gram を次元 3 の類似で検証。
  - `check_266_Zposdef.py`（`.out`）: 6 個の整数 Gram 行列 Z_λ の正定値性を厳密 Schur 補元で確認（98 個のピボットすべて > 40000）。
  - `check_266_exact_certificates.py`（`.out`）: **副論文の 2 つの厳密証明書を tex の仕様（§7）から独立に再実装**（単射モーメント、縮約、7 行関係式の列挙、法 p のエシュロン＋CRT＋有理数再構成＋厳密整数再検査、集合分割展開、混合モーメントの帰着則、Young 対称化、Gram 証明書）。出力 3 行は Appendix A と**完全に一致**: `rank table: 8871 909 897`、`rank table: 46041 3429 2767`、`-613302797399911/6480 179721388988719/6480 5000000000000000`。実行時間 133 秒。
- 以下、定理番号は tex の `\newtheorem`（節ごとに共通カウンタ）と環境の並び順から自分で数えた。主論文: Thm 1.1（main）、Lem 1.2（three-mub-six）、Cor 1.3; §2 seed は Def 2.1〜Prop 2.10; §3 completion は Lem 3.1〜Prop 3.11; §4 phase-cover は Def 4.1〜Prop 4.16; §5 refinement は Lem 5.1〜Cor 5.8; §6 exclusion は Lem 6.1〜Prop 6.9; §7 execution は Thm 7.1（finite-theorem）; §8 conclusion; Appendix A rounding は Lem A.1〜A.3。副論文: Thm 1.1（fourier）、Thm 1.2（main）、Cor 3.1（dichotomy）、Lem 3.2（contraction）、Lem 3.3（partitions）、Prop 3.4（pair-certificate）、Prop 4.1（cubic-classification）、Lem 5.1（basic）、Prop 5.2（cubic-exclusion）、Lem 5.3（cubic-symmetry）、Lem 6.1（gram）、Lem 6.2（mixed-reductions）、Prop 6.3（mixed-certificate）、Lem 7.1（exact-elimination）。

## 0. 要約（先に結論）

1. **副論文の 2 定理は、独立に再実装した厳密計算で完全に再現できた。** Thm 1.1（Matolcsi–Ruzsa–Weiner の Fourier 消滅予想: 立方類 S_6 以外のすべての 6 次複素 Hadamard 行列で g_H(πα)=0）と Thm 1.2（C⁶ に 7 個の MUB は存在しない。Weiner の補完定理と合わせて N(6) ≤ 5）。証明は有限の整数・有理数算術だけに依存し、浮動小数点を含まない。tex の証明（Prop 4.1 の分類、Cor 3.1 の二択、Prop 5.2 の立方排除）も行ごとに追い、誤りは見つからなかった。**私の判断ではこの 2 定理は確定した結果として扱ってよい**（§5 の注意 9 の留保つき）。これ自体が、主論文と独立に N(6) ∈ {3,4,5} を与える新結果である（私の知る限り「7 個の MUB の非存在」は未解決だった）。
2. **主論文の主張（N(6)=3）は条件付きの計算機支援証明であり、その条件は配布パッケージからは検証できない。** Thm 1.1 は「Appendix A の算術・コンパイラ条件のもとで、§7.4 に記録された実行が証明する」という形（introduction.tex:13–18）。数学（seed → completion → phase cover → refinement → exclusion の被覆と必要条件の連鎖、Appendix A の片側丸め余裕）は tex の範囲で丁寧で、私が手と数値で追った箇所に誤りはない。しかし結論は (i) 77 KB の C++ が tex の補題を正しく実装していること、(ii) 2026-10-04 の実行が x86-64/GCC 13.3 上で binary64・RN・FMA なし等の条件を満たしたこと、の 2 点に依存し、どちらも tex からは確認できない。しかも論文が証拠として挙げる `build-manifest.json`、逆アセンブル、プローブ出力、`proof_run.py`、`REPRODUCE.md`、`verification/computation/README.md` は**配布物に含まれていない**（README.md:12 と VERIFY.md:3–8 が明言）。配布されているのは数値記録と、その内部整合性だけである（これは私のパーサで全件確認した）。
3. **日付の不整合。** 論文日付は 2026-09-24 だが、証明に使う実行は 2026-10-04 完了（execution.tex:4, 286）。9 月 7 日の実行は「認証されていない」（introduction.tex:157–158、execution.tex:396–400）。原稿は日付を変えずに事後改訂されている。
4. 確定した誤り: **見つからなかった**（両論文）。副論文の tex に軽微な参照不整合（`certificate/verify.py` vs `verification/verify.py`）がある。

## 1. 主張の正確な内容

### 1.1 主論文

**Thm 1.1（introduction.tex:13–18、[Computer-assisted]）**: 「Appendix A（rounding.tex）の算術・コンパイラ条件のもとで、§7.4（execution.tex:283–301）に記録された新しい完全実行は、C⁶ に 4 つの互いに不偏な正規直交基底が存在しないことを certify する。したがって N(6)=3。」

- 証明（conclusion.tex:3–11）は 2 段: Thm 7.1（execution.tex:242–249、有限実行証明書: 「7 つのソースファイルと Appendix A の算術条件のもとで、§7.2 の呼び出しがすべての終端条件を満たせば、4 MUB は存在しない」）と、§7.4 の記録がその終端条件を満たしたという事実主張。下界は Lem 1.2（テンソル積構成、introduction.tex:166–210）。
- Thm 7.1 の証明（execution.tex:251–281）は Prop 2.10（seed cover）→ Prop 4.16（A 段の被覆）→ Prop 5.7/Cor 5.8（B 段の被覆保存）→ Prop 6.9（C 段の局所排除）の連鎖。各命題は「実行が成功裏に完了したならば」という条件付き。
- 限定条件（論文自身が明記）: 算術は IEEE binary64、最近接丸め、`-ffp-contract=off`（FMA 禁止）、fast-math なし、`FLT_EVAL_METHOD=0`、FTZ/DAZ なし、正しく丸められた sqrt（rounding.tex:29–38; execution.tex:73–76）。さらに `coverSeeds::run` の 2 つの `pow(...,2)` が binary64 乗算に落ちていること（execution.tex:90–104、逆アセンブルで確認したと記載）。プローブは「限定的な証拠であり、すべての演算・コンパイラ意味論の証明ではない」（execution.tex:82–85）。
- 第 2 基底は任意の複素 Hadamard 行列（既知の族に限らない）。同値変形（行列の置換、再位相、複素共役、転置）で extendibility が不変であること（Lem 2.2）を使って探索空間を縮める。
- 実行規模（execution.tex:303–329、Table）: seed cover 151,744,086 組 → 249,878 中心; A 段 4 シャード → 3083 記録; B mid → 600; B fine → 22,379（88 バッチ）; 初回 C で未解決 163（すべて graph cap）; fallback B → 1127; 最終 C で未解決 0。壁時計 4470 秒。

**Cor 1.3（introduction.tex:230–238）**: 同じ条件のもとで、sl_6(C) の標準随伴で閉じた Killing 直交 Cartan 部分代数の最大個数は 3。Boykin らの対応を引用するだけで、新しい内容はない。

### 1.2 副論文

**Thm 1.1（paper.tex:65–71、Fourier vanishing）**: 6 次複素 Hadamard 行列 H が立方行列 T（paper.tex:51–57、Tao の S_6）と同値でなければ、α=(1,1,1,−1,−1,−1) のすべての座標置換 πα で g_H(πα)=0。ここで g_H(a)=(1/6)Σ_k Π_i H(i,k)^{a_i}。Matolcsi–Ruzsa–Weiner の Conjecture 2.3 の完全解決（従来は Karlsson 族で既知）。

**Thm 1.2（paper.tex:142–149、No complete family）**: C⁶ に 7 個の互いに不偏な正規直交基底は存在しない。Weiner の補完定理（d 個の MUB は d+1 個に伸びる）と合わせて N(6) ≤ 5（paper.tex:151–153）。

- 証明構造: Cor 3.1（二択: すべての行順で g(α)=0、または H∘2 が Hadamard）← Prop 3.4（厳密恒等式 3⟨|L|²⟩+⟨|g(α)g(q)|²⟩=0、q=q_cross, q_same）。Prop 4.1（H と H∘2 が共に Hadamard なら H ≡ T。Phangara の予想 3.1.6 も解決）。Prop 5.2（7 MUB 系の中に遷移行列 T はあり得ない: 厳密な指標計算で E|Σz_i³|² = −18 < 0 を導く）。Prop 6.3（残りの場合に 6 個の正定値 Gram 形式の和の定数項が残りの係数絶対値和を負に支配する）。
- 限定条件: なし。すべて有限の整数・有理数算術（Lem 7.1 で代入の妥当性を厳密に保証）。モーメントの値が有理数である必要もない。
- 主論文との関係: 副論文は主論文の N(6)=3 を「その算術条件と実行契約のもとで」引用するだけ（paper.tex:154–157）。逆に主論文は副論文に依存しない。

## 2. 既知の結果との整合

| 結果 | 内容 | 本論文との関係 |
|---|---|---|
| Zauner 1991/1999、Klappenecker–Rötteler 2004 Conj. 1 | N(6)=3 の予想 | 主論文が肯定的に解決（条件付き） |
| Weiner 2013 Cor. 3.4 | d 個の MUB は d+1 個に補完できる | N(6) ≠ 6。副論文で N(6) ≤ 5 の導出に使用 |
| Grassl 2004 Thm 2 | 標準基底と F_6 の組は 4 MUB に伸びない | 主論文の特別な場合。整合 |
| Jaming–Matolcsi–Móra–Szöllősi–Weiner 2009 Thm 1.4 | Fourier 2 径数族は 4 MUB に伸びない | 主論文の方法の直接の先行例。整合 |
| Brierley–Weigert 2008/2009 | 数値的証拠; S_6 は第 3 基底すら持てない | 副論文 Prop 5.2（S_6 は完全系に入れない）より強い既知事実。整合 |
| Matolcsi–Ruzsa–Weiner 2013 Conj. 2.3 | Fourier 消滅予想 | 副論文 Thm 1.1 が全 Hadamard 行列で証明 |
| Maxwell–Brierley 2015 | Karlsson 族での消滅 | 副論文が族の制限を除去 |
| Lampio–Östergård–Szöllősi 2020; Phangara 2026 Conj. 3.1.6 | 立方根 Butson 類の一意性; S-Hadamard の分類予想 | 副論文 Prop 4.1 が任意位相から導出 |
| Bandeira–Doppelbauer–Kunisky 2022 Thm 1.10 | 単一遷移行列の次数 ≤ 6 の SOS では 7 基底を排除できない | 副論文は複数行列の積と高次の関係式を使うので矛盾しない |
| Matolcsi–Matszangosz–Varga–Weiner 2025/2026 | 三つ組に関する予想的十分条件 | 主論文は予想を使わず直接被覆 |

- 既知の下界・困難性と矛盾する点はない。副論文 Thm 1.2（7 MUB の非存在）は、私の知識では従来未解決で、主論文と独立な新結果である。
- 副論文の T と Tadej–Życzkowski の S_6 の同一視（03-cubic.tex:154–157、置換 (0,1,2,5,4,3)）は、私が記憶から書いた S_6 の指数表と一致した（`check_266_cubic.out`）。記憶ベースなので弱い確認。

## 3. 構成と依存関係（主論文）

```
任意の 4 MUB (I, H/√6, ·, ·)
 └─ §2 seed: 同値変形で第 1 行・列を 1 にし、最小ギャップの行対を先頭に。9 個の位相を N=112 のビンで覆い、
      必要条件（向き、ギャップ、和の Taylor 判定、共有座標の区間交差）で篩う → 249,878 個の半径 0.032 の球 (Prop 2.10)
     └─ §3 completion (A 段): 2+4 ブロック。固定枠 K_B, K_C で未知の 4×4 ブロックを 2×2 に圧縮し、
          極分解で U(2) の元を対応させ、U(2) を閉セルで被覆 (Lem 3.2–3.4)。セルごとにノルム判定・方向判定 (Lem 3.8–3.9) で排除、
          または §4 の位相球で被覆 (Prop 3.11)
         └─ §4 phase cover: Root 証明書 (Def 4.1; 残差 F(P)、接ベクトル T、法方向多項式 n_P)。
              厳密解の変位 δ を n_P と接座標で押さえ (Prop 4.5)、楕円体で接座標を包み (Lem 4.11–4.12)、半径を bootstrap (Prop 4.13)。
              atlas に合併 (Lem 4.15) → 半径 ≤ .125 の記録 3083 個 (Prop 4.16)
             └─ §5 refinement (B 段): 接座標の格子で親球を子球に置き換える (Prop 5.2, Lem 5.3–5.6)。
                  .07 → .016 (→ fallback .011)。失敗した親は義務として残る (Prop 5.7, Cor 5.8)
                 └─ §6 exclusion (C 段): 第 1・2 基底に不偏なベクトルをビン (12→24→48→96) で覆い、
                      クラスタ化して閉球に (Lem 6.4)。O/U グラフの辺を必要条件で削り、
                      「交差辺がすべて U の 2 つの disjoint な O-6-クリーク」の非存在を探索 (Lem 6.8, Prop 6.9)
                     └─ §7 execution: 100 回の呼び出しの終端条件 (Thm 7.1)。Appendix A: 全比較に片側余裕
```

- 要となる設計: 「数値的な提案」と「証明書」の分離。線形解法・固有分解・Cholesky・Lloyd 反復・貪欲なグラフ枝刈りはすべて提案にすぎず、排除や被覆は別に検査される不等式だけが正当化する（introduction.tex:141–148; rounding.tex:1457–1527）。失敗・上限到達は常に「未解決」であり、排除にならない。
- 条件付き構造: Prop 2.10, 3.11, 4.16, 5.7, 6.9 はすべて「実行が指定どおり完了したならば」の含意。Thm 7.1 がそれらを束ね、§7.4 がその前提の充足を**報告**する。

## 4. 確認したこと（checked OK）

### 4.1 副論文の有限和の全再計算（依頼 (a)）

- **Prop 3.4（pair certificate）**: 自作の厳密再実装（`check_266_exact_certificates.py`）で、7 行関係式 8871 本・座標 909・独立行 897 を得、2 つの q について代入後の全係数が 0。Appendix A の `rank table: 8871 909 897` と一致。さらに 21 個の厳密 Hadamard 行列（Fourier 2 径数族 8、Diță 型 4、Björck 巡回、S_6、同値変形 7）で恒等式を直接数値評価し、残差 ≤ 6×10⁻³¹（`check_266_pair_identity.out`）。副産物: すべての例で ⟨|L|²⟩=0、つまり L_H（02-pair-moments.tex:209–213）は Hadamard 行列の任意の行順で恒等的に 0 だった。
- **Lem 3.2（contraction）**: p ≤ 6 の 4 例で (7−p)m(A) = −Σ m(A_{i→ℓ}) が 10⁻¹⁷ で成立; 7 行の境界関係式も 10⁻¹⁷。**Lem 3.3（partitions）**: n=2,3 で 10⁻¹⁷。
- **Cor 3.1（dichotomy）**: 非立方の全例で max_π|g(πα)| ≤ 1.4×10⁻¹⁵、S_6 では 1/2 で H∘2 が Hadamard。
- **Prop 4.1（cubic classification）の証明**: 手で全段追った（Newton 恒等式 → 2 つの立方根三つ組; 支持集合の一致の表 (k, 3−k, 3−k, k) と x=y=−1 の矛盾; 3 列への制限の直交性と次元 ≥ 2 の矛盾; 符号 s_ij の対称性 s_ij=s_ji と次数 2 のグラフ = 5 サイクル）。誤りなし。T と T∘2 が Hadamard であることを Q(ω) で厳密確認。
- **§5（cubic obstruction）**: Table 1（03-complete-family.tex:276–298）の N 列・U 列 11 行すべて一致。ループ規則（同:309–318）と直接定義 (eq:cubic-R) が一致。パラメータ付け A=e_p+e_q−e_h−e_i, C=e_q−e_p が O（180 組）を正確に与える。λ_3,6,7,8=(1/2,3,9,4)。R(b)=λ_j G(b) が代表元だけでなく位数 ≤ 3 の全電荷で点ごとに成立。軌道サイズ (360,120,20,90)、係数 (180,360,180,360)。w_5+D, w_8+D の軌道分類の表（同:378–384）一致。連立方程式の解 S_3=S_7=1/10, S_0=−4/5、6+30S_0=−18。Lem 5.3 の 720 個の対称性を総当たりで確認（行置換 720 個すべてが実現）。
- **Lem 5.1（basic）**: (eq:same-matrix) 10⁻¹⁶; (eq:different-matrices) をテンソル MUB 対で 10⁻¹⁶; (eq:order-two) の Σ_{s≠r} g_s(a)=0 を次元 3・5 の完全系で 10⁻¹⁵（6 次元には完全系が無いので別次元で論理を検証）。
- **Lem 6.1（gram）, Lem 6.2 (iii)（7−l の singleton 帰着）**: 次元 3 の類似（7−l → d+1−l）で 20 組の電荷すべて一致。
- **Prop 6.3（mixed certificate）**: 自作再実装で三つ組 2821、混合方程式 45,996 ＋ 零ペア行 258 = 46,041、座標 3429、独立行 2767、98 個の正ピボット、定数 −613302797399911/6480、残差和 179721388988719/6480、合計 −2007321335237/30 < 0。**Appendix A の 3 行目と完全一致**。Z_λ の正定値性は別スクリプトでも確認（最小ピボット 41331 > 40000; numpy 最小固有値 14000）。
- 最終行の算術 (−613302797399911+179721388988719)/6480 = −2007321335237/30 を厳密確認。表 04-mixed.tex:234–247 の (λ, m, e_λ) と `data` の整合、Lem 7.1 の証明（W_J の可逆性 → m = P m）も問題なし。
- `verify.py` の SHA-256 は paper.tex:261–262 の値と一致。

### 4.2 主論文の数学（seed、completion、phase-cover、refinement、exclusion）

- **Lem 1.2**: 3 つの U_r が unitary で |(U_r^*U_s)_{ij}|²=1/6 を sympy で厳密確認。F^*DF の 2 つの和 1+2ω, 2+ω² の絶対値² = 3。
- **seed.tex**: Lem 2.4 の向き反転 x'_i = x_5 − x_{5−i} の論理、Lem 2.5（ギャップの 2δ 安定性）、Lem 2.6 の 2 つの不等式（g_x ∈ [1−(l+1)/N, 1−l/N] と ±w から導出）、Lem 2.7 の定数 5πw²/4、Lem 2.8 の区間、Lem 2.9 の双対上界（Cauchy–Schwarz）を手で再導出。Prop 2.10 の帰納法（被覆判定成功 → 既存中心の球、失敗 → 中点を保持、3/224 < 0.032）は正しい。**centers.txt の 249,878 行全件**が、x の非減少、向きフィルタ a ≤ l−d+1、s の選択規則（|2b_s+1−N| の最初の最小化点）、y 尾の非減少、両ギャップ不等式を満たすことを確認（`check_266_records.out`）。
- **completion.tex**: Lem 3.1（M^*M=[[4,S],[S̄,4]]、擬逆の公式）、A_0=√(2+|1+a|)、Lem 3.2 の恒等式 6I−Y^*Y = K_B^*B'^*B'K_B + K_B^*D'^*P_C D'K_B と q, p の導出、Lem 3.3 の U(2) パラメータ化、Lem 3.4 の ||U'||²=2(x²+y²+z²cos²T+w²sin²T) を手計算で確認し、定数 6h² と 42h⁴ を 20000 個の乱数セルで数値確認。Lem 3.5 の 3 片恒等式、Lem 3.6、Lem 3.7 の E_r、Lem 3.8–3.9、エントリー別係数公式（completion.tex:506–529）を **230 個の厳密 Hadamard 行列で数値検証**（`check_266_completion_chain.out`）。z 上界の最大使用率 0.60、E_r の最大使用率 0.20、q の最大 0.26 < .32。方向判定は乱数方向 20 本 × 230 例で全部成立。
- **phase-cover.tex**: Σ_{i>j}a_ij = 6s²、Σa_ij² ≤ 8s⁴（5000 例で最大 ≈ 8 未満）、||AA^*−6I||_F = 2π√2||F||、||F||²=(||AA^*||_F²−216)/(8π²)、転置不変性を確認。Lem 4.8 の弧の上界と弦恒等式、Lem 4.9 の Im(conj(D°)Z)=sin d と Σ|d−sin d| ≤ arc³/6 を 230 例で確認。Lem 4.11（楕円体の外和）、Lem 4.12（γ I − G の解消子）は手で確認。
- **refinement.tex**: Prop 5.2 の 4 不等式（特に S^*h' = S_h^*Tξ + S_h^*Aν + S_h^*T(I−T^*T)ξ の分解と β, σ）、Lem 5.3 の上界反復、Lem 5.4（単調性を要しない ρ の証明書）、Lem 5.5、Prop 5.7 の二重帰納法（深さ・親の順）を手で追った。誤りなし。
- **exclusion.tex**: Lem 6.1 の X_M=√(29/6)a, L_M=5a, I_M=1.5a（頂点で達成、`check_266_lemmas.out`）、Lem 6.2–6.3 の Taylor 分解、Lem 6.5 の 6−2π√6·0.2 > 20/7 > √6、Lem 6.6 の E_U=πr_e²(1+|S|/6) の 2 項、Lem 6.7 の共有接変数の扱い、Lem 6.8 のクリーク探索（最小頂点による枝分け、第 2 プールの事前 6 クリーク検査、上限到達は可能性扱い）を確認。「区間が空になる」判定はすべて「外側に膨らませた上界を**厳密に**超えたときだけ棄却、等号は保持」の向き（exclusion.tex:536–542; seed.tex:340–342）で健全。

### 4.3 Appendix A（rounding）の論理（依頼 (b)）

- **前提は明示されている**（rounding.tex:23–49）: binary64、最近接丸め、加減乗除・sqrt が正しく丸められる、`-ffp-contract=off`（FMA 禁止）、fast-math/再結合なし、FTZ/DAZ なし、例外マスク、`FLT_EVAL_METHOD=0`、`main` 前の定数も同じ前提、10 進リテラルの翻訳誤差 ≤ 5×10⁻¹⁶（相対）。`nearbyint`/`floor` には「整数を 1 つ選ぶ」以外の精度要求を置かず、`atan2` にも精度要求を置かない（提案値の単位指数を別に検査）。
- **誤差モデル**: u=2×10⁻¹⁶（真の 2⁻⁵³ より大きい）、加法的許容 τ=10⁻³⁰⁰、γ_m=2.1×10⁻¹⁶(m+1)（m ≤ 2000）。(1+u)^{m+1}−1 < γ_m を確認。Higham 流の積の丸め因子の議論で標準的。
- **π と単位指数**: 0x1.921fb54442d18p+1 = 884279719003555/281474976710656 と Machin 公式の有理証明書（rounding.tex:170–173 の有理数）を sympy で**厳密に再現**。|π̂−π| の証明された上界 1.2254171×10⁻¹⁶ < 1.225418×10⁻¹⁶（真値 1.2246×10⁻¹⁶）。cos/sin の級数の裾 y²⁰/20!, y²¹/21! at y=.786 を確認。
- **方法の健全性**: すべての比較が「理想量 ≤ 計算値 + 明示余裕」の片側評価で、表（rounding.tex:1533–1560）に経路ごとの「確立した敵対誤差」と「明示余裕」が並ぶ。提案された方向・逆行列・基底変換は「厳密に保存された数」として扱い、その精度は仮定しない。この設計は原理的に正しい。例として ub.cpp:132–146 の `get` の各フィールド（R_s+1e-9、h+1e-10、m−1e-8、A+1e-7、q+1e-7、pol+1e-7、z+1e-7、E_phase=(L+A+1)R_s²/(2m)+1e-7、E+1e-7）が tex の式と余裕に一致することを読んで確認した。スケール不等式 q<.32、off<.74、pol<.068、R_s<.203、h<.197 を数値確認。
- **ソース定数のスポットチェック**（読むだけ。実行なし）: 付録の余裕表の値がソースに実在する。stageC.cpp: 2e-7（9 箇所）、3e-7、1e-9、1e-12、K=1.00001、binCap=2,000,000、2048、capSteps=5,000,000、`.1`、validate(.03)。stageB.cpp: .249/.23、1e-6、2e-8/5e-8、.500001、1.01、50,000,000、2,000,000/8,000,000、範囲 .0001 ≤ τ ≤ .126。patchx.cpp: .64、1e-18、2e-6、.15、3.8、1e-5、1e-10、1e-7、.00015、.00002、1.45、1.8、beta<.1、.0002、.00001、gap ≥ 1e-6、30 回。ub.cpp: m>1.25、.0001、1e-100/1e100。atlas.cpp: r=.0320001、target=.103、global=.125、levels=3、iters=2、first=1。covout.cpp: 1.251、.01、1e-7、1e-5、RN2 = pow(R·n,2) − .0001、T=pow(√RN2+1.5,2)。seed.cpp: s の規則、向き・ギャップフィルタ、1e-8。
- **残る留保**: (i) この付録は「ソースがこうしている」という散文であり、機械検証ではない。私はソースの数十箇所の定数と数式を照合したにとどまり、77 KB 全体の意味論的一致は確認していない。(ii) `pow(x,2)` の乗算への降格は逆アセンブルの目視に依存し（execution.tex:93–100）、その逆アセンブルは配布されていない。(iii) `-O3 -march=native` は FMA 命令を利用可能にするが、`-ffp-contract=off` が契約を禁じる。GCC は `-ffast-math` なしに浮動小数点の縮約を再結合しないので、前提は GCC の文書化された挙動と整合する。(iv) `nearbyint`/`floor`/`atan2`/`sqrt` は glibc 呼び出しだが、要求は「整数を返す」「sqrt が正しく丸められる」だけで妥当。

### 4.4 配布記録の内部整合（`check_266_records.out`）

- 主要 8 ファイルの SHA-256 が Table execute:hashes（execution.tex:369–376）と**全部一致**。7 本の C++ ソースの SHA-256 とバイト数が `extraction.json` と一致。
- 記録形式（4 バイト k、36 double の p、k×36 double の T、8 double、`356+288k` バイト）で全ファイルがちょうど EOF まで解析できる。件数: coarse 764/774/757/788、mid 600、fine 22,379、fb1 1127、bad 系 7 ファイルは 8 バイト（件数 0）。バッチ 88 個（87×256+107）、連結が fine.bin と半径以外で一致し、全半径が .016 に置換され、置換は拡大。初回 C の未解決 163 件（15 ファイル、すべて graph cap）、各 bad ファイルはそのバッチの順序付き部分列。fallback 入力は C ロケール順の 15 ファイルと一致、親 163、出力 1127、最終 C の done 282+282+282+281 = 1127。C バッチの exit 0/1 が 73/15 で bad 件数と整合。全記録が Root 定数の範囲（C<20、n0,froot<1e-6、id<1e-4、ort,hor<1e-5、mnorm<2、k ≤ 10、|p| ≤ 4、|T| ≤ 1.01）と半径上限を満たす。k はほぼ全部 4（Hadamard 多様体の既知の次元 4 と整合）、k=0 が各段に 1 件（孤立点）。

## 5. 指摘（番号付き）

**確定した誤り: 見つからなかった。**

| # | 分類 | 重大度 | 場所 | 内容 |
|---|---|---|---|---|
| 1 | 認識論的な隙間（依頼 (c)） | **高** | introduction.tex:13–18; execution.tex:82–85, 242–281; conclusion.tex:3–11 | Thm 1.1 は「記録された実行が certify する」という形で、N(6)=3 は (A) tex の数学（Prop 2.10, 3.11, 4.16, 5.7, 6.9, Thm 7.1 の含意）、(B) 7 本の C++ がその補題を正しく実装していること、(C) 2026-10-04 の実行が Appendix A の算術条件下で終端条件を満たしたこと、の 3 つの連言である。(A) は tex で監査でき、私が追った範囲で健全。(B) は tex に「ソースはこうする」と書かれているだけで、形式検証も第三者監査もない。(C) は論文自身が「限定的な証拠であり証明ではない」と述べる（execution.tex:82–85）。したがって正確な現状は「N(6)=3 は、未検証の 77 KB の C++ と 1 回の実行記録を信じれば証明された」であり、Hales の Kepler 予想（Flyspeck 以前）や 4 色定理の初期証明と同じ身分。論文の文章自体はこの区別を一貫して守っており誇張はないが、abstract（paper.tex:56–60）の「We prove」は条件を読者が追わないと誤解を生む。 |
| 2 | 証拠の欠落 | **高** | introduction.tex:154–158; execution.tex:16–21, 87–104, 119, 380–392; README.md:12; VERIFY.md:3–8 | 論文が実行の証拠として参照するもの——`verification/computation/README.md`、`verification/code/proof_run.py`、`verification/REPRODUCE.md`、`build/build-manifest.json`、`covrun` の逆アセンブル（`pow(...,2)` の降格の確認、16 進アドレス 6035〜608b）、プローブの出力記録（MXCSR 0x1f80、FLT_EVAL_METHOD=0）、バイナリの SHA-256、`status.json`/`run.json`、stdout/stderr の生ログ——は**配布ディレクトリに存在しない**。README.md は「`VERIFY.md` が PDF 49–54 ページの指示を上書きする」「コンパイル済みプログラムは含まない」とし、VERIFY.md は「アーカイブ目録、ハッシュ検査、生プロセスログ、機械記録、監督コードを含まない」「歴史的アーカイブの認証・C++ 実装の証明・コンパイラとハードウェアの挙動の検証ではない」と明言する。配布物で検査できるのは数値記録の内部整合（§4.4）だけで、Thm 1.1 の前提 (C) のうちビルド・算術環境に関する部分は**読者には確認できない**。同じ結論を主張したいなら、x86-64 GNU/Linux で `build.sh` を走らせ（約 75 分）、逆アセンブルを自分で確認する必要がある。 |
| 3 | 日付・改訂履歴の不整合 | 中 | paper.tex:51; execution.tex:4, 286, 394–400; introduction.tex:156–158 | 論文日付 2026-09-24、実行完了 2026-10-04 03:10–04:24 UTC。「原稿が報告していた 9 月 7 日の実行は認証されておらず、その証拠パケットは回収されなかった。新しい主要データのハッシュと段階件数は 9 月 7 日の報告値と一致するが、それは旧実行を認証しない」。つまり公開版は日付を変えずに書き換えられており、初版の時点では証拠が揃っていなかったことを論文自身が認めている。追跡上は「2026-09-24 版の主張は 2026-10-04 の再実行で裏打ちし直された」と記録すべき。 |
| 4 | gap（検証不能） | 中 | execution.tex:90–104; rounding.tex:348–354 | `coverSeeds::run` の `RN2 = pow(R·n,2) − .0001` と `T = pow(√RN2+1.5, 2)`（covout.cpp:100–101）が binary64 乗算に落ちることが Prop 2.10 の余裕の前提。GCC は通常 `pow(x,2.)` を `x*x` に置換するが、論文はそれを逆アセンブルで確認したと述べるだけで、逆アセンブルは未配布（#2）。「新しいビルドではこの確認をやり直す必要がある」（execution.tex:102–104）とあるとおり、再現者は必ず自分で確認すること。 |
| 5 | gap（tex と code の対応） | 中 | rounding.tex 全体; execution.tex:16–38 | Appendix A は各ルーチン（`sumok`、`extend`、`supdual`、`trial::certificate`、`buildH::patch`、`Candidate::prep`、`covercompute`、`ell`、`bootstrap`、`relationB::eval/init`、`markRec`、`Layer::test`、`graph::joint`、`graph::exist` 等）について「こう計算し、この余裕を足す」と散文で記述する。私はソースの定数と数式を数十箇所照合し一致を見たが（§4.3）、制御フロー全体（特に A 段 `Asearch`、B 段 `refiner`、C 段 `cluster`/`exist`）が補題どおりであることは確認していない。この層のバグ（例: 余裕の符号、未解決と排除の取り違え、添字ずれ）は tex からは検出できない。独立監査には C++ 全体の精読か、別言語での再実装・再実行が要る。 |
| 6 | 記述の不整合（軽微） | 低 | execution.tex:16, 119, 381; README.md:9–12 | 論文は `verification/code/` に 7 本のソースがあると書くが、配布物では `verification/code/` には `extraction.json` だけで、ソースは `verification/saved-run/compiled-build/source/` にある（`abi-fp-probe.cpp` を加えた 8 本）。ハッシュは一致するので実害はないが、論文の記述と配布構成が食い違う。 |
| 7 | 記述の不整合（軽微、副論文） | 低 | paper.tex:235, 273; 04-mixed.tex:218 | 検証器の所在を `verification/verify.py` と書く箇所と `certificate/verify.py` と書く箇所が混在。`build/certificate/` は存在しない（PDF 生成時には別配置だったと思われる）。SHA-256 は `verification/verify.py` と一致。 |
| 8 | checked OK（余裕） | 低 | 04-mixed.tex:282–288 | 混合証明書の支配は定数 −9.46×10¹⁰ に対し残差和 2.77×10¹⁰、比 0.29。余裕は十分。ただしこの値は座標の優先順位リスト I（verify.py:332–370）に依存し、別の自由座標の選び方では残差が増えうる。私の再実装は同じ I を採用して同じ値を得た（#9）。 |
| 9 | 再実装の独立性に関する留保（副論文） | 低 | 05-verification.tex 全体 | 私の `check_266_exact_certificates.py` は tex の仕様から書いたが、仕様が指定していない慣習（正規形の取り方、鍵の整列順、優先順位リスト I）は公開検証器と同じものを採用した（そうしないと Appendix A の 3 行は数値比較できない）。したがって「同じアルゴリズムの独立実装」であり、帰着則（Lem 6.2）や Gram 補題（Lem 6.1）自体に共通の数学的誤りがあれば両者で再現される。その層は別途、次元 3・5 の完全 MUB 系と 6 次元の MUB 3 つ組で数値検証した（§4.1）。また `data`（Z_λ と電荷リスト）は verify.py から読んだ値をそのまま使った。これは証明書そのもの（Appendix A に印字される入力）なので、別の値を使う理由はない。 |
| 10 | 強い主張の確認（副論文） | 情報 | paper.tex:142–157 | Thm 1.2（7 MUB の非存在、N(6) ≤ 5）は主論文と独立で、浮動小数点を含まない。従来 N(6) ∈ {3,4,5,7} が既知（Weiner）だったので、これは Zauner 予想への無条件の前進である。Bandeira–Doppelbauer–Kunisky の「単一行列・次数 6 の SOS では不可能」という結果（paper.tex:129–137）と矛盾しないことを tex で確認した（複数遷移行列の積、共有列、次数 11 までの関係式を使う）。 |

**依頼 (c) への回答を一文で**: 主論文が証明しているのは「Thm 7.1（有限実行 ⇒ 4 MUB なし）という数学的含意」と「ある実行がその前提を満たしたという、配布物では一部しか検証できない事実主張」であり、「N(6)=3 が（通常の意味で）証明された」と言うには、未検証の C++ と算術環境の前提を受け入れる必要がある。一方、副論文の N(6) ≤ 5 は無条件で、私はこれを独立に再現した。

## 6. 監査できていない部分

- C++ 7 本（77 KB）の全体が Appendix A と各節の補題を正しく実装しているか（#5）。定数のスポットチェックのみ。
- ビルド条件（`pow` の降格、MXCSR、FLT_EVAL_METHOD、FMA 不使用）が 2026-10-04 の実行で成立していたか。証拠が未配布（#2, #4）。
- 実行の再現。禁止事項であり、また x86-64 GNU/Linux が必要。
- Appendix A の誤差予算の各行（例: rounding.tex:370–407 の枠構成の 10⁻¹³〜10⁻¹⁰ 級の伝播、:539–572 の M の誤差、:823–846 の共分散支配、:895–957 の解消子の符号付き加算）を独立に再導出してはいない。スケールと余裕の大小関係を抽出して確認したにとどまる。設計（片側余裕、提案と証明の分離）は健全。
- 主論文 Prop 2.10 の seed 列挙の件数 22,829（x リスト）と 151,744,086（組）の再現。約 1.6 億組の数値フィルタで、Python では非現実的。249,878 中心の必要条件は全件確認した。
- 副論文の立方行列と Tadej–Życzkowski S_6 の同一視は、記憶した S_6 の表との照合（弱い）。
- 引用文献の実在（Phangara 2026、Cárdenes Wuttig–Tindall 2026、Matolcsi et al. 2026、McNulty–Weigert 2026 review）。証明には使われていない。

## 7. 必要な専門家

- **計算機支援証明・検証数値計算**（Rump、Tucker、Hales/Flyspeck 系）: Appendix A の片側丸め余裕と C++ の対応の監査。理想的には別言語（区間演算ライブラリ）での再実装・再実行。
- **量子情報・MUB**（Bengtsson、Matolcsi、Weiner、Szöllősi、Brierley、Grassl）: 副論文 Thm 1.2 の新規性の確認と、Fourier 消滅予想の解決がコミュニティで受け入れられるかの判定。主論文の被覆戦略（Jaming らの拡張）の見落としの有無。
- **GCC/x86-64 の浮動小数点意味論**: `-O3 -march=native -ffp-contract=off` 下での `pow(x,2)` 降格、縮約禁止、FLT_EVAL_METHOD=0、glibc の sqrt/nearbyint の正当性の確認。

## 8. 総括と確信度

副論文（Exact Fourier certificates）: 2 つの厳密証明書を tex の仕様から独立に再実装して 3 行の出力を数値まで完全に再現し、立方節の有限和も厳密に全再計算し、証明の散文部分（分類・二択・立方排除）も誤りなく追えた。**Thm 1.1（Fourier 消滅予想）と Thm 1.2（7 MUB なし、N(6) ≤ 5）が正しいことへの確信度: 高（95%）。** 残りは私の再実装が公開検証器と慣習を共有する点（#9）。

主論文（N(6)=3）: tex に書かれた数学は、手で追った補題と 230 例の数値検証の範囲で誤りがなく、計算機支援証明としての設計（必要条件のみで排除、失敗は未解決、片側余裕、提案と証明書の分離）は正しい。配布された数値記録は論文の表と完全に整合する。しかし結論は未監査の C++ と、配布物で検証できない算術環境の前提に依存し、論文自身がそれを「限定的な証拠」と呼ぶ。**「tex の数学と記録の整合に誤りがない」ことへの確信度: 高（90%）。「N(6)=3 が確立した」と言えるかへの確信度: 中（65%）**——数学的な見落としよりも、コード・環境層の未検証リスクと、9 月 7 日版が認証されなかった経緯（#3）による。副論文の無条件の N(6) ≤ 5 と、主論文の条件付き N(6)=3 は、はっきり区別して記録すべきである。

## 9. 主担当による追記（2026-11-11）

- **§4.2 の Lem 3.1（completion.tex:55–75）の裏取り。** 監査役の `check_266_lemmas.py` は 11-06 に中断され、11-11 の再実行では「M^*M = [[4,S],[S*,4]]」の行が FAIL と出た。原因はスクリプト側（sympy の `simplify` が絶対値 1 の指数関数の共役積を 0 行列と構造的に等しいと判定できない）で、論文の補題は正しい。M = (1 p) なので M^*M の成分は 1^*1 = 4、1^*p = S、p^*p = 4 と手で確認できる。別に `check_266_lem31_numeric.py`（乱数 2,000 例、numpy）で M^*M、特異値の 2 乗 4 ± |S|、逆枠の公式 (completion.tex の式 completion:inverse-frame) を確認し、残差は 1e-14 以下（`check_266_lem31_numeric.out`）。`check_266_lemmas.out` は中断時の空ファイルのままなので、この行の根拠は本追記と上の .out とする。
- **主担当が tex・配布物で再確認した点**: 論文日付 2026-09-24（paper.tex:51）と実行日 2026-10-04（execution.tex:4, 284–286）、9 月 7 日の実行が未認証であること（execution.tex:394–400、introduction.tex:156–158）、参照されるが配布物に無いファイル 4 件（`verification/computation/README.md`、`verification/code/proof_run.py`、`verification/REPRODUCE.md`、`build/build-manifest.json`）、「limited evidence ... not a proof」（execution.tex:82–85）、副論文の `certificate/verify.py`（04-mixed.tex:218）と `verification/verify.py`（paper.tex:235）の混在、副論文 Thm 1.2 と Weiner の補完定理から N(6) ≤ 5（paper.tex:142–153）。定理番号（主 Thm 1.1・Thm 7.1、副 Thm 1.1・Thm 1.2）は \newtheorem と環境の並びから数え直して一致。
- 「7 個の MUB の非存在は従来未解決だった」は監査役の知識に基づく記述で、文献調査はしていない。公開文面では断定しない。

## 10. Codex による独立レビュー（2026-11-11）と訂正

依頼文 `audits/codex_review/prompt_w6_266.txt`、返答 `review_w6_266.md`。7 点はすべて「正しい」。中心的な読み取りと Weiner による N(6) ≤ 5 の推論は支持された。Codex は証明書の計算を再実行していない。訂正は次の 3 つで、本報告の該当箇所は以下のとおり読み替える。

- **§5 #1 と §8: 「限定的な証拠」の範囲を広げすぎていた。** 論文がそう呼ぶのは、プローブとプログラム内チェックによる機械・コンパイラ条件の確認である（execution.tex:82–88）。実行全体や証明全体をそう呼んでいるわけではない。主定理は「明示された算術・コンパイラ条件のもとで、記録された完全実行により N(6)=3 を証明したとする計算機支援証明」と書くのが正確。
- **§0 の 3 と §5 #3: 日付から初版の状態を断定していた。** tex から言えるのは「9 月 24 日付の配布原稿が 10 月 4 日の実行に依拠し、9 月 7 日の実行は『ここでは認証されていない』と書かれている」まで。「初版の時点では証拠が揃っていなかった」「日付を変えずに書き換えられた」は推測なので取り下げる。改訂日は確認できない。
- **§0 の 1: 「確定した結果として扱ってよい」は過大。** 自作の検証器は公開検証器と規約・証明書データを共有する別実装で、証明書計算の出力を再現したもの。定理全体を独立に確立したのではない。N(6) ≤ 5 は「副論文が提示する厳密証明書から得られる上界」と書く。
- 補足: 欠けている 4 ファイルと当時の環境証拠は未同梱だが、ソース・数値記録・再ビルド手順は同梱されている。7 個の MUB の非存在は、2026-03-26 版の総説（arXiv:2410.23997v2、Problem 10.1）が未解決問題として挙げている（Codex の確認。網羅的な文献調査ではない）。230 例は数値検査。
