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
