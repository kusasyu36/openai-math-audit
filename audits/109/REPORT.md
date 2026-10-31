# 監査レポート 109: Integer multiplication below n log n（OpenAI, 2026-09-23）

- 監査日: 2026-10-31（監査役: Claude サブエージェント。報告書は主担当が代筆保存）
- 対象: `data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/`（main.tex + sections 00–11 + figures 2 本、計 5,324 行）。数式はすべて TeX 原文で確認した。
- 自作スクリプト（`audits/109/`、出力は同名 `.out`）:
  - `check_counts.py`: 有限ネットワークの定数、η、τ、指数の連鎖、κ の導出、誤差の連鎖（厳密有理数・sympy）。
  - `check_network.py`: 小さい h での scalar network の 8 行 x 3 段のシミュレーション、ラベル部分空間の入れ子と残差次元（Q 上と F_2 上）、次元総和の記号計算。
  - `check_phase.py`: 位相フレーム C_U と核分解、Z 補正、H_0^{⊗D} の恒等式（sympy 厳密）。
  - `check_packed.py`: 8 操作の恒等式と packed 加算の S_0 / T_0 / 例外集合 B（総当たり＋乱択）。
  - `check_resampling.py`: Gaussian 恒等式、Neumann 評価、chirp、置換の相殺、合成変換、CRT 写像（mpmath 40 桁）。
  - `check_shear.py`: 下三角分解、shear の 3 ステップ、再帰式、chunk-swap の行数条件。
- 結果: 全スクリプトが OK（`check_network.py` は h=6 で 1 件 NG が出たが、論文が h>6 を使う箇所の小 h 特有の現象で、h=7 では OK。`.out` に注記）。
- **リポジトリのコード・Lean は実行していない。** 外部文献の確認には HvdH 2025（arXiv:2503.22848）と Adler–Harvey–Jain–Kleinberg–Rasala Lehman 2006（SODA 版）の PDF 本文を取得して読んだ。AFKL 2019 は要旨のみ。

## 1. 主張の正確な内容

**Thm main**（`00-introduction.tex:21–30`）
- 1 台の決定的 Turing 機械 A。固定有限アルファベット、固定有限本の 1 次元テープ。
- 入力 x#y（x, y ∈ {0,1}^n）。出力 bin_{2n}(val(x)·val(y))。**すべての n ≥ 1 とすべての入力対**で正しい。
- 最悪時間 T_A(n) = O(n (lg n)^{1−κ})、κ = 2^{−182}、lg n = max{⌈log_2 n⌉, 1}。
- 一様性: 機械は 1 台で、n ごとの advice はない。素数・形状パラメータは機械が計算する（`08-assembly.tex:172–191`）。
- 「すべての入力長で」: 巨大な閾値以下は筆算で処理し、定数に吸収する（`08-assembly.tex:1001–1003`）。閾値は天文学的である（§5, §6 参照）。

**付随する主張**
- Cor exact-arithmetic（`00-introduction.tex:194–206`）: 整数除算と平方根も同じ時間。Newton 反復（Brent）による。
- 転置（`10-transposition.tex`）: HvdH 2025 の帰着 Cor 6.2 と組み合わせ、n_1×n_2 の b ビット行列の転置が O(m (lg m)^{1−κ})、m = n_1 n_2 b。特に k×k の 2 値行列が O(k²(lg k)^{1−κ})。
- 計算モデルの注意: テープ本数は「固定」だが、ネットワークの役割ごとに 1 本を使うので W_b ≈ 1.8×10^{20} 本、層では W_c ≈ 1.9×10^{21} 本である（`04-swap.tex:215–217`「Reserve W role tapes」。`05-layers.tex:586` は「fixed number of tapes」とだけ書く）。数学的には有限で固定だが、非構成的な「定数」である。

## 2. 既知の結果との比較（最重要）

### (a) Afshani–Freksen–Kamma–Larsen 2019（ネットワーク符号化予想の下での Ω(n log n)）

- AFKL の主張: 無向 k-pairs 予想（Li–Li 2004、NCC）が真なら、2 つの n ビット整数を掛ける**定数次数のブール回路**は Ω(n lg n) ゲートを要する（arXiv:1902.10935 の要旨で確認。本文は未取得）。
- 本論文との接続は「oblivious か」で決まる。
  - head の移動が入力長だけで決まる（oblivious な）固定テープ TM が時間 T で動くなら、各ステップが O(1) 個の既知セルを読む O(1) サイズの部分回路になり、全体で**大きさ O(T)、fan-in・fan-out とも O(1) のブール回路**になる。
  - oblivious でない TM は Pippenger–Fischer の模倣で O(T log T) の回路しか得られず、n (lg n)^{2−κ} > n log n で矛盾しない。
- 本論文のアルゴリズムは、**読んだ範囲では head の移動が payload に依存しない**（機械的検査ではなく読みによる判断）。payload 依存が隠れうる場所は 3 つで、いずれも TeX で確認した:
  - 例外アドレスの修正と基数ソート: 「アドレスだけに依存する」と明記（`05-layers.tex:187–188`、`342–344`）。鍵は rank（アドレス）。
  - 桁上げ伝播と丸め: 値に依存するが線形走査で head の動きは固定（`08-assembly.tex:783–818`）。制御付きシフトのオフセットはアドレス欄から作る（`02-streams.tex:124–141`）。素数探索・形状・閾値は n だけの関数（`08-assembly.tex:63–99, 172–191`）。
  - 外部の HvdH O(N log N) 乗算器はブラックボックス（`06-transforms.tex:453–456`）。仮に oblivious でなくても、Schönhage–Strassen（oblivious）に置き換えると「packed polynomial products」行が log log 因子だけ増えるが、この行の冪は 1−ε で、論文自身の slack 議論（`08-assembly.tex:995–997`）で log の冪は p^κ に吸収される。
- したがって、**本論文が正しければ、AFKL の定理と合わせて NCC（無向 k-pairs 予想）が偽になる**。論文はこの帰結を一言も述べず、references.bib に AFKL はない。`01-history.tex:110–139` は下界の範囲を論じるが、回路モデルと NCC には触れない。
- HvdH 2025 §1.4（取得した本文で確認）は「ブール回路では NCC の下で Ω(n log n) [AFKL19]」「2 テープ TM の転置は NCC の下で Ω(n² log n) [AHJ+06]」と明記している。
- 分類: **誤りではない。強い警戒信号**（011 の EH 型の帰結と同種）。論文の主張の強さは「NCC の反証」と同等である。

### (b) Paterson–Fischer–Meyer 1974（オンライン乗算の Ω(N log N)）
- 対象は**オンライン**多テープ TM（入出力を 1 桁ずつ）。本論文はオフラインなので矛盾しない。回避している仮定は「オンライン性」である。

### (c) 転置: HvdH 2025、Stoß 1973、AHJ+06
- HvdH 2025 Thm 1.2: 転置が Ω(m lg m) なら乗算も Ω(m lg m)。本論文はその対偶として転置の o(k² log k) を得る（`10-transposition.tex`）。
- Stoß の Ω(n² log n) は「データを動かすだけで計算しない」機械向け。本論文は中間で XOR 結合を作るので、この制限の外にある（`01-history.tex:135–139`、論文自身が指摘）。論文の address interchange は、文字どおり「符号化（coding）が経路制御（routing）に勝つ」構成であり、NCC が否定する現象そのものである。
- **より局所的な衝突**: `04-swap.tex` Lemma chunk-swap を P=G=B=1 で使うと、2^u×2^u ビット行列の転置が O(2^{2u} u^τ) になる。この補題は数値ルーチンを使わず、head 移動はすべてアドレス決定なので oblivious 性に疑いがない。AHJ+06（SODA 版 §6, Thm 6.2–6.3、本文で確認）は NCC ⇒ Floyd の I/O モデルでの oblivious 転置 Ω(p log p) ⇒ oblivious **2 テープ** TM の転置 Ω(n² log n) を示す。本論文の機械は固定 k（≈10^{20}）テープなので直接は当たらない。AHJ+06 の I/O 模倣が固定 k テープへ拡張されるなら（高速メモリ kp で同じ評価になると期待されるが、**未確認**）、Lemma chunk-swap だけで NCC と矛盾する。この場合、監査対象は `03-motifs.tex`（bit network）+ `04-swap.tex` の約 1,100 行に局所化され、§5–8 とは独立になる。

### (d) Morgenstern 1973
- 有界係数の線形算法に対する DFT の Ω(n log n)。本論文はビット計算モデルなので適用外。同伴論文（OpenAIExactDFT2026、`01-history.tex:141–164`）は係数の大きさを無制限にして回避している。本論文の監査対象ではない。

### (e) Ailon 型の下界（2×2 ユニタリゲートのエントロピー）
- 本論文の butterfly 層は、ユニタリ核 C_v に加えて非ユニタリな shear ゲート（係数 ±1, ±1/2）を使い、中間値の増大を guard ビット（Δ = O(d^{20}) ビット、`05-layers.tex:695–824`）で吸収する。ユニタリゲートのみの下界には当たらない。矛盾ではない。

## 3. 構成と依存関係

**log 因子をどこで削るか**（`08-assembly.tex:822–867`、Table 1。コストは V = Tp = Θ(n) で割った p の冪）

| 行 | 正規化コスト | p の冪 | 1 との差 g_i |
|---|---|---|---|
| prefix slot の個別 butterfly | dK | ε(1+c) | > 1/2 |
| **変換レイアウトのチャンク交換** | p K^{τ−1} | 1 − εc(1−τ) | **2^{−181}** |
| 同時 butterfly 層 | ℓ d^{λ'} | 1 − ε(1−λ') | 2^{−129} |
| CRT と軸の並べ替え | d²(1+ℓ^τ) | τ + ε(2−τ) | > 2^{−51} |
| Gaussian 線写像 | d p^{1/2+δ} α | 3/4+δ+3ε/2 | > 1/8 |
| chirp・twist・スカラー積 | d p^δ | ε+δ | > 1/2 |
| packed 多項式積（HvdH 乗算器） | log(rp) | 1−ε | 2^{−75} |

- 素朴な FFT 乗算は正規化コスト ~p（= 変換のレベル数 log T ≈ b = p/6）。
- 2 つの新しい部品が log を削る。
  1. **アドレス欄の高速交換**（`04-swap.tex`）: u ビットのチャンク 2 つの交換を O(V u^τ)、τ = 1 − 2^{−50}。
  2. **同時 butterfly 層**（`05-layers.tex`）: D 本の軸に H_0 を同時に作用させて O(V d^{λ'})、λ' = 1 − 2^{−54}。
- どちらも §3 の**有限線形ネットワーク**（h=100、m=h³=10^6、W ≈ 10^{20}〜10^{21} 本の wire）から来る。役割（wire）ごとに配列を持たせ、配列の「フレーム」（表現）を辺ごとに変え、ゲートは点ごとに作用させる。フレーム変更のコスト（有理階数 / 残差次元の総和 s）が **s < Wm** を満たす（η_b = 339/22587335000000 ≈ 1.5×10^{−11}、η_c = 73/19906842167500 ≈ 3.7×10^{−12}）。これが F_k ≤ (s/W) F_{k−1} + O(1) という再帰を生み、F_k = O(m^{kτ}) になる。
- **κ = 2^{−182} の出所**: 変換レイアウト（`06-transforms.tex:254–269`）で D×q 個の K ビットチャンクを転置するのに Dq−1 回のチャンク交換を使い、各交換が O(V K^τ)。Dq ≈ p/K だから正規化コスト p K^{τ−1} = p^{1−εc(1−τ)}。ε = 2^{−75}、c = 2^{−56}、1−τ = 2^{−50} なので g_2 = 2^{−181}。κ = g_2/2 で、残り半分で log p の冪を吸収する（`08-assembly.tex:849–867, 995–999`）。`check_counts.py` で厳密に一致を確認した。
- **依存関係の鎖**: Lemma scalar-motifs（スカラー網） → Lemma motif-residuals（ラベルの入れ子・残差次元） → Prop bit-motif-interface / complex-motif-interface（s_b, s_c） → Prop power-interchange（再帰）→ Lemma chunk-swap（任意幅）→ Lemma packed-selected-bit-rectangle（packed XOR、chunk-swap を使う）→ Prop simultaneous-layer（層の再帰）→ Lemma synthetic-transform-cost → Prop synthetic-convolution → §8 の組み立て。Gaussian 再標本化（§7）と CRT（§8）は HvdH 2021 と Agarwal–Cooley の既知の恒等式で、置換を変換式に残す形にしている。
- 外部依存: HvdH 2021（Thm 4.2, Lemma 4.6, 4.8–4.12, 2.12, 5.1 の数値ルーチンと評価）、HvdH 2021 Thm 1.1 の乗算器（部分ルーチン）、Rosser–Schoenfeld、Bluestein、Nussbaumer。NCC・AFKL は引用なし。

## 4. 確認したこと（checked OK）

**定数と指数（`check_counts.py`）**
- v, N, I, z_b, z_c, W_b, W_c, L_b/N = 100/539, L_c/N = 101/539, s_b, s_c, η_b, η_c の表示値がすべて一致。隣接数の公式は小 h で総当たり確認。
- η_b, η_c > 20·2^{−50}（約 845 倍、206 倍の余裕）、m^{1−2^{−50}} > m(1−η) を厳密に確認。log m < 20。
- 真の指数 log_m(s/W) = 1 − 1.1×10^{−12}（bit）、1 − 2.7×10^{−13}（complex）。τ = σ = 1 − 2^{−50} は十分ゆるい。
- `05-layers.tex:695–715` の guard 定数: m^{10} < B < m^{11}、ν = 11、C_1 = 20、24W³+4s+4W+4 < E。
- `08-assembly.tex:30–48` の 9 本の不等式すべて。Table 1 の 7 つの冪を独立に導出して一致。min g_i = 2^{−181} = 2κ。
- 誤差の連鎖 `08-assembly.tex:663–739`: 最終誤差 < 448 b 2^{−b+γ}（記号計算で一致）、b ≥ 2^{24} で γ ≤ b/4 かつ 448 b 2^{−3b/4} < 1/2。桁の容量 c_j < 2^{3b}、‖w‖ ≤ 1/16。Neumann 条件 α⁴θ > 3db、(1−1/(2d))^d > 1/2、2.01e^{−πx/2} < 2^{−x}（x > 0.8）。
- L < N/2 が成り立つ最小の h は 39（bit）、40（complex）。h = 100 は余裕がある。

**有限ネットワーク（`check_network.py`）**
- Lemma scalar-motifs: h = 5, 6（bit, F_2）と h = 5（complex, Q(i)）で、任意の scratch 初期値に対し 3 段で (X, Y) ↦ (−Y, X) となり、全 scratch が復元される。
- Lemma motif-residuals の表（`03-motifs.tex:376–396`）: h = 5 で Q 上（形式 I − J/9）に実際の部分空間を作り、19 本の辺すべてで入れ子・非退化・残差次元が表と一致、段の境界でラベルが接続、終端が F と U_a^⊥ になることを確認。P_{U_a^⊥} + P_{U_a} = I、2P_{U_a} の階数 1、rank(P_V − P_U) = dim V − dim U。
- F_2 上（h = 5, 7）: 入れ子・残差次元・残差の非退化かつ非交代（正規直交基底の存在条件）を確認。
- 次元総和: 一般の h で Σ|Δdim| = Wm − 2N + 2L、L = Ich を記号計算で確認（X, Y wire: m−1、side: m、center: m+2h）。

**位相フレーム（`check_phase.py`）**
- `05-layers.tex:60–68` の恒等式すべて（C = HSH = aI+bX、H_0 = bSCS、C² = X、C^{−1} = −iZCZ など）。
- C_U = H^{⊗m} diag(i^{q_U}) H^{⊗m} を m = 3, 4 で厳密に構成し、U ⊂ V に対し C_V C_U^{−1} = Π(aI+bX_v)^{±1} を確認。重み 3 の残差ベクトルでは符号が −1 になる（論文の「ε ∈ {1, −1}」どおり）。
- X_a → Y_a の経路: C_{u^⊥} C_u^{−1} = supp(u) 上 C^{−1}、他は C。Z 補正 i^{wt(u)k} Z^{⊗k}(·)Z^{⊗k} で C^{⊗mk}（k = 1, 2 列）。q_{u^⊥} − q_u = wt(x) + 2[u·x] mod 4。
- H_0^{⊗D} = b^D S_D C^{⊗D} S_D（D = 1, 2, 3）、(C^{−1})^{⊗f} = (−i)^f Z^{⊗f} C^{⊗f} Z^{⊗f}、P_M X_v P_M^{−1} = X_{Mv}、C_U の成分が Z[i]/2^m にあること。

**packed 加算（`check_packed.py`）**
- 8 操作の恒等式 (eq:eight-packed-result) を全パリティ・c ∈ {0,1} で確認。最大変位は 2（論文の 4、8 より小さい）。最初の 4 操作が (00 01 10)、後の 4 操作が (01 10 11) の巡回。
- K = 9, f = 3 の長方形で、S_0 = T_0（B の外）、S_0(B) ⊂ B、逆方向の 8 本で S_0 を反転、|B|/M = 0.28 ≤ 80(f−1)2^{−K} = 0.3125。

**変換の恒等式（`check_resampling.py`、mpmath 40 桁）**
- T P_s F_s = P_t F_t S、‖N − I‖ < 2.01 e^{−πα²θ/2}、eq:no-sort-one-dimensional の 2^{2α²} 正規化、‖S'‖ < 3/4、‖J'‖ < 7/8、‖D'‖ < 1（3 組のパラメータ）。q_j の単調性。
- chirp 恒等式 (eq:permuted-chirp)、a の周期性、置換の相殺 (eq:source-permutation-cancel)、F_s(u*v) = S(F_s u)(F_s v)、F^+F = I/S。
- R_r 上: F_t^+ F_t^− = I/t、指標和、DIF 式 (eq:synthetic-dif)、再帰が出力 j をスロット rev(j) に置くこと、M F^+((F^−f)(F^−g)/r) = f*g/(rM)、twist x ↦ ζy、CRT の (eq:triangular-crt) と畳み込みの保存。

**shear と再帰（`check_shear.py`）**
- Lemma lower-lower を消去手順どおり実装し、40 個の乱択有理行列で A = E_1 Π E_2 を確認。pivot の 3 ステップ、Z/q^b 上の shear 列、3 ステップの交換、Φ_{M'}Φ_M = Φ_{M'+M}。
- 再帰 F_k ≤ (s_b/W_b)F_{k−1} + 1 で F_k/m^{kτ} が有界（a/m^τ = 0.999999999985）。chunk-swap の ρ, R ≥ W^k, ρ < e、幾何和、descriptor の評価。

**手で確認した点**
- `03-motifs.tex` のゲート順（Y-side, Y-center, X-side, X-center, Y-center, Y-side, X-center, X-side）が stage 2 の逆スケジュールでも物理 wire について同じになること。
- 各時刻のゲートに触れる全 wire のラベルが共通であること（t に依存するのは「one gate per target/source」の時刻だけ）。
- フレーム恒等式 (eq:common-frame-identity) の帰納法（辺・点ごとのゲート・終端）。
- 時間の帳簿: 各ノードの非再帰コストが k に依存しない定数（辺 ≈ 5W 本、辺ごとに 2(m²+2m) 回の packed 加算、三角変換 m 回の制御付きシフト）で O(V (eK)^τ) に収まること。parking の O(sV) も定数。
- 層の再帰の解法（`05-layers.tex:631–686`）: 内部和は a/m^λ < 1 で O(e^λ)、葉は e^σ d^{β(1−σ)}、群の数 (m−1)(1+⌊log_m d⌋) を log(2d) で吸収。
- 行の供給（`05-layers.tex:544–561`）: R_row = 2^{q_0 K} ≥ W^{k_0}、padding は 2 倍以内、深さ j で W^{k_0−j} で割り切れる。
- RAM 的な議論の混入は見つからなかった。表引きはなく、ソートはテープ上の基数ソート、descriptor の算術は記録幅 Θ(rp) に吸収されている。

## 5. 懸念事項（分類・重大度・TeX の場所）

1. **【強い警戒信号・最重要】NCC の反証と同等の主張**（§2(a), (c)）。oblivious な固定テープ TM の O(n (lg n)^{1−κ}) は、AFKL 2019 の定数次数ブール回路の Ω(n lg n)（NCC 条件付き）と両立しない。論文はこの帰結を論じず、AFKL を引用しない。場所: `00-introduction.tex:21–30`、`01-history.tex:110–139`。
   - oblivious 性の判断は読みによるもので、機械的検査ではない。payload 依存が隠れうる 3 箇所（HvdH ブラックボックス、桁上げ・丸めの走査、例外レコードのソート鍵）は §2(a) に挙げた。
   - 弱い環: AFKL の「定数次数」の定義（fan-in と fan-out の両方か）。oblivious TM 由来の回路は両方 O(1) なので、どちらの定義でも成立する見込み。
   - 局所版: Lemma chunk-swap（`04-swap.tex:285–389`）単体の転置 o(N log N) は、AHJ+06 の k テープ拡張が成り立てば NCC と矛盾する（未確認）。
2. **【gap・中】有限ネットワークの「鍵」が正しいかは、小 h の構造検証と次元計算でしか確かめていない。** s < Wm の数え方（`03-motifs.tex:351–457`、`492–526`、`568–669`）は構造的に整合し、小 h で部分空間を作っても矛盾は出なかった。しかし「フレームを変えながらゲートを点ごとに作用させた総コストが、各役割に C^{⊗e}（または Φ_I）を直接作用させるより少ない」という主張の正しさは、§2(a) の帰結と表裏である。構造的誤りがあるなら、「label at a gate is common to all touched wires」（`03-motifs.tex:329–342`）か「scratch role も終端作用素を受ける」（`03-motifs.tex:179–186`）のどちらかにある可能性が高いが、特定できなかった。
3. **【gap・中】`05-layers.tex:524–543` の層のノード overhead O(V((eK)^τ+1))。** 辺ごとの基底変換 P_M, P_M^{−1} は m(m−1)+3m ≈ 10^{12} 回の packed 加算で、各々 O((V/W)(L^τ+1))。辺の本数 ≈ 5W。合計は O(10^{13} V (eK/m)^τ) で確かに O(V(eK)^τ) だが、論文は「fixed」とだけ書き、辺数と行加算数の積を明示しない。正当化は可能で、誤りではない。
4. **【gap・小】`05-layers.tex:632–636` の σ の条件 max{0, log_m(s/W)} ≤ σ** は「rational で strict に」と書くが、実際の log_m(s_c/W_c) = 1 − 2.7×10^{−13} < 1 − 2^{−50} = σ で成り立つ（検算済み）。
5. **【記述不明瞭】閾値の大きさ。** d = ⌊b^{ε}⌋ ≥ 2 には b ≥ 2^{2^{75}}、K = ⌊d^c⌋ ≥ 6 には b ≥ 6^{2^{131}} が要る（`08-assembly.tex:63–69`、`05-layers.tex:954–959`）。つまり宇宙の原子数をはるかに超える入力長まで、機械は筆算（または HvdH）で動く。定理の O(·) 表現としては正しいが、「every input length」が実質的に意味するのは「閾値以下は筆算」であることを論文は 1 文で済ませている（`08-assembly.tex:1001–1003`）。
6. **【記述不明瞭】テープ本数 ≈ 10^{21}（`04-swap.tex:215–217`）。** `05-layers.tex:586` は「fixed number of tapes」としか書かず、本数が W であることが読み取りにくい。「fixed finite number of tapes」の文言は正しいが、Thm main の読者はこの規模を想定しない。
7. **【記述不明瞭・軽微】** `04-swap.tex:206` の O(1) は実際には ≈ s ≈ 10^{26}（parking コスト）。`05-layers.tex:288–290` の変位「at most four」は実際には 2 以下（害はない）。`03-motifs.tex:433–435` の「complement of two neighboring triple lines contains a unit」は h > 6 を使うが明示していない（h = 100 なので問題ない）。
8. **確定した誤りは見つからなかった。**

## 6. 監査できていない部分

- HvdH 2021 から引用する数値ルーチンの時間評価（Lemma 4.8–4.12、2.12）と、本論文での record 変換・テープ移動の帳簿（`07-resampling.tex:203–296`）。恒等式と縮小性は数値で確認したが、時間は読んだだけ。
- `02-streams.tex` の基本走査のカウンタの帳簿（O(L+w) など）。標準的なので深追いしていない。
- `04-swap.tex:213–268` の固定テープ depth-first schedule の細部（parking、descriptor stack）。読んで矛盾は見つけていないが、実装レベルの検証はしていない。
- 層の再帰を「配列全体」で実装したときの、役割の分割と子呼び出しの descriptor の整合（`05-layers.tex:563–629`）。
- Cor exact-arithmetic（§9）と転置（§10）は読んだだけ。Newton の誤差式 (*) は手で追って問題なし。
- AFKL 2019 の定理の正確なモデル定義（本文は未取得。要旨のみ）。AHJ+06 の k テープ拡張は未確認。
- 大きな h, m, W での実機シミュレーションは不可能（m = 10^6、W ≈ 10^{21}）。「s < Wm」の数え方は小 h の構造検証と記号計算に頼っている。

## 7. 必要な専門家

- **計算複雑性（下界・ネットワーク符号化）**: §2(a)(c) の判定。まず Lemma chunk-swap（`03-motifs.tex` bit network + `04-swap.tex`、約 1,100 行）を「無向グラフで coding が routing に勝つ」具体例として検討し、AHJ+06 の k テープ拡張の可否を判定する。次に oblivious 性の確認と AFKL の回路モデルへの帰着。これが最優先。
- **整数乗算アルゴリズム（Harvey–van der Hoeven 系統）**: §5–8 のテープ上の帳簿付けと、HvdH 2021 の引用の妥当性。
- **有限体上の双線形形式・符号理論（Frankl–Wilson、index coding）**: §3 の有限ネットワークのラベル構成の独立検証。

## 8. 総括（約 150 語）

- 主張は「固定テープ多テープ TM で乗算が O(n(lg n)^{1−κ})、κ = 2^{−182}」。Schönhage–Strassen 予想の否定である。
- 削減の源は、h=100 の 3 元集合から作る有限線形ネットワークで、フレーム変更の総階数が Wm より 10^{−11} の割合だけ小さいこと。κ は変換レイアウトのチャンク交換の行から出る。
- 定数、指数の連鎖、κ の導出、誤差の連鎖、小 h でのネットワーク、ラベルの入れ子、位相フレーム、packed 加算、変換恒等式をすべて自作コードで検算し、誤りは見つからなかった。
- 最大の懸念は外部との整合である。アルゴリズムは oblivious に見える。そうなら AFKL 2019 と合わせてネットワーク符号化予想が偽になる。chunk-swap 単体でも AHJ+06 と衝突しうる。論文はこれを論じない。
- 確信度: 検算した部分が正しいことは高い。定理全体が正しいことは低い（2 割程度）。NCC 反証と同等の主張であり、§3–4 の有限ネットワークの専門家による検証が不可欠である。
