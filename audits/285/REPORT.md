# 監査レポート 285: A torsion-free counterexample to the Kadison–Kaplansky projection conjecture（OpenAI, 2026-09-23）

- 監査日: 2026-10-31（監査役: Claude サブエージェント）
- 主対象: `data/openai_math/preprints/A-Torsion-Free-Counterexample-to-the-Kadison-Kaplansky-Projection-Conjecture-September-23-2026/build/sections/*.tex`（11 ファイル、約 5,500 行）
- 副対象（序文と構成の確認のみ）: `A-Torsion-Free-Counterexample-to-Reduced-Baum-Connes-Injectivity-September-23-2026/build/sections/00_introduction.tex`
- 第 3 論文（irrational trace）は取得していない。
- 自作スクリプト（`audits/285/`。出力は同名 `.out`）:
  - `check_scalar_compression.py` — 自由積による scalar 射影化（scalar.tex Prop. `scalar:compression`）の下界と形式的恒等式
  - `check_nilpotent_group.py` — 冪零群 E の群法則・コサイクル・次数付け・有限表現 R・Magnus 検出子（nilpotent.tex）
  - `check_parameters.py` — 定数と指数の帳簿 81 項目（setup / nilpotent / spectral / geometry / sampling / expansion / pruning）。**81/81 OK**
  - `check_model_band.py` — モデル作用素のバンド構造（Lemma `spec:model`）と「1 辺 1 電圧 1 重み」サンプリング機構の有限モデル
- **リポジトリのコード・スクリプト・Lean は実行していない。**

## 1. 主張の正確な内容

**Thm 1.1**（introduction.tex:13–21）
- 有限生成・ねじれなしの離散群 G_proj と、射影 e ∈ C_r^*(G_proj) で 0 < τ(e) < 1/2 を満たすものが存在する。
- τ は標準トレース τ(x) = ⟨xδ_1, δ_1⟩（introduction.tex:9）。C*-環は **reduced** 群 C*-環（左正則表現のノルム閉包、scalar.tex:10–13）。
- 「非自明」は e ≠ 0, 1 の意味。これは 0 < τ(e) < 1/2 から従う。トレースの値は「正で 1/2 未満」としか述べられず、具体値はない。監査役の概算（∫_0^δ N t^{n_f}(δ−t)^{n_f} dt からの粗い見積り）では q^{−n_f−1} 程度（q = d^{A_*}）と極端に小さい。論文はこの値を述べていない。
- 群は有限生成。与えられるのは無限グラフ表示（P の全閉路を関係子）に有限生成自由群 F_N を自由積したもの（introduction.tex:25–26、geometry.tex:90–94、conclusion.tex:13）。**有限表示は与えられておらず、主張もされていない**（有限表示可能でないことの証明もない）。可算である。
- G_proj = G_base * F_N。G_base は Thm `setup:matrix-target`（setup.tex:6–17）の群で、自己共役 T ∈ M_N(C[G_base]) が正則表現で孤立したスペクトルバンドをもち、その射影 p_0 ∈ M_N(C_r^*(G_base)) が **非正規化**トレース 0 < τ_{G_base,N}(p_0) < 1/2 をもつ。
- スカラー射影は e = (Up_0)(p_0U^*Up_0)^{−1}(Up_0)^*（scalar.tex:69–73、U = (λ(t_1),…,λ(t_N))）。τ(e) = τ_{G_base,N}(p_0)。

**既知の帰結との整合（論文自身の主張）**
- 伴論文 OpenAIBass2026 により C[G_proj] には非自明な冪等元はない（introduction.tex:60–68）。つまり e は群環には入らず、完備化の中にだけある。本監査ではこの伴論文を確認していない。
- 非整数トレースは、係数なし reduced Baum–Connes 全射性の反例でもある（conclusion.tex:29–32、BGW2016 Prop 10.3 を引用）。

## 2. 既知の結果との比較

| 既知の正の結果 | 対象クラス | 本論文の群が入るか |
|---|---|---|
| Higson–Kasparov（係数付き BC） | a-T-menable（Haagerup）群、可算アメナブル群 | 不明。論文は論じていない（下記） |
| Mineyev–Yu、Lafforgue | 双曲群と、そのねじれなし部分群（MY2002 Thm 21。部分群は有限表示でなくてよい） | 独立には排除できない。無限表示が与えられているだけで、有限表示可能でないことは示されていない。排除されるのは定理の帰結（矛盾するため）としてのみ |
| Lafforgue の Banach KK（RD + 有限 BG など） | property RD と有限次元分類空間をもつ群など | 適用根拠はない。ただし排除もできない（同上） |
| Guentner–Higson–Weinberger | 線型群（単射性のみ） | 無関係（全射性の話） |
| Pimsner–Voiculescu | 自由群 | G_proj は自由積だが因子 G_base が自由でない |

- **矛盾する既知定理は見つからなかった。** ただし次の点は論文が説明していない。
  - Thm 1.1 が正しければ G_base は a-T-menable でない（Higson–Kasparov に反するため）。しかし論文は G_base が Hilbert 空間に粗埋め込みできないことも、Haagerup 性をもたないことも論じていない（`a-T-menable`, `Haagerup`, `coarse embedding`, `property A` の語は本文にない。grep で確認）。
  - Cayley グラフに埋め込まれるのは expander そのものではなく、expander 塔 Y の **冪零群 F による電圧被覆 P** の成分である（geometry.tex:64–74）。アメナブル群による被覆は expander ではないので、「Gromov monster だから粗埋め込み不可」という標準的な論法はそのままでは使えない。
  - したがって「この群が Haagerup 性をもたない理由」は定理の帰結としてしか得られず、独立の裏付けがない。専門家が最初に疑うべき点である。
- 相対双曲性との関係（checked 非矛盾）: Lemma `geo:thin`（geometry.tex:424–430）はコピーを錐で潰した X̂ が細い三角形をもつことを示し、錐の安定化群は冪零群 F に埋め込まれる（geometry.tex:828–849）。もし X̂ が有限個の G_base 軌道のスポークしかもたなければ、Bowditch 型の相対双曲性から「BC を満たす部分群に相対双曲な群は BC を満たす」型の既知結果（出典は専門家の確認要）が効いて矛盾しうる。しかし Prop `geo:incidences`（geometry.tex:759–790）によりスポークの軌道は V(Y) で添字付けられ無限個なので、相対双曲性は従わず、矛盾はない。
- Kadison–Kaplansky と整数性: 論文はねじれなし群に対する「assembly の像上でトレースは整数」（Lück Thm 0.3、BGW Prop 10.3）を正しく引用している（introduction.tex:40–47）。主張は「完備化の K_0 全体では整数でない」であり、既知の整数性と矛盾しない。
- Atiyah 予想との関係: p_0 は核射影ではなく、ギャップで孤立したバンドの射影なので、Atiyah 予想（核の次元の整数性）とは直接関係しない。
- irrational trace: 本論文はトレースの無理数性を主張していない（第 3 論文の主張。未取得）。
- HLS 2002・Willett–Yu との区別は正しく述べられている（係数付きのみ、introduction.tex:70–78）。

## 3. 構成と依存関係

**群の作り方（概要）**
1. 層 i ∈ Z をもつ無限塔 Y^0（tower.tex）: |A_i| = d·2^{H_i}, |B_i| = 2^{H_i}, d = 2^{4ℓ}, D = d², H_i = 20010ℓ + 1 + |i|（setup.tex:72–79、tower.tex:49–52）。隣接層の二部接続は (d, D_{ik})-biregular で、正規化作用素は定数ベクトルを保ち直交補空間で ≤ 1 − χ、χ = 1/800（expansion.tex:530–542）。
2. 同じ頂点集合に「長い語から経路を復元できる」ラベル付け Y^#（tower.tex、pruning.tex）。両者を一致させるため頂点の o(1) 割合を削除し、生存部分 Y を得る（pruning.tex:252–277）。Y は folded、girth > g_i = H_i/(100 log_2 d)、長さ ≥ L_i = 100H_i/log_2 d の非後退路はラベル語で一意。
3. 電圧群 F（nilpotent.tex:255–278）: 上三角冪零群の affine 拡大 E_Z 二つを中心で貼り合わせたもの。可算・冪零・ねじれなし。中心生成元 z_0 の Fourier 変数 t ∈ T でファイバー分解し、ファイバーは twisted 群環。
4. 各辺に電圧 z_e ∈ F と行列重み W_e ∈ M_N(C) を **ランダムに 1 回ずつ** 選ぶ（sampling.tex:607–740）。同時に (a) 被覆 P の girth > C_* H_i/log_2 d（Schwartz–Zippel、sampling.tex:742–839）と (b) ‖T_P − ΠOΠ‖ ≤ η（モーメント法、sampling.tex:841–1046）が正の確率で成り立つことを示す。
5. G_base = ⟨S | P の全閉路のラベル語⟩（geometry.tex:90–94）。S は Y のラベルに重みタグを加えた有限アルファベット、|S| = O(d^{89/32}(log d)^3) = o(dD)（sampling.tex:690–699）。
6. 幾何（geometry.tex、Thm `geo:main`）: P の各成分は Cayley グラフ X に埋め込まれる。コピーに整列順序を入れ、各コピーの「先行コピーとの交わり」を 1 つの有限木球 B(Q) に閉じ込める（|B(Q)| ≤ 2^{C_• H_i}）。割当領域 D_Q は互いに素で、D_Q 内部の辺を除いた残りは森。G_base はねじれなし。
7. スペクトル（spectral.tex）: モデル O(t) の上バンド（固有値 1、固有空間 = 活性 2 層の定数ベクトル ⊗ p(t)）を被覆へ、さらに Cayley 作用素 T へ移す。森の誤差 f_ℓ = O(√(|S|/dD)) → 0。トレースは全コピーの局所上空間への制限写像の単射性で τ(p_0) ≤ S_0 ≤ S ≤ 96μ_ℓ = o(1)。
8. 自由積 F_N を付け加えて scalar 射影へ（scalar.tex）。

**「小さいトレース」の仕組み（nilpotent.tex:728–982）**
- ファイバー射影 p(t) = p_t ⊗ p_{t−δ}[V_t]。p_θ は Rieffel 型の Schwartz 射影の格子版で tr_θ(p_θ) ≤ C|θ|^{n_f}。有限「磁気」表現 R（nil:R）が第 2 因子の実効パラメータを t − δ にずらす。
- 層の窓 h_i(t) は b(t) = κ^{−1}log_2(t/(δ−t)) で決まり、t → 0 で i → −∞、t → δ で i → +∞ に逃げる。したがって射影の場 t ↦ E_P(t) は t = 0, δ で不連続だが、各行列成分は t について連続で、ギャップは t について一様。これが「被覆の C*-環に非整数トレースの射影がある」ことの核心である。
- 重み付き評価: 𝔴_i N t^{n_f}(δ−t)^{n_f} ≤ C d^{−(A_* n_f − C_• h_*)} 2^{−(κ n_f − C_•)|i|}（nilpotent.tex:956–961）。指数の正値性は `check_parameters.py` で確認した。

**外部定理への依存**
- van Kampen 補題、Gross–Tucker の電圧被覆、Gersten 型の等周不等式→細い三角形（自前で証明）、Schwartz–Zippel、Magnus 展開、Rieffel の Heisenberg 表現の Schwartz 係数（nilpotent.tex:401–403）、Echterhoff–Williams の twisted absorption（nilpotent.tex:793–794）、Brown の周期分解（geometry.tex:931–932）、BGW2016 Prop 10.3。
- **OpenAI 自身の別論文への依存**: 主定理の証明は伴論文に依存しない。伴論文 OpenAIBass2026 は「e が群環に入らない」という付随的主張にだけ使われる。Sauers2026 は条件付きの先行結果として言及されるだけで、本論文は独自の Prop `scalar:compression` を使う。
- 副論文（BC 単射性）は別の群 G_inj（Heisenberg 群の電圧、Bott 消去）を作る。本論文とは構成が異なる。共有するのは塔・ラベル付け・電圧被覆の技法である。

## 4. 確認したこと（checked OK）

**手計算で確認した論理**
- Prop `scalar:compression`（scalar.tex:52–183）: 自由積の正規形による good/bad 分解、(scalar:pointwise) |ξ_j(v)|² ≤ a_j m_v、(bad-input) ≤ a‖ξ‖²、(good-output) 等長性、(bad-output) ≤ a‖ξ‖²、c_a = √(1−a) − √a > 0 ⇔ a < 1/2、A ≥ c_a² p、e² = e、τ(e) = a（トレース性）。自由積のねじれなし性の論法（scalar.tex:176–182）も標準どおり。
- spectral.tex の作用素論: Lemma `spec:compression`（a + 20ε）、Lemma `spec:perturbation`（u/(r_0 − u)）、Prop `spec:cover-projection`（η + 20√β_ℓ ≤ α、α/(r_0−α) < 1）、Prop `spec:localization`（4θ^k、6 max 𝔴_j、96μ_ℓ）、Lemma `spec:copy-cut`（HS ノルムによる (spec:ball-HS)）、Lemma `spec:forest-norm`（2w_1√m_1）、Prop `spec:global-gap`（‖T − C‖ ≤ 3α、1 − χ/4 の位置）、Lemma `spec:lost-restrictions`（Cauchy–Schwarz の |B(Q)| を 𝔴_i で支払う）、Lemma `spec:kernel-form`（(1−χ+α) + χS + f_ℓ ≤ 1 − χ + 3α、p_0H ∩ K = 0、極分解による Tr p_0 ≤ Tr P_{K^⊥} ≤ S_0）。
- Lemma `nil:trace` の adjoint-square 恒等式と射影の比較。Lemma `spec:model` の対角 h_i² p/(2|A_i|) と tr E_* = tr p。
- geometry.tex: Lemma `geo:arcs` の Euler 計算（b ≤ (30/7)(m−1) ≤ 6m、m = 1 で J = 0）、Lemma `geo:chain-filling` の 12n、Lemma `geo:thin` の数値（R ≤ 345600）、Prop `geo:balls` の L^* = 4D_0 + 8、|B(Q)| の指数 7500 → 8000 → C_•/2 = 10^4、Prop `geo:cycle-directness`、Prop `geo:forest`、Prop `geo:incidences`、Lemma `geo:stabilizers`、Lemma `geo:acyclic`、Prop `geo:torsion`（自由胞体作用 + 長さ 2 の自由分解 vs 周期分解）。
- sampling.tex: Lemma `sam:boxes`（アメナブル群の Følner 箱による正則ノルムの有限テスト）、Lemma `sam:empirical`（モーメント法）、Lemma `sam:walk-count`、(sam:block-probability) の指数 C_1 − C_2/2 ≤ −4、girth の和集合評価 E_0 > 3C_* + 5、(sam:second-width)・(sam:first-width) の指数、Lemma `nil:word` の多項式が全変数で非零であること（木の辺を単位行列に特殊化しても非零）。
- expansion.tex:544–612 の Cheeger 計算（(1/20)²/2 = 1/800）。

**自作コードで確認した事実**
| スクリプト | 内容 | 結果 |
|---|---|---|
| check_scalar_compression | K = Z/3 * F_2、非対角な p ∈ M_2(C[Z/3])（a = 1/3）。200 本のランダムな ξ = pξ で ‖Uξ‖/‖ξ‖ ≥ c_a = 0.239（実測最小 0.960）。(bad-input)(good-output)(bad-output) も各ベクトルで成立。形式的恒等式 e² = e、τ(e) = τ_N(p) は有限次元の代用で確認 | OK |
| check_nilpotent_group | m = 3 で群法則の結合律・逆元公式・整数多項式性、コサイクル恒等式、c(Δ_yγ,Δ_yη) = y^m c、π_ε の乗法子、R mod q（q = 5）が表現で乗法子 e(−δc)、Λ_q 上で R = 1。Magnus 検出子は m = 4, 5 で 600 語すべて非自明、次数 ≤ m − 1 | OK |
| check_parameters | 81 項目の帳簿。h_* = 5003（H_0/log_2 d = 5002.5 + 1/(4ℓ)）、アルファベット指数 89/32、窓の不等式、スペクトル定数、Euler 計算、R ≤ 345600、|B(Q)| の指数、girth 和集合評価、(sam:C2-order)、Cheeger、pruning の 62.5H_i | 81/81 OK |
| check_model_band | Part A: 4 層の有限塔・ランダム biregular 接続・階数 2 の p で、O(t) の固有値 +1 の重複度が正確に rank p、残りが ≤ 1 − χ、対角が h_i²p/(2|A_i|)、全頂点の対角和 = tr p。Part B: Z/24 の畳み込み射影（階数 5）を「1 辺 1 電圧 1 重み」で標本化。‖T − U⊗K‖ ≈ 2.75/√d（d = 16…1024）、上バンドの重複度は常にモデルどおり 5 | OK |

**パラメータの実感**: C_* ≈ 4×10^13、C_• ≈ 8×10^10、m > L(201C_*) ≈ 6×10^50、n_f ≈ 2×10^101、N = q^{n_f}。「有用な評価は求めない」（introduction.tex:27–28）は正直な記述である。

## 5. 懸念事項

| # | 分類 | 重大度 | 場所 | 内容 |
|---|---|---|---|---|
| 1 | gap（独立の裏付けなし） | 中 | 本文全体（introduction.tex:49–58） | 定理は G_base が a-T-menable でないことを含意するが、論文はそれを独立に示さない。Cayley グラフに入るのは expander ではなくその冪零群被覆なので、粗埋め込み不可の標準論法は使えない。専門家は「G_base が Haagerup 性をもつ理由がないか」をまず検討すべき |
| 2 | 未監査の核心 | 高 | tower.tex, expansion.tex, pruning.tex（約 1,400 行） | Prop `prune:tower-properties` の 5 条件（一様な expansion、削除割合 ε_ℓ → 0、folded、経路復元、girth）を本監査は検証していない。失敗すれば spec:model と geo:arcs が同時に倒れる |
| 3 | 未監査の核心 | 高 | sampling.tex:346–740 | モデル h_i h_k p(t) を「1 辺 1 電圧 1 重み」の期待値で η 以内に近似する部分（テンプレート、Fourier マスク、合同部分群による共通リスト、(sam:edge-mean)）は読んだが細部は検証していない。特に (sam:subgroup-product) の比較（l.545–594）は込み入っている |
| 4 | gap（一様性） | 中 | 全体 | 「十分大きい ℓ」の条件が多数あり、層 i ∈ Z について一様であることが各所で主張される。帳簿は確認したが、すべての o(1) が i に一様かは個別に追っていない |
| 5 | 記述不明瞭 | 低 | spectral.tex:11–20、sampling.tex:38–41 | 「identity trace」の規約（基底頂点 v の和と N 行列指標の非正規化和）が複数箇所で繰り返されるが、ℒ ≅ ℓ²(G)⊗ℓ²(V(Y))⊗C^N の同一視で左作用の等変性を使う点は Lemma `nil:trace` が右作用で述べられており、読者は nilpotent.tex:316–317 の一言に頼る |
| 6 | 記述不明瞭 | 低 | sampling.tex:161–191 | Thm `sam:realization` の「full regular norm」は F の右正則表現のノルム。F がアメナブルだから Følner 箱の有限テストで済む（Lemma `sam:boxes`）。この依存関係（アメナブル性が本質）が明示されていない |
| 7 | 記述不明瞭 | 低 | nilpotent.tex:949–963 | Prop `nil:weighted` は sup_t の評価だが、実際に必要なのは t での積分。sup で足りるのは spec:localization の点ごとの評価のためで、設計上の余裕が大きい |
| 8 | checked OK（注意） | — | scalar.tex:149–150 | c_a > 0 は a < 1/2 でだけ成立。Thm の「< 1/2」はここから来る。本構成では a は極端に小さいので問題ない |

- **確定した誤り: 見つからなかった。**
- 「行番号を引ける gap」も見つからなかった。上の 2, 3 は「読んだが検証できなかった」の意味である。

## 6. 監査できていない部分
- tower.tex 全体（ストリームと slot、padded permutation 規則、matching 変数）。
- expansion.tex:40–529（条件付き Fourier 評価、平均のギャップ、cut の tail）。Cheeger 部分のみ確認。
- pruning.tex:1–251（prefix 復元、短い閉路の数え上げ、削除予算）。Prop の証明の第 4 項のみ読んだ。
- sampling.tex:346–740 の近似誤差の細部（上記 #3）。
- geometry.tex の van Kampen 図式操作（l.140–343）は読んだが、Gruber/Ollivier 型の「source assignment を保つ簡約」の細部は専門家の確認が要る。
- 伴論文 OpenAIBass2026、Sauers2026、FKOP2026、Robert2025 は未確認。
- 副論文（BC 単射性）の本体、第 3 論文。

## 7. 必要な専門家
- (a) 幾何学的群論（graphical small cancellation、Gruber・Osajda の手法）: geometry.tex 全体と、G_base が Haagerup 性をもたないことの独立検証（懸念 #1）。
- (b) 作用素環・K 理論（Baum–Connes、uniform Roe 環、twisted 群環の連続場）: 「被覆の等変 Roe 環に非整数トレースの射影がある」という spectral.tex の仕組み、および Rieffel 型射影の格子版（nilpotent.tex:573–692）。
- (c) 確率的組合せ論・ランダム行列: sampling.tex のモーメント法と girth 評価、expansion.tex の expander 構成。

## 8. 総括
本論文は、有限生成ねじれなし群 G_proj = G_base * F_N の reduced 群 C*-環に 0 < τ(e) < 1/2 の射影があると主張し、Kadison–Kaplansky 予想と係数なし Baum–Connes 全射性を否定する。本監査では、scalar 射影化、冪零群と有限磁気表現、スペクトル移送とトレース評価、幾何的埋め込みとねじれなし性、定数の帳簿（81 項目）を手計算と自作コードで確認し、誤りは見つからなかった。「層の窓が t → 0, δ で無限遠に逃げる」ことで射影場の不連続性を回避する仕組みは整合的で、有限モデルでも再現できた。一方、塔・expansion・pruning（約 1,400 行）と標本化の近似誤差の細部は未検証であり、G_base が a-T-menable でないことの独立の根拠も論文にない。確信度: 確認した範囲の正しさは高い。定理全体については、未監査部分が大きく判定できない。
