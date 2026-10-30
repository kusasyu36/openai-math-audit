# 監査 077: 正曲率曲面に対する Fourier restriction(OpenAI, 2026-09-24)

対象 (A) "Elliptic capacity propagation and Fourier restriction to the sphere"(103頁, tex: build/sections/00–08)
対象 (B) "Diagonal Fourier extension for positively curved surfaces in three dimensions"(25頁)
スクリプト: audits/077/check_exponents.py, check_l4.py, check_examples.py(リポジトリのコードは実行していない)

## 1. 主定理の正確な内容
- **(A) Thm 1.1**: 単位球面 S² の extension 作用素 E は、すべての p>3 で L^∞(S²)→L^p(ℝ³) 有界。定数 C_p は p のみに依存し、端点 p=3 については何も主張しない。これは 3次元の有界データ版 restriction 予想そのものであり、部分範囲ではない。
- **(A) Thm 1.2 = Thm 6.4**(packet 定理): 任意の有限複素係数配列 d と任意の座標マスクに対して 𝒜(T_m d)³ ≤ C m^{10κ+ε} G(d)³ が成り立つ(0<κ<1/10)。G は楕円 capacity K(c)=θ² sup Σ|c|²/W(E)(W=L^{1+κ}w^{1−κ}、速度と位置を独立に平行移動したテスト)を使ったゲージで、𝒜 はそこから作る atomic norm。
- **(A) Cor 1.3**: Kakeya maximal を主張する。‖K_δ f‖_{L³(S²)} ≤ C_ε δ^{−ε}‖f‖_{L³}(管の位置について sup を取る強い形)。Thm 1.1 に Bourgain の factorization と双対性を組み合わせて導く。
- **(A) §9**: 球面の L^q→L^q(q>3)と、混合ノルム L^a→L^b(b>2a', b>3)。
- **(B) Thm 1.1**: 境界付きでもよいコンパクト滑らかな正曲率曲面に対し、すべての 3<p<∞ で L^p(Σ)→L^p(ℝ³)。証明の中身は各グラフについて L³(Σ)→L^p(p>3)を示すもので、Hölder で L^p ⊂ L³ に落としている。「diagonal」は入力と出力が同じ指数 p という意味。系として Schrödinger 局所平滑化(s>2−6/p, p>3)と Bourgain 条件付き振動積分を述べている。
- 定量的主張は ((B) 式(1.3)、(5.x)): C^best ≤ C(p−3)^{−7/p}P(θ,θ)^{1/p}M(θ)^{3/(2p)}(θ=10^{−6}(p−3)²)と、下界 ≥ c(p−3)^{−1/p}。

## 2. 最先端の結果および必要条件との比較
- 先行結果(論文自身の引用による。ウェブでの独立確認はしていない): Guth 2016 は p>13/4(有界データ)、Wang の broom は p>42/13≈3.231(放物面)、Wang–Wu 2024 は diagonal で p>22/7≈3.143(境界込み)。依頼文にあった「Wang–Wu p>3.04」は論文中に見当たらない。
- 必要条件 p>3(定数関数)と p≥2q'(Knapp)を厳密有理数で確認した。結果はすべて範囲内で、**範囲外の指数ペアは一つもない**。
  - (∞, p>3) は範囲内で、端点 p=3 は主張から除外されている。
  - (p,p), p>3 と (3, p>3) は境界線 p=2q' 上の点 (3,3) のすぐ内側にある。
  - 混合ノルムは b>2a' という厳密不等式で述べられている。
  - Schrödinger の閾値 2−6/p は Rogers の線 2n(1/2−1/p)−2/p(n=2)と一致し、臨界 p=2(n+1)/n=3 も一致する。
- 球面の L^4 定数 32π³ は g=1 で**等号**が成り立つことを確認した(∫sin⁴x/x²=π/4 から ‖E1‖₄⁴=32π³)。

## 3. 構成と依存関係
- **幾何部 (A §2–5, tex 約4,700行)**: Theorem 2.1(geo:main)が q* ≥ 1+2s−10κ を与える。道具立ては以下のとおり。
  - 楕円 capacity を使った hairbrush(Lemma 2.4)と heavy box の数え上げ、角度の排除
  - calibrated width path
  - 有限 Shannon エントロピーのプロファイル、接線、peeling、signed envelope(§3)
  - Ren–Wang の平面 Furstenberg 定理 Thm 4.1 から得る射影・radial 不等式と equal-slope tie(§4)
  - aspect 速度の3分岐(停留・増加・減少)による背理(§5)
- **packet 部 (A §6–7)**: 正確な frame 恒等式(Lemma 6.1)と atomic norm、合成則、局所性を整えたうえで、GIOW 流の refined decoupling(Thm 6.6、外部入力は Bourgain–Demeter ℓ² のみ)と離散4乗伝播(Prop 6.7)を使う。続いて二重極値の議論を行う。まず時間次元の下限 s を取り、次にパス数を最小化する。これで条件付き時間法則(Lemma 7.6)を作り、geo:main を適用する(Lemma 7.7–7.8)。最後に、frame から来る2次モーメント、decoupling から来る4次モーメント、capacity の3つを打ち消し合わせて γ≤10κ を得る。
- **大域化 (A §8)**: 疎な球族の上の3乗評価(Gram/Schur)、L^4 による峰の数え上げ、疎被覆、再生核、layer-cake。Tao の ε-removal の型に従う。
- **依存関係のまとめ**
  - (A) Thm 1.1 は Wang–Zahl も OpenAI Kakeya 論文も**使っていない**(Wang–Zahl は「文脈」として触れているだけ)。
  - (A) Cor 1.3 は Thm 1.1 だけから出るので、(A) 単独で 3D Kakeya maximal 予想を主張していることになる。
  - (B) は (A) Thm 6.4 と **OpenAI Kakeya maximal Thm 1.1** の両方に依存しており、未査読の OpenAI 原稿2本を前提とする条件付きの結果である。
  - (A) は OpenAI の Bochner–Riesz Cor 12.4 と Tao 1999 を経由する別ルートにも言及している。
  - **geo:main が破綻すれば、(A) Thm 1.1・Cor 1.3・(B) が同時に倒れる。**

## 4. 確認したこと(すべて checked OK)
- §7 の帳簿付け:
  - terminal capacity の指数 2(1−s)+20κ+2ζ(min(re_j,δ)=δe_j/f_j と、底上げでアスペクト比が減ることを含む)
  - (7.40)/閾値 = (7.41)
  - h^{1/2}min(m^{−2(1−s)},h^{−1}) ≤ m^{−(1−s)}
  - 最終的な増幅 ≤ m^{10κ+ζ}
  - γ_σ≤2σ
  - r₊w₊²ℓ⁴=r₋w₋²、N²wθ=1
  - 時間法則の one-bin 評価の論理、represented roots の数え上げ、2つの mass 評価
- hairbrush の Cauchy–Schwarz(m ≤ AM^{2e}η^{2κ})、heavy-count の指数 −6e−8τ、angular exclusion の指数 −7e−(c_κ+4)d−10τ(sympy で恒等的に一致)。
- §8 の大域化:
  - 疎な3乗評価のスケール N⁶δ⁶=1
  - Gram の減衰 H^{−1/2}
  - 被覆の分離不等式、峰の数 μ^{−7}、切断誤差 ρ^{−4}≤μ⁴
  - layer-cake の積分範囲 p>3+α
  - E1=2sin(2π|x|)/|x|、‖σ*σ‖₂²=32π³
- (B):
  - weighted Kakeya の双対計算、bundle 下界の定数 π/3、capacity 分解(段数 O(log n))
  - 共鳴除去(面積 Cμ^{24}、‖F−F*‖∞≤Cμ^{16})と位相の下界
  - 50Jγ₀=17τ/800、J≤34/τ、βJ>8
  - layer 積分 = 2/τ+12/τ²+48/τ³+96/τ⁴ ≤ 158τ^{−4}
  - L^4 補題の Jacobian と π/c、下界 (p−3)^{−1/p}
  - 相互参照番号(Thm 6.4、(6.8)、Def 6.3)が PDF と一致すること
- Kakeya 転送の指数 12/q−4 と −4(q−3)/(3(q−2))、および混合ノルムの補間恒等式。
- 例による整合性チェック(反例にならないかの確認で、証明の検証ではない):
  - Knapp の単一 packet と平面波では、Thm 1.2 は m^{10κ} を除いて**等号で飽和**する。
  - 集束波では焦点ビンが N^{−6}=G₁³ で、他のビンからは log N が出る。これは端点 p=3 での発散に対応している。
  - 平面に集中する直線族では q* ≤ 2+s−κ となり、主張 1+2s−10κ と矛盾しない。停留 aspect の場合の最終不等式 2+s−κ>1+2s−10κ は、この例と正確に対応している。

## 5. 懸念事項(深刻度つき)
- **確定的な誤り: 見つからなかった。**
- **ギャップ/正当化されていない段階: 行番号を引用できるものは見つからなかった。**
- **最重要リスク(未監査の核心)**: geo:main(A, Thm 2.1, 01-geometric.tex l.107)は、s=1, D=1 の場合だけで 3D Kakeya 集合予想(Minkowski 形、δ^{−10κ} の損失つき)を**単独で導く**。帰着は次のとおり。
  1. 各 M^{−1} 方向セルに直線を1本ずつ置き、T を一様分布とする。このとき K=O(1)、s=1、ε=0 となり仮定を満たす。
  2. 出力範囲は単位円板だけに縮退する。
  3. したがって (M^{−1})³ セルあたりの質量は ≤ M^{−3+10κ+o(1)} となり、管の和集合の体積は ≥ M^{−10κ−o(1)} となる。
  
  s=1 では指数 3 が体積の天井なので、余裕はゼロである。つまり Wang–Zahl 級の内容を、Ren–Wang + エントロピー + hairbrush という新しい方法で約4,700行の中に収めていることになる。ここは本監査ではほぼ未検証で、正しさの判断は全面的にこの部分にかかっている。
- **説明が不明瞭(軽微)**:
  - §8 の Gram 行あたりの項数 Cμ^{−8}R² は、窓の重なりを考えると過大評価である。上界としては無害。
  - §9 の管の包含の説明「axial displacement ≤ 1/2+δ²<3/4」は、2本の管が軸方向 [−3/4,3/4] を覆うという意味だと読めば正しい。
  - Thm 6.6 は SinghSingh2026 の反例(arXiv:2608.28711)を引用しているが、この文献は確認できていない。

## 6. 監査できていない部分
- (A) §2:
  - 確認済みは hairbrush、heavy count、angular exclusion のみ。
  - Bounded complexity、range extension、width band、axis stability、calibrated widths は未検証。
- (A) §3: エントロピー全体(未検証)。
- (A) §4:
  - Ren–Wang の引用形は双対形として妥当であることを確認した。
  - radial、scalar 不等式、product coefficients、equal-slope tie は未検証。
- (A) §5: 停留ケースは流し読みのみ。増加・減少ケースは未検証。
- (A) Thm 6.6: refined decoupling の帰納法の細部は部分的にしか確認していない。
- (A) Prop 6.7 と座標マスクの相性: 主要な点は確認した。
- (B) §6: 曲面への転送と系。
- OpenAI Kakeya maximal 論文本体。

## 7. 必要な専門家
- (a) Furstenberg/Kakeya の専門家: geo:main(§2–5)の全体。上記の s=1 の帰着を入口にする。
- (b) decoupling の専門家: Thm 6.6 と Prop 6.7(マスク付きの多重度 k の扱い)。機械が書いた wave packet 議論で典型的に破綻する箇所である。
- (c) (B) の条件付き性について: OpenAI Kakeya 論文の独立監査が必要。

## 8. 総括(150語以内)
(A) は 3次元球面の有界データ版 restriction 予想を全範囲 p>3 で主張している。(B) は正曲率曲面について L^p→L^p(p>3)を主張している。どちらの指数も必要条件の範囲内で、端点は主張していない。本監査では §6–8、(B) 全体の帳簿付け、hairbrush 系の指数、定数を厳密計算で確認し、誤りは見つからなかった。Knapp・平面波・集束波の例でも主張は破れず、飽和するだけだった。ただし核心の geo:main は、s=1 の場合だけで 3D Kakeya 予想を導くほど強い。これを新手法の約4,700行で示しているが、その部分はほぼ未監査である。(B) はさらに OpenAI Kakeya 論文に依存している。確信度は、監査した帳簿部分に誤りがないことについては高い。定理全体が正しいかどうかは、核心が未監査のため判定できない。
