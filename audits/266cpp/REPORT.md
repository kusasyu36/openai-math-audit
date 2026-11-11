# 監査報告: 成果 266「N(6)=3」— 配布 C++ ソースと TeX の補題の照合

- 監査役: Claude サブエージェント
- 監査日: 2026-11-11
- 前提となる前回の監査: `audits/266/REPORT.md`(tex の数学と保存記録の整合。#5「C++ が補題どおりか」が未監査だった)

## 冒頭

**結論(先に):** 配布された C++ 8 本を全行読み、tex の 5 つの段(seed、completion、phase-cover、refinement、exclusion)と付録 A のルーチン記述に 1 つずつ突き合わせた。**健全性を破るバグ(排除できていない領域を排除済みと記録する誤り)は見つからなかった。** 排除・被覆を決める比較は、読んだ範囲で全部が tex の必要条件と同じ向きで、余裕は外向き、失敗は未解決として次段に渡る。残るのは「読んで確かめた」ことの限界と、付録 A の誤差予算の数値・ビルド・実機の前提で、これらは §7 と §8 に書く。

**対象と読んだ量**

| 対象 | 読んだ行 / 全行 | 備考 |
|---|---|---|
| C++ ソース 8 本(seed 119、covout 143、ub 386、patchx 473、atlas 241、stageB 336、stageC 567、abi-fp-probe 51) | **2,316 / 2,316** | 全行。証明経路に入らない関数(§1.3)も読んだ |
| `build.sh` | 29 / 29 | |
| `VERIFY.md` | 80 / 80 | |
| `execution/commands.json`(2,594 行) | 目では全行を読んでいない | 100 件全部の argv・returncode を自作スクリプト `check_266cpp_handoff.py` で機械的に検査(88 個の初回 C の引数が全部同一で part/parts=0/1、各バッチが自分の f<j>.bin だけを読む、など 16 項目。§3.5) |
| `execution/check_records.py` | 読んでいない・実行していない | |
| tex: seed 490、completion 715、phase-cover 1,008、refinement 656、exclusion 715、execution 400、rounding 1,566 | **5,550 / 5,550**(この 7 本は全行) | 原文(.tex)で読んだ |
| tex: introduction 246、conclusion 11、paper.tex 79 | grep のみ | paper.tex は `\newtheorem` と `\appendix` の位置だけ確認 |

**自作スクリプト(`audits/266cpp/`)**

- `extract_266cpp_literals.py` → `.out`: ソースをテキストとして読み、浮動小数点リテラル(91 種 411 箇所)、7 以上の整数リテラル、数学ライブラリ関数の呼び出し箇所、`long double`/`float`/`#pragma` などを行番号つきで抽出。
- `check_266cpp_handoff.py` → `.out`: `commands.json` の 100 件の argv・終了コード、`results/*.txt` の最終行の件数、データファイルの先頭 8 バイトの件数を突き合わせる 16 項目。全部 OK。
- `check_266cpp_formulas.py` → `.out`: 「tex から導いた式」と「ソースを読んで Python に写した式」の数値比較 41 項目(うち数個は tex のスケール主張の検算)。全部 OK。

**守ったルール:** OpenAI の C++、`build.sh`、`check_records.py`、`lean/` は一切コンパイル・実行していない。読んだだけ。
**ルールからの逸脱 1 件(申告):** 作業の最初に 1 度だけ、`commands.json` の argv を要約する `python3 -I -c` のワンライナーを、カレントディレクトリが `data/…/saved-run` の状態で実行した(OpenAI のコードは import していない。`-I` つき)。以後の Python は指示どおりプロジェクト直下から `uv run … python -I audits/266cpp/<name>.py` で実行した。

**定理番号** は `paper.tex:27–34` の `\newtheorem`(全環境が theorem カウンタを共有、節ごと)と環境の並びから数え、前回の報告書の番号付けと一致した(§2 Def 2.1〜Prop 2.10、§3 Lem 3.1〜Prop 3.11、§4 Def 4.1〜Prop 4.16、§5 Lem 5.1〜Cor 5.8、§6 Lem 6.1〜Prop 6.9、§7 Thm 7.1、付録 A Lem A.1〜A.3)。

---

## §1 プログラム全体の流れ

### 1.1 実行ファイルと include

`build.sh:17–25` が作るのは 4 本(+プローブ 2 本)。フラグは全部 `-O3 -march=native -ffp-contract=off`。

| 実行ファイル | 翻訳単位 | include の連鎖 | 役割 |
|---|---|---|---|
| `covrun` | covout.cpp(`-std` 指定なし) | covout → seed | 9 位相の seed 被覆 |
| `arun` | atlas.cpp(`-std=gnu++17`) | atlas → patchx → ub → seed | A 段: U(2) セル探索と粗い atlas |
| `brun` | stageB.cpp | stageB → atlas → … | B 段: 半径の細分 |
| `crun` | stageC.cpp | stageC → stageB → … | C 段: ベクトル被覆とグラフ排除 |

下位ファイルの `main` は `HEADER`/`UB_HEADER`/`PATCH_HEADER`/`ATLAS_HEADER`/`REFINE_HEADER` マクロで無効化される(seed.cpp:117、ub.cpp:368、patchx.cpp:454、atlas.cpp:216、stageB.cpp:283)。

### 1.2 段と入出力(execution.tex:123–134、VERIFY.md、commands.json と一致)

```
covrun centers.txt
  seeds::build(112)            seed.cpp:83–96   X 5 ビン列・Y 尾 4 ビン列を列挙し必要条件で篩う
  coverSeeds::run(S,.032,2)    covout.cpp:99–133 共有ビンで対にし、既存中心の球に入らなければ中点を保持
  → centers.txt(n, 件数, 各行 X5 個 + s + Y 尾 4 個)  covout.cpp:141

arun centers.txt i 4 coarse_i.bin bad_i.txt   (i=0..3、添字 ≡ i mod 4)
  CenterList::load             atlas.cpp:173–187
  block::init → Asearch::runA(.103)   patchx.cpp:441–451  N_u=16→32→64、以後 512 まで再帰
    testA: ノルム判定 → trial::pass/certificate → local::test(Candidate::prep, covercompute::value, buildH::snap/patch)
  Atlas::add(R, R.rad, 1)      atlas.cpp:145–168  上限 .125 で合併
  → coarse_i.bin(Record 列)、bad_i.txt(未解決の seed 添字)、未解決があれば exit 1

brun .07 .013 mid.bin midbad.bin - coarse0..3.bin
brun .016 .005 fine.bin finebad.bin batches/f mid.bin      (256 件ずつ batches/f<j>.bin、半径は .016 に置換)
  readRecords → validate → 半径降順 stable_sort → 親ごとに refiner::run   stageB.cpp:285–335
  → 子の atlas、失敗した親(bad)、bad があれば exit 1

crun 12 4 1 .042 1 1 5000000 0 1 batches/bad<j>.bin batches/f<j>.bin    (88 バッチ)
  negative(): Layer::test ×4 層(12,24,48,96)→ cluster → graph::candidate(pairtest, joint, exist)
  → 排除できなかった記録を bad<j>.bin に保存、あれば exit 1

brun .011 .0038 fb1.bin fb1bad.bin - (空でない bad<j>.bin 全部)
crun 12 4 1 .042 1 1 5000000 i 4 fb1cbad_i.bin fb1.bin   (i=0..3)
```

### 1.3 証明経路に入らないコード(読んだが、排除にも被覆にも使われない)

呼び出し元を追って確認した。tex も「entry point から呼ばれない補助ルーチンは追加の排除判定ではない」と明記する(rounding.tex:7–9、completion.tex:659–662)。

- seed.cpp: `share`(:53–71)、`seeds::list`(:97–115)、`main`(:118)。
- ub.cpp: `scalar`(:200–210)、`derr`(:212–231)、`Usearch::test`(:310–332)、`runblock`(:333–357)、`main`(:369–385)。`bound` の `zK`, `E1`, `E2` と `z2off`(:138, :142)は `derr` 専用。`block` の `g`, `sumerr`, `qpol` は未使用。
- patchx.cpp: `main`(:455–472)。
- stageC.cpp: `graph::duo`、`startDuo`(:445–482、searchMode=2 専用。実行は searchMode=1)。
- 引数 `linear=1`, `detail=1` なので `Layer::test` と `graph::joint` の反復は 1 回。

---

## §2 tex ↔ C++ の対応表

依頼文が名指ししたルーチンは**全部ソースに実在し、名前も同じ**。

| tex | 内容 | C++(ファイル:行) | 照合 |
|---|---|---|---|
| seed.tex:192–193, 244–252 | s の選択、中点 gap | `seeds::mk` seed.cpp:74–82 | 一致 |
| Lem 2.4, 2.6(seed.tex:186–190, 229–233) | 向き・gap フィルタ | `seeds::build` seed.cpp:83–96、covout.cpp:122 | 一致 |
| Lem 2.7(seed.tex:264–277) | 和の Taylor 判定 | `sumok` seed.cpp:35–51 | 一致 |
| Lem 2.8(seed.tex:306–324) | 共有座標の区間 | `extend` covout.cpp:9–25、交差判定 :123 | 一致 |
| Lem 2.9(seed.tex:383–404) | 中心再利用の支持関数 | `coverSeeds::supdual` covout.cpp:40–72、`testBucket` :73–88 | 一致 |
| Prop 2.10(seed.tex:449–458) | seed 被覆 | `coverSeeds::run` covout.cpp:99–133、`main` :135–143 | 一致 |
| Lem A.2(rounding.tex:187–195) | 単位指数 | `expair` seed.cpp:21–31 | 一致(数値でも確認) |
| Lem 3.1(completion.tex:55–70) | 固定枠 | `frame::init` ub.cpp:37–68 | 一致 |
| Lem 3.2(completion.tex:136–162)、(completion:A-bounds) | q, p, A_* | `block::get` ub.cpp:130–148 | 一致 |
| Lem 3.3(completion.tex:206–223) | U(2) の閉被覆 | `uc::init` ub.cpp:159–169、格子 patchx.cpp:447–448 | 一致 |
| Lem 3.4(completion.tex:251–259) | セルの Taylor 定数 √6h、(√42/2)h² | `block::get` の `6*b.h`, `sqrt(6.*42)/2*b.h*b.h` ub.cpp:140, 143 | 一致 |
| Lem 3.5–3.7(completion.tex:313–447) | z、L、E_r | `block::setup` ub.cpp:95–129、`get` | 一致 |
| completion.tex:506–531 | 線形化の成分公式 | `block::setup`、`makeD` ub.cpp:179–184、`makeJ` :185–199、`uc::initDeriv` :170–174 | 一致(数値でも確認) |
| Lem 3.8(completion.tex:541–548) | ノルム排除 | `Asearch::testA` patchx.cpp:423–424 | 一致 |
| Lem 3.9, Cor 3.10(completion.tex:583–642) | 方向排除 | `trial::certificate` ub.cpp:236–251、提案は `trial::pass` :252–301 | 一致 |
| Prop 3.11(completion.tex:667–715) | 再帰と未解決 | `Asearch::testA`/`runA` patchx.cpp:416–451、atlas.cpp `main` :217–240 | 一致 |
| Def 4.1、Prop 4.5–4.6(phase-cover.tex:43–85, 247–354) | Root 証明書 | `buildH::gen` patchx.cpp:49–72、`buildH::patch` :108–187、`Root` :40–45 | 一致 |
| Lem 4.4(phase-cover.tex:188–245) | 剰余定数 2 種 | patchx.cpp:143(√8)、:144–166(テンソル) | 一致 |
| Lem 4.7(phase-cover.tex:358–373) | 同値の輸送 | `mapH::transform` atlas.cpp:19–28、`inverseTransform` stageB.cpp:4–12 | 一致 |
| Lem 4.8(phase-cover.tex:447–470) | 弧の上界 | `Candidate::prep` patchx.cpp:194–215 | 一致 |
| Lem 4.9–4.10(phase-cover.tex:507–628) | 中心補正、直接距離 | `Candidate::delta` patchx.cpp:216–221、`covercompute::value` :291–309 | 一致 |
| Lem 4.11–4.12(phase-cover.tex:676–756) | 楕円体 | `covercompute::value` patchx.cpp:310–347、`ell` :243–276、`spectralUpper` :224–231(Lem A.3) | 一致 |
| Prop 4.13(phase-cover.tex:852–865) | bootstrap | `bootstrap` patchx.cpp:278–288 | 一致 |
| phase-cover.tex:892–902 | 局所被覆 | `local::test` patchx.cpp:366–410 | 一致 |
| Lem 4.15、Prop 4.16(phase-cover.tex:917–1008) | atlas 合併 | `variants` atlas.cpp:47–88、`canonicalMap` :89–106、`Atlas::find` :120–144、`Atlas::add` :145–168 | 一致 |
| Prop 5.2(refinement.tex:106–135) | 親子の関係 | `relationB::init` stageB.cpp:61–106 | 一致 |
| Lem 5.3–5.4(refinement.tex:229–300) | 上界反復、ρ | `relationB::eval` stageB.cpp:45–60、ρ 探索 :101–105 | 一致 |
| Lem 5.5(refinement.tex:377–387) | 最遠の隅 | `relationB::corner` stageB.cpp:107–110、`markRec` :137–150、`mark` :151–160 | 一致 |
| refinement.tex:327–375 | 接座標の格子 | `refiner::run` stageB.cpp:250–272、`outside` :125–128、`solveCell` :168–202 | 一致 |
| Lem 5.6(refinement.tex:435–444) | 割り当てと半径の成長 | `assign` stageB.cpp:14–35、`setRadius` :129–136、`addNear` :161–167 | 一致 |
| Prop 5.7、Cor 5.8(refinement.tex:502–593) | 被覆の保存、バッチ | `preload` stageB.cpp:213–249、`main` :285–335 | 一致 |
| Lem 6.1(exclusion.tex:93–103) | ビンの変位 | `Layer::init` stageC.cpp:42–52 | 一致 |
| Lem 6.2–6.3(exclusion.tex:131–202, 257–280) | 単一ベクトル判定 | `Layer::test` stageC.cpp:53–149 | 一致 |
| Lem 6.4(exclusion.tex:318–355) | クラスタ化 | `cluster` stageC.cpp:182–226、`diff` :166–178 | 一致 |
| Lem 6.5(exclusion.tex:365–386) | ラベルの相異 | stageC.cpp:220 の `s<.1`、`pairtest` が自己対を作らない(:409) | 一致 |
| Lem 6.6(exclusion.tex:420–443) | 辺の Taylor | `graph::pairtest` stageC.cpp:368–400 | 一致 |
| Lem 6.7(exclusion.tex:487–503) | 共有接変数つき結合判定 | `graph::joint` stageC.cpp:308–367、`single::init`/`support` :233–278 | 一致 |
| Lem 6.8(exclusion.tex:570–629) | クリーク探索 | `graph::clique` stageC.cpp:416–427、`graph::exist` :428–444、`candidate` :483–488 | 一致 |
| Prop 6.9(exclusion.tex:639–664) | 局所排除 | `negative` stageC.cpp:496–531、`main` :534–566 | 一致 |
| Thm 7.1(execution.tex:242–281) | 有限実行 | 4 つの `main` の終了コードと bad ファイル | §3.5 |
| rounding.tex:40–46 | 環境検査 | `environment` seed.cpp:13–16、`checkcpu` atlas.cpp:6–12 | 一致(covrun は MXCSR を見ない。tex も明記) |

---

## §3 健全性の向きの確認結果

記号: **向き**=排除条件が必要条件の否定と同じ向きで余裕が外向きか。**失敗**=判定不能のとき保持・未解決側に倒れるか。**被覆**=ループが全体集合を尽くすか。

### 3.1 Seed 段

| 判定(行) | 捨てる条件 | tex の必要条件 | 向き | 失敗時 |
|---|---|---|---|---|
| seed.cpp:86 | `a>e-d+1` | a ≤ l−d+1(seed.tex:189) | 整数、同じ | — |
| seed.cpp:88 | `p.g>1.-(e-1)*w+1e-8` | g_x,mid ≤ 1−(l−1)/N(seed.tex:231) | 外向き | NaN なら保持 |
| seed.cpp:46(`sumok`) | `l>r+1e-8`(`!(…)` で保持) | (seed:sum-support) seed.tex:274–275 | 外向き | NaN なら保持 |
| covout.cpp:122 | `B.g+1e-8<1.-(l+2)*w` | g_y,mid ≥ 1−(l+2)/N(seed.tex:232) | 外向き | NaN なら保持 |
| covout.cpp:123 | 区間が厳密に離れる | Lem 2.8、端点共有は保持(seed.tex:340–342) | 外向き(分子 +1e-7、端点 ±1e-5) | \|a\|<1e-6 は制約を落とすだけ |
| covout.cpp:82, 85(`testBucket`) | 成功 → 中心を**追加しない** | ‖k‖²+9/4+2k·h の上界 < (NR)²−.0001(Lem 2.9) | 上界 +.01、閾値 −.0001 | 提案が使えないと `best=1e10` のまま → 中点を保持 |

- **被覆:** X は `F(e,n)F(d,e+1)F(c,d+1)F(b,c+1)F(a,b+1)` で a≤b≤c≤d≤e の全組(seed.cpp:85)。Y は先頭 a∈[0,n) × 尾 b≤c≤d≤e の全組(:92)。境界 0 と n−1 を含む。対は共有ビン `Y[A.p.b[A.p.s]]` の全件(covout.cpp:121)。
- **再利用の正当性:** 同じ s の中心だけを試す(covout.cpp:113)。k は「現在 − 既存」で X 5 成分(:116)と Y 4 成分(:77)の符号が揃っている。`supdual` の行の配置(Y の共有列が行 s、Y の尾が行 5–8)は tex の J_*(seed.tex:372–378)と同じであることを数値でも確認(スクリプト 7d)。`int(T)`、バケツ、`s-a>RN2` は試す中心の選別にしか使われず、covered が偽なら必ず中心を追加する(:125–128)。
- **添字:** 既存中心との差は `int8_t`(covout.cpp:28)に入るが、値域は [−111,111](rounding.tex:341–342)。

### 3.2 A 段

| 判定(行) | 捨てる/被覆する条件 | tex | 向き | 失敗時 |
|---|---|---|---|---|
| ub.cpp:41, 61(`frame::init`) | `!(m>1.25)`、`!(r>.45*.45)` で false | 枠の構成の受理条件(completion.tex:104–107) | — | atlas.cpp:229 で bad に記録、exit 1 |
| patchx.cpp:424 | `sum>b.z*b.z+1e-5` でセル排除 | Lem 3.8: Σ(\|D\|−1)² ≤ z² | 外向き | NaN なら保持 |
| ub.cpp:250(`certificate`) | `yc > sup+.0001` でセル排除 | Lem 3.9 の右辺 5 項 | 外向き。`neg` は 0 始まりで max(0, max(−y_D)) | 非有限・範囲外の y は false(:238–239) |
| ub.cpp:136 | `!(q<5.9)` | q<6(completion.tex:155, 198) | 厳しい側 | z=E=1e10 → 排除も被覆もできず未解決へ |
| patchx.cpp:389, 405(`local::test`) | `s>=0 && s+1e-8<=target` で被覆 | Prop 4.13 の上界 ≤ 目標 | 上界側に余裕 | false → 分割、深さ上限で unresolved |
| patchx.cpp:433 | N≥512 で未被覆 | — | — | `unresolved++` |
| patchx.cpp:418 | 追加ノード > 1e6 | — | — | `giveup`、`unresolved++`。以後 :417 で false を返すが unresolved>0 は残る |
| atlas.cpp:160 | `radius+dist+2e-8<maxR` で合併 | Lem 4.15: 半径 r+d | 既存半径を max(…, r+d+1e-8) に**拡大** | 合併できなければ挿入(捨てない) |

- **被覆:** U(2) の基底格子は `F(a,N/2)F(t,N/4)F(c,N)F(s,N)`(patchx.cpp:447)で Lem 3.3 の定義域 α∈[0,½]、t∈[0,¼]、β,γ∈[0,1] をちょうど覆う。子は各座標 2p, 2p+1 の 16 個(:448, :435)。`uc::init` の中点は α=(b0+.5)/N など、合成角は `whole[(b0+b2+1)%N]`、`whole[(b0-b2+N)%N]` で正しい(ub.cpp:164–165)。数値でも tex のパラメータ化と一致(スクリプト 3a, 3b)。
- **線形化:** `fixed[0]`、`fixed[1+k]`、`fixed[5+k]`、`bm`、`cm`、`W` を tex の線形写像 (completion:linear-map) から手で展開して一致を確認し、さらに乱数 200 組で `jac` と行列形の J(s,θ) を比べて差 3.4e-15(スクリプト 4b)。c_k の係数が D でなく D_ker(`d.K`)を使う点も tex どおり(completion.tex:521)。
- **レベルの扱い:** l=0(N=16)は it=0 でノルム判定のみ、l=1(N=32)はノルム+方向、l=2(N=64)以降が final(patchx.cpp:446, :426)。completion.tex:696–702 と一致。
- **Root の一覧は seed ごとにリセット**(patchx.cpp:365)され、seed の終了時に全部 atlas に入る(atlas.cpp:232)。未解決の seed の Root も入るが、その seed は bad に記録される。
- **`Candidate` の使い回し:** `local` のメンバで再利用されるが、`P` と `X` に書き込む位置の集合は毎回同じ(種 9 箇所と D ブロック 16 箇所)で、他は 0 初期化のまま。前の候補の値は残らない(patchx.cpp:192–206)。

### 3.3 B 段

| 判定(行) | 条件 | tex | 向き | 失敗時 |
|---|---|---|---|---|
| stageB.cpp:127(`outside`) | Σmax(0,\|u\|−h−1e-8)² > E² | 箱と原点の最小距離 > E_P(refinement.tex:353–361) | 下界を小さく、E を大きく(+1e-7) | NaN なら箱は残る |
| stageB.cpp:147(`markRec`) | 部分和 ≤ (ρ−5e-8)² の箱だけ印 | Lem 5.5 | 各軸 +1e-8、最後 +2e-8 でも ρ 未満 | 印が付かなければ通常の走査に残る |
| stageB.cpp:173 | `t < b.rho-1e-8` | Lem 5.5(t_B ≤ ρ) | 外向き | 次の手段へ |
| stageB.cpp:185 | `d<target-1e-7`(d=eval(t)+2e-8) | 直接割り当て(refinement.tex:431–433) | 外向き | 分割、深さ 6 か k=0 で false |
| stageB.cpp:102–105 | ρ 探索 | Lem 5.4 | 最後に成功した試行点 lo だけを使う | ρ≤1e-7 は不使用 |
| stageB.cpp:132–134(`setRadius`) | t≤ρ、d≤target | refinement.tex:314–323 | — | **abort**(黙って進まない) |
| stageB.cpp:27 | `radius+ℓ+2e-8<maxR` で合併 | Lem 5.6 | 半径拡大 | 挿入 |
| stageB.cpp:254, 193 | 格子 > 5e7、深さ上限 | — | — | 親は bad |

- **被覆:** `L=⌈E/w+.500001⌉`、中心 (b−L)w、半幅 w/2 で [−(L+½)w, (L+½)w] ⊇ [−E,E](乱数 2 万例でも違反 0、スクリプト 8f)。状態 0 の箱は全部 `solveCell` に渡り(stageB.cpp:266–270)、false が 1 つでもあれば親は bad。分割は 2^k 個の子 u±h/2 を全部訪れ、1 つ失敗で親セルも失敗(:197–200)。
- **添字:** `markRec` の `at+b*stride[axis]` と `run` の復号 `(idx/stride[a])%naxis` は同じ混合基数。b は [0,naxis−1] に切り詰め済み(:155–156)。
- **blocker の帰納:** `done` に入るのは `ref.run()` が成功した親だけ(stageB.cpp:318)。照会半径は元の入力半径 `rec.radius`(:218)。失敗した親は blocker にならず `bad` に入って保存される(:320, :328)。Prop 5.7 の二重帰納と同じ構造。
- **stageB.cpp:188**(`addNear` の戻り値を捨てる)は、被覆が :186 の `assign` で確立済みなので無害。tex も明記(refinement.tex:476–478)。
- **写像した親** V=ΦP への関係を同じ格子添字で使えること((ΦT)^*(Φδ)=T^*δ、refinement.tex:480–487)、`inverseTransform` が `transform` の逆であることを式で確認。

### 3.4 C 段

| 判定(行) | 捨てる条件 | tex | 向き | 失敗時 |
|---|---|---|---|---|
| stageC.cpp:63 | \|S_i\|² < (1−τ)²−2e-7(τ<1 のみ)、> (1+τ)²+2e-7 | (exclude:direct-tests) | 外向き | NaN なら保持 |
| stageC.cpp:70 | Σ(\|S_i\|−1)² > (2π(KX+r))²+2e-7 | 同上 | 外向き | 同上 |
| stageC.cpp:85 | \|c_i\| > aΣ_{j>0}\|J_ij\|+r‖J_i‖+min(ρ₁,ρ₂)+2e-7 | (exclude:scalar-linear) | 外向き | 同上 |
| stageC.cpp:138 | \|y·c\| > bJ+bB+bE+2e-7 | (exclude:bin-support) | 外向き | 提案が非有限なら break → 保持 |
| stageC.cpp:147 | 保持ビン ≥ 2,000,000 | — | — | overflow → 記録は bad |
| stageC.cpp:204, 219–220 | 球 > 2048、\|v\|>2、s∉[0,.1) | Lem 6.4, 6.5 の前提 | — | clusterFail → bad |
| stageC.cpp:378, 389, 397 | 辺を張らない | (exclude:edge-distance)、(orth-support)、(unbiased-support) | 外向き(r+1e-10、+2e-7、ε+2e-8) | — |
| stageC.cpp:364(`joint`) | \|yc\| > sup+3e-7 で辺を張らない | (exclude:joint-inequality) | 外向き | 提案が不良なら true(辺を残す) |
| stageC.cpp:520 | ある層で保持ビンが空 → **排除** | Prop 6.9 の経路 1 | — | — |
| stageC.cpp:527 | `candidate()` が false → **排除** | Prop 6.9 の経路 2 | — | cap は `gaveup` で true |

- **被覆:** `gen0` は N^5 全ビン(第 0 座標は 0 固定、stageC.cpp:151)。`gen1` は 32 個の子全部(:154–157)。`cluster` は全ビンに球を割り当てる(近い球がなければ中点を新しい中心に、:203–206)。
- **最終半径は動かさない中心から測る:** 最終反復(`it+1==iters`)では中心を更新せず(:215–217)、半径は「その反復で距離を測った `balls[j].v`」に対する √dst+xnorm+1e-9 の最大(:209)。Lem 6.4 と一致。最終反復で新設された球は中点が中心で dst=0。
- **`exist` の正当性(Lem 6.8 の証明をコードで追った):** 配置 (A,B) があり a=min(A∪B)∈A とする。根では試した頂点を o と u の両方から消す(`n==6` のときだけ、:433)ので、i=a の時点で A∖{a} と B は全部残る。以後の再帰では o だけから消し、pool は u∩N_U(i) で縮む。枝を捨てる 3 条件(`pool.count()<6`、`v.count()<n-1`、`!clique(pool,6)`)は、B⊆pool かつ A の残りが v にある限り成立しない。`clique` は最小頂点で枝分けする全列挙。`tick` の超過は true を返して `gaveup` を立て、`candidate` は `q || gaveup`。
- **`joint` の符号:** z=x₀−x₁、Q=e(u₀−u₁)。x₀ の係数は J⁰ᵀy₀+Dᵀy_e、x₁ は J¹ᵀy₁−Dᵀy_e(:354)。`support` を先に呼ぶので G_y は辺の係数と相殺されない(:353、exclusion.tex:525–532)。zB は 2 球ぶんの −Bᵀy を累積してからノルム(:361–362)= 共有する 1 つの η。J、B、D_O、D_U の符号と T の添字の転置(`TT[k][j][i]`)は有限差分で確認(スクリプト 6e, 6g, 6h, 6i)。
- **生の `R.p`, `R.T` を使うこと:** C 段は記録の p と T をそのまま使い(stageC.cpp:500–501)、A/B 段の距離は dephase して wrap した差の射影で測る。両者をつなぐのは整列式 (phase:alignment)(phase-cover.tex:403–422):d=Π(D₀(P)−D₀(Q)−N) なら e(Q+d) は e(P) と行・列の位相倍で同値。記録の意味は「e(P+δ)、δ 水平、‖δ‖≤r、F(P+δ)=0」(Def 4.14)で、C 段の仮定 (exclude:root-decomposition) と同じ。

### 3.5 段の間の受け渡し

| 受け渡し | 書く側 | 読む側 | 欠落の可能性 |
|---|---|---|---|
| centers.txt | covout.cpp:141(件数と行は同じ vector から。書き込みエラーは abort) | atlas.cpp:173–187(n==112、件数 1..2e6、各値の範囲、fscanf 失敗は abort) | 途中で切れたファイルは abort。A の 4 分割は `i=part; i+=parts`(atlas.cpp:227)で添字を分割 |
| coarse_i.bin → mid.bin → fine.bin / fb1.bin | `saveRecords` atlas.cpp:191–202(fwrite 失敗・fclose 失敗は abort) | `readRecords` :203–215(件数 ≤2e6、k∈[0,10]、fread 不足は abort)+ `validate` stageB.cpp:274–282 | 読み込み後に余りバイトがあっても検出しない(§6 #9)。件数より短ければ abort |
| バッチ f<j>.bin | stageB.cpp:300–311(挿入順に 256 件、半径を target に置換、tmp→rename) | crun | 置換は拡大(半径 ≤ target を :327 で全件再検査)。バッチの網羅は `fdone.txt` と外部の検査に依存(C++ 自体は確認しない) |
| bad<j>.bin | stageC.cpp:551, 564(読んだ記録をそのまま保存) | brun fallback | 15 ファイルを全部渡すことは実行手順側の義務 |
| 最終 C の 4 分割 | — | stageC.cpp:550 `(offset+i)%parts==part` | 分割は添字の剰余で漏れなし |

**記録された呼び出しの検査(`check_266cpp_handoff.out`、16 項目全部 OK):** 引数が 1 つ違うだけで記録が黙って飛ばされる(例: あるバッチだけ `parts` が 2)ので、100 件全部を機械的に確認した。

- 初回 C の 88 件は argv[0:10] が全部同一(`crun 12 4 1 .042 1 1 5000000 0 1`)。バッチ j は入力が `batches/f<j>.bin` ただ 1 つ、出力が `batches/bad<j>.bin`。
- 各バッチで、ログの `done` = `f<j>.bin` の件数、ログの `bad` = `bad<j>.bin` の件数、失敗 4 分類の和 = bad、終了コードは bad>0 のとき 1。
- `fdone.txt` = (88, 22379) = バッチ件数の和 = `fine.bin` の件数。初回 C の bad は合計 163。
- fallback の入力は、件数が 0 でない 15 個の bad ファイル全部で、C ロケール順。親 163、出力 1127 = `fb1.bin` の件数、bad 0。
- 最終 C は part i / 4(i=0..3)で `fb1.bin` を読み、`done` は各剰余類の大きさ(282, 282, 282, 281)、和 1127、bad 0。
- A は part i / 4 で `centers.txt` を読み、`done` は剰余類の大きさ(和 249,878)、`failedBalls0`、`global` = `coarse_i.bin` の件数(和 3083)、`bad_i.txt` は 0 バイト。
- Bmid は親 3083 → 600、Bfine は親 600 → 22,379、どちらも bad 0・終了コード 0。引数は execution.tex:127–129 のとおり。

`results/*.txt` は配布者が整形したログの抜粋で、生の stdout ではない(VERIFY.md:3–8)。バッチの連結が `fine.bin` と半径以外で一致すること、bad が入力の順序付き部分列であることは、前回の監査(`check_266_records.out`)に依拠し、今回はやり直していない。

### 3.6 NaN/∞ と min/max

tex は「min/max は未検査の NaN を小さい証明書に変えるためには使われていない」と述べる(rounding.tex:1519–1527)。排除・被覆を決める比較を全部見て、NaN のときは保持・失敗側に倒れることを確認した(§3.1–3.4 の表の「失敗時」)。例外的な形が 2 つある。

- `bootstrap`(patchx.cpp:284)と `eval`(stageB.cpp:56)の `if(z>=s)break; s=z;` は、z が NaN だと s が NaN になる。ただし呼び出し側の比較(`s<=target` など)が偽になり被覆は成立しない。`setRadius` は abort。
- `markRec`(stageB.cpp:147)の `if(q>… || rho<5e-8)continue;` は q が NaN だと印を付ける側に倒れる。q は有限の `eta`, `w`, `L` からしか作られず、`eta` が NaN なら `nu0` が NaN → `eval` が NaN → ρ=0 → `mark` が :152 で即 return するので、**到達不能**。

---

## §4 浮動小数点とビルド条件

### 4.1 フラグと付録 A の前提

| 付録 A の前提(rounding.tex:23–49) | build.sh / ソース | 判定 |
|---|---|---|
| binary64、中間結果も binary64 | x86-64 の SSE2 が既定。`static_assert(FLT_EVAL_METHOD==0)` はプローブだけ(abi-fp-probe.cpp:28)。本体は `is_iec559 && digits==53`(seed.cpp:14) | 整合。本体は FLT_EVAL_METHOD を直接は検査しない |
| 最近接丸め | `fesetround(FE_TONEAREST)` と確認(seed.cpp:15)、A/B/C は MXCSR の丸めビットも検査(atlas.cpp:10–11) | 整合 |
| FMA・縮約なし | 4 本とも `-ffp-contract=off`(build.sh:22–25)。ソースに `fma` の呼び出しなし(抽出スクリプトで 0 件) | 整合。GNU 方言の既定は `-ffp-contract=fast` なので、この明示は必須で、実際に付いている |
| fast-math・再結合なし | `-ffast-math`/`-Ofast` なし。プローブは `__FAST_MATH__` を `#error`(abi-fp-probe.cpp:18–20) | 整合 |
| FTZ/DAZ なし、例外マスク | A/B/C は `checkcpu` で検査。**covrun は検査しない**(covout.cpp:136 は `environment()` のみ) | tex も明記(rounding.tex:43–45)。§6 #3 |
| `main` 前の定数も最近接 | `const double sqrt6=sqrt(6.)`(ub.cpp:24)、`Basis bas`(patchx.cpp:4–10)、`inv=1./sqrt(6.)`(stageC.cpp:5)が `main` の前に評価される | tex が前提として明記(rounding.tex:45–46)。`checkcpu` はその後に走る。§6 #4 |
| sqrt が正しく丸められる | プローブが √2 の 1 例を検査(abi-fp-probe.cpp:44–46) | 前提 |

- **`-O3 -march=native`:** `-march=native` は FMA や AVX の命令を使用可能にするが、縮約は `-ffp-contract=off` が禁じる。GCC は `-ffast-math`(`-fassociative-math`)なしに浮動小数点の和の順序を変えるベクトル化をしない、というのが文書化された挙動。**バイナリは配布されておらず、逆アセンブルも見ていない**ので、これは GCC の仕様に基づく判断で、実物の確認ではない。
- **`long double`、`float`:** ソースに無い(抽出スクリプトで 0 件)。x87 は使われない。
- **`#pragma STDC FENV_ACCESS ON`**(seed.cpp:9): GCC はこの pragma に未対応と文書化されている。影響があるのはコンパイル時の定数畳み込みの丸めだが、GCC の畳み込みは最近接で、実行時も最近接を検査しているので食い違わない。tex に記述なし(§6 #7)。
- **covrun だけ `-std` 指定がない**(build.sh:22)。execution.tex:72 が「この コンパイラの既定の GNU C++17」と明記。

### 4.2 数学ライブラリ関数の全使用箇所(`extract_266cpp_literals.out`)

| 関数 | 箇所 | 証明での役割 | tex の扱い |
|---|---|---|---|
| `pow` | covout.cpp:100, 101 の 2 箇所だけ | :100 `RN2=pow(R*n,2)-.0001` は被覆判定の閾値。:101 の `T` は試行の選別のみ | 乗算に落ちることを前提(rounding.tex:348–354、execution.tex:90–104)。§6 #2 |
| `sqrt` | 105 箇所 | 上界・ノルム全般 | 正しい丸めを前提(rounding.tex:29–31) |
| `nearbyint` | seed.cpp:23、patchx.cpp:38, 195, 218、atlas.cpp:13、stageC.cpp:171 | seed.cpp:23 は `expair` の象限還元。他は位相の整数リフト | 「整数を 1 つ選ぶ」(rounding.tex:51–53)。ただし Lem A.2 の証明は最近接を使う(§6 #5) |
| `floor` | atlas.cpp:14(`turn`、並べ替えとバケツのキー)、stageB.cpp:155(印を試す範囲) | 選別のみ | 同上 |
| `ceil` | atlas.cpp:122、stageB.cpp:156, 221(選別)、**stageB.cpp:251**(L) | :251 は格子の被覆に効く(ceil(x)≥x が必要) | rounding.tex:1131–1135 |
| `atan2` | patchx.cpp:189 の 1 箇所 | 位相の**提案**。直後に `\|ex(P)−D/m\|<1e-7` で検査(:202) | 精度の前提なし(rounding.tex:54–55) |
| `copysign` | patchx.cpp:25 | Jacobi 法(提案のみ) | — |
| `isfinite` | ub.cpp:238、patchx.cpp:83, 119、stageC.cpp:113, 341 | 提案の検査 | — |
| `sin`/`cos`/`tan`/`exp`/`log`/`hypot`/`fma` | **0 件** | cos/sin はソース内の Taylor 級数 `expair`(seed.cpp:21–31)、ノルムは `hp`=`sqrt(x*x+y*y)`(seed.cpp:32) | Lem A.2 |
| `atof` | stageB.cpp:288(target, w)、stageC.cpp:537(RB)。ub/patchx の main は経路外 | 下記 | 「10 進リテラルは保存された厳密なパラメータとして使う」(rounding.tex:55–57) |
| `atoi`/`strtoull` | atlas.cpp:220、stageC.cpp:536–539 | 分割の添字、層数、上限 | — |

- **`expair` の誤差:** ソースを Python の binary64 演算に写し、mpmath(40 桁)の e(x) と 4,369 点で比べて最大 3.4e-16。tex の上界 2e-14(rounding.tex:187–189)の内側(スクリプト 1)。これは標本であって証明ではない。
- **`atof` と実行時の閾値:** B 段の target は `Atlas A(target,…)` の `maxR` として保存され、以後の比較は全部その保存値に対して行われ、バッチの半径にもその double がそのまま書かれる(stageB.cpp:298, :305)。C 段はファイルの double を読む。10 進表記と別の表現を突き合わせて被覆を結論する箇所は見つからなかった。したがって `atof` の変換精度は健全性に効かない。A 段の `.0320001` はコンパイル時リテラルで、.032+1e-14 以上であれば足りる。

### 4.3 未定義動作の探索

読んだ範囲で、証明経路上に未定義動作は見つからなかった。確認した点:

- **配列の範囲:** `trig[8]` の添字は N=16..512 で 0..5(ub.cpp:162。N≥512 では再帰しない、patchx.cpp:433)。`picks[8]` は `topK=4` まで。`counts[8]` は l<3。`BS::data[32]` は球 ≤2048 個(stageC.cpp:204 が 2048 個目の追加前に止める、:404 で words>32 は abort)。`int low[maxK]`。バケツ添字は `min(bins-1,max(0,…))` で切り詰め(atlas.cpp:116)。
- **未初期化:** LDL/solve は下三角だけ読む(ub.cpp:25–31、stageC.cpp:19–31)ので、下三角だけ埋めた行列で足りる。`covercompute` の作業配列、`dataU`、`G[maxK][maxK]`(patchx.cpp:324)は読む前に k×k 部分が書かれる。`Root`、`mapH best{}`、`Candidate` は既定メンバ初期化子で 0。`Usearch::bestres` は未初期化だが経路外の `Usearch::test` でしか読まれない。
- **`mapH best{}` が全 0 のまま返る可能性:** `canonicalMap` は少なくとも 1 回コールバックされる(tol≥0 なので最大・最小が必ず採用される、phase-cover.tex:955–964)。入力の p が有限であること(`gen` の `\|p\|≤4`、`validate`)が前提で、満たされている。
- **符号付き整数:** 距離の二乗和は ≤ 9·111²。格子総数は乗算前に 5e7 で検査(stageB.cpp:254)。`hash` は uint64_t。double→int は範囲が保証された値だけ(seed.cpp:23 は `\|x\|≤16` の後、stageB.cpp:251 は L≤1262)。
- **型のエイリアシング:** 本体には無い(`reinterpret_cast` はプローブの 1 箇所、`unsigned char*` 経由で適法)。
- **`abs`:** `using namespace std` と `<cmath>` のもとで double 引数には `std::abs(double)` が選ばれる、と読んだ(コンパイルして確かめてはいない)。

---

## §5 定数の照合表

抽出スクリプトが拾った**浮動小数点リテラル 91 種**を、読んだ文脈で種類ごとに分類した。411 箇所それぞれの行番号は `extract_266cpp_literals.out` にあり、下の表は種類ごとに主な行を挙げる(同じ値が P と H の両方で使われる場合は両方の表に出る)。**表に行番号を挙げていない出現箇所もある**ので、「411 箇所を 1 行ずつ表にした」わけではない。分類: **P**=証明に効く余裕・閾値(tex の行つき)、**G**=失敗・拒否側のガード(通らなければ未解決)、**H**=提案・選別だけ(値は健全性に無関係)、**X**=証明経路外。`0.`、`.5`、`1.`、`2.`、`4.`、`6.` などの構造的な値は式の一部として §2–§3 で照合済み。

### 5.1 証明に効く定数(P)— tex と全部一致

| 値 | ソース | 役割 | tex |
|---|---|---|---|
| `0x1.921fb54442d18p+1` | seed.cpp:12 | π̂ | rounding.tex:151–155 |
| 16. / .25 / 次数 9 | seed.cpp:22, 24, 26 | `expair` の範囲、象限、18/19 次 | rounding.tex:188, 198–213 |
| 112、.032、.0320001 | covout.cpp:138, 140、atlas.cpp:175, 221 | N、seed 半径 | seed.tex:173, 427, 457 |
| 5.(×πw²/4)、1e-8 | seed.cpp:41, 46, 88、covout.cpp:122 | Taylor 定数、gap と sumok の余裕 | seed.tex:269、rounding.tex:276–290 |
| 5.(×πw/4)、1e-6、1e-7、1e-5 | covout.cpp:12, 16, 19, 20 | E、除数の下限、分子と端点の余裕 | seed.tex:302, 338、rounding.tex:309–316 |
| 1.251、1e4、1e10、.01 | covout.cpp:47, 63, 49, 71 | 双対の係数、提案の上限、初期値、支持の余裕 | seed.tex:422–424、rounding.tex:319–339 |
| 2.25、.0001 | covout.cpp:82, 85, 100 | 9/4、RN2 の引き下げ | seed.tex:391、rounding.tex:343–350 |
| 1.25、.45 | ub.cpp:41, 61 | 枠の受理 | rounding.tex:374, 398 |
| 1e-9、1e-10、1e-8、1e-7 の列 | ub.cpp:132–146 | `get` の各フィールドの余裕 | rounding.tex:454–494 |
| 5.9 | ub.cpp:136 | q の受理(tex は q<6) | completion.tex:155, 198(厳しい側。§6 #8) |
| 6·42 | ub.cpp:143 | √252/2 | completion.tex:439 |
| 1e-5 | patchx.cpp:424 | A のノルム判定 | rounding.tex:451–452 |
| .0001、1e-100/1e100 | ub.cpp:250, 239 | A の方向判定、方向の範囲 | rounding.tex:497–514 |
| 4.、1e-18 | patchx.cpp:51, 110 | \|p\|≤4、val | rounding.tex:522, 536 |
| 2.(因子)、2.(mnorm)、8.、1e-8×3、1e-7 | patchx.cpp:130–143 | 因子の範囲、mnorm、C=π√8·mnorm | rounding.tex:540–572、phase-cover.tex:315 |
| .1、.0002、.00001、1e-8 | patchx.cpp:161, 164 | テンソル改良 | rounding.tex:596–628 |
| 2e-6、1e-8×2、1e-9×2 | patchx.cpp:182–183 | id、ort、hor、froot、n0 | rounding.tex:630–667 |
| 1e-4、1e-5×2、20 | patchx.cpp:184 | Root の受理 | rounding.tex:668–671 |
| .15、.51、1e-7 | patchx.cpp:200, 202 | 大きさ、位相の範囲、弦 | rounding.tex:264, 717–727 |
| 1e-6×2、1e-9 | patchx.cpp:205 | mi、invmag | phase-cover.tex:498–500 |
| 1e-5、1e-10、1e-8、3.8、1e-7、1e-8×3 | patchx.cpp:207–213 | slack、q、ratio、arc、weighted、d2 | phase-cover.tex:491–502 |
| 1e-10、2e-7 | patchx.cpp:230 | spectralUpper | rounding.tex:859–877 |
| 10.、100.、1e-8、1e-9、1e-6、1e8、1e-10×3、.000001×2 | patchx.cpp:244–275 | `ell` のガードと余裕 | rounding.tex:884–968 |
| 10.、3、20、.00002、1e-6、30、1e-8 | patchx.cpp:279–283 | bootstrap | rounding.tex:971–987 |
| .0001、.00001、1e-7 | patchx.cpp:307–308 | s0 | phase-cover.tex:640–649 |
| 1.45、.45、1.8、.8、1e-7 | patchx.cpp:336, 343 | 楕円体の β | phase-cover.tex:717–741、rounding.tex:823–846 |
| .00015、.00002 | patchx.cpp:338, 346 | 接上界の外部誤差 | phase-cover.tex:838–845 |
| 1e-8 | patchx.cpp:389–390, 405–406 | 局所半径 | phase-cover.tex:894–895 |
| .103、.125、levels=3、iters=2、first=1 | atlas.cpp:221–222 | A の目標・上限・レベル | execution.tex:113–115 |
| 1e-7、1e-8、2e-8 | atlas.cpp:140, 160–161、stageB.cpp:27–28 | 距離と合併 | phase-cover.tex:966–974、rounding.tex:703–709 |
| 2e-8、1e-9、1e-8×6 | stageB.cpp:47–59 | `eval` | refinement.tex:254–257、rounding.tex:1073–1109 |
| 1e-8×4、1e-7×2 | stageB.cpp:77–81, 98 | `init` の各上界、sig | rounding.tex:1021–1071 |
| 1.01、1e-7、25、2e-7、1e-8 | stageB.cpp:101–105 | ρ 探索 | refinement.tex:276–288 |
| 1e-8(corner、outside、markRec)、1e-7(E)、2e-8、5e-8 | stageB.cpp:109, 126, 124, 140, 145, 147 | 格子と印 | refinement.tex:333–337, 358–360, 398–408 |
| .500001 | stageB.cpp:251 | L | refinement.tex:337 |
| 1e-8、2e-8、1e-7 | stageB.cpp:133, 173, 184–185 | 半径更新、隅、直接割り当て | refinement.tex:431–433、rounding.tex:1187–1189 |
| .126、20、1e-6、1e-4、1e-5、2、4、1.01 | stageB.cpp:274–281 | `validate` | refinement.tex:603–618 |
| .0001、.126、.05 | stageB.cpp:289 | 引数の範囲 | refinement.tex:603–607 |
| 5e7、6、2e6、8e6、256 | stageB.cpp:254, 117, 31, 293, 302 | 有限の上限 | refinement.tex:623–626 |
| 1e-12×4、5−1/6、5、1.5 | stageC.cpp:43–44 | a、X_M、L_M、I_M | exclusion.tex:98–100, 674–676 |
| 2e-7×5 | stageC.cpp:63, 70, 85, 138 | 単一ベクトル判定 | exclusion.tex:682–688 |
| 2e6、2048、2.、.1、1e-9 | stageC.cpp:35, 204, 219, 220, 209 | 上限、中心と半径の検査、半径の余裕 | exclusion.tex:304, 341–354 |
| 1e-10、2e-7×5、2e-8×2、3e-7 | stageC.cpp:371, 378–397, 388, 396, 364 | 辺と結合判定 | exclusion.tex:690–697 |
| 1e-100/1e100 | stageC.cpp:114, 342 | 方向の範囲 | exclusion.tex:291–294 |
| 1.00001、.03、1e-9×2 | stageC.cpp:502, 497, 504–505 | K、r≤.03、r_H、n | exclusion.tex:30, 77, 677–680 |
| 5e6、32 語 | stageC.cpp:303, 282 | グラフの上限、ビット集合 | exclusion.tex:631–637 |
| 1e6、512 | patchx.cpp:418, 433 | A の追加ノード上限、最深格子 | completion.tex:702–706 |

### 5.2 失敗側のガード(G)

`.092`(patchx.cpp:315、`toff>.092²` で失敗の番兵 10. を返す。**tex に記述なし**)、`.145`(patchx.cpp:293、rounding.tex:762)、`.249`/`.23`(stageB.cpp:64、refinement.tex:638–639)、`1e-6`(stageB.cpp:73)、`2.`(stageB.cpp:95)、`1e-7`/`2e-7`/`3e-7`(atlas.cpp:140, 150–151、stageB.cpp:18–19, 152, 205 の margin と ρ の下限)、`1e-14`(patchx.cpp:238、Cholesky の pivot)、`.012`/`.09`(stageC.cpp:541、RB の範囲)、`1e10`(ub.cpp:136、不正な q の番兵)、`10.`(stageB.cpp:46、`eval` の範囲外の番兵。patchx.cpp:279, 293, 315 も同じ番兵)、`1024`(stageC.cpp:510、M の範囲)、`6`/`24`/`10`/`100`/`1e9`(stageC.cpp:540–543、引数の範囲)。どれも通らなければ被覆不成立・未解決・abort になる向き。

### 5.3 提案・選別だけの定数(H)— tex に値の記述がないものを含む

`.2`、`.02`、`.03`(covout.cpp:48, 67–68、supdual の重み)、`.51`、`1.5`(covout.cpp:101–102、T と K)、`1e-5`、`.2`、`.02`、`.01`、`.001`、`1e-30`、`1e80`、`1e-100`(ub.cpp:254–298、`trial::pass` の重み)、`1e-13`、`1e-12`、40(patchx.cpp:19–34、Jacobi)、`1e-20`、`.12`、`.25`、`.02`、24、5(patchx.cpp:81–106、`snap`)、`.64`(patchx.cpp:122、固有値の分割。phase-cover.tex:313 に記載)、`4.35`(patchx.cpp:144、テンソル改良を試す閾値。**tex に記述なし**)、`.03`、`2e-6`、`.28`、`.58`、`1.`(patchx.cpp:249–250、γ の提案)、`.145`、`.15`(patchx.cpp:377, 379、近い Root の選別)、`1e-8`(atlas.cpp:92–100、署名の比較)、`.08`、`.061`、120(atlas.cpp:153、stageB.cpp:21, 240, 317、写像の候補)、`.00001`、`2e-6`(atlas.cpp:122、stageB.cpp:221、バケツの到達範囲)、`.03`(stageB.cpp:298、バケツ数 16/8)、`1e-6`、`1e-5`、`1e-300`(stageC.cpp:108, 139–141、重み)、`1.11`、`1e-10`、`.042`(stageC.cpp:183、クラスタの候補閾値。exclusion.tex:357–360 に記載)、`.000001`、`1e-6`(stageC.cpp:336–363、結合判定の重み)、`1e-18`(patchx.cpp:106、`snap` の戻り値の判定。証明書側の `val<1e-18` は :110 で別に検査)、`1e10`(patchx.cpp:196、最小値の初期値)。

**式の一部である値(P でも H でもない):** `.25`(ub.cpp:48、基底ベクトル p−S/4)、`1e-9`(patchx.cpp:309、k=0 のとき `bootstrap` に渡す正の接入力。数学的には 0 でよく、正の値は上界を大きくする側。phase-cover.tex:881–882)。

### 5.4 証明経路外(X)

`1e-10`×2、`1e-6`、`5.`(seed.cpp:58–65、`share`)、`2.0001`(seed.cpp:103、`list`)、`1e-12`(ub.cpp:230、`derr`)、`1e-7`×2(ub.cpp:138, 142、`zK`)、`.22`(stageC.cpp:455、`duo`)、既定値で実行時に上書きされるもの: `.04`(stageC.cpp:227、RB)、`.126`(atlas.cpp:114、Atlas の既定上限。A は `.125`、B は target を渡す)。プローブの 3 つの 16 進定数(abi-fp-probe.cpp:45–46)。

**照合の結果:** P の定数はすべて、tex が述べる値・位置・符号と一致した。食い違いは 0 件。tex に値の記述がない定数(`.092`、`4.35`、`5.9`、および H の多く)は、失敗側か提案専用で、健全性に効かない。

---

## §6 番号付きの指摘

**健全性を破るバグ: 見つからなかった。**

| # | 分類 | 重大度 | 場所 | 内容 |
|---|---|---|---|---|
| 1 | checked OK | — | 全段 | §3 の全判定(seed 6、A 8、B 8、C 10)が tex の必要条件と同じ向き、余裕は外向き、失敗は未解決。ループは tex の全体集合を尽くす。段の受け渡しで領域を落とす経路は見つからなかった。 |
| 2 | 健全性に影響しうる未確認点 | 低〜中 | covout.cpp:100、rounding.tex:348–354、execution.tex:90–104 | `RN2=pow(R*n,2)-.0001` が乗算に落ちることを tex は前提にし、逆アセンブルで確認したと述べるが、逆アセンブルは配布されていない(前回 #4)。**今回の補足(私の評価):** この比較で必要なのは「RN2 が真の (R·n)² を上回らない」ことで、引き下げ幅は 1e-4、真値との差は乗算なら 1e-12 以下。`pow` が乗算に落ちなくても、その誤差が 1e-4 未満なら不等式は保たれる。したがってこの前提への依存は、tex の書き方が示唆するほど脆くはない。ただし「ライブラリの `pow` の誤差が 1e-4 未満」も証明ではなく仕様上の期待である。:101 の `T` は試行の選別だけ。 |
| 3 | 健全性に影響しうる未確認点 | 低 | covout.cpp:136、rounding.tex:43–45 | covrun は `environment()`(丸めモード)だけを呼び、MXCSR の FTZ/DAZ・例外マスクを検査しない。tex は「同じ環境が必要だが entry point は検査しない」と明記。seed 段の量は 1e-6 以上のスケールで、アンダーフローが効く演算は読んだ範囲に見当たらないが、前提であって検査ではない。 |
| 4 | 健全性に影響しうる未確認点 | 低 | ub.cpp:24、patchx.cpp:4–10、stageC.cpp:5、rounding.tex:45–46 | `sqrt6`、`Basis bas`(Helmert 基底、`sqrt` と除算)、`inv` は `main` の前に初期化され、`checkcpu()` はその後。tex は前提として明記。`bas` は Root の証明書の計算(`gen` の J、`patch` の hvec)に入るが、id の検査は 36 座標の完全な恒等式を `PI` で直接評価するので(patchx.cpp:168–176)、`bas` の誤差は提案の質にしか効かない、と読んだ。 |
| 5 | tex に明記のない前提(無害と判断) | 低 | seed.cpp:23、rounding.tex:51–53 対 198–199 | 前提の節は `nearbyint`/`floor` について「有限の引数で整数を 1 つ選ぶ。境界でどの整数でもよい」とだけ述べるが、Lem A.2 の証明は「k は 4x の最近接整数、ゆえに \|x−k/4\|≤1/8」を使う。つまり `expair` については `nearbyint` が最近接整数を返すことが必要。同様に stageB.cpp:251 の `ceil` は引数以上の整数を返すことが必要。どちらも IEEE で厳密に定義された演算で、実装が外すとは考えにくいが、前提の一覧(加減乗除と sqrt)には入っていない。 |
| 6 | 健全性に影響しうる未確認点 | 中 | rounding.tex 全体 | 付録 A の誤差予算(「この経路の敵対誤差は X 未満」)の**数値そのもの**は再導出していない。今回確認したのは「ソースの余裕の値・位置・符号が tex の記述どおり」と、スケールの主張のいくつか(R_s<.203、q<.32、hor·√(R_s²+arc²)<8.81e-5、K>√(1+…) など、スクリプト 8a–8d)。余裕 ÷ 主張された誤差の比を表(rounding.tex:1529–1566)から計算すると、大半は 30 倍〜10⁶ 倍だが、小さい行がある: **Candidate の接上界 1.25 倍**(.00012 対 .00015 rad)、**A のノルム判定 10 倍**(1e-6 対 1e-5)、**第 2 の G の成分誤差 10 倍**(k≤10 個 × 1e-9 対 1e-7·I)、A の方向判定 25 倍、C の結合判定 30 倍、Root の恒等式 40 倍。接上界については、表の .00015 rad とは別に、ソースは 2π で割った後に .00002 サイクル(≈1.26e-4 rad)を足し(patchx.cpp:338, 346)、`bootstrap` がさらに .00002 サイクルを足す(:280)。tex はこの 2 つを特定の誤差項に割り当てておらず(phase-cover.tex:844–845、rounding.tex:819–821, 975–976)、追加の余裕になっている。合計すると (1.5e-4+2×1.26e-4)/1.2e-4 ≈ 3.3 倍。また .00012 の内訳として tex が挙げる項は 8.81e-5+1.27e-5+3.4e-6 ≈ 1.04e-4 で、最大項 8.81e-5 は q<3.8、hor<1e-5 から再計算して確認した(スクリプト 8c)。予算が 1 桁ずれると破れうるのは、この接上界と A のノルム判定・第 2 の G の 3 行。 |
| 7 | tex に記述のない処理(無害) | 低 | seed.cpp:9 | `#pragma STDC FENV_ACCESS ON`。GCC は未対応と文書化している。実行時に最近接を検査しているので、定数畳み込みとの食い違いは生じない。 |
| 8 | tex との食い違い(無害) | 低 | ub.cpp:136、patchx.cpp:315, 144 | ソースにあって tex に値の記述がない定数: q の受理 `5.9`(tex は q<6。厳しい側)、`toff>.092²` の拒否(失敗側)、テンソル改良を試す閾値 `C>4.35`(試さなくても第 1 の C が有効)。他の H 分類の定数(§5.3)も同様。 |
| 9 | tex との食い違い(無害) | 低 | atlas.cpp:203–215 | `readRecords` は宣言件数ぶん読んだ後、余りのバイトがあっても検出しない。件数より短いファイルは `fread` 不足で abort するので、領域が落ちる向きにはならない。EOF までの一致は `check_records.py` の役目(VERIFY.md。今回は読んでいない)で、前回の自作検査が確認済み。 |
| 10 | tex に記述のない処理(無害) | 低 | atlas.cpp:140 | `Atlas::find` の成功条件 `margin>3e-7 && s<(margin-2e-7)²` は tex に式として出てこない。合併の機会を減らすだけで、合併の正当性は `add` の `radius+dist+2e-8<maxR` が別に検査する。 |
| 11 | checked OK(到達不能を確認) | 低 | stageB.cpp:147 | §3.6。NaN が印を付ける側に倒れる唯一の形だが、到達不能。 |
| 12 | checked OK | 低 | stageB.cpp:188 | `addNear` の戻り値を無視。被覆は直前の `assign` で確立済み。tex も明記。 |
| 13 | checked OK | 低 | stageC.cpp:433 | 根でだけ `u.del(i)`。「和集合の最小頂点が第 1 クリークに属する」という規約を課すだけ(Lem 6.8 の証明と一致)。 |
| 14 | checked OK | 低 | stageC.cpp:207–221 | 最終反復の半径は動かさない中心から測る(Lem 6.4)。 |
| 15 | checked OK | 低 | stageB.cpp:288、stageC.cpp:537 | `atof` で読む閾値は保存された double として一貫して使われ、10 進表記との突き合わせで被覆を結論する箇所はない(§4.2)。 |
| 16 | tex との食い違い(無害、既報) | 低 | execution.tex:16 | ソースの所在を `verification/code/` と書くが、配布物では `verification/saved-run/compiled-build/source/`(前回 #6)。 |
| 17 | 範囲外の義務(C++ の外) | 情報 | execution.tex:117–122, 135–136 | バッチの網羅(`fdone.txt`)、bad ファイルを全部 fallback に渡すこと、最終 C の 4 分割を全部走らせることは、C++ ではなく実行手順(未配布の `proof_run.py`)と事後検査の義務。C++ 側は「未解決があれば exit 1 と bad の保存」までを正しく行う。 |

---

## §7 監査できていない部分

- **コンパイルした結果。** ソースを読んだだけで、オーバーロード解決(`abs`)、テンプレートの実体化、最適化後の命令列は確かめていない。バイナリも逆アセンブルも配布されていない。`pow` の降格(#2)、縮約なし、ベクトル化で和の順序が変わらないことは、GCC の文書化された挙動に基づく判断である。
- **付録 A の誤差予算の数値**(#6)。依頼の範囲外とした。
- **実行の再現。** 禁止事項。したがって「この C++ を走らせると配布された記録が出る」ことは確認していない。確認したのは「C++ が tex の補題どおりに書かれている」ことだけ。
- **`check_records.py`、`lean/`。** 読んでいない。`commands.json` は自作スクリプトで全 100 件を検査したが、目で全行を読んではいない。バッチの連結と bad の部分列性は前回の自作検査に依拠。
- **introduction.tex、conclusion.tex。** grep のみ。ルーチンの記述は無かった。
- **提案アルゴリズムの品質。** `snap`、`jacobi`、`chol`、`trial::pass` の重み更新、Lloyd 反復は「提案にすぎない」ことを確認しただけで、収束や精度は見ていない(健全性には不要)。
- **私の読みの見落とし。** 2,316 行を 1 人で 1 回読んだ。主要な式(A の線形化、Root の残差と微分、C の J・B・辺の微分、seed の J_*)は数値で裏を取ったが、制御フローは目で追っただけで、別の読み手や別言語での再実装による照合はしていない。数値比較は「tex の式」と「私がソースから写した式」の一致を見るもので、写し間違いがソースの誤りを隠す可能性は残る(写しは行ごとに元の式の形を保った)。

---

## §8 確信度

**「配布された C++ が、付録 A の算術前提のもとで、tex の補題を健全に実装している(排除できていない領域を排除済みと記録する経路がソースレベルで存在しない)」: 約 90%。**

根拠:

- 8 本 2,316 行を全行読み、tex の 7 本 5,550 行を全行読んで、名指しされた全ルーチンと全判定を 1 対 1 で突き合わせた。食い違いは無害なもの(§6 #7–#10, #16)だけだった。
- 健全性に効く判定(§3.1–3.4 の表の 32 行)すべてで、不等号・余裕の符号・NaN の倒れ方・失敗時の行き先を個別に確認した。記録された 100 件の呼び出しの引数と件数も機械的に確認した(§3.5)。
- 符号や添字を取り違えやすい式(線形化の 13 列、Helmert 基底での Jacobian、C 段の J・B・辺の微分、seed の J_* の行配置、U(2) の微分)を、tex から独立に組んだ式と数値比較し、41 項目が全部 OK だった。
- 証明に効く定数は tex と全部一致した(§5.1)。
- 設計が「提案と証明書の分離」で一貫しており、複雑な数値計算(固有分解、線形解法、クラスタ化)は誤っても未解決が増えるだけ、という構造をコード上で確認できた。

残る 10% の内訳(主観):

- 私の読み落とし(制御フローの稀な経路、C++ の意味論の誤解): 5〜6%。
- 数値比較で使った「写し」の誤りがソースの誤りと偶然一致する: 1% 未満。
- 付録 A の誤差予算のどこかが実際の誤差を覆っていない(#6。余裕が 10 倍以下なのは、A のノルム判定、第 2 の G、接上界(追加の余裕を含めて約 3.3 倍)の 3 行): 3〜4%。これは「C++ が tex どおりか」ではなく「tex の予算が正しいか」の問題だが、健全性には同じく効く。

**これは N(6)=3 の確信度ではない。** 上の数字は前回の報告の (B)「C++ が補題を実装している」に当たる。(C)「2026-10-04 の実行が付録 A の算術条件のもとで行われた」(ビルド、`pow` の降格、MXCSR、ハードウェア)は今回も配布物からは確認できず、前回の指摘 #2・#4 がそのまま残る。前回 #5 の「未監査の最大部分」は、今回の精読で「1 人の読み手が全行を読み、誤りを見つけなかった」状態に進んだ。形式検証や別言語での再実装・再実行に代わるものではない。

## 9. Codex による独立レビュー（2026-11-11）と訂正

依頼文 `audits/codex_review/prompt_w7_266cpp.txt`（「見落としたバグを探せ」という形で依頼）、返答 `review_w7_266cpp.md`。Codex は C++ 8 本と対応する tex を読み、保存データを独自の読取処理で確認した（コンパイル・実行はしていない）。結論は「指定された実行条件と有効な Record を前提に、健全性を破る具体的なバグは見つからなかった」。網羅性（seed、U(2)、B 格子、C 格子）、stageB.cpp:188 と stageC.cpp:433 の無害性、定数 16 個の照合は支持された。**ただし本報告の次の 3 つの記述は強すぎたので訂正する。** 3 か所とも主担当がソースで確認した。

1. **「NaN・失敗はすべて未解決側」は誤り。排除側に倒れる NaN の経路が 3 つある。**
   - `stageC.cpp:385` `if(disto<=r+2e-7){`: `disto` か上限が NaN だと比較が偽になり、直交の辺を張らない。辺の欠落はグラフの排除につながり得る。
   - `stageC.cpp:209` `next[ix].s=max(next[ix].s,sqrt(dst)+layer.xnorm+1e-9);`: 候補が NaN だと `max` が既存の小さい半径を残し得る。`:220` の最終検査はこの NaN を検出できない。
   - `stageB.cpp:147`: 部分和が NaN だと被覆済みの印を付ける向き。
   - Codex の評価: 指定経路では C の位相と半径が検査され、指数関数の引数は範囲内、平方和は有限なので、**有効な入力からこれらの経路に到達する例は見つからない**。論文も有限性に依拠している（rounding.tex:1500–1517）。正確な言い方は「提案計算の NaN を拒否し、残る証明計算の有限性は入力範囲から保証する設計」。
2. **「段間の読込で領域が落ちる向きにはならない」（§6 #9）は無条件には誤り。** `atlas.cpp:203–205` の `readRecords` は先頭の件数だけを信じ、件数ぶん読んだ後の余りを検査しない（ソース全体に EOF 検査は無い）。件数ヘッダが過小なファイル（例: 件数 0 の 8 バイト + 実際の Record 列）を渡すと、C は処理 0 件・bad 0 件で正常終了できる。ただし Codex は独自の読取処理で、配布された全 190 個のバイナリが宣言件数どおりで余剰なく EOF に達すること、fine の 22,379 件と 88 バッチの一致、未解決 163 件、最終 1,127 件の 4 分割（282、282、282、281）、最終 bad 0 を確認した。正確な言い方は「保存データの整合は確認した。C++ の読込関数単独では過小なヘッダを検出しない」。`commands.json` は、元の監督処理が終了コードを実際に検査したことまでは証明しない。
3. **未定義動作は「指定パラメータでは見つからなかった」に限定する。** `atlas.cpp:220, 227`: `part=1, parts=2147483647` は引数検査を通り、`i+=parts` が符号付き整数のオーバーフローを起こす。保存された実行は parts=4 なので影響しない。
- 浮動小数点: フラグは算術の前提と整合するが、実際の命令列とライブラリ内部は未確認。「long double が無いから x87 不使用」という説明は不十分（x86-64 では SSE が既定、というのが正しい根拠）。`pow` の「脆さは小さい」は感度の評価としては妥当だが、論文の乗算化の確認条件（execution.tex:90–104）を不要にする根拠ではない。
- tex に値の記載が無い定数（.092、4.35、5.9）は「成功する機会を制限する値で、指定経路の誤った排除を作らない」と書くのが正確。
- **確信度「約 90%」は公開文面から外す。** 誤差予算と生成コードが未確認のまま第三者が裏付けられる数値ではない。結論は「静的な読解と保存データの照合で、指定条件下で健全性を破るバグを発見しなかった」。
