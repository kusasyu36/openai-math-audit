**公開前に修正が必要です。全体をobliviousとする読みには反例となる記述があります。ただし、NCCへの条件付きの懸念自体は残ります。** コードは実行せず、TeXと原論文を確認しました。

1. **判定：誤り――記述どおりの機械全体がoblivious、とは言えません。**

   決定的なのは[07-resampling.tex:215](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/07-resampling.tex:215)です。数値を通常の符号・絶対値表現へ変換し、**先頭ゼロを除いて出力する**と明記しています。出力長が数値に依存し、238行では出力ヘッドがそのフィールドの直後に来ると説明しています。したがってヘッド位置も値に依存します。

   報告が挙げた箇所では、例外ソートの鍵は確かにアドレスです（[05-layers.tex:395](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/05-layers.tex:395)）。桁上げも固定幅の走査で実装できます（[08-assembly.tex:783](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/08-assembly.tex:783)）。しかし外部乗算器のoblivious性は保証されていません（[06-transforms.tex:453](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/06-transforms.tex:453)）。

   固定幅化や別の乗算器への置換で修復できる可能性と、**元のアルゴリズムが既にobliviousであること**は別です。

2. **判定：正しい。**

   固定有限アルファベット・固定本数のテープというモデルは[02-streams.tex:4](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/02-streams.tex:4)のとおりです。

   ヘッド位置が既知なら、各ステップを定数サイズの回路に展開できます。同じセルを再読する場合も、その都度記号を引き継ぐ節点を置けば、fan-outも定数に保てます。合計サイズは \(O(T)\)。約 \(10^{21}\) 本でも、入力長に依存しない固定数なら問題ありません。定数が巨大になるだけです。

3. **判定：一部正しい。**

   **AFKLとのモデルの対応は正しいです。** [AFKL 2019の定理1と§2](https://drops.dagstuhl.de/storage/00lipics/lipics-vol132-icalp2019/LIPIcs.ICALP.2019.10/LIPIcs.ICALP.2019.10.pdf)は、任意のBooleanゲートを許し、入力次数・出力次数の両方を定数に制限しています。(2)の回路は対象内で、\(O(n(\log n)^{1-\kappa})=o(n\log n)\)ならNCCと矛盾します。

   ただし、(1)が成立しないので、**報告の推論経路はそのままでは未完成**です。また「NCCの反証と同等」は不適切で、成立しても一方向の含意です。一方、構成全体が正しければ、次の局所的な経路からNCCの否定が従います。

4. **判定：正しい。**

   指定build内の全`.tex`と[references.bib](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/references.bib:1)を検索しました。`obliv`、AFKL各著者名、`network coding`、`k-pairs`に該当する記載・引用はありません。これは確認できますが、不記載自体は数学的誤りの証拠ではありません。

5. **判定：正しい――補題が正しい場合の、条件付きの衝突です。**

   [04-swap.tex:285](./data/openai_math/preprints/Integer-multiplication-below-n-log-n-September-23-2026/build/sections/04-swap.tex:285)で \(P=G=B=1\) とすれば、辺長 \(m=2^u\) のビット行列を
   \[
   O(m^2(\log m)^\tau),\qquad \tau<1
   \]
   で転置します。この部分の分割・再帰・走査は形状で決まり、全体の数値ルーチンにあった可変長変換を使いません。

   [AHJ+06 §6、定理6.2・6.3](https://www.cs.cornell.edu/~rdk/papers/AHJKRL05.pdf)のTM模倣定理は、確かに**2テープ**についての記述です。ただし固定 \(k\) 本への拡張は導けます。ページ長を \(m\)、各テープのブロック長を十分小さい \(\Theta(m/k)\) とすれば、同じ長さの時間区間で触る全テープの内容が一定数のページに収まります。それらを集め、区間を模倣し、書き戻す操作は \(O_k(1)\) 回のI/Oです。位置はobliviousな予定から決まります。

   よって時間 \(\beta\) の計算を \(O_k(\beta/m+m)\) 回のoblivious I/Oで模倣でき、上の補題なら \(o(m\log m)\) 回になります。これは定理6.2とNCCに反します。**この拡張部分は私の導出であり、原論文の定理の引用そのものではありません。**

公開文面への修正案として、[READMEの109行](./public_repo/README.md:55)のoblivious性に関する部分は、次の内容へ差し替えるのが適切です。

> アルゴリズム全体のoblivious性は確認できない。数値表現の変換には入力値に依存するヘッド移動がある。同じ時間上界のoblivious実装が得られれば、AFKL (2019) により無向k-pairsネットワーク符号化予想は否定される。また、chunk-swap補題には独立した条件付きの衝突がある。補題・構成の正しさは専門家による検証を要する。

総合すると、**報告は全体のoblivious性を読み込みすぎています。AFKLの適用モデルは合っており、局所補題への警戒は維持できますが、これを論文の誤りが確定したという扱いにはできません。**