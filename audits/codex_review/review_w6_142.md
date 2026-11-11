**依存関係と \(q=2\) の帰結は確認できました。主な修正点は、「未証明」と断定する表現と、補助素数を与えた場合の計算量の説明です。** 以下、142・029は各論文の `build/` を指します。ファイルの変更・コード実行はしていません。

1. **判定：正しい。**  
   根拠：142の [10-analytic.tex:4–18](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/10-analytic.tex:4) は、029の定理1.2を証明なしで引用すると明記しています。[同:103–111、160–169](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/10-analytic.tex:103) では、それを補助素数の上界に使い、その上界を定理1.1へ代入しています。本文に別の**無条件の、一様な多項式上界**は見当たりません。小標数分岐には不要ですが、一般入力に対する計算量保証には必要です。

2. **判定：一部正しい。**  
   **Hecke \(L\) 関数のGRHから従う、という数学的含意は正しい**です。非自明零点の実部が \(1/2\) なら、実部 \(>1-10^{-6}\) に零点はありません。通常のゼータ関数のRHだけを仮定する話ではありません。142の [10-analytic.tex:172–178](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/10-analytic.tex:172) も、この意味を明記しています。  
   また、引用定理を \(K=\mathbb Q(\mu_{12})\)、自明指標に適用すると、\(\zeta_K(s)\) が同じ半平面で零点を持ちません。\(\zeta_K\) は \(\zeta(s)\) と非主Dirichlet \(L\) 関数の積なので、\(\zeta(s)\) に対する準リーマン予想が従います。適用条件は [同:9–17](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/10-analytic.tex:9)。  
   ただし「文献では未証明」は断定が強すぎます。029自身が証明を提示しています（[029・06-zero-free.tex:740–764](./data/openai_math/preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/06-zero-free.tex:740)）。正確には「既存の標準的な定理ではなく、029が証明を主張しているが、その証明は監査で検証できていない」です。既存の零点なし領域との違いは、[Thorner–Zamanの定理3.1](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-p.pdf)でも確認できます。

3. **判定：正しい。**  
   条件と上界は [142・00-introduction.tex:47–55](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/00-introduction.tex:47)、[10-analytic.tex:91–98](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/10-analytic.tex:91) に一致します。  
   \(q=2\) では、Eulerの規準と相互法則から
   \[
   -1=\left(\frac p\ell\right)
      =(-1)^{((p-1)/2)((\ell-1)/2)}
        \left(\frac\ell p\right)
      =\left(\frac\ell p\right).
   \]
   \(\ell\equiv1\pmod{24}\) により符号は常に \(+1\)。したがって、**任意の奇素数 \(p\) について、条件を満たす \(\ell\) があれば素数平方非剰余です**。ただし \(p=3\) では、そのような \(\ell\) 自体が存在しません。  
   \(n=2\) なら \(B=23+3\lceil\log_2p\rceil\)。指数関数は固定次数の多項式より速く増えるので、\(p>B^{200000}\) は十分大きい全素数で成立します。帰結は空虚ではなく、最小素数平方非剰余は
   \[
   \ll(\log p)^{20000000}
   \]
   となります。残る有限個の奇素数は定数に吸収できます。Burgess型評価との比較も妥当です。[2026年の査読論文の序論](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/distributions-associated-to-small-quadratic-nonresidues-and-their-partitions/33B7B8E394F7EB99757D42F21B381B61)も、既知の一般的上界を \(p^{1/(4\sqrt e)+\varepsilon}\) としています。  
   なお、補助素数そのものが「最小」である必要はありません。最小値への**上界を与える**、が正確です。

4. **判定：正しい。**  
   [142・main.tex:17–26、60–70](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/main.tex:17) の共有カウンタと入力順を確認しました。導入部の二定理が1.1・1.2、最後の節の定理・補題・命題が11.1・11.2・11.3です。報告書の `main.tex` 行番号にはずれがありますが、定理番号は合っています。

5. **判定：正しい。**  
   029から独立した帰着であることは [142・00-introduction.tex:61–80](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/00-introduction.tex:61) と、[09-complexity.tex:304–308](./data/openai_math/preprints/Deterministic-Polynomial-Factorization-over-Prime-Fields-October-4-2026/build/sections/09-complexity.tex:304) の証明構造から確認できます。HeckeのGRHで補助素数を小さくできれば、多項式時間となり、Evdokimovの準多項式評価を改善します。比較対象の評価は [Altmanの論文・序論](https://arxiv.org/html/2509.12705v1#S1)でも確認できます。  
   **ただし計算量は \(O((B+E)^C)\) で、\(E\) は補助素数の数値です。** 任意の補助素数を二進入力で与えれば、その入力長に関して多項式時間になる、とは述べていません。また、この依存関係の確認は、代数的証明全体の正しさの認定ではありません。

6. **判定：正しい。**  
   [029・01-introduction.tex:38–53](./data/openai_math/preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/01-introduction.tex:38) は、同じ体の族、全有限位数Hecke指標、同じ幅、導手・高さ制限なし、主指標の極は許容、と明記しています。142の定理11.1との取り違えはありません。途中の評価定数が体や指標に依存することと、零点なし幅が共通であることも区別されています。

公開文面への修正案として、[README.md:57](./public_repo/README.md:57) は次の点を直すべきです。

- 「全面依存」「確信度は029が上限」→「**提示された無条件証明は、029の定理1.2の正しさに依存する**」。142の結論の真偽と、採用した証明の依存関係を区別する。
- 「文献では未証明」→「**未査読の029が証明を主張しているが、監査ではその核心を検証できていない**」。「GRH」は「有限位数Hecke \(L\) 関数のGRH」と明記する。
- 非剰余の帰結は「**命題11.3が正しければ、十分大きい全素数 \(p\) に対して**」とする。
- 「補助素数を与えれば多項式時間」→「**\(B\) と補助素数の数値上限 \(E\) に関して多項式時間。\(E\le B^{O(1)}\) なら元の入力長に関して多項式時間**」。
- 「自己完結」は「029の解析定理に依存しない」とする。16件の検算はClaudeの報告として帰属させ、今回の再確認が全証明の検証まで含むようには書かない。

全体として、報告の主要な依存関係と数学的帰結は支持できます。相互法則の符号や大標数条件による反証はありません。公開時には、未検証の証明をどう位置づけるかと、計算量の入力尺度を明確にしてください。