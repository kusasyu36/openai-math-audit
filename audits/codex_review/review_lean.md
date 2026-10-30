- **定理1.1・`Separated`：判定＝論文どおり。**  
  [Lean:28–41](./public_repo/lean/Koeru/Hindman164.lean:28) と [論文・序論:14–29](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/01_introduction.tex:14) は一致します。量化順序は `∀ r m R D χ, ∃ a c` で、列はこれらすべてに依存できます。`a 0` には `R < a 0` だけ、それ以降には前項すべての和・積を使う厳密不等式を課しています。指数は実数、自然数の和・積から実数へのキャストも適切です。`StrictMono` は原文の昇順・相異性そのものです。**修正不要。**

- **予想：判定＝論文どおり。**  
  [Lean:20–24](./public_repo/lean/Koeru/Hindman164.lean:20) は、任意の彩色と指定された有限サイズに対する存在主張です。[論文・序論:65–70](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/01_introduction.tex:65) に対応し、無限列や、全サイズに共通する色を要求していません。`m=0` の追加は自明な場合で、意味を変えません。**修正不要。**

- **系1.2・`FS`・`FP`：判定＝論文どおり。**  
  [Lean:43–57](./public_repo/lean/Koeru/Hindman164.lean:43) は [論文・序論:4–11](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/01_introduction.tex:4)、[同:39–46](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/01_introduction.tex:39) に一致します。空集合を除き、単元集合を含み、各元の重複使用はありません。集合像による重複除去、二つの濃度、交わりの等式も正確です。「特に」の合併の濃度は明記されていませんが、型から直ちに従います。**修正不要。**

- **補題2.5：判定＝ずれあり（入力の定義域を一般化）。**  
  [論文・§2:259–262](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/02_framework.tex:259) は **正整数**の列ですが、[Lean:37–41](./public_repo/lean/Koeru/Lemma25.lean:37) の `x : ℕ → ℕ` には正値条件がなく、0 も許します。この場合は `χ 0` も使われ得ます。原文を弱める差ではありませんが、入力条件は一致していません。  
  **修正案：** `x` の後に `(∀ j < n, 0 < x j) →` を追加する。  
  その他は一致します。`n` は `m,r` のみに依存し、鎖・色は `χ,x` に依存できます。`range n` は `[n]` の添字変更です。[`IsChain`:26–34](./public_repo/lean/Koeru/Lemma25.lean:26) は [論文の(2.2):29–34](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/02_framework.tex:29) を正確に表し、値自体の相異性は要求していません。

- **有限和・Ramsey の補助命題：判定＝論文どおり。**  
  [`FSGood`:21–23、`finite_sums_finite`:51](./public_repo/lean/Koeru/FiniteSums.lean:21) は [論文・§2:241–251](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/02_framework.tex:241) に対応します。`F` は彩色に依存せず、正の部分和すべてが `F` 以下です。補助列に相異性がない点も原文どおりです。全自然数上の彩色を使う形は、有限区間上の彩色の形と同値です。  
  [`RBound`・`ramsey`:20–25](./public_repo/lean/Koeru/Ramsey.lean:20)、[`Homog`・`ramsey_simultaneous`:97–103](./public_repo/lean/Koeru/Ramsey.lean:97) も [論文・§2:266–273](./data/openai_math/preprints/Monochromatic-finite-sums-and-products-in-the-positive-integers-September-23-2026/build/sections/02_framework.tex:266) の選択を表しています。空部分集合も含める差は自明で無害です。**修正不要。**

正整数と `r=0` の扱いは、主結果では問題ありません。正の項の非空和・積しか彩色しないため、`χ : ℕ → Fin r` の 0 での値は結論に影響しません。`r=0` ではそのような関数自体が存在せず、予想・系の該当ケースは適切に空虚です。正の `r,m` について主張全体を自明・空虚にする落とし穴は見当たりません。

対象4ファイルの明示的な **`sorry` は [`theorem_1_1`:325–326](./public_repo/lean/Koeru/Hindman164.lean:325) のみ**です。補題2.5と補助定理にはありません。`conjecture_of_theorem` と `corollary_of_theorem` は定理1.1を仮定した条件付き定理です。主定理はコメントにも「未証明」と明記され、証明済みとの過大主張はありません。ただし冒頭7行の「証明は `sorry`」は「主定理の証明は `sorry`」と限定すると正確です。

**総合判定：主結果の命題は忠実。補題2.5には正整数条件の省略という一般化があります。定理1.1の無条件の形式証明は未完成です。**