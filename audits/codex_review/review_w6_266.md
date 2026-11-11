実行・ビルドはせず、TeX、配布構成、監査側のソースと保存出力を照合しました。「主」「副」は、それぞれ指定の主論文・副論文を指します。

1. **判定：正しい。**  
   根拠：主 [introduction.tex:13–18](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/introduction.tex:13) は算術・コンパイラ条件を明記しています。[execution.tex:242–249](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/execution.tex:242) が有限実行から非存在への含意、[conclusion.tex:3–10](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/conclusion.tex:3) が記録された実行への適用です。加えて、3個の構成から等号を得ます。  
   **公開文面への修正案：**「明示された算術・コンパイラ条件のもとで、記録された完全実行により N(6)=3 を証明したとする計算機支援証明」。単に未解決の数学的仮定に依存する、という印象は避けてください。

2. **判定：正しい。ただし引用対象の限定が重要。**  
   根拠：主 [execution.tex:82–88](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/execution.tex:82)。限定的な証拠と呼ばれているのは、プローブとプログラム内チェックによる**機械・コンパイラ条件の確認**です。実行全体や数学的証明全体をそう呼んでいるわけではありません。REPORT §5項目1と§8は、この範囲を広げすぎています。  
   **公開文面への修正案：** 現行行の「platform probes」という限定を維持し、「論文自身が実行を証明ではないと認める」とは書かない。

3. **判定：正しい。**  
   根拠：指定4パスはいずれも存在しませんでした。参照箇所は主 execution.tex:87、119、380–390。[README.md:9–12](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/README.md:9) は配布手順の置き換えを明記し、[VERIFY.md:18–31](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/verification/saved-run/VERIFY.md:18) は内部整合性検査と、実行履歴・実装・機械動作の検証を区別しています。ただし、ソース、数値記録、再ビルド手順は配布されています。  
   **公開文面への修正案：**「参照された4ファイルと当時の環境証拠は未同梱。ソース・数値記録・再計算手順は同梱」と補足する。

4. **判定：正しい。改訂履歴についての推測は別。**  
   根拠：主 [paper.tex:51](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/paper.tex:51) は9月24日、[execution.tex:286–301](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/execution.tex:286) は10月4日の実行を記載。[introduction.tex:154–158](./data/openai_math/preprints/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026/build/sections/introduction.tex:154) は9月7日の実行が「ここでは認証されていない」と述べます。ただし、これだけでは初版の公開内容や証拠の有無までは断定できません。  
   **公開文面への修正案：**「9月24日付の配布原稿は10月4日の実行に依拠する。改訂日は確認できない」。REPORT の「初版の時点では証拠が揃っていなかった」は削除する。

5. **判定：正しい〔定理の記載・依存関係・推論〕。**  
   根拠：副 [paper.tex:142–157](./data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/build/paper.tex:142)、[05-verification.tex:4–15](./data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/build/sections/05-verification.tex:4)。主論文の浮動小数点計算を使わない、整数・有理数による証明です。[Weiner 原論文、Corollary 3.4](https://www.renyi.hu/~mweiner/lmsmub.pdf) も、任意の d 個の MUB を d+1 個に補完できると述べています。したがって7個が存在しなければ6個も存在せず、N(6) ≤ 5 は正しい推論です。  
   従来の状況については、少なくとも[2026年3月26日版の総説、Problem 10.1](https://arxiv.org/html/2410.23997v2#S10.SS3) が「7個の非存在」を未解決問題として掲げています。発表直前までの先行研究を網羅確認したわけではありません。  
   **公開文面への修正案：** 推論は維持。「確定した新結果」との断定は避け、「副論文が提示する厳密証明書から得られる上界」とする。本レビューでは証明書計算を再実行していません。

6. **判定：正しい。**  
   根拠：副 [paper.tex:235、273](./data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/build/paper.tex:235) と [04-mixed.tex:215–218](./data/openai_math/preprints/Exact-Fourier-certificates-for-complex-Hadamard-matrices-of-order-six-September-24-2026/build/sections/04-mixed.tex:215)。実際の配布先は `verification/verify.py` です。  
   **公開文面への修正案：**「検証器のパス表記に不一致がある」という軽微な指摘のままでよいです。

7. **判定：正しい〔現行の公開行は概ね公平〕。**  
   根拠：[public_repo/README.md:61](./public_repo/README.md:61) は、主論文の条件と未確認部分を区別し、副論文の再実装についても共有部分を明記しています。**「reproduced independently」という無限定な表現は、現行行にはありません。**  
   REPORT §5項目9の説明は、[自作検証器:4–7](./audits/266/check_266_exact_certificates.py:4)、[同:491–502](./audits/266/check_266_exact_certificates.py:491) と整合します。同じ証明書データを使うこと自体は、別実装による検算を無効にはしません。ただし、同じ数学的帰着や規約にある誤りを共有する可能性は残ります。  
   **公開文面への修正案：**「公開検証器と規約・証明書データを共有する別実装で、証明書計算の出力を再現した」とする。「独立に証明した」「定理全体を独立に確立した」は過大です。また230例の検査は「数値検査」と明記すると明確です。

総評：中心的な読み取りと Weiner による推論は妥当です。修正すべきなのは、プローブへの留保を証明全体へ広げる表現、日付から初版の状態を断定する表現、そして再実装の一致を定理全体の確定と同一視する表現です。現行の公開行は、REPORT 本文の一部より慎重で、公平です。