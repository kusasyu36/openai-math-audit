# 監査レポート 059: Zariski 重複度予想の反例（重複度 2 と 3）

- 対象: `data/openai_math/preprints/Ambiently-homeomorphic-isolated-hypersurfaces-of-multiplicities-two-and-three-September-24-2026/`。姉妹論文 `...-germs-in-C4-with-multiplicities-four-and-five-September-27-2026/` も対象。
- 監査日: 2026-10-29（監査役: Claude サブエージェント。報告書は主担当が代筆保存）
- **リポジトリの .py は一度も実行していない。** 計算はすべて自作の `a1〜a4.py` で行い、出力は `a*_output.txt` に保存した。

## 結論

**はっきりした誤りは見つからなかった。** 本論文の主定理が正しい確率は、監査役の見積もりで約 80%。

## 確認済み（OK）

- **補題 3.1 の積関係**（d = 3^11）
  - 849 行すべてで、素因数の指数が完全に打ち消し合う。
  - E > 0（140 桁）。
- **付録 A の階数証明書**（1118/1100 → 896/867、F_1009 上で階数 868）
- **補題 3.2 の鎖構成**: 3396 本すべてで a_i ≥ 2、重みが 1/2 未満、次数が 3 以上。
- **変数の数と不変量**
  - 変数の数がつり合っている。
  - Hessian の階数は 16E と 0。
  - Milnor 数が一致する。
  - N ≡ 0 (mod 8)。
  - 6 個の立方の符号数表 (20,2)/(6,15)/(6,15)。
  - 三倍角恒等式。
- **スペクトルの mod 2 一致**: 小さいモデル（d = 9, 27）で実際の鎖多項式を作り、F_p 上ですべての指標を厳密に照合した。不一致は 0。命題 3.3 の Case 1〜3 を端から端まで確認した。
- **§4〜§6 の論理**（一行ずつ追った）: いずれも仮定を満たしていて、論理は正しいと判断した。
  - 補題 4.1、命題 4.3（奇素数）、命題 4.4（2 進での一意性）、補題 4.5〜4.6（Hasse–Minkowski）
  - 補題 5.2（行列式の式を自分で導出し直した）、補題 5.3（数値）、命題 5.4（SU 群の強近似）
  - §6（Levine の分類、ホモトピー球面、錐への延長）
- **Hertling–Mase**: arXiv:2009.08053 で、Thm 1.3(d) が鎖型の Thom–Sebastiani 和をカバーすることを確認した。

## 懸念（重大度つき）

1. **（中・専門家確認）定理 2.4（Némethi / Balnojan–Hertling）の適用。**
   - 「スペクトルが mod 2 で一致すれば、実 Seifert 形式も同型」という主張が、実（∞）成分への唯一の入力になっている。
   - 原典の正確な主張（固有値 1 の部分、符号の規約）を照合できていない。
   - **論文の正しさを左右する最大の確認点。**
2. **（軽微）** Hertling–Mase の §10 で、x² や x^a を鎖型に含めてよいか。
3. **（軽微）** Levine 1970 Thm 3 の Seifert 形式の規約。Durfee 1974 で代用できる。
4. **（低、結論には影響しない）** 補題 5.1 の複素共役の符号の議論が非常に簡潔。ただし C⊗C を使うので、符号が ± どちらでも消える。
5. **（無害）** 付属の explicit_relation.py は立方の個数を「1 または 2」としているが、論文では常に 2 個。s_0 が偶数なので結果は変わらない。
6. **未精査**: 系 6.4〜6.6（右同値、初期形式の Milnor ファイバー）。主定理には影響しない。
7. **参考**: N は約 10^144 と天文学的だが、有限なので主張に問題はない。

## 姉妹論文（C⁴、重複度 4 と 5）

- 依頼文には「C⁴ では高次元の道具が使えない」と書いたが、これは誤りだった。Kato と Durfee の分類は n ≥ 3 で成り立ち、C⁴ はちょうど n = 3 の境界にあたる。
- 初期形式 T·yz·(yz−T²) と yz·X² は、手計算で確認した（OK）。
- トポロジーの部分は Kato Thm B（n = 3）と Łojasiewicz の錐定理を使っていて、仮定は満たしている。
- **最大のリスク**: 実際の整数 Seifert 形式を計算する「tube」定理（§4〜8）。Mayer–Vietoris、Bockstein、mod 4 での貼り合わせの計算は未検証。

## 再現方法

```
cd ~/chiken/openaikoeru/audits/059
J=../../data/openai_math/preprints/Ambiently-homeomorphic-isolated-hypersurfaces-of-multiplicities-two-and-three-September-24-2026/verification/support/explicit-relation.json
python3 -I a1_relation.py $J        # 約 1 秒
python3 -I a2_spectrum_model.py $J  # 約 2 分
python3 -I a3_cubes.py              # 約 3 秒
python3 -I a4_rank.py               # 約 30 秒
```
