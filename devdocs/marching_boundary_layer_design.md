# 変動外縁マーチング境界層ソルバー 設計メモ

状態: **提案（未実装）**。基準コミットは `6b16c0d`、記載日は 2026-11-06。

最終目標は、円錐–円柱のような軸対称物体に沿う `Cf(s)` と境界層厚さの分布を、
上流の履歴を含めて求めることです。このメモは第1段階（平面2次元の変動外縁）の
API と検証計画を定め、第2段階（軸対称）への接続を示します。

## 1. 決定事項

| 項目 | 決定 | 理由 |
| --- | --- | --- |
| 置き場所 | `aerophysics` 内の新モジュール | 依存が NumPy/SciPy で足り、物性、例外、検証基盤を共有できる |
| 公開関数 | `flat_plate_boundary_layer` とは別の関数 | 既存 API の配列は「独立した条件」、本機能の配列は「同じ流れの下流分布」で意味が異なる |
| 外部流 | 範囲外。外縁条件を入力で受け取る | 形状依存の処理と外部流ソルバーは [PROJECT_NOTES](../PROJECT_NOTES.md) の方針どおり別扱い |
| 既存 SA スクリプト | 移動せず凍結 | 文書上の「独立した検証器」という位置づけを保つ |
| 実装言語 | Python と NumPy/SciPy | 8節の計測を参照 |

## 2. 第1段階の範囲

**対象:** 鋭い前縁から始まる平面2次元の付着境界層、完全気体（`PerfectGas`）、
一定の `Pr` と `Prt`、層流または SA-noft2、指定壁温または断熱壁、
流れ方向に変化する等エントロピーの外縁条件。

**対象外:** 軸対称（第2段階）、外部流との連成、指定遷移と遷移予測、粗さ、
吸込み・吹出し、壁面曲率、エントロピー層（鈍頭物体）、高温気体、逆流、GUI。

## 3. 支配方程式

### 変数と規約

- 単位は SI。`s` は前縁から下流へ測る壁面距離、`n` は壁面法線距離。
- 変換は Levy–Lees 型（`j=0` が平面、`j=1` が第2段階の軸対称）。

  ```text
  xi(s)  = int_0^s rho_e mu_e U_e r_w^(2j) ds
  eta    = U_e r_w^j / sqrt(2 xi) * int_0^n rho dn
  f' = u/U_e,  g = T/T_e,  n_sa = nu_tilde/nu_e,  Y = int_0^eta g d(eta)
  C  = rho mu / (rho_e mu_e),  chi = n_sa / (C g^2),  eps = chi fv1
  ```

- `beta > 0` は順圧力勾配（加速）を表す。
- `Cf` は既存 API と同じく外縁動圧基準。熱流束は壁へ入る向きを正とする。

### 方程式

`x` を一定外縁の [`build_turbulent_sa_reference.py`](../docs/scripts/build_turbulent_sa_reference.py)
の式から一般化します。変更点は式の後にまとめます。

```text
beta   = 2 xi d(ln U_e)/d(xi)
Lambda = 2 xi d(ln nu_e)/d(xi)
R      = sqrt(2 xi) / (r_w^j mu_e)      （平板では sqrt(2 Re_x)）

(C(1+eps) f'')' + f f'' + beta (g - f'^2)            = 2 xi (f' f'_xi - f'' f_xi)

(C/Pr (1 + eps Pr/Prt) g')' + f g'
    + (gamma-1) M_e(xi)^2 C (1+eps) f''^2            = 2 xi (f' g_xi - g' f_xi)

((C g + n/g) n')'/(sigma g) + cb2 n'^2/(sigma g^2) + f n'
    + cb1 St n - cw1 fw n^2/Y^2 - Lambda f' n        = 2 xi (f' n_xi - n' f_xi)

St = R |f''|/g + n fv2/(kappa^2 Y^2)
y  = sqrt(2 xi) Y / (U_e r_w^j rho_e)
Cf = 2 C_w f''_w / R
```

変更点は次のとおりです。

- **運動量式:** `beta (g - f'^2)` を追加。
- **エネルギー式:** 陽な `beta` 項は現れません。外縁で等エントロピーなら、圧力仕事と
  `h_e` の変化が打ち消し合います。粒子が等エントロピー的に膨張しても `T/T_e` は
  変わらない、という見方と整合します。ただし `M_e`、`T_e`、`g_w` は `xi` に依存します。
- **SA 式:** `- Lambda f' n` を追加し、渦度スケールを `R(xi)` に置き換えます。
  拡散項と破壊項の形は変わりません。

### 外縁の SA 境界条件

外縁では `n_e` を固定せず、**有次元の `nu_tilde` を流線に沿って保存**します
（`n_e ∝ 1/nu_e`）。`n_e` を固定すると `Lambda` 項と釣り合わず、人工的な外縁勾配が
生じます。試作では、外層の `n` の変動が固定で約50%、保存で約0.2%でした。
`Cf` への影響は 0.001% 未満なので、主な理由は方程式との整合です。

### 外縁条件の整合性

上の式は次の2条件を前提とします。

```text
dp_e = -rho_e U_e dU_e      （外縁の Euler 式）
dh_e = -U_e dU_e            （外縁の全エンタルピー一定）
```

したがって独立な関数は1つです。`p_e(s)` と外縁の全圧・全温を基本入力にして残りを
等エントロピー関係から計算し、不整合な入力が入らないようにします。鋭い先端の
円錐–円柱では、物体表面の流線は同じ円錐衝撃波を通過した後に等エントロピーなので、
この仮定を満たします。

## 4. 公開 API 案

`src/aerophysics/marching_boundary_layer.py` に公開 API を置き、数値計算の本体は
`_laminar_similarity.py` に倣って非公開モジュール `_boundary_layer_march.py` に分けます。

```python
@dataclass(frozen=True, slots=True)
class BoundaryLayerEdge:
    distance: FloatArray          # m, s >= 0, 狭義単調増加
    pressure: FloatArray          # Pa
    temperature: FloatArray       # K
    velocity: FloatArray          # m/s
    gas: PerfectGas
    viscosity_model: DynamicViscosityModel

    @classmethod
    def from_pressure(
        cls, distance, pressure, *, total_pressure, total_temperature,
        gas=AIR, viscosity_model=AIR_VISCOSITY,
    ) -> "BoundaryLayerEdge": ...


def marching_boundary_layer(
    edge: BoundaryLayerEdge,
    *,
    regime: BoundaryLayerRegime,          # LAMINAR または TURBULENT
    wall_temperature: ArrayLike | None = None,   # None は断熱壁、スカラーまたは s の配列
    prandtl_number: float = 0.72,
    turbulent_prandtl_number: float = 0.9,
    edge_sa_ratio: float = 3.0,           # 始点の nu_tilde/nu_e
) -> MarchingBoundaryLayerResult: ...
```

- **既存規約の踏襲:** 引数名は `flat_plate_boundary_layer` に揃えます。結果は凍結
  dataclass で、配列は呼び出し元と共有しません（`tests/test_array_ownership.py` の規約）。
- **`TRANSITIONAL`:** 第1段階では `ModelRangeError` とします。
- **許容誤差:** 0.2 の設計決定どおり内部に固定し、引数にしません。
- **格子:** 層流の `eta` は、初回の比 1.015、最大間隔 0.05、`eta_max=16`（519節点）です。
  比 1.04 では基準の 0.05% を満たせませんでした。SA は既存の検証で確認済みの幾何格子
  （比 1.04、`eta_max=150`）を初期値に再評価します。流れ方向は `ln p_e` の PCHIP 補間と、`ln xi`・`ln U_e` の刻み上限で内部的に
  細分し、結果は利用者の `distance` で返します。上限値は検証で決めます。

```python
@dataclass(frozen=True, slots=True)
class MarchingBoundaryLayerResult:
    distance: FloatArray
    reynolds_number: FloatArray
    boundary_layer_thickness: FloatArray
    displacement_thickness: FloatArray
    momentum_thickness: FloatArray
    local_skin_friction_coefficient: FloatArray   # 外縁動圧基準
    wall_shear_stress: FloatArray                 # 一様流基準などへの換算用
    wall_temperature: FloatArray
    wall_heat_flux: FloatArray                    # 壁へ入る向きが正
    pressure_gradient_parameter: FloatArray       # beta
    eta: FloatArray                               # (N_eta,)
    wall_distance: FloatArray                     # (N_s, N_eta)
    velocity_ratio: FloatArray                    # (N_s, N_eta)
    temperature_ratio: FloatArray                 # (N_s, N_eta)
    nu_tilde_ratio: FloatArray | None             # (N_s, N_eta)
```

`wall_heat_flux` は既存平板 API の Reynolds アナロジーとは異なり、温度解の壁面勾配から
求める値です。文書で区別します。

### 失敗と警告

| 状況 | 通知 |
| --- | --- |
| 入力が物理的に不正（負の圧力、単調でない `s` など） | `ValueError` |
| 壁面せん断応力が0へ向かう（剥離） | `BoundaryLayerSeparationError(ModelRangeError)`。剥離位置と上流までの結果を属性に持つ |
| Newton と Picard がともに収束しない | `BoundaryLayerConvergenceError(RuntimeError)` |
| 圧力の変化長が `delta` に対して短い（鋭い肩など） | `ApplicabilityWarning`（`beta` の閾値は検証で決定） |
| 壁第1点の `y+ > 1`、強い粘性干渉 | `ApplicabilityWarning` |

## 5. 数値解法

既存スクリプトの構成を踏襲します。流れ方向は可変刻み BDF2（始点のみ BDF1）で `xi` に
ついて離散化します。未知数は節点あたり `(f, u, g, n, Y)` の5個で、減衰 Newton と
帯行列 LAPACK で解き、失敗時に緩和 Picard へ切り替えます。

新規に必要なものは次のとおりです。

1. 上記の `beta`、`Lambda`、`R(xi)` と、`xi` 依存の `M_e`、`T_e`、`g_w`、`n_e`
2. 断熱壁の Neumann 条件 `g'(0)=0`。既存スクリプトの `AD` は回復温度を指定した等温壁です
3. `beta` と `Lambda` を含む Picard フォールバック
4. 任意の粘性モデルに対する `C` と `nu_e` の評価（既存スクリプトは Sutherland 固定）
5. 外縁速度、壁面せん断、壁面温度勾配の変化を見る流れ方向の刻み制御

**最大のリスクは SA の正値性です。** 強い膨張では `n=nu_tilde/nu_e` が下がり、乱流の
外縁フロントで Newton の正値性ダンピングが詰まります。試作では `beta` と `Lambda` を
含む Picard で通過できましたが、フォールバックは最悪で約450ステーション中95回でした。
まず実績のある Newton + Picard で実装し、回数や所要時間が問題になれば、負の `nu_tilde` を
滑らかに扱える SA-neg（Allmaras, Johnson and Spalart 2012）を別途評価します。

## 6. 検証計画

既存スクリプトとの一致は整合確認であり、独立検証には数えません。層流は解析解、
乱流は積分収支と外部 CFD で確認します。

「試作」欄はリポジトリ外の使い捨てスクリプトでの結果です。既存スクリプトのクラスを
継承して項を追加したものなので、正式実装では全項目をテストとして再度確立します。

| ID | 検証 | 基準 | 試作の結果 |
| --- | --- | --- | --- |
| V1 | ZPG 極限と `turbulent_sa_boundary_layer.csv` | 既存の格子不確かさ 0.1% 以内 | 継承のため差0（検証としては自明） |
| V2 | Falkner–Skan（`beta = -0.15, 0.1, 0.5`）と独立 BVP | `f''(0)` が 0.05% 以内 | 最大 0.017% |
| V3 | Howarth の線形減速流の剥離位置 | Hartree (ARC R&M 2426) の `x_s = 0.9589`、すなわち `x/L = 0.1199` と 0.5% 以内 | `x/L = 0.1195–0.1196` |
| V4 | Illingworth–Stewartson 変換（`M0=3`、`Pr=1`、`C=1`、断熱） | V3 の位置と一致。`g` が厳密解 `1+m(1-u^2)` と一致 | V3 と4桁一致、`max abs(dg) = 4.4e-4` |
| V5 | Cohen–Reshotko の相似解（壁面加熱あり） | 独立 BVP と 0.05% 以内。NACA TR 1293 Table II は `beta >= 0` の行で 0.2% 以内（表の精度の限界のため） | 未実施 |
| V6 | 圧縮性層流の ZPG 極限と `_laminar_similarity` | 0.05% 以内 | 未実施（既存スクリプトの記録では 0.01% 以内） |
| V7 | 運動量積分式（SA、変動外縁） | 中央値 0.05% 以内。上限は実装で設定 | 600分割の中央値 0.004–0.015% |
| V8 | 全エンタルピー積分式（SA、壁面冷却） | 刻み細分で頭打ちの原因を説明できること | 中央値 0.15–0.23%、下記参照 |
| V9 | 流れ方向の格子収束 | 全ての刻み上限を同時に細分し、`Cf` の差 0.1% 以内、熱流束の差 0.5% 以内。圧力変化と壁温変化の両方で確認 | `Cf` のみ、最大 0.21%（M2→3、幅 0.1） |
| V10 | 壁法線方向の格子収束 | 比 1.04→1.02 の `Cf` 差 0.3% 以内 | 未実施 |
| V11 | 圧力勾配つき SA の外部比較 | SU2 の超音速膨張ケース、NASA TMR の 2D bump（`M=0.2`、CFD どうしの比較）、ARC R&M 3633（`M=0.6–2.8` の局所 `Cf` 実験、精度は約 ±0.0002） | 未実施。データ確認後に基準を設定 |

**層流の正式実装（[#81](https://github.com/pandorobo11/aerophysics/pull/81)）の結果:** V2 は
最大 0.03%、V3 は −0.07%、V4 は 0.01% 未満、V5 は BVP と最大 0.01%、Table II と最大
0.07%、V6 は最大 0.02%、V7 と V8 は M2→3 の中央値がそれぞれ 0.02% と 0.07% でした。
V9 は熱流束が最大 0.40%、`Cf` が最大 0.12% です。V10 は2次精度で収束し、Richardson 推定の
誤差は壁面せん断で約 0.01%、熱流束で約 0.08% です。

**V9 の注意点:** 当初の V9 は `relative_step` だけを変えて `Cf` を比べていましたが、これでは
収束を確かめられませんでした。この上限は多くの条件で刻みを縛っていないからです。その結果、
100 K の壁温変化を 1 cm で与えると熱流束が 22% ずれる問題を見逃していました（`Cf` への影響は
0.002%）。そこで壁面温度勾配の変化にも上限を設け、下限値を入れて熱流束のゼロ交差に対応しました。

**V3 の基準値について:** 二次資料には `x/L = 0.125` とするものがありますが、根拠が
確認できないので使いません。

**V11 の注意点:** TMR は [tmbwg.github.io/turbmodels](https://tmbwg.github.io/turbmodels/index.html)
へ移転しました。2D bump の `Cf` は無次元の `x` に対する値で、最大壁面傾斜が 12.8° あるので
`s` と `x` を区別します。`Cf` の正規化は要確認です。円錐–円柱の肩を挟む `Cf(s)` の計測データは
今回の調査では見つかりませんでした。R&M 3633 の物体は、円錐状の先端を持つ細長い軸対称物体です。

**V8 の注意点:** エネルギー収支の中央値は流れ方向に細分してもほとんど減りません
（例: 0.227%→0.152%、0.167%→0.186%）。壁法線方向の格子、または非保存形の離散化が
原因と考えられますが、未確定です。積分収支の評価自体が中央差分であり、鋭い変化の中では
誤差が目立ちます（幅 0.02 の最大偏差は運動量 3.1%→0.8%、エネルギー 8.4%→2.1% と
細分で減少）。正式実装で原因を切り分けます。

**現時点で言えないこと:** 乱流の圧力勾配の結果は、数値的な整合性を確認しただけです。
実験や他コードに対する精度は未確認です。SA と同じモデルを使う CFD と一致しても、
実験に対する精度の保証にはなりません。

### テストの配置

- 通常の pytest: API、例外、配列の非共有、V2〜V10。SA の V1 と V7〜V10 は高速版
- `scripts/check-generated.sh`: 文書用の検証表と、計算量の大きい SA の V1、V11
- フォールバックと失敗経路も分岐カバレッジ 95% の対象とする

## 7. 実装の順序と第2段階

1. **層流の PR:** `BoundaryLayerEdge`、公開関数、圧力勾配、断熱壁、V2〜V6、文書
2. **SA の PR:** `Lambda`、外縁の `nu_tilde` 保存、Picard フォールバック、V1、V7〜V11
3. **第2段階の PR（軸対称）:** 任意引数 `wall_radius` を追加

第2段階では `xi`、`R`、`y` に `r_w` が入るだけで、薄い層の近似のもとでは
変換後の方程式の形は変わりません。解析的な検証として、層流の鋭い円錐の `Cf` が同じ `s`
の平板の `sqrt(3)` 倍になる Mangler の関係を使えます（局所値の比で、物体全体の平均 `Cf`
の比は `2/sqrt(3)` です）。長い円柱では `delta/r_w` が
大きくなり得るので、横曲率の無視には `ApplicabilityWarning` を設けます。NASA CR-4668 も
細長い物体では曲率項が効く可能性を指摘しています。

円錐–円柱の形状と外縁条件を組み立てる処理はパッケージに入れず、docs の例として
示します。円錐部は既存の Taylor–Maccoll 解で求められますが、肩より下流の外縁条件は
`flow_sequence` だけでは決まらず、外部の Euler 解などが必要です。

## 8. 実装言語の判断

| 計測（1ケース、単一プロセス） | 結果 |
| --- | --- |
| 既存 SA、ZPG、約300ステーション、314節点 | 1.2–1.5 秒 |
| 試作、変動外縁、約450〜900ステーション | 1.6–6.6 秒 |
| 層流の正式実装（519節点） | 0.4–2 秒 |
| 実行時間の内訳（M6、cProfile） | 残差とヤコビアン約60%、`solve_banded` 約19%、Picard 約17%（一部重複） |

帯行列の解法を LAPACK のまま使うなら、組み立てを移植しても速度向上の上限は約5倍です。
1万ケースでも単一プロセスで数時間から半日程度で、ケース単位で容易に並列化できます。よって
コンパイル言語は導入しません。解析ヤコビアンと SA-neg のほうが先に試す価値があります。

大規模なスイープや化学反応への拡張で必要になった場合に、プロファイルを取ったうえで
数値計算の本体だけを Rust（PyO3 + maturin）へ移すことを検討します。

## 9. 未決事項と推奨

| 項目 | 推奨 |
| --- | --- |
| `delta*` と `theta` の積分範囲 | 両方とも数値的な外縁まで積分する。既存スクリプトの `delta*` は `delta99` までなので、比較時に揃える |
| 剥離時の戻り方 | 例外に上流までの結果を持たせる |
| 指定遷移 | 第2段階の後に検討する。現在の平板 API の仮想原点は変動外縁へそのまま適用できない |

## 10. 参考文献の状態

| 状態 | 文献 |
| --- | --- |
| 内容を確認済み | Hartree (ARC R&M 2426)、Curle (ARC R&M 3164)、Cohen and Reshotko (NACA TR 1293)、Harris and Blanchard (NASA TM-83207)、Iyer (NASA CR-4668)、Tauber (NASA TP-2914)、Seiff (NACA TN 3284)、Moore (NACA TN 2521) |
| 書誌情報のみ確認（本文は未入手） | Howarth (1938)、Stewartson (1949)、Illingworth (1949)、Mangler (1948)、Smith and Clutter (1963) |

- Howarth の値は Hartree と Curle の引用で確認しました。Mangler の関係は Tauber、Seiff、Moore で
  確認しました。原論文との直接の照合は残っています。
- BL2D は方程式と検証ケースの参考にとどめ、名称は使いません。NASA のカタログではコードは
  U.S. Release Only ですが、報告書 CR-4668 は NTRS で公開されています。
