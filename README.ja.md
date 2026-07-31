# aerophysics

`aerophysics` は、標準大気や圧縮性流れなどの基礎物理モデルを、
出典と適用範囲を追跡できる Python API として提供する科学計算パッケージです。

熱量的完全気体とNASA多項式による熱的完全気体、U.S. Standard Atmosphere
1976、等エントロピー流れ、垂直・斜め衝撃波、Prandtl–Meyer 膨張、飛行条件、
航空慣用単位の明示的変換を対象とします。滑面平板の層流・乱流・指定遷移
境界層と圧縮性補正も利用できます。公開計算 API の単位は SI、角度は radian
を基本とします。

## インストール

Python 3.12 以上が必要です。

```console
python -m pip install aerophysics
```

### ローカルGUI

GUI用の追加依存関係を含めてインストールすると、ブラウザ上で標準大気、
飛行条件、斜め衝撃波、平板境界層を計算・プロットできます。

```console
python -m pip install "aerophysics[gui]"
aerophysics-gui
```

起動後にローカルURLがブラウザで開きます。入力値は選択した表示単位から
SIへ明示変換され、結果CSVと再現用の設定JSONをダウンロードできます。

## 使用例

```python
from aerophysics import FlightCondition, standard_atmosphere

sea_level = standard_atmosphere(0.0)
print(sea_level.temperature)       # 288.15 K
print(sea_level.speed_of_sound)    # 340.294... m/s

condition = FlightCondition.from_mach(
    geometric_altitude=10_000.0,
    mach=0.8,
    characteristic_length=2.0,
)
print(condition.dynamic_pressure)  # Pa
print(condition.reynolds_number)
```

計算 API は SI 単位を使用します。航空慣用単位は `aerophysics.units` の
明示的な変換関数を使用してください。

## NASA多項式による熱的完全空気

`AIR_NASA7` と `AIR_NASA9` は、N₂/O₂/Ar/CO₂のモル分率を固定した乾燥空気
です。各成分のNASA多項式をモル分率で混合し、

$$
c_p(T)=\frac{\sum_i x_i\bar c_{p,i}(T)}{M},\qquad
c_v(T)=c_p(T)-R,\qquad
\gamma(T)=\frac{c_p(T)}{c_v(T)}
$$

を計算します。

```python
from aerophysics import AIR_NASA7, AIR_NASA9

temperature = [300.0, 1000.0, 2000.0, 6000.0]

print(AIR_NASA9.cp(temperature))                  # J/(kg K)
print(AIR_NASA9.cv(temperature))                  # J/(kg K)
print(AIR_NASA9.heat_capacity_ratio(temperature))
print(AIR_NASA9.speed_of_sound(temperature))      # m/s

# NASA標準エンタルピーは生成エンタルピーを含みます。
print(AIR_NASA9.standard_enthalpy(1000.0))         # J/kg

# 温度変化に使う顕熱は基準温度を明示して求めます。
print(AIR_NASA9.sensible_enthalpy(
    1000.0,
    reference_temperature=298.15,
))
```

乾燥空気の計算例は次のとおりです。今回採用したNASA7/9データは、4成分
すべてについて200–6000 Kが公称範囲なので、表中に外挿値はありません。

| T [K] | 形式 | cp [J/(kg K)] | cv [J/(kg K)] | γ |
|---:|:---:|---:|---:|---:|
| 300 | NASA7 | 1004.844 | 717.789 | 1.399914 |
| 300 | NASA9 | 1004.829 | 717.774 | 1.399923 |
| 1000 | NASA7 | 1140.675 | 853.621 | 1.336279 |
| 1000 | NASA9 | 1141.033 | 853.978 | 1.336138 |
| 2000 | NASA7 | 1251.921 | 964.866 | 1.297507 |
| 2000 | NASA9 | 1250.334 | 963.280 | 1.297997 |
| 3500 | NASA7 | 1308.475 | 1021.420 | 1.281035 |
| 3500 | NASA9 | 1309.321 | 1022.266 | 1.280802 |
| 5000 | NASA7 | 1342.190 | 1055.136 | 1.272054 |
| 5000 | NASA9 | 1340.968 | 1053.914 | 1.272370 |
| 6000 | NASA7 | 1357.248 | 1070.194 | 1.268227 |
| 6000 | NASA9 | 1360.666 | 1073.611 | 1.267373 |

任意の係数データを公称範囲外で評価する場合、通常は `ModelRangeError` に
なります。`allow_extrapolation=True` を明示すると端の多項式区間で外挿し、
`ApplicabilityWarning` が通知されます。

6000 Kまで数値評価できても、ここでの組成は凍結されています。実際の高温
空気で生じる解離・電離・化学平衡・熱的非平衡は扱いません。また、
`AIR_NASA7`/`AIR_NASA9` の局所的な $\gamma(T)$ を、定数 $\gamma$ を仮定
する既存の衝撃波・等エントロピー・Prandtl–Meyer解析式へ代入することは
できません。

境界層内の単独突起について、自由流中の抗力係数と前面積上の有効動圧積分
から直接抗力を推定できます。既定の乱流 1/7 乗速度分布、任意の速度・密度
プロファイル、および Walz 温度関係による圧縮性近似に対応します。

## 開発

```console
uv sync --all-groups --all-extras
uv run pytest
uv run ruff check .
uv run mypy
uv run sphinx-build -W -b html docs docs/_build/html
```

モデルの数式、出典、仮定、適用範囲は英語の API 文書に記載します。

## リリース

`pyproject.toml` のバージョンを更新してコミットし、一致する `vX.Y.Z`
タグを push します。

```console
git tag v0.3.0
git push origin v0.3.0
```

リリースワークフローはすべての検証を実行し、wheel、sdist、
`aerophysics-docs-X.Y.Z.zip` を非公開リポジトリの GitHub Release に
添付します。ドキュメントZIPを展開し、
`aerophysics-docs-X.Y.Z/index.html` をブラウザで開いてください。
Releaseの取得には、このリポジトリの読み取り権限が必要です。
