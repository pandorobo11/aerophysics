# aerophysics

`aerophysics` は、標準大気や圧縮性流れなどの基礎物理モデルを、
出典と適用範囲を追跡できる Python API として提供する科学計算パッケージです。

0.3 では完全気体、U.S. Standard Atmosphere 1976、等エントロピー流れ、
垂直・斜め衝撃波、Prandtl–Meyer 膨張、飛行条件、航空慣用単位の明示的変換を
対象とします。滑面平板の層流・乱流・指定遷移境界層と圧縮性補正も利用できます。
公開計算 API の単位は SI、角度は radian を基本とします。

## インストール

Python 3.12 以上が必要です。

```console
python -m pip install aerophysics
```

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

## 開発

```console
uv sync --all-groups
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
