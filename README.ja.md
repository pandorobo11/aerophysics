# aerophysics

`aerophysics` は、大気・空気力学の基礎物理モデルを、出典と適用範囲を
追跡できるベクトル化 Python API として提供する科学計算パッケージです。
標準大気と気体物性、飛行条件、圧縮性流れ、衝撃波・膨張波、滑面平板
境界層を扱います。各モデルの仮定、適用範囲、出典、検証は英語マニュアルを
正本として記録します。

公開計算 API は SI 単位と radian を使用します。航空機全体の形状設計や
CFD 解析は対象外です。

## インストール

Python 3.12 以上が必要です。

```console
python -m pip install aerophysics
```

ローカル GUI を利用する場合は、追加依存関係を含めてインストールして
起動します。

```console
python -m pip install "aerophysics[gui]"
aerophysics-gui
```

## 最小例

```python
from aerophysics import FlightCondition, standard_atmosphere

sea_level = standard_atmosphere(0.0)
print(sea_level.temperature)       # K
print(sea_level.speed_of_sound)    # m/s

condition = FlightCondition.from_mach(
    geometric_altitude=10_000.0,
    mach=0.8,
    characteristic_length=2.0,
)
print(condition.dynamic_pressure)  # Pa
print(condition.reynolds_number)
```

## ドキュメント

- [英語マニュアルのソースと目次](docs/index.rst)
- [Python クイックスタート](docs/getting_started/quickstart.rst)
- [ローカル GUI ガイド](docs/guides/gui.rst)
- [検証記録](docs/verification/index.rst)
- [API リファレンス](docs/api/index.rst)
- [HTML マニュアルとリリース配布物](https://github.com/pandorobo11/aerophysics/releases/latest)
- [開発・リリース手順](DEVELOPMENT.md)

## ライセンス

[MIT](LICENSE)
