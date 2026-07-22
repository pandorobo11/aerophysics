# aerophysics

`aerophysics` は、標準大気や圧縮性流れなどの基礎物理モデルを、
出典と適用範囲を追跡できる Python API として提供する科学計算パッケージです。

0.1 では完全気体、U.S. Standard Atmosphere 1976、等エントロピー流れ、
飛行条件、航空慣用単位の明示的変換を対象とします。公開計算 API の単位は
SI を基本とします。

## 開発

```console
uv sync --all-groups
uv run pytest
uv run ruff check .
uv run mypy
```

