# 調和振動子・Beattie–Bridgemanモデル文献

このworktreeでは、モデルの理論原典と、空気定数・風洞計算法・検証値を
区別して記録する。既存のNACA Report 1135は親作業ツリーの文献集に収録済み
なので重複保存しない。`literature/` はリポジトリ設定によりGit追跡対象外。

| 区分 | 文献 | 状態 | ローカルPDF・識別子 | 用途・メモ |
|---|---|---|---|---|
| 理論原典 | E. H. Kennard, *Kinetic Theory of Gases*, McGraw-Hill, 1938 | 未取得 | LCCN 38006323; OCLC 537197 | 調和振動子による振動比熱の理論。市販書籍のため書誌のみ。 |
| 理論原典 | J. A. Beattie and O. C. Bridgeman, *A New Equation of State for Fluids*, Proc. American Academy of Arts and Sciences 63, 229–306, 1928 | 未取得 | DOI `10.2307/20026205` | Beattie–Bridgeman状態方程式。非OAのため書誌・DOIのみ。 |
| 空気定数・式 | R. E. Randall, *Thermodynamic Properties of Air: Tables and Graphs Derived from the Beattie–Bridgeman Equation of State Assuming Variable Specific Heats*, AEDC-TR-57-8, AD 135331, 1957 | 未取得 | AD 135331 | JAXA参考文献記載を確認。DTIC/NTIS公式PDFを再検索したが公開ダウンロードを確認できなかった。 |
| 空気定数・式 | R. E. Randall, *Thermodynamic Properties of Gases: Equations Derived from the Beattie–Bridgeman Equation of State Assuming Variable Specific Heats*, AEDC-TR-57-10, AD 135332, 1957 | 未取得 | AD 135332 | JAXA参考文献記載を確認。DTIC/NTIS公式PDFを再検索したが公開ダウンロードを確認できなかった。 |
| 風洞実装 | M. Watari, *Air Models Used in Flow Calculations for the JAXA Hypersonic Wind Tunnel*, JAXA-RR-06-011, 2007 | 取得済み | [PDF](2007_Watari_Air-Models-JAXA-Hypersonic-Wind-Tunnel.pdf) | 50ページ。定数、統一式、風洞計算手順、FORTRANコードを収録。モデルの発明元ではない。 |
| 風洞実装 | M. Watari, *Air Models Used in the JAXA Hypersonic Wind Tunnels*, JAXA-SP-06-020, 23–26, 2007 | 取得済み | [PDF](2007_Watari_Air-Models-JAXA-Hypersonic-Wind-Tunnels-Conference.pdf) | 4ページの会議論文。完全版との関係確認用。 |
| 独立検証 | D. W. Witte and K. E. Tatum, *Computer Code for Determination of Thermally Perfect Gas Properties*, NASA TP-3447, 1994 | 取得済み | [PDF](1994_Witte_Computer-Code-Thermally-Perfect-Gas-Properties.pdf) | 80ページ。熱的完全気体流れの独立検証資料であり、調和振動子モデルの原典ではない。 |

JAXA資料は、0.5 m極超音速風洞と0.44 m衝撃風洞でThermally Perfectモデル、
1.27 m極超音速風洞でBeattie–Bridgemanモデルを採用した再実装資料である。
物理モデル名にはJAXAを含めず、コードの空気プリセット名も
`AIR_HARMONIC_OSCILLATOR` と `AIR_BEATTIE_BRIDGEMAN` とする。

最終確認日: 2026-08-05
