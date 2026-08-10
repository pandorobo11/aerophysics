# 文献収集メモ

## 指示

- `docs` で参照されている文献を収集する。
- ダウンロードしたPDFは、内容が分かりやすい名前へ変更する。
- 見つかった文献と見つからなかった文献を表として `literature/` 内に記録する。
- この指示自体も `literature/` 内にメモとして残す。

## 収集方針

- 対象は、2026-08-09時点の作業ツリーにある `docs/references.rst` の全45件。
- 記載済みURL、DOI、正式題名の順で調査し、公式機関、出版社、著者、大学リポジトリを優先する。
- 一般公開された第三者サイトも調査対象とするが、認証、購読、アクセス制限の回避は行わない。
- 実際に開けるPDFを保存し、`pdfinfo` でページ数を取得できた場合のみ「取得済み」とする。
- ファイル名は `年_筆頭著者または機関_短縮題名.pdf` とし、ASCII文字、ハイフン区切り、空白なしを基本とする。
- 出版社版を取得できず、正規に公開された著者版・プレプリントを取得した場合はメモ欄に版を記載する。

## Git worktreeへの引き継ぎ

このリポジトリで次の設定を一度実行すると、`git worktree add` で作成した
worktreeへ、メインworktreeの `literature/*.pdf` が自動的にコピーされる。
既に同名のファイルがある場合は上書きしない。

```console
git config core.hooksPath "$(git rev-parse --show-toplevel)/.githooks"
```

## 収集結果

| No. | 分野 | 年 | 著者・機関 | 題名 | 状態 | ローカルファイル名 | 入手元URL | メモ |
|---:|---|---:|---|---|---|---|---|---|
| 1 | Atmosphere | 1976 | NOAA, NASA, USAF | U.S. Standard Atmosphere, 1976 | 取得済み | `1976_NOAA-NASA-USAF_US-Standard-Atmosphere.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf) | 243ページ。NASA-TM-X-74335。 |
| 2 | Compressible flow | 1953 | Ames Research Staff | Equations, Tables, and Charts for Compressible Flow | 取得済み | `1953_Ames_Equations-Tables-Charts-Compressible-Flow.pdf` | [NASA](https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf) | 70ページ。NACA Report 1135。 |
| 3 | Boundary layers | 2000 | Schlichting, H.; Gersten, K. | Boundary-Layer Theory, 8th ed. | 未取得 | — | [Springer](https://link.springer.com/book/10.1007/978-3-642-85829-1) | 出版社では購読コンテンツ。一般公開された8版の全文PDFを確認できなかった。 |
| 4 | Boundary layers | 2006 | White, F. M. | Viscous Fluid Flow, 3rd ed. | 未取得 | — | [McGraw Hill](https://www.mheducation.com/highered/product/Viscous-Fluid-Flow-White.html) | 市販教科書。一般公開された3版の全文PDFを確認できなかった。 |
| 5 | Boundary layers | 1988 | Glass, C. E.; Hunt, L. R. | Aerothermal Tests of Quilted Dome Models on a Flat Plate at a Mach Number of 6.5 | 取得済み | `1988_Glass_Aerothermal-Tests-Quilted-Dome-Models.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19880012941/downloads/19880012941.pdf) | 71ページ。NASA TP-2804。 |
| 6 | Boundary layers | 1956 | Lee, D. B.; Faget, M. A. | Charts Adapted from Van Driest's Turbulent Flat-plate Theory | 取得済み | `1956_Lee_Charts-Adapted-Van-Driest-Turbulent-Flat-Plate.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19930084604/downloads/19930084604.pdf) | 18ページ。NACA TN-3811。 |
| 7 | Boundary layers | 2011 | Gnoffo, P. A.; Berry, S. A.; Van Norman, J. W. | Uncertainty Assessments of 2D and Axisymmetric Hypersonic Shock Wave–Turbulent Boundary Layer Interaction Simulations at Compression Corners | 取得済み | `2011_Gnoffo_Uncertainty-Hypersonic-SWBLI-Compression-Corners.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/20110013216/downloads/20110013216.pdf) | 44ページ。AIAA Paper 2011-3142。 |
| 8 | Boundary layers | 1971 | Hopkins, E. J.; Inouye, M. | An Evaluation of Theories for Predicting Turbulent Skin Friction and Heat Transfer on Flat Plates at Supersonic and Hypersonic Mach Numbers | 未取得 | — | [NASA NTRS](https://ntrs.nasa.gov/citations/19710049171)<br>[DOI](https://doi.org/10.2514/3.6323) | NTRSレコードにダウンロードなし。出版社版は非OAで、公開リポジトリ版も確認できなかった。 |
| 9 | Boundary layers | 1972 | Hopkins, E. J. | Charts for Predicting Turbulent Skin Friction from the Van Driest Method (II) | 取得済み | `1972_Hopkins_Charts-Turbulent-Skin-Friction-Van-Driest-II.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19730001588/downloads/19730001588.pdf) | 49ページ。NASA TN D-6945。 |
| 10 | Boundary layers | 2013 | Willems, S.; Gülhan, A. | Experiments on Shock Induced Laminar-Turbulent Transition on a Flat Plate at Mach 6 | 取得済み | `2013_Willems_Shock-Induced-Laminar-Turbulent-Transition-Mach-6.pdf` | [DLR](https://elib.dlr.de/89126/1/2013_Willems_EUCASS_5.pdf) | 12ページ。DLR公開版。 |
| 11 | Mean-velocity transformations | 1956 | Coles, D. | The Law of the Wake in the Turbulent Boundary Layer | 取得済み | `1956_Coles_Law-of-the-Wake-Turbulent-Boundary-Layer.pdf` | [CaltechAUTHORS](https://authors.library.caltech.edu/records/7qymy-0pf51) | 36ページ。Caltech機関リポジトリの公開版。 |
| 12 | Mean-velocity transformations | 1951 | Van Driest, E. R. | Turbulent Boundary Layer in Compressible Fluids | 未取得 | — | [DOI](https://doi.org/10.2514/8.1895) | 出版社版は非OAで、公開リポジトリ版を確認できなかった。 |
| 13 | Mean-velocity transformations | 1961 | Spalding, D. B. | A Single Formula for the Law of the Wall | 未取得 | — | [DOI](https://doi.org/10.1115/1.3641728) | 出版社版は非OAで、公開リポジトリ版を確認できなかった。 |
| 14 | Mean-velocity transformations | 2014 | Zhang, Y.-S.; Bi, W.-T.; Hussain, F.; She, Z.-S. | A Generalized Reynolds Analogy for Compressible Wall-Bounded Turbulent Flows | 未取得 | — | [DOI](https://doi.org/10.1017/jfm.2013.620) | 出版社版は非OA。検索結果に本文表示はあったが、有効な公開PDFを保存できなかった。 |
| 15 | Mean-velocity transformations | 2016 | Trettel, A.; Larsson, J. | Mean Velocity Scaling for Compressible Wall Turbulence with Heat Transfer | 未取得 | — | [DOI](https://doi.org/10.1063/1.4942022) | 出版社版は非OAで、公開リポジトリ版を確認できなかった。 |
| 16 | Mean-velocity transformations | 2020 | Volpiani, P. S.; Iyer, P. S.; Pirozzoli, S.; Larsson, J. | Data-Driven Compressibility Transformation for Turbulent Wall Layers | 未取得 | — | [DOI](https://doi.org/10.1103/PhysRevFluids.5.052602) | 出版社版は非OA。機関リポジトリではファイルが管理者限定だった。 |
| 17 | Mean-velocity transformations | 2021 | Griffin, K. P.; Fu, L.; Moin, P. | Velocity Transformation for Compressible Wall-Bounded Turbulent Flows with and without Heat Transfer | 取得済み | `2021_Griffin_Velocity-Transformation-Compressible-Wall-Flows.pdf` | [arXiv](https://arxiv.org/pdf/2108.07397) | 8ページ。著者プレプリント。出版版DOIは `10.1073/pnas.2111144118`。 |
| 18 | Mean-velocity transformations | 2022 | Bai, T.; Griffin, K. P.; Fu, L. | Compressible Velocity Transformations for Various Noncanonical Wall-Bounded Turbulent Flows | 取得済み | `2022_Bai_Compressible-Velocity-Transformations-Noncanonical-Flows.pdf` | [arXiv](https://arxiv.org/pdf/2204.00874) | 26ページ。著者プレプリント。出版版DOIは `10.2514/1.J061554`。 |
| 19 | Mean-velocity transformations | 2023 | Griffin, K. P.; Fu, L.; Moin, P. | Near-Wall Model for Compressible Turbulent Boundary Layers Based on an Inverse Velocity Transformation | 取得済み | `2023_Griffin_Near-Wall-Model-Inverse-Velocity-Transformation.pdf` | [Cambridge University Press](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/98624970B28D611D978347620F87D30D/S0022112023006274a.pdf/nearwall_model_for_compressible_turbulent_boundary_layers_based_on_an_inverse_velocity_transformation.pdf) | 22ページ。出版社公開PDF。 |
| 20 | Mean-velocity transformations | 2024 | Danis, M. E.; Durbin, P. | On the Accuracy of Compressibility Transformations | 取得済み | `2024_Danis_Accuracy-Compressibility-Transformations.pdf` | [arXiv](https://arxiv.org/pdf/2410.01915) | 17ページ。著者プレプリント。出版版はOA。 |
| 21 | Boundary-layer protrusions | 1981 | Young, A. D.; Paterson, J. H. | Aircraft Excrescence Drag | 未取得 | — | [NTIS](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA106030.xhtml) | NTISに176ページの公開レコードとPDF導線はあるが、DTICのPDF取得先が403を返し保存できなかった。 |
| 22 | Boundary-layer protrusions | 1965 | Hoerner, S. F. | Fluid-Dynamic Drag | 取得済み | `1965_Hoerner_Fluid-Dynamic-Drag.pdf` | [Western Norway University of Applied Sciences](https://home.hvl.no/ansatte/gste/ftp/MarinLab_files/Litteratur/Hoerner_1965_Fluid-dynamic_drag.pdf) | 455ページ。大学サイトで公開されているスキャン。 |
| 23 | Boundary-layer protrusions | 1971 | Johnson, D. F.; Mitchell, G. A. | Experimental Investigation of Two Methods for Generating an Artificially Thickened Boundary Layer | 取得済み | `1971_Johnson_Artificially-Thickened-Boundary-Layer.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19710013205/downloads/19710013205.pdf) | 40ページ。NASA TM X-2238。 |
| 24 | Boundary-layer protrusions | 1973 | Stallings, R. L., Jr.; Lamb, M.; Howell, D. T. | Drag Characteristics of Circular Cylinders in a Laminar Boundary Layer at Supersonic Free-stream Velocities | 取得済み | `1973_Stallings_Drag-Circular-Cylinders-Laminar-Boundary-Layer.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19740003973/downloads/19740003973.pdf) | 45ページ。NASA TN D-7369。 |
| 25 | Units | 2008 | National Institute of Standards and Technology | Guide for the Use of the International System of Units (SI) | 取得済み | `2008_NIST_Guide-Use-International-System-Units.pdf` | [NIST](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf) | 90ページ。NIST SP 811、2008年版。`docs` の書誌には年の記載なし。 |
| 26 | Thermochemistry | 1993 | McBride, B. J.; Gordon, S.; Reno, M. A. | Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species | 取得済み | `1993_McBride_Coefficients-Thermodynamic-Transport-Properties.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19940013151/downloads/19940013151.pdf) | 94ページ。NASA TM-4513。`AIR_NASA7` の7係数形式とデータの原典。 |
| 27 | Thermochemistry | 2002 | McBride, B. J.; Zehe, M. J.; Gordon, S. | NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species | 取得済み | `2002_McBride_NASA-Glenn-Coefficients-Thermodynamic-Properties.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf) | 295ページ。NASA/TP-2002-211556。`AIR_NASA9` の9係数形式とデータベースの原典。 |
| 28 | Transport properties | 1893 | Sutherland, W. | LII. The Viscosity of Gases and Molecular Force | 取得済み | `1893_Sutherland_Viscosity-Gases-Molecular-Force.pdf` | [Zenodo](https://zenodo.org/records/1602745) | 26ページ。Zenodo公開版。DOI `10.1080/14786449308620508`。 |
| 29 | Transport properties | 1951 | Keyes, F. G. | A Summary of Viscosity and Heat-Conduction Data for He, A, H2, O2, CO, CO2, H2O, and Air | 未取得 | — | [ASME / DOI](https://doi.org/10.1115/1.4016339) | 出版社レコードを確認。一般公開された全文PDFを確認できなかった。 |
| 30 | Transport properties | 2010 | Bova, S. W.; Bond, R. B.; Kirk, B. S. | Stabilized Finite Element Scheme for High Speed Flows with Chemical Non-Equilibrium | 取得済み | `2010_Bova_Stabilized-Finite-Element-High-Speed-Flows.pdf` | [OSTI](https://www.osti.gov/servlets/purl/1124439) | 16ページ。Sandia Report SAND2010-1168C。乾燥空気用Keyes式の再録元。 |
| 31 | Transport properties | 1971 | Blottner, F. G.; Johnson, M.; Ellis, M. | Chemically Reacting Viscous Flow Program for Multi-Component Gas Mixtures | 取得済み | `1971_Blottner_Chemically-Reacting-Viscous-Flow.pdf` | [UNT Digital Library](https://digital.library.unt.edu/ark:/67531/metadc1030568/) | 320ページ。Sandia Report SC-RR-70-754、DOI `10.2172/4658539`。Blottner係数と混合則の出典。 |
| 32 | Transport properties | 2010 | Doraiswamy, S. | Computational Study of Nonequilibrium Chemistry in High Temperature Flows | 取得済み | `2010_Doraiswamy_Nonequilibrium-Chemistry-High-Temperature-Flows.pdf` | [University of Minnesota](https://hdl.handle.net/11299/99626) | 111ページ。Ph.D. dissertation。Appendix AのAr/CO2係数表を参照。 |
| 33 | Transport properties | 1950 | Wilke, C. R. | A Viscosity Equation for Gas Mixtures | 未取得 | — | [AIP Publishing](https://pubs.aip.org/aip/jcp/article-pdf/18/4/517/18796913/517_1_online.pdf) | 出版社では購入対象。一般公開された全文PDFを確認できなかった。DOI `10.1063/1.1747673`。 |
| 34 | High-temperature gas models | 1938 | Kennard, E. H. | Kinetic Theory of Gases: With an Introduction to Statistical Mechanics | 未取得 | — | LCCN `38006323`、OCLC `537197` | 調和振動子による振動比熱の理論原典。市販書籍のため書誌のみ。 |
| 35 | Real-gas equation of state | 1928 | Beattie, J. A.; Bridgeman, O. C. | A New Equation of State for Fluids | 未取得 | — | [DOI](https://doi.org/10.2307/20026205) | Beattie–Bridgeman状態方程式の原典。非OAのため書誌・DOIのみ。 |
| 36 | Real-gas air properties | 1957 | Randall, R. E. | Thermodynamic Properties of Air: Tables and Graphs Derived from the Beattie–Bridgeman Equation of State Assuming Variable Specific Heats | 未取得 | — | AD 135331 | AEDC-TR-57-8。JAXA参考文献記載を確認。DTIC/NTIS公式PDFを再検索したが公開ダウンロードを確認できなかった。 |
| 37 | Real-gas air properties | 1957 | Randall, R. E. | Thermodynamic Properties of Gases: Equations Derived from the Beattie–Bridgeman Equation of State Assuming Variable Specific Heats | 未取得 | — | AD 135332 | AEDC-TR-57-10。JAXA参考文献記載を確認。DTIC/NTIS公式PDFを再検索したが公開ダウンロードを確認できなかった。 |
| 38 | Wind-tunnel gas models | 2007 | Watari, M. | Air Models Used in Flow Calculations for the JAXA Hypersonic Wind Tunnel | 取得済み | `2007_Watari_Air-Models-JAXA-Hypersonic-Wind-Tunnel.pdf` | [JAXA](https://jaxa.repo.nii.ac.jp/records/2224) | 50ページ。JAXA-RR-06-011。空気定数、統一式、風洞計算手順、FORTRANコードを収録。モデルの発明元ではない。 |
| 39 | Wind-tunnel gas models | 2007 | Watari, M. | Air Models Used in the JAXA Hypersonic Wind Tunnels | 取得済み | `2007_Watari_Air-Models-JAXA-Hypersonic-Wind-Tunnels-Conference.pdf` | [JAXA](https://jaxa.repo.nii.ac.jp/records/5839) | 4ページ。JAXA-SP-06-020会議論文。完全版との関係確認用。 |
| 40 | Thermally perfect flow | 1994 | Witte, D. W.; Tatum, K. E. | Computer Code for Determination of Thermally Perfect Gas Properties | 取得済み | `1994_Witte_Computer-Code-Thermally-Perfect-Gas-Properties.pdf` | [NASA NTRS](https://ntrs.nasa.gov/citations/19950005582) | 80ページ。NASA TP-3447。熱的完全気体流れの独立検証資料であり、調和振動子モデルの原典ではない。 |
| 41 | Conical flow | 1964 | Sims, J. L. | Tables for Supersonic Flow Around Right Circular Cones at Zero Angle of Attack | 取得済み | `1964_Sims_Tables-Supersonic-Flow-Right-Circular-Cones.pdf` | [NASA NTRS](https://ntrs.nasa.gov/api/citations/19640009035/downloads/19640009035.pdf) | 430ページ。NASA SP-3004。γ=1.4のTaylor–Maccoll表を円錐衝撃波検証へ使用。 |
| 42 | Detached shocks | 1962 | Ambrosio, A.; Wortman, A. | Stagnation-Point Shock-Detachment Distance for Flow around Spheres and Cylinders | 未取得 | — | [DOI](https://doi.org/10.2514/8.5988) | ARS Journal 32(2), 281。球・二次元円柱の離脱距離相関の原典。出版社レコードと式を確認したが、公開PDFは保存していない。 |
| 43 | Detached shocks | 1967 | Billig, F. S. | Shock-Wave Shapes around Spherical- and Cylindrical-Nosed Bodies | 未取得 | — | [DOI](https://doi.org/10.2514/3.28969) | 頂点曲率半径と双曲線衝撃波形状の原典。出版社レコードと式を確認したが、公開PDFは保存していない。 |
| 44 | Detached shocks | 1964 | Seiff, A. | Recent Information on Hypersonic Flow Fields | 未取得 | — | [NASA SP-24 / Google Books](https://books.google.com/books?id=fygCAAAAIAAJ) | 密度比による球の離脱距離相関の原典。書籍表示で書誌と式を確認したが、全文PDFは保存していない。 |
| 45 | Detached shocks | 1965 | Inouye, M. | Blunt Body Solutions for Spheres and Ellipsoids in Equilibrium Gas Mixtures | 未取得 | — | [NASA NTRS](https://ntrs.nasa.gov/citations/19650012766) | NASA TN D-2780。Seiff式を0.04<ρ1/ρ2<0.16で独立検証した資料。公式レコードと本文を確認したが、PDFは保存していない。 |

## 集計

- 対象: 45件
- 取得済み: 26件
- 未取得: 19件
- 最終確認日: 2026-08-09
