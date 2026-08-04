References
==========

Atmosphere
----------

NOAA, NASA, and USAF. *U.S. Standard Atmosphere, 1976*. NOAA-S/T 76-1562,
NASA-TM-X-74335, 1976.

https://ntrs.nasa.gov/citations/19770009539

Transport properties
--------------------

Sutherland, W. *LII. The Viscosity of Gases and Molecular Force*.
Philosophical Magazine Series 5, 36(223), 507--531, 1893.

https://doi.org/10.1080/14786449308620508

Keyes, F. G. *A Summary of Viscosity and Heat-Conduction Data for He, A, H2,
O2, CO, CO2, H2O, and Air*. Transactions of the ASME, 73, 589--596, 1951.

https://doi.org/10.1115/1.4016339

The dry-air Keyes equation and SI coefficients used here are also reproduced
by Bova, S. W., Bond, R. B., and Kirk, B. S. in *Stabilized Finite Element
Scheme for High Speed Flows with Chemical Non-Equilibrium*, Sandia Report
SAND2010-1168C, 2010.

https://www.osti.gov/servlets/purl/1124439

Blottner, F. G., Johnson, M., and Ellis, M. *Chemically Reacting Viscous Flow
Program for Multi-Component Gas Mixtures*. Sandia Report SC-RR-70-754, 1971.

https://digital.library.unt.edu/ark:/67531/metadc1030568/

Doraiswamy, S. *Computational Study of Nonequilibrium Chemistry in High
Temperature Flows*. Ph.D. dissertation, University of Minnesota, 2010.
Appendix A records the Ar and CO2 Blottner coefficients used by the frozen
dry-air preset.

https://conservancy.umn.edu/handle/11299/99626

Wilke, C. R. *A Viscosity Equation for Gas Mixtures*. Journal of Chemical
Physics, 18(4), 517--519, 1950.

https://doi.org/10.1063/1.1747673

Thermochemistry
---------------

McBride, B. J., Gordon, S., and Reno, M. A. *Coefficients for Calculating
Thermodynamic and Transport Properties of Individual Species*. NASA TM-4513,
1993. This is the source of the NASA seven-coefficient data.

https://ntrs.nasa.gov/citations/19940013151

McBride, B. J., Zehe, M. J., and Gordon, S. *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*.
NASA/TP-2002-211556, 2002. This documents the NASA nine-coefficient form and
database.

https://ntrs.nasa.gov/citations/20020085330

NASA Glenn Research Center. *Chemical Equilibrium with Applications (CEA)*,
``data/thermo.inp``, commit
``2c33ae8ab74635beb85b580c37caccd61778554d``.

https://github.com/nasa/cea

Cantera Developers. ``data/nasa_gas.yaml``, commit
``39cfc1436347cfe515cfe499c131ea1220743f16``. The file distributes the
NASA TM-4513 species coefficients used by ``AIR_NASA7``.

https://github.com/Cantera/cantera

Compressible flow
-----------------

Ames Research Staff. *Equations, Tables, and Charts for Compressible Flow*.
NACA Report 1135, 1953.

https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf

Boundary layers
---------------

Schlichting, H. and Gersten, K. *Boundary-Layer Theory*, 8th ed. Springer,
2000.

White, F. M. *Viscous Fluid Flow*, 3rd ed. McGraw-Hill, 2006.

Glass, C. E. and Hunt, L. R. *Aerothermal Tests of Quilted Dome Models on a
Flat Plate at a Mach Number of 6.5*. NASA TP-2804, 1988. The report records
the Eckert reference-temperature expression used here.

https://ntrs.nasa.gov/citations/19880012941

Lee, D. B. and Faget, M. A. *Charts Adapted from Van Driest's Turbulent
Flat-plate Theory for Determining Values of Turbulent Aerodynamic Friction and
Heat-transfer Coefficients*. NACA TN-3811, 1956.

https://ntrs.nasa.gov/citations/19930084604

Gnoffo, P. A., Berry, S. A., and Van Norman, J. W. *Uncertainty Assessments of
2D and Axisymmetric Hypersonic Shock Wave--Turbulent Boundary Layer
Interaction Simulations at Compression Corners*, 2011. Appendix A states the
Van Driest II transformation factors used here.

https://ntrs.nasa.gov/citations/20110013216

Hopkins, E. J. and Inouye, M. *An Evaluation of Theories for Predicting
Turbulent Skin Friction and Heat Transfer on Flat Plates at Supersonic and
Hypersonic Mach Numbers*. AIAA Journal, 9(6), 993--1003, 1971.

https://doi.org/10.2514/3.6323

Hopkins, E. J. *Charts for Predicting Turbulent Skin Friction from the Van
Driest Method (II)*. NASA TN D-6945, 1972.

https://ntrs.nasa.gov/citations/19730001588

Willems, S. and Gülhan, A. *Experiments on Shock Induced Laminar-Turbulent
Transition on a Flat Plate at Mach 6*. 5th European Conference for Aeronautics
and Space Sciences (EUCASS), 2013. Equation (7) gives the local implicit
skin-friction relation used here.

https://elib.dlr.de/89126/

Compressible turbulent mean-velocity transformations
-----------------------------------------------------

.. _ref-coles-1956:

Coles, D. *The Law of the Wake in the Turbulent Boundary Layer*. Journal of
Fluid Mechanics, 1(2), 191--226, 1956.

https://doi.org/10.1017/S0022112056000135

.. _ref-van-driest-1951:

Van Driest, E. R. *Turbulent Boundary Layer in Compressible Fluids*.
Journal of the Aeronautical Sciences, 18(3), 145--160 and 216, 1951.

https://doi.org/10.2514/8.1895

.. _ref-spalding-1961:

Spalding, D. B. *A Single Formula for the Law of the Wall*. Journal of Applied
Mechanics, 28(3), 455--458, 1961.

https://doi.org/10.1115/1.3641728

.. _ref-zhang-bi-hussain-she-2014:

Zhang, Y.-S., Bi, W.-T., Hussain, F., and She, Z.-S. *A Generalized Reynolds
Analogy for Compressible Wall-Bounded Turbulent Flows*. Journal of Fluid
Mechanics, 739, 392--420, 2014.

https://doi.org/10.1017/jfm.2013.620

.. _ref-trettel-larsson-2016:

Trettel, A. and Larsson, J. *Mean Velocity Scaling for Compressible Wall
Turbulence with Heat Transfer*. Physics of Fluids, 28, 026102, 2016.

https://doi.org/10.1063/1.4942022

.. _ref-volpiani-2020:

Volpiani, P. S., Iyer, P. S., Pirozzoli, S., and Larsson, J. *Data-Driven
Compressibility Transformation for Turbulent Wall Layers*. Physical Review
Fluids, 5, 052602(R), 2020.

https://doi.org/10.1103/PhysRevFluids.5.052602

.. _ref-griffin-fu-moin-2021:

Griffin, K. P., Fu, L., and Moin, P. *Velocity Transformation for Compressible
Wall-Bounded Turbulent Flows with and without Heat Transfer*. Proceedings of
the National Academy of Sciences, 118(34), e2111144118, 2021.

https://doi.org/10.1073/pnas.2111144118

.. _ref-bai-griffin-fu-2022:

Bai, T., Griffin, K. P., and Fu, L. *Compressible Velocity Transformations for
Various Noncanonical Wall-Bounded Turbulent Flows*. AIAA Journal, 60(7),
4325--4337, 2022.

https://doi.org/10.2514/1.J061554

.. _ref-griffin-fu-moin-2023:

Griffin, K. P., Fu, L., and Moin, P. *Near-Wall Model for Compressible
Turbulent Boundary Layers Based on an Inverse Velocity Transformation*.
Journal of Fluid Mechanics, 970, A36, 2023.

https://doi.org/10.1017/jfm.2023.627

.. _ref-danis-durbin-2024:

Danis, M. E. and Durbin, P. *On the Accuracy of Compressibility
Transformations*. Physics of Fluids, 36, 126119, 2024.

https://doi.org/10.1063/5.0242189

Boundary-layer protrusions
--------------------------

Young, A. D. and Paterson, J. H. *Aircraft Excrescence Drag*. AGARD-AG-264,
1981.

https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA106030.xhtml

Hoerner, S. F. *Fluid-Dynamic Drag*. Published by the author, 1965,
Chapter 5.

https://home.hvl.no/ansatte/gste/ftp/MarinLab_files/Litteratur/Hoerner_1965_Fluid-dynamic_drag.pdf

Johnson, D. F. and Mitchell, G. A. *Experimental Investigation of Two Methods
for Generating an Artificially Thickened Boundary Layer*. NASA TM X-2238,
1971.

https://ntrs.nasa.gov/citations/19710013205

Stallings, R. L., Jr., Lamb, M., and Howell, D. T. *Drag Characteristics of
Circular Cylinders in a Laminar Boundary Layer at Supersonic Free-stream
Velocities*. NASA TN D-7369, 1973.

https://ntrs.nasa.gov/citations/19740003973

Units
-----

National Institute of Standards and Technology. *Guide for the Use of the
International System of Units (SI)*, NIST Special Publication 811, Appendix B.

https://www.nist.gov/pml/special-publication-811
