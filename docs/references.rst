References
==========

Atmosphere
----------

.. _ref-us-standard-atmosphere-1976:

NOAA, NASA, and USAF. *U.S. Standard Atmosphere, 1976*. NOAA-S/T 76-1562,
NASA-TM-X-74335, 1976.

https://ntrs.nasa.gov/citations/19770009539

Transport properties
--------------------

.. _ref-sutherland-1893:

Sutherland, W. *LII. The Viscosity of Gases and Molecular Force*.
Philosophical Magazine Series 5, 36(223), 507--531, 1893.

https://doi.org/10.1080/14786449308620508

.. _ref-keyes-1951:

Keyes, F. G. *A Summary of Viscosity and Heat-Conduction Data for He, A, H2,
O2, CO, CO2, H2O, and Air*. Transactions of the ASME, 73, 589--596, 1951.

https://doi.org/10.1115/1.4016339

.. _ref-bova-bond-kirk-2010:

The dry-air Keyes equation and SI coefficients used here are also reproduced
by Bova, S. W., Bond, R. B., and Kirk, B. S. in *Stabilized Finite Element
Scheme for High Speed Flows with Chemical Non-Equilibrium*, Sandia Report
SAND2010-1168C, 2010.

https://www.osti.gov/servlets/purl/1124439

.. _ref-blottner-johnson-ellis-1971:

Blottner, F. G., Johnson, M., and Ellis, M. *Chemically Reacting Viscous Flow
Program for Multi-Component Gas Mixtures*. Sandia Report SC-RR-70-754, 1971.

https://digital.library.unt.edu/ark:/67531/metadc1030568/

.. _ref-doraiswamy-2010:

Doraiswamy, S. *Computational Study of Nonequilibrium Chemistry in High
Temperature Flows*. Ph.D. dissertation, University of Minnesota, 2010.
Appendix A records the Ar and CO2 Blottner coefficients used by the frozen
dry-air preset.

https://conservancy.umn.edu/handle/11299/99626

.. _ref-wilke-1950:

Wilke, C. R. *A Viscosity Equation for Gas Mixtures*. Journal of Chemical
Physics, 18(4), 517--519, 1950.

https://doi.org/10.1063/1.1747673

Thermochemistry
---------------

.. _ref-mcbride-gordon-reno-1993:

McBride, B. J., Gordon, S., and Reno, M. A. *Coefficients for Calculating
Thermodynamic and Transport Properties of Individual Species*. NASA TM-4513,
1993. This is the source of the NASA seven-coefficient data.

https://ntrs.nasa.gov/citations/19940013151

.. _ref-mcbride-zehe-gordon-2002:

McBride, B. J., Zehe, M. J., and Gordon, S. *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*.
NASA/TP-2002-211556, 2002. This documents the NASA nine-coefficient form and
database.

https://ntrs.nasa.gov/citations/20020085330

.. _ref-nasa-cea-data:

NASA Glenn Research Center. *Chemical Equilibrium with Applications (CEA)*,
``data/thermo.inp``, commit
``2c33ae8ab74635beb85b580c37caccd61778554d``.

https://github.com/nasa/cea

.. _ref-cantera-nasa-gas-data:

Cantera Developers. ``data/nasa_gas.yaml``, commit
``39cfc1436347cfe515cfe499c131ea1220743f16``. The file distributes the
NASA TM-4513 species coefficients used by ``AIR_NASA7``.

https://github.com/Cantera/cantera

Compressible flow
-----------------

.. _ref-naca-report-1135:

Ames Research Staff. *Equations, Tables, and Charts for Compressible Flow*.
NACA Report 1135, 1953.

https://www.nasa.gov/wp-content/uploads/2023/03/equations-tables-charts-compressibleflow-report-1135.pdf

.. _ref-kennard-1938:

Kennard, E. H. *Kinetic Theory of Gases: With an Introduction to Statistical
Mechanics*. McGraw-Hill, 1938. This is the theoretical source used for the
harmonic-vibrational heat-capacity model.

.. _ref-beattie-bridgeman-1928:

Beattie, J. A. and Bridgeman, O. C. *A New Equation of State for Fluids*.
Proceedings of the American Academy of Arts and Sciences, 63, 229--306, 1928.

https://doi.org/10.2307/20026205

.. _ref-randall-1957-air:

Randall, R. E. *Thermodynamic Properties of Air: Tables and Graphs Derived
from the Beattie--Bridgeman Equation of State Assuming Variable Specific
Heats*. AEDC-TR-57-8, AD 135331, 1957.

.. _ref-randall-1957-gases:

Randall, R. E. *Thermodynamic Properties of Gases: Equations Derived from the
Beattie--Bridgeman Equation of State Assuming Variable Specific Heats*.
AEDC-TR-57-10, AD 135332, 1957.

.. _ref-watari-2007-report:

Watari, M. *Air Models Used in Flow Calculations for the JAXA Hypersonic Wind
Tunnel*. JAXA-RR-06-011, 2007. This is an implementation and wind-tunnel
application source, not the origin of either physical model.

https://jaxa.repo.nii.ac.jp/records/2224

.. _ref-watari-2007-proceedings:

Watari, M. *Air Models Used in the JAXA Hypersonic Wind Tunnels*. Proceedings
of the Wind Tunnel Technology Association 76th Meeting, JAXA-SP-06-020,
23--26, 2007.

https://jaxa.repo.nii.ac.jp/records/5839

.. _ref-witte-tatum-1994:

Witte, D. W. and Tatum, K. E. *Computer Code for Determination of Thermally
Perfect Gas Properties*. NASA TP-3447, 1994. This is an independent
thermally-perfect-flow verification source, not the harmonic-oscillator
model's origin.

https://ntrs.nasa.gov/citations/19950005582

.. _ref-sims-1964:

Sims, J. L. *Tables for Supersonic Flow Around Right Circular Cones at Zero
Angle of Attack*. NASA SP-3004, 1964. The Taylor--Maccoll reference values use
an ideal-gas heat-capacity ratio of 1.4.

https://ntrs.nasa.gov/citations/19640009035

.. _ref-ambrosio-wortman-1962:

Ambrosio, A. and Wortman, A. *Stagnation-Point Shock-Detachment Distance for
Flow around Spheres and Cylinders*. ARS Journal, 32(2), 281, 1962.

https://doi.org/10.2514/8.5988

.. _ref-billig-1967:

Billig, F. S. *Shock-Wave Shapes around Spherical- and Cylindrical-Nosed
Bodies*. Journal of Spacecraft and Rockets, 4(6), 822--823, 1967.

https://doi.org/10.2514/3.28969

.. _ref-seiff-1964:

Seiff, A. *Recent Information on Hypersonic Flow Fields*. In *The High
Temperature Aspects of Hypersonic Flow*, NASA SP-24, 1964. This is the source
of the density-ratio sphere standoff correlation.

https://books.google.com/books?id=fygCAAAAIAAJ

.. _ref-inouye-1965:

Inouye, M. *Blunt Body Solutions for Spheres and Ellipsoids in Equilibrium
Gas Mixtures*. NASA TN D-2780, 1965. The report independently compares the
Seiff relation over :math:`0.04<\rho_1/\rho_2<0.16`.  Detached-shock
verification also uses the independently calculated air sphere entry in
Table I, printed page 10 (PDF page 12): :math:`M_\infty=8.949`,
:math:`\rho_\infty/\rho_2=0.1253`, and :math:`\Delta/R_b=0.0994`.

https://ntrs.nasa.gov/citations/19650012766

Boundary layers
---------------

.. _ref-schlichting-gersten-2000:

Schlichting, H. and Gersten, K. *Boundary-Layer Theory*, 8th ed. Springer,
2000.

.. _ref-white-2006:

White, F. M. *Viscous Fluid Flow*, 3rd ed. McGraw-Hill, 2006.

.. _ref-glass-hunt-1988:

Glass, C. E. and Hunt, L. R. *Aerothermal Tests of Quilted Dome Models on a
Flat Plate at a Mach Number of 6.5*. NASA TP-2804, 1988. The report records
the Eckert reference-temperature expression used here.

https://ntrs.nasa.gov/citations/19880012941

.. _ref-lee-faget-1956:

Lee, D. B. and Faget, M. A. *Charts Adapted from Van Driest's Turbulent
Flat-plate Theory for Determining Values of Turbulent Aerodynamic Friction and
Heat-transfer Coefficients*. NACA TN-3811, 1956.

https://ntrs.nasa.gov/citations/19930084604

.. _ref-gnoffo-berry-van-norman-2011:

Gnoffo, P. A., Berry, S. A., and Van Norman, J. W. *Uncertainty Assessments of
2D and Axisymmetric Hypersonic Shock Wave--Turbulent Boundary Layer
Interaction Simulations at Compression Corners*, 2011. Appendix A states the
Van Driest II transformation factors used here.

https://ntrs.nasa.gov/citations/20110013216

.. _ref-hopkins-inouye-1971:

Hopkins, E. J. and Inouye, M. *An Evaluation of Theories for Predicting
Turbulent Skin Friction and Heat Transfer on Flat Plates at Supersonic and
Hypersonic Mach Numbers*. AIAA Journal, 9(6), 993--1003, 1971.

https://doi.org/10.2514/3.6323

.. _ref-hopkins-1972:

Hopkins, E. J. *Charts for Predicting Turbulent Skin Friction from the Van
Driest Method (II)*. NASA TN D-6945, 1972.

https://ntrs.nasa.gov/citations/19730001588

.. _ref-willems-gulhan-2013:

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

.. _ref-young-paterson-1981:

Young, A. D. and Paterson, J. H. *Aircraft Excrescence Drag*. AGARD-AG-264,
1981.

https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA106030.xhtml

.. _ref-hoerner-1965:

Hoerner, S. F. *Fluid-Dynamic Drag*. Published by the author, 1965,
Chapter 5.

https://home.hvl.no/ansatte/gste/ftp/MarinLab_files/Litteratur/Hoerner_1965_Fluid-dynamic_drag.pdf

.. _ref-johnson-mitchell-1971:

Johnson, D. F. and Mitchell, G. A. *Experimental Investigation of Two Methods
for Generating an Artificially Thickened Boundary Layer*. NASA TM X-2238,
1971.

https://ntrs.nasa.gov/citations/19710013205

.. _ref-stallings-lamb-howell-1973:

Stallings, R. L., Jr., Lamb, M., and Howell, D. T. *Drag Characteristics of
Circular Cylinders in a Laminar Boundary Layer at Supersonic Free-stream
Velocities*. NASA TN D-7369, 1973.

https://ntrs.nasa.gov/citations/19740003973

Units
-----

.. _ref-nist-sp-811:

National Institute of Standards and Technology. *Guide for the Use of the
International System of Units (SI)*, NIST Special Publication 811, Appendix B.

https://www.nist.gov/pml/special-publication-811
