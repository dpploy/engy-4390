# Tritium Calculation for MSR, No Reaction

UMass Lowell Fall 2024 <br>
Dept. of Chemical Engineering, Nuclear Program <br>
Engy-4390: Nuclear Systems Design and Analysis

View the project on `NBViewer`: [![NBViewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/dpploy/engy-4390/blob/main/projects/2024/tritium/report.ipynb)

Run the project on `Binder`: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dpploy/engy-4390/HEAD?filepath=projects%2F2024%2Ftritium%2Freport.ipynb)

 >[James R. Frye](https://github.com/JamesFrye03), and [Prof. Subash L. Sharma](https://github.com/SubashSharma1008) <br>
 >[Dept. of Chemical Engineering (Nuclear Energy)](xxx) <br>
 >University of Massachusetts Lowell, USA <br>


# Overview
Tritium (T), a radioactive isotope of hydrogen, plays a crucial role in fusion reactions. The U.S. currently produces tritium using lithium control rods at the Watts Bar Nuclear Power Plant, after which it is extracted at the Savannah River Site. This is a batch process that requires tritium to be shipped between facilities for purification, costing between $40-60k per gram. Tritium decays at an annual rate of 5.5%, with existing demand primarily met by recycling tritium from dismantled nuclear weapons. However, the demand for tritium is projected to rise due to its use in fusion reactors and the need to replenish the U.S. government’s stockpile. Tritium is also generated as a byproduct in Molten Salt Reactors (MSRs), which use lithium in the molten salt. In these reactors, tritium is a waste product that must be removed from the primary loop before the steam generator to prevent its release into the environment. One method for tritium removal involves permeation through the primary piping under vacuum. This project aims to model this removal method in a 1D framework and explore different materials to identify the ideal pipe material for efficient tritium separation, enabling its use in fusion applications.

|  |
|:---:|
| <img width="380" src="pics/problem.png" title="Problem domain"> |
| <p style="text-align:center;"><b>Problem domain sketch.</b></p> |

References:

- [1] V. F. de Almeida, [*Engy-5330: Computational Continuum Transport Phenomena*](https://github.com/dpploy/engy-5330),  University of Massachusetts Lowell, Dept. of Chemical Engineering (Nuclear Energy Program).
- [2] V. F. de Almeida, [*Engy-4390: Nuclear System Design and Analysis*](https://github.com/dpploy/engy-4390),  University of Massachusetts Lowell, Dept. of Chemical Engineering (Nuclear Energy Program).
- [3] Multiphysics Object-Oriented Simulation Environment [(MOOSE)](https://mooseframework.org)
- [4] Stempien D. John, “Tritium Transport, Corrosion, and Fuel Performance Modeling in Fluoride Salt-Cooled High-Temperature Reactor (FHR)”. Massachusetts Institute of Technology. PDF. June 2017
- [5] R. Serrano-Lópeza, J. Fraderaa, S. Cuesta-Lópeza. “Molten salts database for energy applications”. PDF. September 2014.
