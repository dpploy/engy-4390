# Uranium Mill Network Dynamics Simulation

UMass Lowell Fall 2024 <br>
Dept. of Chemical Engineering, Nuclear Program <br>
Engy-4390: Nuclear Systems Design and Analysis

View the project on `NBViewer`: [![NBViewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/dpploy/engy-4390/blob/main/projects/2024/u-mill/project-report.ipynb)

Run the project on `Binder`: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dpploy/engy-4390/HEAD?filepath=projects%2Fu-mill%2Fproject-report.ipynb)

Run usecases on `Binder`: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dpploy/engy-4390/HEAD?filepath=projects%2Fu-mill%2Fusecases)

 >**Students:** Alyssa Puglisi, Cole Connerty, Gui Lima  <br>
 >**Mentor/Collaborator:** [Prof. Valmor F. de Almeida](https://github.com/dealmeidavf) <br>
 >[Dept. of Chemical Engineering (Nuclear Energy Program)](https://www.uml.edu/Engineering/Chemical/faculty/de-Almeida-Valmor.aspx) <br>
 >University of Massachusetts Lowell, USA <br>

Conventional uranium mills are chemical plants that extract uranium through the following proces:
1. Uranium ore is delivered to the mill where is it crushed and refined before being leached (extracted by sulfuric acid or alkaline solutions)
1. Extracted uranium is concentrated to a substance known as "yellowcake"
1. Yellowcake is then transported from the mill to a uranium conversion facility where it is processed through the stages of the nuclear fuel cycle


This reposititory is oriented around modeling the uranium milling processes that take place at the White Mesa Mill. The White Mesa Mill is the only fully-licensed and operating conventional mill in the United States. The mill is central to the highest-grade uranium mines and deposits in the U.S. 

Our goal is to provide user-friendly White Mesa Mill case scenarios, accessed and run through Jupyter Notebook, of the uranium leaching, decantation/filtration, solvent extraction, precipitation, and evaporation/calcination processes. Cortix is an open-source Python library for enabling development of simulation of network modules and will be used to create a dynamic network of each separation process. 

|  |
|:---:|
| <img width="1200" src="pics/network-0.gv.png" title="Plant Layout"> |
| <p style="text-align:center;"><b>Mill Network Graph (Based on White Mesa, Utah, USA).</b></p> |

#### Installation Instructions
+ Windows OS activate Symlinks
- Enable developer mode in Windows Settings > System > For developers
- Find .gitconfig file in User section of File Explorer
- Edit .gitconfig file using Notepad application by adding the following section:
[core] <br>
&nbsp;&nbsp;&nbsp;&nbsp;symlinks = true
- Reboot your machine
- Clone engy-4390 into the desktop Github application to run any nb-cases

References:

 + [Cortix](https://cortix.org/)
 + [Eng-4350: Nuclear Reactor Engineering: course notes](https://github.com/dpploy/engy-4350)
+ [NRC: Conventional Uranium Mills](https://www.nrc.gov/materials/uranium-recovery/extraction-methods/conventional-mills.html)
+ [Energy Fuels: White Mesa Mill](https://www.energyfuels.com/white-mesa-mill)
