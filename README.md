# Towards-Data-Driven-Modelling-of-Turbulent-Two-Phase-Flows-with-Large-Interfaces

We used Direct numerical simulation (DNS) and Reynolds-averaged Navier–Stokes (RANS) simulations to generate data. In addition, we used the DNS data (which is considered accurate) to correct the RANS model.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

## Installation

```bash

# Basilisk (http://basilisk.fr/src/INSTALL)
# Openfoam 9 (https://openfoam.org/version/9/)
# Python and related packages (numpy,matplotlib,sklearn,pandas, imageio)

```


## **Usage**

It is divided into three main sections. 

Basilisk folder
```bash
It is used for DNS. Open the folder ':\Basilisk\DNS\apps\channelflow'. and use the makefile to compile it. The Basilisk software package is needed.
```

Openfoam folder
```bash
It is used for RANS simulation. Open the folder ':\Openfoam\twophaseopemfoam9'. The Openfoam9 software package is needed.
```

Postprocess folder
```bash
It is used for data analysis.
Four sections are included:
1. laminar Vx.ipynb. The results of the Laminar case are shown here. Open ':\laminar Vx.ipynb' for more details.
2. Vis.ipynb. Visualization of the simulation results was presented. Open ':\Result Visualization\Vis.ipynb' for more details.
3. Wave.ipynb. A gif file of the wave was generated here. Open ':\Wave.ipynb' for more details.
4. RF folder. The random forests model was used here to predict the eddy viscosity. Open ':\RF\Random forest.ipynb' for more details.
```
