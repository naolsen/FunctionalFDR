# Functional FDR

This repository contains the R code used for the article:

N.L. Olsen, A. Pini, S. Vantini: *False discovery rate for functional data*. TEST **30**, 784–809 (2021)
 
https://doi.org/10.1007/s11749-020-00751-x

The code is free to use, modify and replicate, but please cite the above article if doing so. 

## 1D simulation
```Fmax.R```: Implementation of the Fmax method.

```scen.R```: Simulation file.

```plot.R```: Plot simulation results. 

## 2D simulation
```simulation.R```: Runs simulation.

```simulation_results.R```: Extracts results from simulation.

```plot_sim_example```, ```ggplot```: Plot simulation results, cf. figures 3 and 4 of the article. 

## Application: analysis of climate data
```analysis.R```: Analysis file

```plots.R```: Plot results

```Data.RData```: Data file. Contains monthly and year means for each 1x1 degree on Earth for the years 1983-2007. 

These data were obtained from the NASA Langley Research Center Atmospheric Science Data Center Surface meteorological and Solar Energy (SSE) web portal supported by the NASA LaRC POWER Project. 
Data are freely available at: NASA Surface Meteorology and Solar Energy, A Renewable Energy Resource web site (release 6.0): http://eosweb.larc.nasa.gov
