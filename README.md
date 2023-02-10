# Fleet Structure Simulations

## Purpose
To investigate best practices for the treatment of fishery fleet structure and selectivity processes under differential shifts in fishery dynamics. Specifically, we aim to address how considerations of fleet structure and selectivity parameterizations depend on: 
1)	The rate at which fishery fleet structure changes
2)	Differences in selectivity forms among fleets


### Repository Structure
| Folder  | Items |
| --------| --------|
|R_scripts| Contains R scripts for running simulations, utility functions, etc. |
|docs| Auxillary documents (contains some pseudo-code) |
|figs| Contains general figures and model outputs |
|input| Excel files for different life history types (EBS Walleye pollock and Sablefish) |
|src| Contains source code for TMB model, which is compiled through R|
|test| R scripts testing the development of the TMB model. To see examples of running the OM and EM see [here](https://github.com/chengmatt/State_Space_Selex_Sim/blob/master/test/tmb_estimation_test.R)|
