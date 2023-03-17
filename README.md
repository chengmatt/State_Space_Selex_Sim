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


### Operating Model Options
The OM is able to be specified for multiple fishery and survey fleets, as well as multiple sexes. 
| OM Component  | Options |
| --------| --------|
|Recruitment| Beverton Holt Recruitment, Mean Recruitment |
|Fishing mortality pattern| Contrast, Constant, Increase, Decrease, Increase_Plat, Const_Inc, Const_Ramp_Const, Contrast_Const  |
|Selectivity| Logistic, Gamma, Uniform, Exponential Logistic, Double Logistic, Double Normal |
|Compositional Likelihoods| Multinomial, Dirichlet, Dirichlet-Multinomial |

### Estimation Model Options
The EM is able to be specified for multiple fishery and survey fleets, as well as multiple sexes. 
| EM Component  | Options |
| --------| --------|
|Recruitment| Beverton Holt Recruitment, Mean Recruitment |
|Selectivity| Logistic, Gamma, Exponential Logistic, Double Logistic |
|Time-varying selectivity| Random walk, 1DAR1_Year, Time-blocks during any period|
|Compositional Likelihoods| Multinomial, Dirichlet-Multinomial |

