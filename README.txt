[ Main_Data Folder ]
- MainCode.R file replicates the data analysis in Section 5 of the main paper.
- MainFunction.R file contains functions used in MainCode.R file.
- Summary.RData file contains the IPW/nonparametric doubly robust estimates, their standard errors, and the corresponding Wald statistics across 2-dimensional grids of the counterfactual parameters.
- Clean_Data_SelfAttendance.csv is the cleaned data obtained from the raw data which is availble from https://www.openicpsr.org/openicpsr/project/113783/version/V1/view
- Clean_Data_SelfAttendance_Randomization.csv shows which stratum each individual belongs to in the randomization. 
- Estimate folder contains csv files generated from the nonparametric estimation of the outcome regression across 100 sample splitting.

[ Supplementary Folder ]
- Section16.R file replicates the simulations in Section 1.6 of the main paper.
- Function16.R file contains functions used in Section16.R file.
- Section15.R file replicates the simulations in Section 1.5 of the supplementary material.
- Function15.R file contains functions used in Section15.R file.
- Results folder contains csv files generated from the simulation.