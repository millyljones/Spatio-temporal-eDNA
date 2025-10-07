# Spatio-temporal eDNA from qPCR data

This is the repository associated with the work in [More than presence-absence; modelling (e)DNA concentration across time and space from qPCR survey data](https://link.springer.com/article/10.1007/s42519-025-00477-9). The paper looks at inferring concentrations of target species DNA in the environment from qPCR data. The method links inferred DNA concentrations both spatially and temporally, whilst accounting for lab contamination/inhibition and CT heteroscedasticity (increasing variability in CT values as DNA concentrations decrease). See paper for full model details.

There are three models considered in the paper:
1. Model 1, the full model which accounts for contamination/inhibition and CT heteroscedasticity. This is also called 'Full_Code' throughout the code.
2. Model 2, as Model 1 but does not account for CT heteroscedasticity and assumes that variation in CT values is constant across replicates on a plate (but may vary between plates). This is also called 'Sigmay_Code' throughout the code.
3. Model 3, as Model 1 but does not account for contamination/inhibition. This is also called 'Nocont_Code' throughout the code.

The models are written in NIMBLE (version 1.3.0) which must be installed prior to using this code. All other functions are written in R (4.2).

# Files

The repository contains the following files:

1. Manuscript (pre-publication) and Supplementary Material.
2. Model Codes.R : Contains the NIMBLE models for Models 1, 2, and 3 from the Manuscript.
3. Simulate Data.R : Contains the function to simulate an example qPCR data set.
4. Simulate-Data-Walkthrough : A walkthrough to simulating an example qPCR data set with the functions from `Simulate Data.R'

### Simulating data

See Simulate-Data-Walkthrough.

### Analysing qPCR survey data

See Analyse-qPCR-survey-data. 
