# Covid Study - National Lockdowns in England
Data and code to reproduce the analysis from "National lockdowns in England: The same restrictions for all, but how do the reductions in Covid-19 mortality risks differ geographically?".

Description of files:  

"coviddata.csv" - contains the weekly number of deaths, expected deaths, and estimated risks for each LAD in the English mainland. 

"W_nbhd_matrix.rds" - contains the neighbourhood matrix that is used to fit the model.

"IMD_index.csv" - contains the average rank of IMD values in the LSOAs within each LAD, used to rank the LADs according to deprivation.

"R - code for paper.R" - R file which will load required libraries, fit the models, read in the data, and produce most plots (except for maps) and numerical summaries presented in the 'Results' section of the paper.  

# Supplementary materials 
The supplementary materials are provided in the file "Supplementary Materials.docx". They contain mathematical formulas for the weekly expected number of deaths by LAD, trace plots for the MCMC simulation, a map that shows the regions of England, and plots for within-cluster sum of squares and average silhouette width that were used for the clustering analysis. 
