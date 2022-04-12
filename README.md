# Covid Study - National Lockdowns in England
Data and code to reproduce the analysis from "National lockdowns in England: The same restrictions for all, but how do the reductions in Covid-19 mortality risks differ geographically?".

Description of files:  

"coviddata.csv" - contains the weekly number of deaths, expected deaths, and estimated risks for each LAD in the English mainland. 

"IMD_index.csv" - contains the average rank of IMD values in the LSOAs within each LAD, used to rank the LADs according to deprivation.

"R - data section.R" - R file which will load required libraries, read in the data, and produce the plots, numerical summaries, and tests presented in the 'Materials and Methods' section of the paper.  

"R - results.R" - R file which will load required libraries, read in the data, and produce the plots and numerical summaries presented in the 'Results' section of the paper.  

"R - MCMC.R" - R file which will load required libraries, read in the data, and run the functions that were used to obtain the fitted values. The code also produces trace plots for the AR(2) model which were used to check the convergence of the algorithm.  

