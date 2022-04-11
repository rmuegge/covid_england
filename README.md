# Covid Study - National Lockdowns in England
Data and code to reproduce the analysis from "National lockdowns in England: The same restrictions for all, but how do the reductions in Covid-19 mortality risks differ geographically?".

Description of files:
"coviddata.csv" - contains the weekly number of deaths, expected deaths, and estimated risks for each LAD in the English mainland. 
"R - results.R" - R file which will load required libraries, read in the data, and produce the plots and numerical summaries presented in the paper. 
"R - MCMC.R" - R file which will load required libraries, read in the data, and run the functions that were used to obtain the fitted values. The code also produces trace plots which were used to check the convergence of the algorithm.  
