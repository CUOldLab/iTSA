# iTSA
Implementing and testing iTSA single-temperature thermal profiling

This code repository supports the publication at DOI: XXXX/XXXXX for comparing samples of drug-treated and control proteins at a single temperature.

The input file (proteinGroups.txt) needed to reproduce the results is too large to store in GitHub, however you can access it at ftp://massive.ucsd.edu/MSV000083640/search/

Sourcing the file Rscript_EMpiricalBayes_TMT_expandedEdition.R with the proteinGroups.txt file located in the same directory should reproduce the results shown in the publication. This script requires ggvolcano.R and fdrfunctions.R to be present in the same working directory. 

Sub-sampling analysis can be replicated by sourcing the SubSampling.R script, with the file 20190204_Staurosporine_iTSA52m_KB-EBayesAnalysis.csv located in the same directory.
