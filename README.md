# iTSA

An isothermal shift assay for proteome scale drug-target identification

Existing approaches for drug target identification employ complex workflows with low throughput. The isothermal shift assay (iTSA) is a proteome-wide method for identification of drug targets within lysates or living cells. iTSA uses a simplified experimental design with improved statistical power to detect small molecule binding targets across the proteome.
This code repository supports the publication  Ball, K.A., Webb, K.J., Coleman, S.J., Cozzolino, K.A., Jacobsen, J.,Jones, K.R., Stowell, M.H.B., Old, W.M. "An isothermal shift assay for proteome scale drug-target identification" Communications Biology (2019/2020).

The input file (proteinGroups.txt) needed to reproduce the results is too large to store in GitHub, however you can access it at ftp://massive.ucsd.edu/MSV000083640/search/

Sourcing the file Rscript_EmpiricalBayes_TMT_V20190402KB.R with the proteinGroups.txt file located in the same directory should reproduce the results shown in the publication. This script requires ggvolcano.R and fdrfunctions.R to be present in the same working directory. 

Sub-sampling analysis can be replicated by sourcing the SubSampling.R script, with the file 20190204_Staurosporine_iTSA52m_KB-EBayesAnalysis.csv located in the same directory. Alternatively, if you first run the above Empirical Bayes script yourself, use the file ending in "RankedEmpiricalBayesAnalysis.csv" as the input file located in the same directory. 

<a href="https://zenodo.org/badge/latestdoi/178905716"><img src="https://zenodo.org/badge/178905716.svg" alt="DOI"></a>
