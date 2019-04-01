#Load up the results & sub-sample 

library(utils) #For combinatorics
library(dplyr)

setwd("C:/path/to/working/directory/")

FullResults <- read.table(file = "20190204_Staurosporine_iTSA52m_KB-EBayesAnalysis.csv", 
                         header = TRUE, 	
                         sep = ",",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))

#The vectors of vehicle & drug column headers from which we'll sub-sample
VehicleCols <- c("V_log2.i._TMT_1","V_log2.i._TMT_3","V_log2.i._TMT_5","V_log2.i._TMT_7","V_log2.i._TMT_9")
DrugCols <- c("D_log2.i._TMT_2","D_log2.i._TMT_4","D_log2.i._TMT_6","D_log2.i._TMT_8","D_log2.i._TMT_10")

#Sub-sample sets of 4 measurements from the groups of 5
# (This yields 5x5=25 pseudo-experiments)

#Create a column in the FullResults dataframe where we can store the frequency of
# a particular protein's 'discovery' by significance for each psuedo-experiment
FullResults$Rep4 = rep(0,nrow(FullResults))

print("Sub-sampling with 4 replicates")
for(i in 1:5){
  #Columns to subsample in the Vehicle
  SubVcols <- combn(VehicleCols,4)[,i]

  for(j in 1:5){
    #Columns to subsample in the Control
    SubDcols <- combn(DrugCols,4)[,j]
    
    #Our pseudo-experiment dataset
    subsample <- select(FullResults,append(SubVcols,SubDcols))
    design <- model.matrix(~ 0+factor(c(0,0,0,0,1,1,1,1)))
    colnames(design) <- c("V","D")		
    fit <- lmFit(subsample[append(SubVcols,SubDcols)], design)
    fit <- eBayes(fit)
    contrast <- makeContrasts("D-V", levels=design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2, trend = TRUE)
    ptopf <- topTableF(fit2, adjust="BH",genelist = FullResults[,"ID"],number=Inf)
    ptopf$rank <- 1:length(ptopf$F)
    ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
    
    metanms <- c("ID")
    ftest.tab <- topTableF(fit2, adjust="BH",genelist = FullResults[,metanms],number=Inf,p.value=1)
    tab.t1 <- topTable(fit2, adjust="BH",genelist = FullResults[,metanms], coef=1,p.value=1, number=Inf)

#    print(sum(tab.t1$adj.P.Val <= 0.05))
#    print(sum(tab.t1$adj.P.Val <= 0.01))
    
    #Grab the significant proteins, increment their Rep4 column
    #sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.05])
    sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.001])
    FullResults$Rep4 = FullResults$Rep4 + (FullResults[,"ID"] %in% sig_list)

  }

}

#Sub-sample sets of 3 measurements from the groups of 5
# (This yields 10x10=100 pseudo-experiments)

#Create a column in the FullResults dataframe where we can store the frequency of
# a particular protein's 'discovery' by significance for each psuedo-experiment
FullResults$Rep3 = rep(0,nrow(FullResults))
print("Sub-sampling with 3 replicates")

for(i in 1:10){
  #Columns to subsample in the Vehicle
  SubVcols <- combn(VehicleCols,3)[,i]
  
  for(j in 1:10){
    #Columns to subsample in the Control
    SubDcols <- combn(DrugCols,3)[,j]
    
    #Our pseudo-experiment dataset
    subsample <- select(FullResults,append(SubVcols,SubDcols))
    design <- model.matrix(~ 0+factor(c(0,0,0,1,1,1)))
    colnames(design) <- c("V","D")		
    fit <- lmFit(subsample[append(SubVcols,SubDcols)], design)
    fit <- eBayes(fit)
    contrast <- makeContrasts("D-V", levels=design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2, trend = TRUE)
    ptopf <- topTableF(fit2, adjust="BH",genelist = FullResults[,"ID"],number=Inf)
    ptopf$rank <- 1:length(ptopf$F)
    ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
    
    metanms <- c("ID")
    ftest.tab <- topTableF(fit2, adjust="BH",genelist = FullResults[,metanms],number=Inf,p.value=1)
    tab.t1 <- topTable(fit2, adjust="BH",genelist = FullResults[,metanms], coef=1,p.value=1, number=Inf)
    
#    print(sum(tab.t1$adj.P.Val <= 0.05))
    #    print(sum(tab.t1$adj.P.Val <= 0.01))
    
    #Grab the significant proteins, increment their Rep4 column
    #sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.05])
    sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.001])
    FullResults$Rep3 = FullResults$Rep3 + (FullResults[,"ID"] %in% sig_list)
        
  }
  
}

#Create a column in the FullResults dataframe where we can store the frequency of
# a particular protein's 'discovery' by significance for each psuedo-experiment
FullResults$Rep2 = rep(0,nrow(FullResults))
print("Sub-sampling with 2 replicates")

for(i in 1:10){
  #Columns to subsample in the Vehicle
  SubVcols <- combn(VehicleCols,2)[,i]
  
  for(j in 1:10){
    #Columns to subsample in the Control
    SubDcols <- combn(DrugCols,2)[,j]
    
    #Our pseudo-experiment dataset
    subsample <- select(FullResults,append(SubVcols,SubDcols))
    design <- model.matrix(~ 0+factor(c(0,0,1,1)))
    colnames(design) <- c("V","D")		
    fit <- lmFit(subsample[append(SubVcols,SubDcols)], design)
    fit <- eBayes(fit)
    contrast <- makeContrasts("D-V", levels=design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2, trend = TRUE)
    ptopf <- topTableF(fit2, adjust="BH",genelist = FullResults[,"ID"],number=Inf)
    ptopf$rank <- 1:length(ptopf$F)
    ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
    
    metanms <- c("ID")
    ftest.tab <- topTableF(fit2, adjust="BH",genelist = FullResults[,metanms],number=Inf,p.value=1)
    tab.t1 <- topTable(fit2, adjust="BH",genelist = FullResults[,metanms], coef=1,p.value=1, number=Inf)
    
#    print(sum(tab.t1$adj.P.Val <= 0.05))
    #    print(sum(tab.t1$adj.P.Val <= 0.01))
    
    #Grab the significant proteins, increment their Rep4 column
    #sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.05])
    sig_list <- as.character(tab.t1$ID[tab.t1$adj.P.Val <= 0.001])
    FullResults$Rep2 = FullResults$Rep2 + (FullResults[,"ID"] %in% sig_list)
    
  }
  
}


#Write this out (& don't over-write the original file)
write.csv(FullResults,file="20190201_postAnnotate_subsample.csv",row.names = F) ## looks good


#Make some Venn diagrams & UpSet plots to illustrate the information
library("UpSetR")
library("VennDiagram")

x <- list(
#"2 Replicates" = as.character(FullResults$ID[FullResults$Rep2 > 80]),
"3 Replicates" = as.character(FullResults$ID[FullResults$Rep3 > 80]),
"4 Replicates" = as.character(FullResults$ID[FullResults$Rep4 > 20]),
"5 Replicates" = as.character(FullResults$ID[FullResults$Sig_5FDR.FC == 1])
)

print(x)                    
overlap = calculate.overlap(x)
print(overlap)
venn.diagram(x,filename="venn.tiff",fill=c('red','blue','green'),category.names=c("","",""))

x <- list(
  "2 Replicates" = as.character(FullResults$ID[FullResults$Rep2 > 80]),
  "3 Replicates" = as.character(FullResults$ID[FullResults$Rep3 > 80]),
  "4 Replicates" = as.character(FullResults$ID[FullResults$Rep4 > 20]),
  "5 Replicates" = as.character(FullResults$ID[FullResults$Sig_5FDR.FC == 1])
)
pdf("Upset_freqfirst.pdf")
upset(fromList(x),order.by = 'freq')
dev.off()



