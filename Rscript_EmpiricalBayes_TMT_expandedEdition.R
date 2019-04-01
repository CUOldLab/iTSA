############################
## lesson 1: A VERY brief introduction to R
############################
#
source("https://bioconductor.org/biocLite.R")
biocLite("BBmisc")
biocLite("limma")
#
# File&addresses: windows uses "\" but R only recognizes "/"
#
##Quick commands
#	getwd()
#	setwd()
#	library()
#rows <- [1,]
#columns <-[,1]
#rm(list=ls())	#clear all variables from workspace
#ls() #list all variables in workspace


#################################################################
## Emperical Bayes Statistical Anlalysis - script by Kerri Ball
################################################################

library(limma)
library(dplyr)
library(readxl)
library(lattice)
library(ggplot2)
library(plotly)
library(listviewer)
library(tidyverse)
library(BBmisc)
library(stringr)

############################################################################
##import data

#setwd is not needed if proteinGroups.txt is in this working directory
#setwd("C:/path/to/working/directory/") 

proteinGroups <- read.table(file = "proteinGroups.txt", 
		header = TRUE, 	
		sep = "\t",
		quote = "\"'",
           	dec = ".",
		numerals = c("warn.loss"),
		row.names = NULL,
		na.strings = c("NA","NaN","Infinite"))

##Filter data rows
data <- proteinGroups
colnames(data)
data <- subset(data, Only.identified.by.site != "+")
data <- subset(data, Reverse != "+")
data <- subset(data, Potential.contaminant != "+")
data <- subset(data, MS.MS.count > 2)
#data <- subset(data, Reporter.intensity.count.3 > 2)
#add unique filter + razor > 1
#data[ ,c(2:6,8:20,31:51,53:64)] <- list(NULL) #keep: 1=Protein.IDs; 7=Gene.names; 21-30=Reportr.intensity.corrected.x; 52=MS.MS.Count
#colnames(data)

#save data in subsets
data48 <- data
data52 <- data
data56 <- data
dataTC <- data
data.48C <- grep(".48C", colnames(data48))
data.52C <- grep(".52C", colnames(data52))
data.56C <- grep(".56C", colnames(data56))
data.TC <- grep(".TC", colnames(dataTC))


dataTMP <- data48[c("Protein.IDs", "Gene.names", "MS.MS.count")]
colnames(data48)
colnames(data48[,data.48C])
dataTMP48 <- cbind(dataTMP, data[,data.48C])
colnames(dataTMP48)
dataTMP48 <- subset(dataTMP48, Reporter.intensity.count.3.48C > 2)
data.48C <- grep("Reporter.intensity.corrected", colnames(dataTMP48))
colnames(dataTMP48[,data.48C])
data48 <- dataTMP48[c("Protein.IDs", "Gene.names", "MS.MS.count")]
data48 <- cbind(data48, dataTMP48[,data.48C])
colnames(data48)


data <- data48
#Define TMT data columns
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))       # identify column numbers that contain intensity data
colnames(data[,data.TMT])
#log2 transform Reporter.intensity
data[,data.TMT] <- log(data[,data.TMT],base=2)
View(data)

#Create Unique ID column
data.metanms <- grep("Protein.IDs", colnames(data)) #use specific columns names to identify protein id columns
colnames(data)[data.metanms] <-  "ID"	#change column names
data$Protein.IDs <- data$ID
data$last <- (regexpr(';', as.character(data$Protein.IDs))-1) # if no ";" then = -1
data$ID.i <- ifelse(data$last > 0, (substr(as.character(data$Protein.IDs),1,data$last)), (as.character(data$Protein.IDs)))
data$last <- (regexpr('-', as.character(data$ID.i))-1) # if no ";" then = -1
data$ID <- ifelse(data$last > 0, (substr(as.character(data$ID.i),1,data$last)), (as.character(data$ID.i)))
data$ID.i <- NULL
data$last <- NULL
data <- data[order(-data$MS.MS.count),] #sort by MS.MS. count (largest to smallest)
data$dup <- duplicated(data[,1]) #first occurance is False, second and subsequent say True
data <- subset(data, data$dup != TRUE)
data$dup <- NULL

#Look at data
#pdf("BoxPlot_BeforeQuantileNormalization.pdf") #remove # signs for pdf to print in file
boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE, 
	frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
								"gray75","gray75","gray75","gray75",
								"gray70","gray70","gray70","gray70",
								"gray65","gray65","gray65","gray65",
								"gray60","gray60","gray60","gray60",
								"gray55","gray55","gray55","gray55",
								"gray50","gray50","gray50","gray50",
								"gray45","gray45","gray45","gray45",
								"gray40","gray40","gray40","gray40",
								"gray35","gray35","gray35","gray35"))
# Close the pdf file
#dev.off() 

data.norm <- data
colnames(data.norm)
##############################################################
###Remove any channels that should not be included in the normalization and analysis
#data.norm$Reporter.intensity.corrected.11 <- NULL
data.TMT <- grep("Reporter.intensity.corrected", colnames(data.norm))       # identify column numbers that contain intensity data
colnames(data.norm[,data.TMT])
data.norm[,data.TMT] <-  normalizeBetweenArrays(as.matrix(data.norm[,data.TMT]),method="quantile" )

boxplot((data.norm[,data.TMT]), main = "After Quantile Normalization", names = TRUE, 
	frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
								"gray75","gray75","gray75","gray75",
								"gray70","gray70","gray70","gray70",
								"gray65","gray65","gray65","gray65",
								"gray60","gray60","gray60","gray60",
								"gray55","gray55","gray55","gray55",
								"gray50","gray50","gray50","gray50",
								"gray45","gray45","gray45","gray45",
								"gray40","gray40","gray40","gray40",
								"gray35","gray35","gray35","gray35"))


View(data.norm)
colnames(data.norm)
############################################################################################
###Rename columns according to experimental design and set up model.matrix appropriately
data.norm <- setNames(data.norm, c("ID","Gene.names", "MS.MS.count",
                                   "V_log2.i._TMT_1", 
                                   "D_log2.i._TMT_2",  
                                   "V_log2.i._TMT_3",  
                                   "D_log2.i._TMT_4", 
                                   "V_log2.i._TMT_5",  
                                   "D_log2.i._TMT_6",  
                                   "V_log2.i._TMT_7", 
                                   "D_log2.i._TMT_8",  
                                   "V_log2.i._TMT_9",  
                                   "D_log2.i._TMT_10",
                                   
                                    "Protein.IDs" ))


################################################################################################################
###change file name
#write.csv(data.norm,file="20190201_quantileNorm_Staurosporine_LyTSA52_noEbayes.csv",row.names = F)


design <- model.matrix(~ 0+factor(c(0,1,0,1,0,1,0,1,0,1)))
colnames(design) <- c("V","D")		

fit <- lmFit(data.norm[,data.TMT], design)
fit <- eBayes(fit) 	##Apply empirical Bayes smoothing to the standard errors.
contrast <- makeContrasts("D-V", levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
ptopf <- topTableF(fit2, adjust="BH",genelist = data[,"ID"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val

###################################################################################################
###coef=1 is the first contrast in cont.wt <- makeContrasts #coef=0 dosn't exist but will not error
###adjust number of tables to number of contrasts
metanms <- c("ID")
ftest.tab <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1)
tab.t1 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=1,p.value=1, number=Inf)
#tab.t2 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=2,p.value=1, number=Inf)
#tab.t3 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=3,p.value=1, number=Inf)

tab.t1[5:6] <- -log(tab.t1[5:6],base=10)
#tab.t2[5:6] <- -log(tab.t2[5:6],base=10)
#tab.t3[5:6] <- -log(tab.t3[5:6],base=10)


##P.Val volcano plots
tab1 <- tab.t1
Xdata <- tab1$logFC     			# x-axis = logFC 
Ydata <- tab1$P.Val         			# y-axis = P.Val  
plot(Xdata, Ydata,            	# plot the variables 
 	main = "",
	xlab="logFC",        		# x-axis label 
	ylab="-Log10(P.Val)")         # y-axis label

sum(tab.t1$adj.P.Val >= -log10(0.05))
sum(tab.t1$adj.P.Val >= -log10(0.01))

###Change column names and merge EB results with nomalized data
colnames(tab.t1)
tab.t1[ ,c(3,4,7)] <- list(NULL)
tab.t1 <- setNames(tab.t1, c("ID","Log2FC","Log10P.Value","Log10adj.P.Val")) 
data.analysis <- merge(data.norm,	tab.t1, by = "ID" ,all = T, sort=F)
data.analysis <- data.analysis[order(-data.analysis$Log10adj.P.Val),]
View(data.analysis)

#write.csv(data.analysis,file="20190201_preFDR.csv",row.names = F) ## looks good


#Create FDR and FC columns
data.analysis$FDR0.001Sig <- ifelse(data.analysis$Log10adj.P.Val>=3, 1, 0)
data.analysis$FDR0.01Sig <- ifelse(data.analysis$Log10adj.P.Val>=2, 1, 0)
data.analysis$FDR0.05Sig <- ifelse(data.analysis$Log10adj.P.Val>=1.3, 1, 0)
data.analysis$absFC <- abs(data.analysis$Log2FC)
plot(data.analysis$absFC, data.analysis$Log10adj.P.Val)
quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile <- quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile80 <- data.frame(percentile[1])
percentile90 <- data.frame(percentile[2])
data.analysis$FC0.10Sig <- ifelse(data.analysis$absFC>= (percentile90[1,1]), 1, 0)
data.analysis$FC0.20Sig <- ifelse(data.analysis$absFC>= (percentile80[1,1]), 1, 0)

data.analysis$Sig_0.1FDR.FC <- ifelse(data.analysis$FDR0.001Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_1FDR.FC <- ifelse(data.analysis$FDR0.01Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_5FDR.FC <- ifelse(data.analysis$FDR0.05Sig >= 1, (ifelse(data.analysis$FC0.20Sig >= 1, 1, 0)), 0)
data.analysis$Sig_FDR.FC <- data.analysis$Sig_1FDR.FC + data.analysis$Sig_5FDR.FC + data.analysis$Sig_0.1FDR.FC

data.analysis$Direction <- ifelse(data.analysis$Log2FC >= 0, (1*data.analysis$Sig_FDR.FC), (-1*data.analysis$Sig_FDR.FC))
data.analysis <- data.analysis[order(-data.analysis$Sig_FDR.FC),]
data.analysis$rank <- 1:length(data.analysis$Sig_FDR.FC)
plot(data.analysis$rank, data.analysis$Log10adj.P.Val)

#write.csv(data.analysis,file="20190201_postFDR_preGeneName.csv",row.names = F) ## looks good


#Create gene.name column
data.analysis$last <- (regexpr(';', as.character(data.analysis$Gene.names))-1) # if no ";" then = -1
data.analysis$Gene.name <- ifelse(data.analysis$last > 0, (substr(as.character(data.analysis$Gene.names),1,data.analysis$last)), (as.character(data.analysis$Gene.names)))
#if no gene name is present then use ID
data.analysis$last <- nchar(data.analysis$Gene.name)
data.analysis$name <- data.analysis$Gene.name
data.analysis$Gene.name <- ifelse (data.analysis$last == 0, as.character(data.analysis$ID), data.analysis$name)
data.analysis$last <- NULL
data.analysis$name <- NULL
colnames(data.analysis)

#write.csv(data.analysis,file="20190201_posteneName_preAnnotate.csv",row.names = F) ## looks good



#annotate kinases
##import annotation file
annotation <- read.table(file = "UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt", 
                         header = TRUE, 	
                         sep = "\t",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))

test2 <- merge(data.analysis,	annotation, by = "ID" ,all = T, sort=F)
test3 <- subset(test2, Sig_FDR.FC >= 0)

write.csv(test3,file="20190201_postAnnotate_preSort.csv",row.names = F) ## looks good
data.analysis <- test3

#reorganize columns
V_TMT <- grep("V_log2.i._TMT_", colnames(data.analysis)) 
D_TMT <- grep("D_log2.i._TMT_", colnames(data.analysis)) 

data.sort <- data.analysis[c("ID", "Gene.name", "rank", "Direction", "Sig_FDR.FC" , "Is_kinase", "Log2FC", "Log10P.Value", "Log10adj.P.Val", "Protein.IDs", "Gene.names", "MS.MS.count")]

test2 <- cbind(data.sort, data.analysis[,V_TMT], data.analysis[,D_TMT])
colnames(test2)
data.analysis <- test2



####prepare volcano plot viewer html files
results <- data.analysis
#colnames(results)
results$adj.P.Val <- 10^(-1*results$Log10adj.P.Val)
sum(results$adj.P.Val <= 0.05)
results <- results[order(results$rank),] #sort by rank (smallest to largest)
View(results)

##within the ggvolcano plot there is a directory to the fdrfunction.R (check that this address is correct)
source('./fdrfunctions.R')
source('ggvolcano.R')

myggma <- results %>% ggmaplot(xvar = Log2FC, yvar = Log10P.Value, adjpvalvar = adj.P.Val  )
ss.gg <- myggma$maplot.ggplot
ss.plotly <- myggma$maplot.plotly

##this is the code that opens the viewer - give it a minute!
ss.plotly
##
##

data.analysis$Gene.name <- paste(data.analysis$Gene.name, ";", sep="")
data.analysis$Gene.names <- paste(data.analysis$Gene.names, ";", sep="")
################################################################################################################
###change file name
write.csv(data.analysis,file="20190204_Staurosporine_iTSA48m_KB-EBayesAnalysis.csv",row.names = F)
data.analysis48 <- data.analysis

###############################################################################################################
colnames(data52)
colnames(data52[,data.52C])
dataTMP <- data52[c("Protein.IDs", "Gene.names", "MS.MS.count")]

dataTMP52 <- cbind(dataTMP, data52[,data.52C])
colnames(dataTMP52)
dataTMP52 <- subset(dataTMP52, Reporter.intensity.count.3.52C > 2)
data.52C <- grep("Reporter.intensity.corrected", colnames(dataTMP52))
colnames(dataTMP52[,data.52C])
data52 <- dataTMP52[c("Protein.IDs", "Gene.names", "MS.MS.count")]
data52 <- cbind(data52, dataTMP52[,data.52C])
colnames(data52)


data <- data52
#Define TMT data columns
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))       # identify column numbers that contain intensity data
colnames(data[,data.TMT])
#log2 transform Reporter.intensity
data[,data.TMT] <- log(data[,data.TMT],base=2)
#View(data)

#Create Unique ID column
data.metanms <- grep("Protein.IDs", colnames(data)) #use specific columns names to identify protein id columns
colnames(data)[data.metanms] <-  "ID"	#change column names
data$Protein.IDs <- data$ID
data$last <- (regexpr(';', as.character(data$Protein.IDs))-1) # if no ";" then = -1
data$ID.i <- ifelse(data$last > 0, (substr(as.character(data$Protein.IDs),1,data$last)), (as.character(data$Protein.IDs)))
data$last <- (regexpr('-', as.character(data$ID.i))-1) # if no ";" then = -1
data$ID <- ifelse(data$last > 0, (substr(as.character(data$ID.i),1,data$last)), (as.character(data$ID.i)))
data$ID.i <- NULL
data$last <- NULL
data <- data[order(-data$MS.MS.count),] #sort by MS.MS. count (largest to smallest)
data$dup <- duplicated(data[,1]) #first occurance is False, second and subsequent say True
data <- subset(data, data$dup != TRUE)
data$dup <- NULL

#Look at data
#pdf("BoxPlot_BeforeQuantileNormalization.pdf") #remove # signs for pdf to print in file
boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE, 
        frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
                                              "gray75","gray75","gray75","gray75",
                                              "gray70","gray70","gray70","gray70",
                                              "gray65","gray65","gray65","gray65",
                                              "gray60","gray60","gray60","gray60",
                                              "gray55","gray55","gray55","gray55",
                                              "gray50","gray50","gray50","gray50",
                                              "gray45","gray45","gray45","gray45",
                                              "gray40","gray40","gray40","gray40",
                                              "gray35","gray35","gray35","gray35"))
# Close the pdf file
#dev.off() 

data.norm <- data
colnames(data.norm)
##############################################################
###Remove any channels that should not be included in the normalization and analysis
#data.norm$Reporter.intensity.corrected.11 <- NULL
data.TMT <- grep("Reporter.intensity.corrected", colnames(data.norm))       # identify column numbers that contain intensity data
colnames(data.norm[,data.TMT])
data.norm[,data.TMT] <-  normalizeBetweenArrays(as.matrix(data.norm[,data.TMT]),method="quantile" )

boxplot((data.norm[,data.TMT]), main = "After Quantile Normalization", names = TRUE, 
        frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
                                              "gray75","gray75","gray75","gray75",
                                              "gray70","gray70","gray70","gray70",
                                              "gray65","gray65","gray65","gray65",
                                              "gray60","gray60","gray60","gray60",
                                              "gray55","gray55","gray55","gray55",
                                              "gray50","gray50","gray50","gray50",
                                              "gray45","gray45","gray45","gray45",
                                              "gray40","gray40","gray40","gray40",
                                              "gray35","gray35","gray35","gray35"))


View(data.norm)
colnames(data.norm)
############################################################################################
###Rename columns according to experimental design and set up model.matrix appropriately
data.norm <- setNames(data.norm, c("ID","Gene.names", "MS.MS.count",
                                   "V_log2.i._TMT_1", 
                                   "D_log2.i._TMT_2",  
                                   "V_log2.i._TMT_3",  
                                   "D_log2.i._TMT_4", 
                                   "V_log2.i._TMT_5",  
                                   "D_log2.i._TMT_6",  
                                   "V_log2.i._TMT_7", 
                                   "D_log2.i._TMT_8",  
                                   "V_log2.i._TMT_9",  
                                   "D_log2.i._TMT_10",
                                   
                                   "Protein.IDs" ))


design <- model.matrix(~ 0+factor(c(0,1,0,1,0,1,0,1,0,1)))
colnames(design) <- c("V","D")		

fit <- lmFit(data.norm[,data.TMT], design)
fit <- eBayes(fit) 	##Apply empirical Bayes smoothing to the standard errors.
contrast <- makeContrasts("D-V", levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
ptopf <- topTableF(fit2, adjust="BH",genelist = data[,"ID"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val

###################################################################################################
###coef=1 is the first contrast in cont.wt <- makeContrasts #coef=0 dosn't exist but will not error
###adjust number of tables to number of contrasts
metanms <- c("ID")
ftest.tab <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1)
tab.t1 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=1,p.value=1, number=Inf)
#tab.t2 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=2,p.value=1, number=Inf)
#tab.t3 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=3,p.value=1, number=Inf)

tab.t1[5:6] <- -log(tab.t1[5:6],base=10)
#tab.t2[5:6] <- -log(tab.t2[5:6],base=10)
#tab.t3[5:6] <- -log(tab.t3[5:6],base=10)


##P.Val volcano plots
tab1 <- tab.t1
Xdata <- tab1$logFC     			# x-axis = logFC 
Ydata <- tab1$P.Val         			# y-axis = P.Val  
plot(Xdata, Ydata,            	# plot the variables 
     main = "",
     xlab="logFC",        		# x-axis label 
     ylab="-Log10(P.Val)")         # y-axis label

sum(tab.t1$adj.P.Val >= -log10(0.05))
sum(tab.t1$adj.P.Val >= -log10(0.01))

###Change column names and merge EB results with nomalized data
colnames(tab.t1)
tab.t1[ ,c(3,4,7)] <- list(NULL)
tab.t1 <- setNames(tab.t1, c("ID","Log2FC","Log10P.Value","Log10adj.P.Val")) 
data.analysis <- merge(data.norm,	tab.t1, by = "ID" ,all = T, sort=F)
data.analysis <- data.analysis[order(-data.analysis$Log10adj.P.Val),]
View(data.analysis)

#write.csv(data.analysis,file="20190201_preFDR.csv",row.names = F) ## looks good


#Create FDR and FC columns
data.analysis$FDR0.001Sig <- ifelse(data.analysis$Log10adj.P.Val>=3, 1, 0)
data.analysis$FDR0.01Sig <- ifelse(data.analysis$Log10adj.P.Val>=2, 1, 0)
data.analysis$FDR0.05Sig <- ifelse(data.analysis$Log10adj.P.Val>=1.3, 1, 0)
data.analysis$absFC <- abs(data.analysis$Log2FC)
plot(data.analysis$absFC, data.analysis$Log10adj.P.Val)
quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile <- quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile80 <- data.frame(percentile[1])
percentile90 <- data.frame(percentile[2])
data.analysis$FC0.10Sig <- ifelse(data.analysis$absFC>= (percentile90[1,1]), 1, 0)
data.analysis$FC0.20Sig <- ifelse(data.analysis$absFC>= (percentile80[1,1]), 1, 0)

data.analysis$Sig_0.1FDR.FC <- ifelse(data.analysis$FDR0.001Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_1FDR.FC <- ifelse(data.analysis$FDR0.01Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_5FDR.FC <- ifelse(data.analysis$FDR0.05Sig >= 1, (ifelse(data.analysis$FC0.20Sig >= 1, 1, 0)), 0)
data.analysis$Sig_FDR.FC <- data.analysis$Sig_1FDR.FC + data.analysis$Sig_5FDR.FC + data.analysis$Sig_0.1FDR.FC

data.analysis$Direction <- ifelse(data.analysis$Log2FC >= 0, (1*data.analysis$Sig_FDR.FC), (-1*data.analysis$Sig_FDR.FC))
data.analysis <- data.analysis[order(-data.analysis$Sig_FDR.FC),]
data.analysis$rank <- 1:length(data.analysis$Sig_FDR.FC)
plot(data.analysis$rank, data.analysis$Log10adj.P.Val)

#write.csv(data.analysis,file="20190201_postFDR_preGeneName.csv",row.names = F) ## looks good


#Create gene.name column
data.analysis$last <- (regexpr(';', as.character(data.analysis$Gene.names))-1) # if no ";" then = -1
data.analysis$Gene.name <- ifelse(data.analysis$last > 0, (substr(as.character(data.analysis$Gene.names),1,data.analysis$last)), (as.character(data.analysis$Gene.names)))
#if no gene name is present then use ID
data.analysis$last <- nchar(data.analysis$Gene.name)
data.analysis$name <- data.analysis$Gene.name
data.analysis$Gene.name <- ifelse (data.analysis$last == 0, as.character(data.analysis$ID), data.analysis$name)
data.analysis$last <- NULL
data.analysis$name <- NULL
colnames(data.analysis)




#annotate kinases
##import annotation file
annotation <- read.table(file = "UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt", 
                         header = TRUE, 	
                         sep = "\t",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))

test2 <- merge(data.analysis,	annotation, by = "ID" ,all = T, sort=F)
test3 <- subset(test2, Sig_FDR.FC >= 0)

write.csv(test3,file="20190201_postAnnotate_preSort.csv",row.names = F) ## looks good
data.analysis <- test3

#reorganize columns
V_TMT <- grep("V_log2.i._TMT_", colnames(data.analysis)) 
D_TMT <- grep("D_log2.i._TMT_", colnames(data.analysis)) 

data.sort <- data.analysis[c("ID", "Gene.name", "rank", "Direction", "Sig_FDR.FC" , "Is_kinase", "Log2FC", "Log10P.Value", "Log10adj.P.Val", "Protein.IDs", "Gene.names", "MS.MS.count")]

test2 <- cbind(data.sort, data.analysis[,V_TMT], data.analysis[,D_TMT])
colnames(test2)
data.analysis <- test2



####prepare volcano plot viewer html files
results <- data.analysis
#colnames(results)
results$adj.P.Val <- 10^(-1*results$Log10adj.P.Val)
sum(results$adj.P.Val <= 0.05)
results <- results[order(results$rank),] #sort by rank (smallest to largest)
View(results)

##within the ggvolcano plot there is a directory to the fdrfunction.R (check that this address is correct)
source('./fdrfunctions.R')
source('ggvolcano.R')

myggma <- results %>% ggmaplot(xvar = Log2FC, yvar = Log10P.Value, adjpvalvar = adj.P.Val  )
ss.gg <- myggma$maplot.ggplot
ss.plotly <- myggma$maplot.plotly

##this is the code that opens the viewer - give it a minute!
ss.plotly
##
##

data.analysis$Gene.name <- paste(data.analysis$Gene.name, ";", sep="")
data.analysis$Gene.names <- paste(data.analysis$Gene.names, ";", sep="")
################################################################################################################
###change file name
write.csv(data.analysis,file="20190204_Staurosporine_iTSA52m_KB-EBayesAnalysis.csv",row.names = F)
data.analysis52 <- data.analysis

###############################################################################################################
colnames(data56)
colnames(data56[,data.56C])
dataTMP <- data56[c("Protein.IDs", "Gene.names", "MS.MS.count")]

dataTMP56 <- cbind(dataTMP, data56[,data.56C])
colnames(dataTMP56)
dataTMP56 <- subset(dataTMP56, Reporter.intensity.count.3.56C > 2)
data.56C <- grep("Reporter.intensity.corrected", colnames(dataTMP56))
colnames(dataTMP56[,data.56C])
data56 <- dataTMP56[c("Protein.IDs", "Gene.names", "MS.MS.count")]
data56 <- cbind(data56, dataTMP56[,data.56C])
colnames(data56)


data <- data56
#Define TMT data columns
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))       # identify column numbers that contain intensity data
colnames(data[,data.TMT])
#log2 transform Reporter.intensity
data[,data.TMT] <- log(data[,data.TMT],base=2)
#View(data)

#Create Unique ID column
data.metanms <- grep("Protein.IDs", colnames(data)) #use specific columns names to identify protein id columns
colnames(data)[data.metanms] <-  "ID"	#change column names
data$Protein.IDs <- data$ID
data$last <- (regexpr(';', as.character(data$Protein.IDs))-1) # if no ";" then = -1
data$ID.i <- ifelse(data$last > 0, (substr(as.character(data$Protein.IDs),1,data$last)), (as.character(data$Protein.IDs)))
data$last <- (regexpr('-', as.character(data$ID.i))-1) # if no ";" then = -1
data$ID <- ifelse(data$last > 0, (substr(as.character(data$ID.i),1,data$last)), (as.character(data$ID.i)))
data$ID.i <- NULL
data$last <- NULL
data <- data[order(-data$MS.MS.count),] #sort by MS.MS. count (largest to smallest)
data$dup <- duplicated(data[,1]) #first occurance is False, second and subsequent say True
data <- subset(data, data$dup != TRUE)
data$dup <- NULL

#Look at data
#pdf("BoxPlot_BeforeQuantileNormalization.pdf") #remove # signs for pdf to print in file
boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE, 
        frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
                                              "gray75","gray75","gray75","gray75",
                                              "gray70","gray70","gray70","gray70",
                                              "gray65","gray65","gray65","gray65",
                                              "gray60","gray60","gray60","gray60",
                                              "gray55","gray55","gray55","gray55",
                                              "gray50","gray50","gray50","gray50",
                                              "gray45","gray45","gray45","gray45",
                                              "gray40","gray40","gray40","gray40",
                                              "gray35","gray35","gray35","gray35"))
# Close the pdf file
#dev.off() 

data.norm <- data
colnames(data.norm)
##############################################################
###Remove any channels that should not be included in the normalization and analysis
#data.norm$Reporter.intensity.corrected.11 <- NULL
data.TMT <- grep("Reporter.intensity.corrected", colnames(data.norm))       # identify column numbers that contain intensity data
colnames(data.norm[,data.TMT])
data.norm[,data.TMT] <-  normalizeBetweenArrays(as.matrix(data.norm[,data.TMT]),method="quantile" )

boxplot((data.norm[,data.TMT]), main = "After Quantile Normalization", names = TRUE, 
        frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
                                              "gray75","gray75","gray75","gray75",
                                              "gray70","gray70","gray70","gray70",
                                              "gray65","gray65","gray65","gray65",
                                              "gray60","gray60","gray60","gray60",
                                              "gray55","gray55","gray55","gray55",
                                              "gray50","gray50","gray50","gray50",
                                              "gray45","gray45","gray45","gray45",
                                              "gray40","gray40","gray40","gray40",
                                              "gray35","gray35","gray35","gray35"))


View(data.norm)
colnames(data.norm)
############################################################################################
###Rename columns according to experimental design and set up model.matrix appropriately
data.norm <- setNames(data.norm, c("ID","Gene.names", "MS.MS.count",
                                   "V_log2.i._TMT_1", 
                                   "D_log2.i._TMT_2",  
                                   "V_log2.i._TMT_3",  
                                   "D_log2.i._TMT_4", 
                                   "V_log2.i._TMT_5",  
                                   "D_log2.i._TMT_6",  
                                   "V_log2.i._TMT_7", 
                                   "D_log2.i._TMT_8",  
                                   "V_log2.i._TMT_9",  
                                   "D_log2.i._TMT_10",
                                   
                                   "Protein.IDs" ))


design <- model.matrix(~ 0+factor(c(0,1,0,1,0,1,0,1,0,1)))
colnames(design) <- c("V","D")		

fit <- lmFit(data.norm[,data.TMT], design)
fit <- eBayes(fit) 	##Apply empirical Bayes smoothing to the standard errors.
contrast <- makeContrasts("D-V", levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
ptopf <- topTableF(fit2, adjust="BH",genelist = data[,"ID"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val

###################################################################################################
###coef=1 is the first contrast in cont.wt <- makeContrasts #coef=0 dosn't exist but will not error
###adjust number of tables to number of contrasts
metanms <- c("ID")
ftest.tab <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1)
tab.t1 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=1,p.value=1, number=Inf)
#tab.t2 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=2,p.value=1, number=Inf)
#tab.t3 <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=3,p.value=1, number=Inf)

tab.t1[5:6] <- -log(tab.t1[5:6],base=10)
#tab.t2[5:6] <- -log(tab.t2[5:6],base=10)
#tab.t3[5:6] <- -log(tab.t3[5:6],base=10)


##P.Val volcano plots
tab1 <- tab.t1
Xdata <- tab1$logFC     			# x-axis = logFC 
Ydata <- tab1$P.Val         			# y-axis = P.Val  
plot(Xdata, Ydata,            	# plot the variables 
     main = "",
     xlab="logFC",        		# x-axis label 
     ylab="-Log10(P.Val)")         # y-axis label

sum(tab.t1$adj.P.Val >= -log10(0.05))
sum(tab.t1$adj.P.Val >= -log10(0.01))

###Change column names and merge EB results with nomalized data
colnames(tab.t1)
tab.t1[ ,c(3,4,7)] <- list(NULL)
tab.t1 <- setNames(tab.t1, c("ID","Log2FC","Log10P.Value","Log10adj.P.Val")) 
data.analysis <- merge(data.norm,	tab.t1, by = "ID" ,all = T, sort=F)
data.analysis <- data.analysis[order(-data.analysis$Log10adj.P.Val),]
View(data.analysis)

#Create FDR and FC columns
data.analysis$FDR0.001Sig <- ifelse(data.analysis$Log10adj.P.Val>=3, 1, 0)
data.analysis$FDR0.01Sig <- ifelse(data.analysis$Log10adj.P.Val>=2, 1, 0)
data.analysis$FDR0.05Sig <- ifelse(data.analysis$Log10adj.P.Val>=1.3, 1, 0)
data.analysis$absFC <- abs(data.analysis$Log2FC)
plot(data.analysis$absFC, data.analysis$Log10adj.P.Val)
quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile <- quantile(data.analysis$absFC, probs = c(80,90)/100)
percentile80 <- data.frame(percentile[1])
percentile90 <- data.frame(percentile[2])
data.analysis$FC0.10Sig <- ifelse(data.analysis$absFC>= (percentile90[1,1]), 1, 0)
data.analysis$FC0.20Sig <- ifelse(data.analysis$absFC>= (percentile80[1,1]), 1, 0)

data.analysis$Sig_0.1FDR.FC <- ifelse(data.analysis$FDR0.001Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_1FDR.FC <- ifelse(data.analysis$FDR0.01Sig >= 1, (ifelse(data.analysis$FC0.10Sig >= 1, 1, 0)), 0)
data.analysis$Sig_5FDR.FC <- ifelse(data.analysis$FDR0.05Sig >= 1, (ifelse(data.analysis$FC0.20Sig >= 1, 1, 0)), 0)
data.analysis$Sig_FDR.FC <- data.analysis$Sig_1FDR.FC + data.analysis$Sig_5FDR.FC + data.analysis$Sig_0.1FDR.FC

data.analysis$Direction <- ifelse(data.analysis$Log2FC >= 0, (1*data.analysis$Sig_FDR.FC), (-1*data.analysis$Sig_FDR.FC))
data.analysis <- data.analysis[order(-data.analysis$Sig_FDR.FC),]
data.analysis$rank <- 1:length(data.analysis$Sig_FDR.FC)
plot(data.analysis$rank, data.analysis$Log10adj.P.Val)

#Create gene.name column
data.analysis$last <- (regexpr(';', as.character(data.analysis$Gene.names))-1) # if no ";" then = -1
data.analysis$Gene.name <- ifelse(data.analysis$last > 0, (substr(as.character(data.analysis$Gene.names),1,data.analysis$last)), (as.character(data.analysis$Gene.names)))
#if no gene name is present then use ID
data.analysis$last <- nchar(data.analysis$Gene.name)
data.analysis$name <- data.analysis$Gene.name
data.analysis$Gene.name <- ifelse (data.analysis$last == 0, as.character(data.analysis$ID), data.analysis$name)
data.analysis$last <- NULL
data.analysis$name <- NULL
colnames(data.analysis)




#annotate kinases
##import annotation file
annotation <- read.table(file = "UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt", 
                         header = TRUE, 	
                         sep = "\t",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))

test2 <- merge(data.analysis,	annotation, by = "ID" ,all = T, sort=F)
test3 <- subset(test2, Sig_FDR.FC >= 0)

write.csv(test3,file="20190201_postAnnotate_preSort.csv",row.names = F) ## looks good
data.analysis <- test3

#reorganize columns
V_TMT <- grep("V_log2.i._TMT_", colnames(data.analysis)) 
D_TMT <- grep("D_log2.i._TMT_", colnames(data.analysis)) 

data.sort <- data.analysis[c("ID", "Gene.name", "rank", "Direction", "Sig_FDR.FC" , "Is_kinase", "Log2FC", "Log10P.Value", "Log10adj.P.Val", "Protein.IDs", "Gene.names", "MS.MS.count")]

test2 <- cbind(data.sort, data.analysis[,V_TMT], data.analysis[,D_TMT])
colnames(test2)
data.analysis <- test2



####prepare volcano plot viewer html files
results <- data.analysis
#colnames(results)
results$adj.P.Val <- 10^(-1*results$Log10adj.P.Val)
sum(results$adj.P.Val <= 0.05)
results <- results[order(results$rank),] #sort by rank (smallest to largest)
View(results)

##within the ggvolcano plot there is a directory to the fdrfunction.R (check that this address is correct)
source('./fdrfunctions.R')
source('ggvolcano.R')

myggma <- results %>% ggmaplot(xvar = Log2FC, yvar = Log10P.Value, adjpvalvar = adj.P.Val  )
ss.gg <- myggma$maplot.ggplot
ss.plotly <- myggma$maplot.plotly

##this is the code that opens the viewer - give it a minute!
ss.plotly
##
##

data.analysis$Gene.name <- paste(data.analysis$Gene.name, ";", sep="")
data.analysis$Gene.names <- paste(data.analysis$Gene.names, ";", sep="")
################################################################################################################
###change file name
write.csv(data.analysis,file="20190204_Staurosporine_iTSA56m_KB-EBayesAnalysis.csv",row.names = F)
data.analysis56 <- data.analysis
###############################################################################################################



dataTMP <- dataTC[c("Protein.IDs", "Gene.names", "MS.MS.count")]
dataTMPTC <- cbind(dataTMP, data[,data.TC])
colnames(dataTMPTC)
dataTMPTC <- subset(dataTMPTC, Reporter.intensity.count.3.TC > 2)
data.TC <- grep("Reporter.intensity.corrected", colnames(dataTMPTC))
colnames(dataTMPTC[,data.TC])
dataTC <- dataTMPTC[c("Protein.IDs", "Gene.names", "MS.MS.count")]
dataTC <- cbind(dataTC, dataTMPTC[,data.TC])
colnames(dataTC)


data <- dataTC
#Define TMT data columns
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))       # identify column numbers that contain intensity data
colnames(data[,data.TMT])
#log2 transform Reporter.intensity
data[,data.TMT] <- log(data[,data.TMT],base=2)
View(data)

#Create Unique ID column
data.metanms <- grep("Protein.IDs", colnames(data)) #use specific columns names to identify protein id columns
colnames(data)[data.metanms] <-  "ID"	#change column names
data$Protein.IDs <- data$ID
data$last <- (regexpr(';', as.character(data$Protein.IDs))-1) # if no ";" then = -1
data$ID.i <- ifelse(data$last > 0, (substr(as.character(data$Protein.IDs),1,data$last)), (as.character(data$Protein.IDs)))
data$last <- (regexpr('-', as.character(data$ID.i))-1) # if no ";" then = -1
data$ID <- ifelse(data$last > 0, (substr(as.character(data$ID.i),1,data$last)), (as.character(data$ID.i)))
data$ID.i <- NULL
data$last <- NULL
data <- data[order(-data$MS.MS.count),] #sort by MS.MS. count (largest to smallest)
data$dup <- duplicated(data[,1]) #first occurance is False, second and subsequent say True
data <- subset(data, data$dup != TRUE)
data$dup <- NULL

#Look at data
#pdf("BoxPlot_BeforeQuantileNormalization.pdf") #remove # signs for pdf to print in file
boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE, 
        frame = FALSE, notch = TRUE, col = c(	"gray80","gray80","gray80","gray80",
                                              "gray75","gray75","gray75","gray75",
                                              "gray70","gray70","gray70","gray70",
                                              "gray65","gray65","gray65","gray65",
                                              "gray60","gray60","gray60","gray60",
                                              "gray55","gray55","gray55","gray55",
                                              "gray50","gray50","gray50","gray50",
                                              "gray45","gray45","gray45","gray45",
                                              "gray40","gray40","gray40","gray40",
                                              "gray35","gray35","gray35","gray35"))
# Close the pdf file
#dev.off() 

data.analysis <- data
colnames(data.analysis)


#Create gene.name column
data.analysis$last <- (regexpr(';', as.character(data.analysis$Gene.names))-1) # if no ";" then = -1
data.analysis$Gene.name <- ifelse(data.analysis$last > 0, (substr(as.character(data.analysis$Gene.names),1,data.analysis$last)), (as.character(data.analysis$Gene.names)))
#if no gene name is present then use ID
data.analysis$last <- nchar(data.analysis$Gene.name)
data.analysis$name <- data.analysis$Gene.name
data.analysis$Gene.name <- ifelse (data.analysis$last == 0, as.character(data.analysis$ID), data.analysis$name)
data.analysis$last <- NULL
data.analysis$name <- NULL
colnames(data.analysis)

View(data.analysis)
data.analysis$Gene.name <- paste(data.analysis$Gene.name, ";", sep="")
data.analysis$Gene.names <- paste(data.analysis$Gene.names, ";", sep="")

############################################################################################
###Rename columns according to experimental design and set up model.matrix appropriately
data.analysis <- setNames(data.analysis, c("ID","Gene.names", "MS.MS.count",
                                   "V_log2.i._TMT_1", 
                                   "V_log2.i._TMT_2",  
                                   "V_log2.i._TMT_3",  
                                   "V_log2.i._TMT_4", 
                                   "V_log2.i._TMT_5",  
                                   "V_log2.i._TMT_6",  
                                   "V_log2.i._TMT_7", 
                                   "V_log2.i._TMT_8",  
                                   "V_log2.i._TMT_9",  
                                   "v_log2.i._TMT_10",
                                   
                                   "Protein.IDs", "Gene.name"  ))

write.csv(data.analysis,file="20190207_K562_TC_filteredData.csv",row.names = F) ## looks good
