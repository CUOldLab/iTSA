#################################################################
## empirical Bayes Statistical Analysis - script by Kerri Ball and Stephen Coleman
#################################################################

### User-defined variables are set below
  # The output filename
prefix <- "date_description"
drawplots <- TRUE

### Define the V and D columnnames to be used for the experiment's multiplexing pattern
columnnames <-c("ID","Gene.names", "MS.MS.count",
                "V_log2.i._TMT_1", "D_log2.i._TMT_2",
                "V_log2.i._TMT_3", "D_log2.i._TMT_4",
                "V_log2.i._TMT_5", "D_log2.i._TMT_6",
                "V_log2.i._TMT_7", "D_log2.i._TMT_8",
                "V_log2.i._TMT_9", "D_log2.i._TMT_10","D_log2.i._TMT_11",
                "Protein.IDs" )

# set up model.matrix according to label strategy/experimental design
design <- model.matrix(~ 0+factor(c(0,1,0,1,0,1,0,1,0,1,1))) #0=Vehicle; 1=Drug
colnames(design) <- c("V","D") #0=Vehicle; 1=Drug


#################################################################
library(limma)
library(stringr)
this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)

### import data
proteinGroups <- read.table(file = "proteinGroups.txt",
                            header = TRUE,
                            sep = "\t",
                            quote = "\"'",
                            dec = ".",
                            numerals = c("warn.loss"),
                            row.names = NULL,
                            na.strings = c("NA","NaN","Infinite"))
### Filter data rows
data <- proteinGroups
colnames(data)
data <- subset(data, Only.identified.by.site != "+")
data <- subset(data, Reverse != "+")
data <- subset(data, Potential.contaminant != "+")
data <- subset(data, Unique.peptides > 0)
data <- subset(data, MS.MS.count > 2)
data <- subset(data, Reporter.intensity.count.3 > 2)

### Filter data columns & define TMT data columns
  # identifies column numbers that contain Label corrected intensity data (TMT lot-specific correction)
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))
data.sort <- data[c("Protein.IDs", "Gene.names", "MS.MS.count")]
data.sort <- cbind(data.sort, data[,data.TMT])
data <- data.sort

### identify column numbers that contain intensity data post column filtering
data.TMT <- grep("Reporter.intensity.corrected", colnames(data))
data[,data.TMT] <- log(data[,data.TMT],base=2) #log2 transform Reporter.intensity

### Create Unique ID column & remove isoforms subjected to losing end of
  # MaxQuant's parsimonious distribution of peptides!
  # use specific columns names to identify protein id columns
data.metanms <- grep("Protein.IDs", colnames(data))
colnames(data)[data.metanms] <- "ID" #change column names
data$Protein.IDs <- data$ID
data$last <- (regexpr(';', as.character(data$Protein.IDs))-1) # if no ";" then = -1
data$ID.i <- ifelse(data$last > 0, (substr(as.character(data$Protein.IDs),1,data$last)),
                    (as.character(data$Protein.IDs)))
data$last <- (regexpr('-', as.character(data$ID.i))-1) # if no ";" then = -1
data$ID <- ifelse(data$last > 0, (substr(as.character(data$ID.i),1,data$last)), (as.character(data$ID.i)))
data$ID.i <- NULL
data$last <- NULL
data <- data[order(-data$MS.MS.count),] #sort by MS.MS. count (largest to smallest)
data$dup <- duplicated(data[,1]) #first occurrence is False, second and subsequent say True
data <- subset(data, data$dup != TRUE)
data$dup <- NULL

### Look at data prior to normalization:
  # NOTE: if data doesn't look pretty good stop here - normalization should not be used to "fix" bad data
if(drawplots){
  boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE,
          frame = FALSE, notch = TRUE,
          col =
            c("gray80","gray75","gray70","gray65","gray60","gray55","gray50","gray45","gray40","gray35"))
  pdf("BoxPlot_BeforeQuantileNormalization.pdf")
  boxplot(data[,data.TMT], main = "Before Normalization", names = TRUE,
          frame = FALSE, notch = TRUE,
          col =
            c("gray80","gray75","gray70","gray65","gray60","gray55","gray50","gray45","gray40","gray35"))
  # Close the pdf file
  dev.off()
}

### Normalize data
data.norm <- data
data.TMT <- grep("Reporter.intensity.corrected", colnames(data.norm)) # identify column numbers that contain intensity data
colnames(data.norm[,data.TMT])
data.norm[,data.TMT] <- normalizeBetweenArrays(as.matrix(data.norm[,data.TMT]),method="quantile" )
if(drawplots){
  boxplot((data.norm[,data.TMT]), main = "After Quantile Normalization", names = TRUE,
          frame = FALSE, notch = TRUE,
          col =
            c("gray80","gray75","gray70","gray65","gray60","gray55","gray50","gray45","gray40","gray35"))
  pdf("BoxPlot_AfterQuantileNormalization.pdf")
  boxplot((data.norm[,data.TMT]), main = "After Quantile Normalization", names = TRUE,
          frame = FALSE, notch = TRUE,
          col =
            c("gray80","gray75","gray70","gray65","gray60","gray55","gray50","gray45","gray40","gray35"))
  # Close the pdf file
  dev.off()
}


### Rename columns according to experimental design and set up model.matrix appropriately
data.norm <- setNames(data.norm, columnnames)
#colnames(data.norm)

### Use limma for empirical bayes


fit <- lmFit(data.norm[,data.TMT], design)
fit <- eBayes(fit) ##Apply empirical Bayes smoothing to the standard errors.
contrast <- makeContrasts("D-V", levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
ptopf <- topTableF(fit2, adjust="BH",genelist = data[,"ID"],number=Inf)
ptopf$rank <- 1:length(ptopf$F)
ptopf$numfp <- ptopf$rank * ptopf$adj.P.Val
metanms <- c("ID")
#f.stat <- topTableF(fit2, adjust="BH",genelist = data[,metanms],number=Inf,p.value=1) #this table reports the moderated F-statistic
#View(f.stat)

t.stat <- topTable(fit2, adjust="BH",genelist = data[,metanms], coef=1,p.value=1, number=Inf) #this table reports the moderated t-statistic
nms <- grep("logFC", colnames(t.stat)) 
colnames(t.stat)[nms] <-  "Log2FC"	#change column names
t.stat$Log10P.Value <- -log(t.stat$P.Value,base=10)
t.stat$Log10adj.P.Val <- -log(t.stat$adj.P.Val,base=10)
t.statNMS <-t.stat[c("ID","Log2FC","Log10P.Value","Log10adj.P.Val","P.Value","adj.P.Val","AveExpr","t","B") ]   
t.stat <- t.statNMS

sum(t.stat$adj.P.Val <= 0.05)
sum(t.stat$adj.P.Val <= 0.01)
sum(t.stat$adj.P.Val <= 0.001)


### preview data
tab1 <- t.stat
if(drawplots){
  Xdata <- tab1$Log2FC
  Ydata <- tab1$Log10adj.P.Val
  plot(Xdata, Ydata, # plot the variables
       main = "",
       xlab="log2FC", # x-axis label
       ylab="-Log10(q.Val)") # y-axis label
}

### Change column names and merge EB results with normalized data
data.analysis <- merge(t.stat, data.norm, by = "ID" ,all = T, sort=F)
#colnames(data.analysis)
data.analysis <- data.analysis[order(-data.analysis$Log10adj.P.Val),]


### Create gene.name column
data.analysis$last <- (regexpr(';', as.character(data.analysis$Gene.names))-1) # if no ";" then = -1
data.analysis$Gene.name <- ifelse(data.analysis$last > 0,(substr(as.character(data.analysis$Gene.names),1,data.analysis$last)),(as.character(data.analysis$Gene.names)))
  # if no gene name is present then use ID
data.analysis$last <- nchar(data.analysis$Gene.name)
data.analysis$name <- data.analysis$Gene.name
data.analysis$Gene.name <- ifelse (data.analysis$last == 0, as.character(data.analysis$ID),data.analysis$name)
  # Clean-up matrix for export: reorganize columns
V_TMT <- grep("V_log2.i._TMT_", colnames(data.analysis))
D_TMT <- grep("D_log2.i._TMT_", colnames(data.analysis))
#colnames(data.analysis) 
data.sort <- data.analysis[c("ID", "Protein.IDs","Gene.names","Gene.name", 
                             "Log2FC","Log10P.Value","Log10adj.P.Val","P.Value","adj.P.Val",
                             "AveExpr","t","B","MS.MS.count")]
test2 <- cbind(data.sort, data.analysis[,V_TMT], data.analysis[,D_TMT])
colnames(test2)
data.analysis <- test2

##CSV files opened in excel will change Sept# and MAR# to dates;
## this addition tells excel that this is text and prevents change!
## (use find & replace in excel if desired)
preData <- test2
preData$Gene.name <- paste(preData$Gene.name, ";", sep="")
preData$Gene.names <- paste(preData$Gene.names, ";", sep="")
if(drawplots){ View(preData) }
write.csv(preData,file=paste(prefix,"_Simple_empiricalBayesAnalysis.csv",sep=""),
          row.names = F)



#Create FDR and FC columns
data.analysis$FDR0.001Sig <- ifelse(data.analysis$Log10adj.P.Val>=3, 1, 0)
data.analysis$FDR0.01Sig <- ifelse(data.analysis$Log10adj.P.Val>=2, 1, 0)
data.analysis$FDR0.05Sig <- ifelse(data.analysis$Log10adj.P.Val>=1.30103, 1, 0)
data.analysis$absFC <- abs(data.analysis$Log2FC)
if(drawplots){ plot(data.analysis$absFC, data.analysis$Log10adj.P.Val) }
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
data.analysis$D <- ifelse(data.analysis$Log2FC >= 0, (1*data.analysis$Sig_FDR.FC), (-1*data.analysis$Sig_FDR.FC))
data.analysis$Direction <- (data.analysis$D/data.analysis$Sig_FDR.FC)
data.analysis <- data.analysis[order(-data.analysis$Sig_FDR.FC),]
data.analysis$rank <- 1:length(data.analysis$Sig_FDR.FC)
#colnames(data.analysis)
data.data <- data.analysis[c("ID","Protein.IDs","Gene.names","Gene.name","rank","Sig_FDR.FC","Direction",
                             "Log2FC","Log10P.Value","Log10adj.P.Val","P.Value","adj.P.Val","AveExpr","t","B","MS.MS.count",
                             "V_log2.i._TMT_1","V_log2.i._TMT_3","V_log2.i._TMT_5","V_log2.i._TMT_7","V_log2.i._TMT_9",
                             "D_log2.i._TMT_2","D_log2.i._TMT_4","D_log2.i._TMT_6","D_log2.i._TMT_8","D_log2.i._TMT_10")]
data.analysis <- data.data
data.results <- data.data
if(drawplots){ plot(data.analysis$rank, data.analysis$Log10P.Value) }

### CSV files opened in excel will change Sept# and MAR# to dates;
  # this addition tells excel that this is text and prevents change!
  # (use find & replace in excel if desired)
data.data$Gene.name <- paste(data.data$Gene.name, ";", sep="")
data.data$Gene.names <- paste(data.data$Gene.names, ";", sep="")
if(drawplots){ View(data.data) }
write.csv(data.data,file=paste(prefix,"_RankedEmpiricalBayesAnalysis.csv",sep=""),row.names = F)

### Annotate data
  # import annotation file:
  # (create annotation file with https://www.uniprot.org/uniprot/ and edit in excel.
  # change "Entry" header to "ID" and I like yes/NA "annotation" columns)
annotation <- read.table(file = "UniprotID_ATPbinding_organism_isHuman_statis-isReviewed.txt",
                         header = TRUE,
                         sep = "\t",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))
test2 <- merge(data.analysis, annotation, by = "ID" ,all = T, sort=F)
test3 <- subset(test2, Sig_FDR.FC >= 0)
data.annotated <- test3
annotation <- read.table(file = "UniprotID_family-isKinase_organism_isHuman_statis-isReviewed.txt",
                         header = TRUE,
                         sep = "\t",
                         quote = "\"'",
                         dec = ".",
                         numerals = c("warn.loss"),
                         row.names = NULL,
                         na.strings = c("NA","NaN","Infinite"))
test2 <- merge(data.annotated, annotation, by = "ID" ,all = T, sort=F)
test3 <- subset(test2, Sig_FDR.FC >= 0)
data.annotated <- test3

##CSV files opened in excel will change Sept# and MAR# to dates;
## this addition tells excel that this is text and prevents change!
## (use find & replace in excel if desired)
data.annotated$Gene.name <- paste(data.annotated$Gene.name, ";", sep="")
data.annotated$Gene.names <- paste(data.annotated$Gene.names, ";", sep="")
data.annotated <- data.annotated[order(data.annotated$rank),]
if(drawplots){ View(data.annotated)}
write.csv(data.annotated,file=paste(prefix,"_AnnotatedEmpiricalBayesAnalysis.csv",sep=""), row.names = F)

summary(fit2$df.prior)
summary(fit2$df.residual)
summary(fit2$df.total)

### plotly for html volcano plot viewer & files
results <- data.results
sum(results$adj.P.Val <= 0.05)
sum(results$adj.P.Val <= 0.001)

library(ggplot2)
library(plotly)
library(tidyverse)
source('./fdrfunctions.R')
source('ggvolcano.R')

myggma <- results %>% ggmaplot(xvar = Log2FC, yvar = Log10P.Value, adjpvalvar = adj.P.Val  )
ss.gg <- myggma$maplot.ggplot
ss.plotly <- myggma$maplot.plotly

### this is the code that opens the viewer - give it a minute!
ss.plotly
## Save plot as html file if desired
##
###############################################################################################################
citation() 
R.Version()