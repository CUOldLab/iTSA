fdrcut <- function(adjpval, fdr.range.cut =  c(0,0.02,0.05,1)) {
  
  adjpval.ranges <- cut(adjpval, fdr.range.cut, include.lowest = T)  
  fdr.range.cut <- sprintf("%.2f", fdr.range.cut)
  if(length(fdr.range.cut) == 4) {
    adjpval.ranges <- factor(adjpval.ranges,levels = rev(levels(adjpval.ranges)),ordered = T,
                             labels = rev( c( paste(fdr.range.cut[1],'< q <=',fdr.range.cut[2],sep = ''),
                                              paste(fdr.range.cut[2],'< q <=',fdr.range.cut[3],sep = ''),
                                              paste(fdr.range.cut[3],'< q <=',fdr.range.cut[4],sep = ''))))
  } else if(length(fdr.range.cut) == 5) {
    adjpval.ranges <- factor(adjpval.ranges,levels = rev(levels(adjpval.ranges)),ordered = T,
                             labels = rev( c( paste(fdr.range.cut[1],'< q <=',fdr.range.cut[2],sep = ''),
                                              paste(fdr.range.cut[2],'< q <=',fdr.range.cut[3],sep = ''),
                                              paste(fdr.range.cut[3],'< q <=',fdr.range.cut[4],sep = ''),
                                              paste(fdr.range.cut[4],'< q <=',fdr.range.cut[5],sep = ''))))
    
  } else {
    stop("fdr.range.cut vector not of length 4")
  }
  return(adjpval.ranges)
}

get.toptable <- function(fitobj, coef, glist, suffix) {
  ttemp <- topTable(fitobj,coef = coef,genelist = glist,number = Inf,p.value = 1)
  datinds <-  which(names(ttemp) != "ID")
  names(ttemp)[datinds] <- paste(names(ttemp)[datinds], suffix, sep="")
  return(ttemp) 
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor,...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex,...) 
  text(.8, .8, Signif, cex=cex, col=2,...) 
}
#pairs(iris[1:4], lower.panel=panel.smooth, upper.panel=panel.cor)

