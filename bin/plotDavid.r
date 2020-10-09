
if (!exists('pars')) pars <- commandArgs(trailingOnly=TRUE)

inFile <- pars[1]
if(file.exists(inFile)==FALSE) {
	write("Input file does not exist.", stderr())
	stop("Quitting.......")
} else {
	davidTable <- read.table(inFile, sep="\t", quote='', header=TRUE)
}

if(nrow(davidTable)==0) {
	message("No entries in the DAVID output")
	q()
}

outDir <- pars[2]
nItems <- as.integer(pars[3])
value <- pars[4]
annotation <- pars[5]
outFile <- paste(outDir, '/', sub('\\.[^\\.]+$//', '', basename(inFile)),".png", sep='')
#Options: pvalue, benjamini, bonferroni, fdr 

if(tolower(value)=="pvalue") {
	value="Pvalue"
} else if(tolower(value)=="benjamini") {
	value="Benjamini"
} else if(tolower(value)=="bonferroni") {
	value="Bonferroni"
} else if(tolower(value)=="fdr") {
	value="FDR"
}

#Split the term column to extract the annotation and add 2 new columns
davidTable$Term <- sub(".*[~:]", "", as.character(davidTable$Term))
head(davidTable[1:5])
#Create a new table with only the 2 required columns. The value is log transformed here

newDavidTable <- data.frame(Term=davidTable$Term, Value=-log10(davidTable[[value]]))
#Sort the table in the increasing order of log transformed value
newDavidSortTable <- newDavidTable[order(newDavidTable$Value, decreasing=TRUE),]

#Extract n number of rows to plot the graph.
#newDavidGraphTable <- head(newDavidSortTable, n=nItems)
if (nItems > nrow(newDavidSortTable)) nItems = nrow(newDavidSortTable)
newDavidGraphTable <- newDavidSortTable[1:nItems, 1:ncol(newDavidSortTable)]
#newDavidGraphTable2 <- newDavidGraphTable[!is.na(newDavidGraphTable$Term),]
Value2 <- newDavidGraphTable$Value
#print(newDavidGraphTable$Term)
#print(max(Value2))
#Draw graph
newDavidGraphTable$Term <- sub("(.{50}).*", "\\1...", perl=T, as.character(newDavidGraphTable$Term))
len = max(nchar(as.character(newDavidGraphTable$Term)))
library(plotrix)
png(outFile, height=120+nItems*20, width=300+len*6)
par(mar=c(5,5+len/2.3,3,1))
xlab <- paste("-log10(",value,")",sep='')
title <- paste("Top",nItems,annotation, "based on", value, sep=' ')
barplot(rev(newDavidGraphTable$Value), las=1, horiz=TRUE, names.arg=rev(newDavidGraphTable$Term), col="slateblue",space=0.5, xlim=c(0, max(Value2)*1.2), xlab=xlab, cex.lab=1.2, cex.axis=1.2, cex.names=1.2)
#axis(3, las=1, at=seq(0, max(newDavidGraphTable$Value)+1))
#mtext(xlab, side=3, line=2, las=0)
mtext(bquote(bold(.(title))), side=3, at=getFigCtr()[1], cex=1.4)

dev.off()
