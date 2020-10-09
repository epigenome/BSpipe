#!/usr/bin/env Rscript
if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

# pars = c('org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'LNCap CCR,LNCap OCR,PrEC CCR,PrEC OCR', 'yyy', '/ddn/work/cjh/single/GSE94361/x1.ccr', '/ddn/work/cjh/single/GSE94361/x1.ocr', '/ddn/work/cjh/single/GSE94361/x2.ccr', 'ddn/work/cjh/single/GSE94361/x2.ocr')

PdfFlag = 1

if (!exists('InFile')) {
	if (length(pars) < 6) {
		message("USAGE: Rscript annotate org_db tx_db names out_prefix input_file [input_file ...]\n")
		quit(save='no', status=-1)
	}

#	OrgDB <- "org.Hs.eg.db"
	OrgDB = pars[1]
	TxDB = pars[2]
	Names = pars[3]
	OutFile = pars[4]
	InFiles = pars[5:length(pars)]
	isdir = file.exists(OutFile) && file.info(OutFile)[2] == T
	last = substr(OutFile, nchar(OutFile), nchar(OutFile))
	if (isdir) {
		if (last != '/') OutFile = paste(OutFile, '/', sep='')
	} else {
		if (last != '.') OutFile = paste(OutFile, '.', sep='')
	}
}

openGraphics <- function(path, width=800, height=700)
{
	flag = F
	if (exists('PdfFlag')) { if (PdfFlag) flag = T }
	if (!flag) { flag = any(grep("\\.pdf$", path)) }

	if (flag)
	{
		if (!any(grep("\\.pdf$", path))) path <- paste(path, ".pdf", sep="")
		pdf(path, w=width/100, h=height/100)
	}
	else
	{
		if (!any(grep("\\.png$", path))) path <- paste(path, ".png", sep="")
		png(path, w=width, h=height)
	}
}

closeGraphics <- function()
{
	graphics.off()
}

loadFile <- function(path, comment.char='', ...)
{
	print(paste("Loading", path, ".."))
	data <- readPeakFile(path, ...)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

saveFile <- function(dat, path, rname=F, cname=T)
{
	write.table(dat, quote=F, sep="\t", row.names=rname, col.names=cname, file=path)
}

draw <- function(fun, x, path=NA, width=800, height=700, close=T, par=NA, ...)
{
	if (!is.na(path)) openGraphics(path, width=width, height=height)
	if (!is.na(par)) par()
	fun(x, ...)
	if (!is.na(path) && close) closeGraphics()
}

ggdraw <- function(fun, x, path=NA, ...)
{
	p <- fun(x, ...)
	draw(print, p, path)
}


suppressMessages(library(ChIPseeker))
#TxDB <- 'TxDb.Hsapiens.UCSC.hg19.knownGene'
suppressMessages(library(TxDB, character.only=T))
txdb <- eval(parse(text=TxDB))
suppressMessages(library(ReactomePA))
suppressMessages(library(clusterProfiler))


cat("Start: ", date(), "\n")


peaks = lapply(InFiles, loadFile)
names(peaks) = gsub(' ', '_' , unlist(strsplit(Names, ',')))
#covplot(peak, weightCol="V5")
#for (i in 1:length(peaks)) {
#	draw(covplot, peaks[[i]], path=paste(OutFile, 'cov.', names(peaks)[i], sep=''), weightCol="V5")
#}

#peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
for (i in 1:length(peaks)) {
	draw(plotAnnoPie, peakAnnoList[[i]], path=paste(OutFile, 'pie.', names(peaks)[i], sep=''))
}
#plotAnnoBar(peakAnnoList)
ggdraw(plotAnnoBar, peakAnnoList, path=paste(OutFile, 'bar', sep=''))
#for (i in 1:length(peaks)) {
#	ggdraw(plotAnnoBar, peakAnnoList[[i]], path=paste(OutFile, 'bar.', names(peaks)[i], sep=''))
#}
#vennpie(peakAnno)
for (i in 1:length(peaks)) {
	draw(vennpie, peakAnnoList[[i]], path=paste(OutFile, 'vennpie.', names(peaks)[i], sep=''))
}
#upsetplot(peakAnno)
for (i in 1:length(peaks)) {
	ggdraw(upsetplot, peakAnnoList[[i]], path=paste(OutFile, 'upset.', names(peaks)[i], sep=''))
}
#upsetplot(peakAnno, vennpie=TRUE)
for (i in 1:length(peaks)) {
	ggdraw(upsetplot, peakAnnoList[[i]], path=paste(OutFile, 'upsetvennpie.', names(peaks)[i], sep=''))
}
#plotDistToTSS(peakAnnoList)
ggdraw(plotDistToTSS, peakAnnoList, path=paste(OutFile, 'tss', sep=''))

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
ggdraw(plotAvgProf, tagMatrixList, path=paste(OutFile, 'tag', sep=''), xlim=c(-3000, 3000))
#plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
ggdraw(plotAvgProf, tagMatrixList, path=paste(OutFile, 'tagsep', sep=''), xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
#tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
draw(tagHeatmap, tagMatrixList, path=paste(OutFile, 'tagheat', sep=''), xlim=c(-3000, 3000), color=NULL)

for (i in 1:length(peaks)) {
	pathway1 <- enrichPathway(as.data.frame(peakAnnoList[[i]])$geneId)
	head(pathway1, 2)
	saveFile(pathway1, path=paste(OutFile, 'pathway1.', names(peaks)[i], '.txt', sep=''))
	gene <- seq2gene(peaks[[i]], tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
	pathway2 <- enrichPathway(gene)
	head(pathway2, 2)
	saveFile(pathway2, path=paste(OutFile, 'pathway2.', names(peaks)[i], '.txt', sep=''))
	ggdraw(dotplot, pathway2, path=paste(OutFile, 'dot.', names(peaks)[i], sep=''))
}

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster = genes, fun = "enrichKEGG", pvalueCutoff = 0.05, pAdjustMethod = "BH")
ggdraw(dotplot, compKEGG, , path=paste(OutFile, 'kegg', sep=''), showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
draw(vennplot, genes, path=paste(OutFile, 'venngene', sep=''))

cat("End: ", date(), "\n")

