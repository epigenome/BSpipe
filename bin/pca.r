#!/bin/env Rscript

suppressMessages(library(getopt))

if (!exists('pars')) pars=commandArgs(TRUE)
cat('Arguments:', pars, '\n')

MinCpgNum <- 5
col.pch = NULL
col.col = NULL
grp.colors = c('green3', 'blue3', 'magenta', 'brown', 'black', 'purple', 'pink', 'darkcyan', 'red', 'cyan', 'gold', 'green', 'darkred')
#grp.colors <- c('green3', 'gold', 'purple', 'orange', 'magenta', 'cyan', 'blue3', 'red', 'pink', 'black')

spec = matrix(c(
		'format',   'b', 1, 'character', 'input format',
		'input',		'i', 1, 'character', 'input file',
		'output',	'o', 1, 'character', 'output file',
		'group',		'R', 1, 'character', 'file for group',
		'only'    ,	'O', 0, 'logical',   'use samples in the group file',
		'pch' ,		'P', 1, 'integer',   'point type [1]',
		'cpg' ,     'c', 1, 'integer',   'minium number of sites in a window',
		'width' ,	'w', 1, 'integer',   'image width',
		'height',	'h', 1, 'integer',   'image height',
		'legend',	'l', 1, 'character', 'legend position: e.g., topright, centerleft',
		'scale',		's', 0, 'logical',   'scale data',
		'pdf',		'Y', 0, 'logical',   'output file in PDF'
), ncol=5, byrow=T)

usage <- function(spec)
{
	message("USAGE: pca.r option")
	for (i in 1:nrow(spec)) message(paste('    ', '-', spec[i,2], ' or ', '--', spec[i,1], ' ', spec[i,4], ' : ', spec[i, 5], sep=''))
	message("\n")
	quit(save='no', status=-1)
}

if (!exists('InFile')) {
	if (length(pars) > 0) {
		opt = getopt(opt=pars, spec=spec[,1:4])
	} else {
		usage(spec)
		opt = list()
	}

	Format = ifelse(is.null(opt$format) || opt$format == 'auto', 0, ifelse(opt$format == 'sitemethy', 1, ifelse(opt$format == 'sitenread', 2, ifelse(opt$format == 'winmethy', 3, 4))))
	if (!is.null(opt$cpg   )) MinCpgNum = as.integer(opt$cpg)
	InFile = opt$input
	OutFile = sub('\\.gz$', '', opt$output)
	OutFile = paste(sub('\\.pca$', '', OutFile), 'pca', sep='.')

	if (!exists('InFile') || is.null(InFile) || is.null(OutFile)) {
		usage(spec)
	}

	PdfFlag = !is.null(opt$pdf) && opt$pdf
	Pch = ifelse(is.null(opt$pch), 1, opt$pch)
	Width = ifelse(is.null(opt$width), 800, opt$width)
	Height = ifelse(is.null(opt$height), 600, opt$height)
}

suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
suppressMessages(library(corrgram))


openGraphics <- function(path, width, height)
{
	flag = F
	if (PdfFlag) flag = T
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

openFile <- function(x) if (x == '-') file('stdin') else x

loadFile <- function(file, header=T, comment.char='', ...)
{
	print(paste("Loading", file, ".."))
	var <- read.table(openFile(file), sep="\t", header=header, check.names=F, comment.char=comment.char, ...)
	print(head(var))
	print(tail(var))
	cat("dim=", dim(var), "\n")
	return(var)
}

draw <- function(fun, x, file=NA, width=Width, height=Height, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	fun(x, ...)

	if (!is.na(file)) closeGraphics()
}

plotPca <- function(pca, grp.pch=NULL, grp.col=NULL, file=NA, width=Width, height=Height, label=F, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	if (is.null(grp.col)) grp.col='black'
	if (is.null(grp.pch)) grp.pch=Pch

	plot(pca$x[,1], pca$x[,2], xlab=paste("PCA 1 (", round(pca$sdev[1]^2, 2), ")", sep=''), ylab=paste("PCA 2 (", round(pca$sdev[2]^2, 2), ")", sep=''), type="p", pch=grp.pch, xlim=c(min(pca$x[,1:2]), max(pca$x[,1:2])), ylim=c(min(pca$x[,1:2]), max(pca$x[,1:2])), col=grp.col, ...)
#	arrows(0,0,pca$rotation[1:3,1]*1,pca$rotation[1:3,2]*1, length=0.1, angle=20, col="red")
#	text(pca$rotation[,1]*10*1.2,pca$rotation[,2]*10*1.2, rownames(pca$rotation), col="red", cex=0.7)
#	text(pca$x[,1],pca$x[,2], rownames(pca$x), col="blue", cex=0.7)
	if (length(grp.col) > 1)
	{
		names = levels(grp$group)
		cols = grp.colors[1:length(names)]
		if (any(grp.col == 'gray'))
		{
			names = c('others', names)
			cols  = c('gray', cols)
		}
		legend(ifelse(is.null(opt$legend), 'topright', opt$legend), names, pch=min(grp.pch):max(grp.pch), col=cols)
	}
	
	if (label) text(pca$x[,1], pca$x[,2], rownames(pca$x), pos=4, pch=grp.pch, col=grp.col)

	if (!is.na(file)) closeGraphics()
}


cat("Start: ", date(), "\n")

DatCol = 4

dat <- loadFile(InFile)
dat.row = nrow(dat)
colnames(dat)[c(1,2,3)] <- c('Ref', 'Start', 'End')
if (colnames(dat)[DatCol] == paste('Key', DatCol, sep='')) {
	colnames(dat)[DatCol] <- 'ID'
	DatCol = DatCol+1
}

dat = na.omit(dat)
cat('Removed features with NA: ', dat.row-nrow(dat), '\n')
col.tab = table(colnames(dat)[DatCol:length(dat)])

if (Format == 1 && any(col.tab != 1)) stop('File format is not the output of bspipe join')
if (Format == 2 && any(col.tab != 1)) stop('File format is not the output of bspipe join')
if (Format >  2 && any(col.tab != 3)) stop('File format is not the output of bspipe win')

if (Format == 0) {
	if (all(col.tab == 1)) {
		Format = 1
	} else if (all(col.tab == 2)) {
		Format = 2
	} else if (all(col.tab == 3)) {
		Format = -1
	} else {
		stop('File format is neither site, methy, nor nread')
	}
}

if (Format > 2 || Format == -1) {
	icpg <- seq(DatCol+0, length(dat), 3)
	imet <- seq(DatCol+1, length(dat), 3)
	itot <- seq(DatCol+2, length(dat), 3)
}

if (Format == -1) Format = ifelse(max(dat[imet[1]], na.rm=T) > 1, 4, 3)
if (Format == 2) {
	dat = dat[apply(dat[seq(DatCol+1, length(dat), 2)], 1, function(x) all(x!=0)), ]
	cat('Removed features with 0: ', dat.row-nrow(dat), '\n')
	print(dim(dat))
	dat = cbind(dat[1:(DatCol-1)], round(data.frame(t(apply(dat[DatCol:length(dat)], 1, function(x) x[seq(1,length(x),2)]/x[seq(2,length(x),2)]))), 3))
} else if (Format > 2) {
	dat.min = apply(dat[icpg], 1, function(x) ifelse(any(is.na(x)), 0, min(x)))
	cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(dat.min < MinCpgNum), '\n')
	dat     = dat    [dat.min >= MinCpgNum,]
	dat.min = dat.min[dat.min >= MinCpgNum ]
	if (colnames(dat)[DatCol-1] == 'ID') dat$ID = factor(dat$ID)
	print(dim(dat))

	# window.methy or window.nread format
	if (Format == 4) {
		dat = cbind(dat[1:(DatCol-1)], round(dat[imet]/dat[itot], 3))
	} else {
		dat = cbind(dat[1:(DatCol-1)], dat[imet])
	}
}

print(head(dat))
print(tail(dat))
print(dim(dat))

if (!is.null(opt$group)) {
	grp <- loadFile(opt$group, header=F, comment.char='#')
	if (length(grp) == 4) {
		colnames(grp) = c('file', 'flag', 'sample', 'group')
		grp <- unique(grp[3:4])
	} else if (length(grp) == 2) {
		colnames(grp) = c('sample', 'group')
		grp <- unique(grp[1:2])
	} else {
		message("Unknown format for a group file")
		q()
	}

	if (!is.null(opt$only)) dat = cbind(dat[1:(DatCol-1)], dat[colnames(dat) %in% as.character(grp$sample)])

	if (length(levels(grp$group)) > length(grp.colors)) grp.colors = rainbow(length(levels(grp$group)))
	col.col = apply(as.matrix(colnames(dat)[DatCol:length(dat)]), 1, function (x) ifelse(sum(x==grp$sample) > 0, grp.colors[grp$group[x==grp$sample]], 'gray'))
	if (is.null(opt$pch)) col.pch = apply(as.matrix(colnames(dat)[DatCol:length(dat)]), 1, function (x) ifelse(sum(x==grp$sample) > 0, grp$group[x==grp$sample], 0))
}

dat.ind = DatCol:length(dat)
pca <- prcomp(t(dat[dat.ind]), retx=T, center=T, scale=!is.null(opt$scale))
plotPca(pca, col.pch, col.col, file=paste(OutFile, 'scatter', sep='.'))
plotPca(pca, col.pch, col.col, label=T, file=paste(OutFile, 'scatter', 'label', sep='.'))
draw(plot, pca, main=NULL, xlab="principle component", file=paste(OutFile, 'scree10', sep='.'))
draw(plot, log(pca$sdev^2), xlab="principle component", ylab="log(variance)", type="b", pch=16, file=paste(OutFile, 'scree', sep='.'))
write.table(pca$x, file=paste(OutFile, 'pos', 'csv', sep='.'), quote = F, sep = "\t", row.names = T, col.names=T)


cat("End: ", date(), "\n")
