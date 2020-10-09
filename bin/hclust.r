#!/bin/env Rscript

suppressMessages(library(getopt))
suppressMessages(library(dynamicTreeCut))

if (!exists('pars')) pars=commandArgs(TRUE)

grp.colors = c('green3', 'gold', 'magenta', 'cyan', 'pink')
spec = matrix(c(
		'input',		'i', 1, 'character', 'input file',
		'output',	'o', 1, 'character', 'output file',
		'column' ,	'c', 1, 'integer',   'starting column for data',
		'end' ,	   'd', 1, 'integer',   'end column for data, default: the same as the staring column',
		'group' ,	'r', 1, 'integer',   'column for groups',
		'bycolumn',	'Y', 0, 'logical',   'output for each column',
		'header' ,	'H', 0, 'logical',   'first line is header',
		'name',	   'N', 2, 'character', 'column names',
		'color',	   'C', 2, 'character', 'color',
		'outline',	'O', 0, 'logical',   'no outline',
		'las',	   'L', 0, 'logical',   'vertical lables',
		'horizontal',	'Z', 0, 'logical',   'horizontal plot',
		'boxwidth',	'W', 1, 'numeric',   'box width',
		'varwidth',	'V', 0, 'logical',   'proportional to the square-roots of the number of observations in the groups',
		'title',	   't', 2, 'character', 'title',
		'xlabel',	'x', 2, 'character', 'X label',
		'ylabel',	'y', 2, 'character', 'Y label',
		'width' ,	'w', 1, 'integer',   'image width',
		'height',	'h', 1, 'integer',   'image height',
		'gct',		'G', 0, 'logical',   'input file in in GCT format',
		'pdf',		'p', 0, 'logical',   'output file in PDF'
), ncol=5, byrow=T)

usage <- function(spec)
{
	message("USAGE: hclust option")
	for (i in 1:nrow(spec)) message(paste('    ', '-', spec[i,2], ' or ', '--', spec[i,1], ' ', spec[i,4], ' : ', spec[i, 5], sep=''))
	message("\n")
	quit(save='no', status=-1)
}

if (!exists('InFile')) {
	if (length(pars) > 0) {
		opt = getopt(opt=pars, spec=spec[,1:4])
	} else {
		opt = list()
	}

#	InFile = opt$input
#	OutFile = opt$output
#	GctFlag = !is.null(opt$gct)
#	DatCol = ifelse(GctFlag, 3, ifelse(is.null(opt$column), -1, opt$column))

#	if (!exists('InFile') || is.null(InFile) || is.null(OutFile) || DatCol == -1) {
#		usage(spec)
#	}

	PdfFlag = !is.null(opt$pdf)
	Width = ifelse(is.null(opt$width), 800, opt$width)
	Height = ifelse(is.null(opt$height), 600, opt$height)
	Outline = ifelse(is.null(opt$outline), T, F)
	Las = ifelse(is.null(opt$las), 1, 2)
	BoxWidth = if (!is.null(opt$boxwidth)) opt$boxwidth else 0.8
	VarWidth = if (!is.null(opt$varwidth)) T else F
}

if (!is.null(opt$name)) opt$name = strsplit(opt$name , ',')[[1]]
if (!is.null(opt$color)) {
	if (sub('rainbow|heat.colors|terrain.colors|topo.colors|cm.colors', '', opt$color) != opt$color) {
		opt$color = eval(parse(text=opt$color))
	} else {
		opt$color = strsplit(opt$color, ',')[[1]]
	}
}

openGraphics <- function(path, width, height)
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


#message("Start: ", date(), "\n")

mat <- read.table('stdin', sep='\t', skip=0, header=T, check.names=F, comment.char='', row.names=NULL)
#dim(mat)
#mat

#as.dist(mat)
hc <- hclust(as.dist(mat), 'complete')
#str(hc)

dc <- cutreeDynamic(hc, method="tree", minClusterSize=2, deepSplit=T)
write.table(t(dc), stdout(), row.names=F, col.names=F)

#message("End: ", date(), "\n")
