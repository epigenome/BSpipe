#!/bin/env Rscript
suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile')) {
	if (length(pars) < 2) {
		message("USAGE: Rscript metcov.r in_file out_prefix [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	if (length(pars) > 2) opt = getopt(matrix(ncol=4, byrow=T, c(
		'title' ,  't', 1, 'character',
#		'key'   ,  'k', 1, 'character',
		'width' ,  'w', 1, 'integer',
		'height',  'h', 1, 'integer',
		'pdf'   ,  'p', 0, 'logical'
		)), pars[3:length(pars)])
	else opt = list()

#	Key = ifelse(!is.null(opt$key), opt$key, 'top')
	Width = ifelse(!is.null(opt$width), opt$width, 600)
	Height = ifelse(!is.null(opt$height), opt$height, 450)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
	Title = ifelse(!is.null(opt$title), opt$title, basename(OutPre))
}


suppressMessages(library(geneplotter))
suppressMessages(library(gplots))


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

loadFile <- function(file)
{
	print(paste("Loading", file, ".."))
	data <- read.table(file, sep="\t", header=F, comment.char='', check.names=F)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

draw <- function(fun, x, file=NA, width=Width, height=Height, par=NULL, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)
	if (!is.null(par)) par(par)

	fun(x, ...)

	if (!is.na(file)) closeGraphics()
}


cat("Start: ", date(), "\n")

PdfFlag = F

dat <- loadFile(InFile)
dat$met = cut(dat[[5]], seq(0, 1, 0.1), include.lowest=T)
dat$cov = cut(dat[[4]], c(0, 1, 5, 10, 20, Inf), include.lowest=T)
met.name = apply(as.matrix(seq(0, 0.9, 0.1)), 1, function(x) paste(x, ifelse(x==0.9, '1', ''), sep='~'))

draw(barplot, table(dat$cov, dat$met), main=Title, xlab='Methylation', ylab='Number of sites', names.arg=met.name, col=c('white', 'yellow', 'green3', 'blue', 'red'), file=paste(OutPre, '.hist.cov', sep=''))

openGraphics(paste(OutPre, '.hist.met', sep=''), Width, Height)
par(mar=c(5, 4, 4, 4) + 0.1)
cov.tab = table(dat$met, cut(dat[[4]], c(seq(0, 100), Inf)))
y = barplot(cov.tab, main=Title, xlab='Read coverage', ylab='Number of sites', space=0, col=rainbow(11), legend.text=levels(dat$met), xaxt='n')
axis(side=1, at=y[1+seq(0,100, 10)], labels=seq(0,100,10))
par(new = T, mar=c(5, 4, 4, 4) + 0.1)
z=rev(cumsum(rev(cov.tab)))
plot(z/z[1], type="l", col="black", axes=F, xlab=NA, ylab=NA)
axis(side=4)
mtext(side=4, line=3, "Cumulative fraction of sites")
closeGraphics()

draw(boxplot, dat[[4]] ~ dat$met, col='green3', names=met.name, main=Title, xlab='Methylation', ylab='Read coverage',                 file=paste(OutPre, '.box.all', sep=''))
draw(boxplot, dat[[4]] ~ dat$met, col='green3', names=met.name, main=Title, xlab='Methylation', ylab='Read coverage', ylim=c(0, 100), file=paste(OutPre, '.box.100', sep=''))

draw(smoothScatter, dat[5:4], colram=topo.colors, main=Title, xlab='Methylation', ylab='Read coverage',                file=paste(OutPre, '.scatter.all', sep=''))
draw(smoothScatter, dat[5:4], colram=topo.colors, main=Title, xlab='Methylation', ylab='Read coverage', ylim=c(0,100), file=paste(OutPre, '.scatter.100', sep=''))

cat("End: ", date(), "\n")
