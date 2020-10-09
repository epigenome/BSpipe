#!/usr/bin/env Rscript
if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile')) {
	if (length(pars) < 1) {
		message("USAGE: Rscript figure.R prefix\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
}


suppressMessages(library(geneplotter))
suppressMessages(library(reshape))


openGraphics <- function(path, width=600, height=500)
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
	data <- read.table(file, sep="\t", fill=T, as.is=T, comment.char='', quote='', na='.', blank.lines.skip=F, row.names=NULL, header=F)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

draw <- function(fun, x, file=NA, width=600, height=500, close=T, par=NA, ...)
{
	if (!is.na(file)) openGraphics(file, width=width, height=height)
	if (!is.na(par)) par()
	fun(x, ...)
	if (!is.na(file) && close) closeGraphics()
}

chart <- function(b, e, j, col='blue3', pie.draw=T, sort=F, encode=T)
{
	name2 = toupper(name)
	name2 = sub('\\.(.+)', ' (\\L\\1\\E)', perl=T, name2)
	k = which(dat[j:e, 1] == '')
	k = ifelse(length(k) == 0, nrow(dat), k[1] + j - 2)
	x = as.integer(dat[j:k,2])
	l = dat[j:k, 1]
	print(dat[j:k, ])
	if (sum(x) == 0) return(invisible())

	if (encode) l = apply(as.matrix(l), 1,
			function (x) if (nchar(x) < 20) x else gsub("([Cei53pud])(DS|xon|ntron|'UTR|romoter|pstream|ownstream)", '\\U\\1\\E', perl=T, x))
	if (sort) {
		l = l[order(x, decreasing=T)]
		x = x[order(x, decreasing=T)]
	}

	if (pie.draw) draw(pie, x, labels=l[1:min(10, length(l))], clockwise=T, col=rainbow(length(l)), main=name2, file=paste(OutFile, name, 'pie', sep='.'), width=400, height=400, cex=1.2)

	x = x[min(10, length(x)):1]
	l = l[min(10, length(l)):1]
	rmargin <- 1.5 + max(nchar(l))/1.5
	draw(barplot, x/sum(x)*100, horiz=T, names.arg=l, las=2, col=col, xlab='(%)', main=name2, file=paste(OutFile, name, 'bar', sep='.'), height=150+length(l)*20, width=400, cex.axis=1.2, cex.names=1.2, par=par(mar=c(4.1, rmargin, 4.1, 2.1)))
}

cat("Start: ", date(), "\n")

PdfFlag = 0
total = 0

OutFile = sub('\\.table', '', InFile)
print(OutFile)
dat = loadFile(InFile)

beg = grep('#####', dat$V1)
end = c(beg[2:length(beg)]-1, nrow(dat))

for (i in 1:length(beg)) {
	name = dat[beg[i],2]
	combination = length(grep('combination', name))>0
	print(paste(name, ':', beg[i], '-', end[i], sep=' '))

	ind = grep('pie chart', dat[beg[i]:end[i], 1]) + beg[i] - 1
	if (length(ind) > 0) chart(beg[i], end[i], ind+1, sort=combination, encode=combination)

	ind = grep('feature plot', dat[beg[i]:end[i], 1]) + beg[i] - 1
	if (length(ind) > 0) chart(beg[i], end[i], ind+1, pie.draw=F, col='green3')
}

cat("End: ", date(), "\n")
