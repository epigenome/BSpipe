#!/usr/bin/env Rscript
if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile')) {
	if (length(pars) < 1) {
		message("USAGE: Rscript figure.R ann_file col_number [out_prefix]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	Idx = as.integer(pars[2])
	OutFile = if (length(pars) > 2) pars[3] else InFile
	isdir = file.exists(OutFile) && file.info(OutFile)[2] == T
	last = substr(OutFile, nchar(OutFile), nchar(OutFile))
	if (isdir) {
		if (last != '/') OutFile = paste(OutFile, '/', sep='')
	} else {
		if (last != '.') OutFile = paste(OutFile, '.', sep='')
	}
}


suppressMessages(library(geneplotter))
suppressMessages(library(reshape))


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

loadFile <- function(file, comment.char='', ...)
{
	print(paste("Loading", file, ".."))
	data <- read.table(file, sep="\t", fill=T, as.is=T, check.names=F, comment.char=comment.char, quote='', ...)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

draw <- function(fun, x, file=NA, width=800, height=700, close=T, par=NA, ...)
{
	if (!is.na(file)) openGraphics(file, width=width, height=height)
	if (!is.na(par)) par()
	fun(x, ...)
	if (!is.na(file) && close) closeGraphics()
}


cat("Start: ", date(), "\n")

PdfFlag = 0


dat = loadFile(pipe(paste("perl -ane 'map {print join(\"\\t\", @F[0..", Idx-1, "], $_, @F[", Idx+1, "..$#F]), \"\\n\"} split(/,/, $F[", Idx, "])' ", InFile, sep='')), comment.char='#', header=F, skip=0)
header = loadFile(InFile, header=F, nrow=1)
colnames(dat)[1:length(header)] = header
dat[Idx] = sub('\\?', '', dat[[Idx]])

cls = unique(unlist(strsplit(dat[[Idx]], ',')))
feat = unique(sub('\\(.*\\)', '', cls))
print(feat)

#dat = do.call(rbind, apply(dat, 1, function(x) mdply(strsplit(x[Idx], ',')[[1]], function(y) c(x[1:(Idx-1)], y))))
#print(dim)

n = round(sqrt(length(feat)))
colors = c('blue3', 'red', 'green3')
for (i in 4:(Idx-1)) {
	openGraphics(paste(OutFile, colnames(dat)[i], sep=''), 50+n*150, 50+n*150)
	par(mfrow=c(n, n))
	for (f in feat) {
		cn = c('upstream','body','downstream')
		sel = paste(f, '(', cn, ')', sep='')
		dat.sel = list(up=NA, body=NA, down=NA)
		if (any(sel[1] == dat[[Idx]])) dat.sel$up   = dat[dat[Idx] == sel[1], i]
		if (any(sel[2] == dat[[Idx]])) dat.sel$body = dat[dat[Idx] == sel[2], i]
		if (any(sel[3] == dat[[Idx]])) dat.sel$down = dat[dat[Idx] == sel[3], i]
		if (length(na.omit(dat.sel$up))+length(na.omit(dat.sel$body))+length(na.omit(dat.sel$down)) == 0) next
		boxplot(dat.sel, main=f, col=colors, ylim=c(0,1), las=1, outline=F)
	}
	closeGraphics()
}

for (f in feat) {
	openGraphics(paste(OutFile, f, sep=''), 50+(Idx-4)*30, 500)
	cn = c('upstream','body','downstream')
	par(mfrow=c(3, 1))
	for (i in 1:length(cn)) {
		try(boxplot(dat[dat[Idx] == paste(f, '(', cn[i], ')', sep=''), 4:(Idx-1)], main=cn[i], col=colors[i], ylim=c(0,1), las=2, outline=F))
	}
	closeGraphics()
}

cat("End: ", date(), "\n")

