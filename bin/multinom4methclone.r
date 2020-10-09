FdrCutoff <- 0.05
PvalCutoff <- 0.01
MaxRowForHeatmap <- 10000
PdfFlag = 0
DatCol = 15


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre')) {
	if (length(pars) < 2) {
		message("USAGE: Rscript multinom.r in_file out_prefix [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	if (length(pars) > 2) opt = getopt(matrix(ncol=4, byrow=T, c(
		'fdr',      'f', 1, 'numeric',
		'pval',     'v', 1, 'numeric',
		'dist',     'D', 1, 'character',
		'linkage',  'L', 1, 'character',
		'thread',   'm', 1, 'integer',
		'width' ,   'W', 1, 'integer',
		'height',   'H', 1, 'integer',
		'pdf'     , 'Y', 0, 'logical'
		)), pars[3:length(pars)])
	else opt = list()

	if (!is.null(opt$fdr     )) FdrCutoff    = as.numeric(opt$fdr)
	if (!is.null(opt$pval    )) PvalCutoff   = as.numeric(opt$pval)
	DistMethod = ifelse(is.null(opt$dist), 'euclidean', opt$dist)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}


suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
suppressMessages(library(nnet))
if (!is.null(opt$thread) && opt$thread > 1) {
	message('Running in parallel mode with ', opt$thread, ' threads')
	suppressMessages(library(parallel))
}


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

loadFile <- function(file, header=T, comment.char='', ...)
{
	print(paste("Loading", file, ".."))
	data <- read.table(file, sep="\t", header=header, check.names=F, comment.char=comment.char, as.is=T, ...)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

drawHeatmap <- function(x, col=bluered(100), file=NA, ...)
{
	Width = ifelse(!is.null(opt$width), opt$width, 100 + ncol(x)*ifelse(ncol(x) < 100, 8, 7))
	Height = ifelse(!is.null(opt$height), opt$height, 300 + sqrt(nrow(x))*10)
	if (!is.na(file)) openGraphics(file, Width, Height)

	hclust2 <- function(x, method=Linkage, ...) hclust(x, method=method, ...)
	dist2 <- function(x, ...) if (DistMethod == 'pearson' | DistMethod == 'spearman' | DistMethod == 'kendall') as.dist(1-cor(t(x), method=DistMethod), ...) else dist(x, method=DistMethod, ...)

	if      (Height <  400) hei=4.0
	else if (Height <  500) hei=3.0
	else if (Height <  600) hei=2.0
	else if (Height <  900) hei=1.5
	else if (Height < 1200) hei=1.0
	else                    hei=0.8

	if      (Width <  500) wid=2.5
	else if (Width <  700) wid=2.0
	else if (Width <  900) wid=1.5
	else if (Width < 1200) wid=1.0
	else if (Width < 1500) wid=0.8
	else                   wid=0.7

	hv <- heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, labCol=NA, cexCol=1.3, cexRow=1, margin=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ColSideColors=c(rep('green3', ncol(x)/2), rep('gold', ncol(x)/2)), ...)

	if (!is.na(file)) closeGraphics()
	return(hv)
}

draw <- function(fun, x, file=NA, width=600, height=500, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	fun(x, ...)

	if (!is.na(file)) closeGraphics()
}

drawVolcano <- function(mq, title, cutoff, prefix, path=NA)
{
	mq.sig  <- sum(mq$q < cutoff)
	mq <- mq[mq$q != 0 & is.finite(mq$q),]
	if (!is.na(path)) openGraphics(path, 600, 500)

	min.q <- min(mq$q)
	max.q <- min(-log10(min.q), 10)

	smoothScatter(mq$m, -log10(mq$q), nbin=512, nrpoints=0, xlab="Difference", ylab=paste("-log10(", prefix, "-value)", sep=''), main=title, ylim=c(0,max.q), cex=ifelse(PdfFlag, 0.5, 1.5), cex.lab=1.2, cex.axis=1.2)
	abline(h=-log10(cutoff), col="darkgrey", lt=2)

	if (min.q < cutoff)
	{
		sig <- mq[mq$q < cutoff,]
		points(sig$m, -log10(sig$q), pch=20, cex=0.3, col='red')
		text(max(mq$m), max.q, adj=c(1,1), mq.sig, col='blue', cex=1.2)
	}

	if (!is.na(path)) closeGraphics()
}

mlogit <- function(var)
{
	x <- var[ 1:16]
	y <- var[17:32]
	z <- data.frame(rbind(cbind(1, seq(1, length(x)), x), cbind(2, seq(1, length(y)), y)))
	colnames(z) = c('group', 'pattern', 'value')

	res1 <- multinom(pattern ~ factor(group), weights=value, data=z, trace=F)
	res2 <- multinom(pattern ~ 1            , weights=value, data=z, trace=F)
	delta <- deviance(res2)- deviance(res1)
	pvalue <- dchisq(delta,15,0)

	return(c(delta, pvalue))
}

analysis <- function(info, data, suffix)
{
	if (nrow(data) <= 1) {
		return(invisible())
	} else if (nrow(data) > MaxRowForHeatmap) {
		data = head(data[order(info$pvalue),], MaxRowForHeatmap)
		info = head(info[order(info$pvalue),], MaxRowForHeatmap)
		suffix = paste(suffix, 'max', MaxRowForHeatmap, sep='')
	}

	hv <- drawHeatmap(data, Rowv=T, Colv=F, dendrogram='row', file=paste(OutPre, suffix, 'heatmap', sep='.'))
	write.table(cbind(info[rev(hv$rowInd), ], data[rev(hv$rowInd), ]), file=paste(OutPre, suffix, 'heatmap.csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
}

output <-function(pct, stat, ext)
{
	pext = paste(ext, PvalCutoff, 'pval', sep='.')
	qext = paste(ext, FdrCutoff, 'qval', sep='.')
	colnames(stat) = c('delta', 'pvalue')
	stat$pvalue[is.na(stat$pvalue)] = 1
	stat$delta[is.na(stat$delta)] = 0
	stat$padjust = p.adjust(stat$pvalue, 'fdr')
	out <- cbind(pct[1:(DatCol-1)], stat)
	write.table(out, file=paste(OutPre, ext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)

	max.delta = max(abs(stat$delta))
	draw(hist, stat$delta  , xlab='Deviance difference', ylab='Number of features', file=paste(OutPre, ext, 'md.hist', sep='.'))
	draw(hist, stat$pvalue , breaks=seq( 0,1,0.01), xlab='P-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'pval.hist', sep='.'))
	draw(hist, stat$padjust, breaks=seq( 0,1,0.01), xlab='Q-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'qval.hist', sep='.'))
	try(drawVolcano(data.frame(m=out$delta, q=out$pvalue ), '', PvalCutoff, 'p', path=paste(OutPre, pext, 'vol', sep='.')))
	try(drawVolcano(data.frame(m=out$delta, q=out$padjust), '',  FdrCutoff, 'q', path=paste(OutPre, qext, 'vol', sep='.')))

	log = out$pvalue<PvalCutoff
	if (sum(log) > 0) write.table(cbind(pct[log,1:(DatCol-1)], stat[log,]), file=paste(OutPre, pext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	try(analysis(cbind(pct[log,1:(DatCol-1)], stat[log,]), pct[log,DatCol:length(pct)], pext))

	log = out$padjust<FdrCutoff
	if (sum(log) > 0) write.table(cbind(pct[log,1:(DatCol-1)], stat[log,]), file=paste(OutPre, qext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	try(analysis(cbind(pct[log,1:(DatCol-1)], stat[log,]), pct[log,DatCol:length(pct)], qext))
}

stat.test <- function(x, fun)
{
	if (is.null(opt$thread) || opt$thread == 1) {
		return(data.frame(t(apply(x, 1, fun))))
	} else {
		return(data.frame(do.call(rbind, mclapply(split(x, 1:nrow(x)), function(y) fun(unlist(y)), mc.cores=opt$thread))))
	}
}


cat("Start: ", date(), "\n")

dat <- loadFile(InFile)
if (nrow(dat) == 0) q()
dat <- dat[1:(DatCol+31)]
#colnames(dat)[1:5] <- c('chr', 'start', 'end', 'strand', 'count')
#if (colnames(dat)[DatCol] == paste('Key', DatCol, sep='')) {
#	colnames(dat)[DatCol] <- 'ID'
#	DatCol = DatCol+1
#}

dat[(DatCol+ 0):(DatCol+15)] = round(dat[(DatCol+ 0):(DatCol+15)] * dat[[ 9]], 0)
dat[(DatCol+16):(DatCol+31)] = round(dat[(DatCol+16):(DatCol+31)] * dat[[10]], 0)

com <- stat.test(dat[DatCol:(DatCol+31)], mlogit)
output(dat, com, 'mlogit')

cat("End: ", date(), "\n")
