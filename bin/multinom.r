MaxRowForHeatmap <- 10000
DatCol = 6

# input file format
#Ref  Start End   Strand   Ncpgs Sam1	Sam1	Sam2	Sam2
#chr19 197220   197285   +  6  28	101110=1:110111=2:001110=3:111111=19	25	110111=10:111110=7:010111=1:111011=2:011111=2:111111=6

# running
# 1. all pairwise comparisons in input file: 
# 2. a given pair comparison in input file: --sample1 -- sample2
# 3. given pair comparisons in pair file: --pair
# 4. all pairwise comparisions in given groups: --group1 --group2 --group
# 5. all pairwise comparisons between group pairs in pair file: --pair --group

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
		'depth',    'd', 1, 'integer',
		'cpg',      'c', 1, 'integer',
		'fdr',      'f', 1, 'numeric',
		'pval',     'v', 1, 'numeric',
		'pair',     'p', 1, 'character',
		'group',    'q', 1, 'character',
		'sample1' , '1', 0, 'logical',
		'sample2' , '2', 0, 'logical',
		'dist',     'D', 1, 'character',
		'linkage',  'L', 1, 'character',
		'thread',   'm', 1, 'integer',
		'width' ,   'W', 1, 'integer',
		'height',   'H', 1, 'integer',
		'pdf'     , 'Y', 0, 'logical'
		)), pars[3:length(pars)])
	else opt = list()

	PvalCutoff <- ifelse(!is.null(opt$pval), opt$pval, 0.01)
	FdrCutoff <- ifelse(!is.null(opt$fdr), opt$fdr, 0.05)
	MinDepth <- ifelse(!is.null(opt$depth), opt$depth, 20)
	MinCpgNum <- ifelse(!is.null(opt$cpg), opt$cpg, 5)
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

#101110=1:110111=2:001110=3:111111=19
mlogit <- function(var)
{
	col.names = c('pattern', 'value')
	x <- data.frame(group=1, matrix(unlist(lapply(strsplit(var[1], ':'), function(x) strsplit(x, '='))), ncol=2, byrow=T, dimnames=list(NULL, col.names)), stringsAsFactors=F)
	y <- data.frame(group=2, matrix(unlist(lapply(strsplit(var[2], ':'), function(x) strsplit(x, '='))), ncol=2, byrow=T, dimnames=list(NULL, col.names)), stringsAsFactors=F)
	y.x <- setdiff(y[2], x[2])
	x.y <- setdiff(x[2], y[2])
	x2 <- if (length(y.x) == 0) x else rbind(x, data.frame(group=1, pattern=y.x, value="0"))
	y2 <- if (length(x.y) == 0) y else rbind(y, data.frame(group=2, pattern=x.y, value="0"))
	if (nrow(x2) == 1 || nrow(y2) == 1) return(c(NA, NA))

	z <- rbind(x2, y2)
	z$value = as.numeric(z$value)

	res1 <- try(multinom(factor(pattern) ~ factor(group), weights=value, data=z, trace=F))
	res2 <- try(multinom(factor(pattern) ~ 1            , weights=value, data=z, trace=F))
	if (class(res1) == "try-error" || class(res1) == "try-error") return(c(NA, NA))
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

	var = cbind(
				t(apply(data[1], 1, function(x) as.numeric(strsplit(x, ':')[[1]]))),
				t(apply(data[2], 1, function(x) as.numeric(strsplit(x, ':')[[1]]))) )
	hv <- drawHeatmap(var, Rowv=T, Colv=F, dendrogram='row', file=paste(OutPre, suffix, 'heatmap', sep='.'))
	write.table(cbind(info[rev(hv$rowInd), ], data[rev(hv$rowInd), ]), file=paste(OutPre, suffix, 'heatmap.csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
}

output <-function(info, pct, stat, title, ext)
{
	pext = paste(ext, PvalCutoff, 'pval', sep='.')
	qext = paste(ext, FdrCutoff, 'qval', sep='.')
	colnames(stat) = c('delta', 'pvalue')
	stat$pvalue[is.na(stat$pvalue)] = 1
	stat$delta[is.na(stat$delta)] = 0
	stat$padjust = p.adjust(stat$pvalue, 'fdr')
	out <- cbind(info, stat)
	write.table(out, file=paste(OutPre, ext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)

	max.delta = max(abs(stat$delta))
	draw(hist, stat$delta  , xlab='Deviance difference', ylab='Number of features', file=paste(OutPre, ext, 'md.hist', sep='.'))
	draw(hist, stat$pvalue , breaks=seq( 0,1,0.01), xlab='P-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'pval.hist', sep='.'))
	draw(hist, stat$padjust, breaks=seq( 0,1,0.01), xlab='Q-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'qval.hist', sep='.'))
	try(drawVolcano(data.frame(m=out$delta, q=out$pvalue ), title, PvalCutoff, 'p', path=paste(OutPre, pext, 'vol', sep='.')))
	try(drawVolcano(data.frame(m=out$delta, q=out$padjust), title,  FdrCutoff, 'q', path=paste(OutPre, qext, 'vol', sep='.')))

	log = out$pvalue<PvalCutoff
	if (sum(log) > 0) write.table(cbind(info[log,], stat[log,]), file=paste(OutPre, pext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
#	try(analysis(cbind(info[log,], stat[log,]), pct[log,], pext))

	log = out$padjust<FdrCutoff
	if (sum(log) > 0) write.table(cbind(info[log,], stat[log,]), file=paste(OutPre, qext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
#	try(analysis(cbind(inf[log,], stat[log,]), pct[log,], qext))
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

org <- loadFile(InFile)
bool <- org[[5]] >= MinCpgNum
cat('Removed features with #CpGs < ', MinCpgNum, ': ', nrow(org)-sum(bool), '\n')
org <- org[bool,]
colnames(org)[1:5] <- c('chr', 'start', 'end', 'strand', 'count')
info <- org[1:(DatCol-1)]
org <- org[DatCol:length(org)]
colnames(org)[seq(2, length(org), 2)] = colnames(org)[seq(1, length(org), 2)]

if (!is.null(opt$group)) {
	grp <- loadFile(opt$group, header=F, comment.char='#')
	if (length(grp) == 4) {
		colnames(grp) = c('file', 'flag', 'sample', 'group')
		grp <- unique(grp[3:4])
	} else if (length(grp) == 2) {
		colnames(grp) = c('sample', 'group')
		grp <- unique(grp[1:2])
	} else {
		stop("Unknown format for a group file")
	}
}

if (!is.null(opt$sample1) || !is.null(opt$sample2)) { #2
	stopifnot(!is.null(opt$sample1) && !is.null(opt$sample2))
	pairs = as.matrix(t(c(opt$sample1, opt$sample2)))
} else if (!is.null(opt$group1) || !is.null(opt$group2)) { #4
	stopifnot(!is.null(opt$group1) && !is.null(opt$group2))
	stopifnot(!is.null(opt$group))
	pairs = as.matrix(grp$sample[grp$group==opt$group1], grp$sample[grp$group==opt$group2])
} else if (is.null(opt$pair)) { #1
	pairs = t(combn(colnames(org)[seq(1, length(org), 2)], 2))
} else {
	pairs.tmp <- loadFile(opt$pair, header=F, comment.char='#')
	if (is.null(opt$group)) { #3
		pairs <- pairs.tmp
	} else { #5
		pairs <- do.call(rbind, apply(pairs.tmp, 1, function(x) expand.grid(grp$sample[grp$group==x[1]], grp$sample[grp$group==x[2]])))
	}
}


for (p in 1:nrow(pairs)) {
	Sam1 = as.character(pairs[p,1])
	Sam2 = as.character(pairs[p,2])
	cat(Sam1, ' vs ', Sam2, '\n')
	name = paste(Sam1, '-', Sam2, sep='')
	ext = paste('mlogit.', tolower(name), sep='')

	dat = org[colnames(org) == Sam1 | colnames(org) == Sam2]
	stopifnot(length(dat) == 4)
	colnames(dat)[c(2,4)] = colnames(dat)[c(1,3)]
	bool <- apply(dat[c(1,3)], 1, function(x) all(!is.na(x)) & min(x) > MinDepth)
	cat('Removed features with depth < ', MinDepth, ': ', nrow(dat)-sum(bool), '\n')
	dat = dat[bool,]
	head(dat)
	dim(dat)
	if (nrow(dat) == 0) next

	com <- stat.test(dat[c(2,4)], mlogit)
	output(info[bool,], dat, com, name, ext)
}

cat("End: ", date(), "\n")
