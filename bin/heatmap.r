FdrCutoff <- 0.05
PvalCutoff <- 0.01
MinCpgNum <- 5
MetDiff <- 0
MaxRowForHeatmap <- 10000
Top <- 500
PdfFlag = 0
DatCol = 4


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre') || !exists('GrpFile')) {
	if (length(pars) < 6) {
		message("USAGE: Rscript groupTest.r in_file out_prefix group_file [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	GrpFile = pars[3]
	if (length(pars) > 4) opt = getopt(matrix(ncol=4, byrow=T, c(
		'format',   'b', 1, 'character',
		'diff',     'd', 1, 'numeric',
		'fdr',      'f', 1, 'numeric',
		'pval',     'v', 1, 'numeric',
		'cpg',      'c', 1, 'integer',
		'top',      'T', 1, 'integer',
		'ttest',    't', 0, 'logical',
		'wilcoxon', 'w', 0, 'logical',
		'raoscott', 'r', 0, 'logical',
		'kstest',   'k', 0, 'logical',
		'vartest',  's', 0, 'logical',
		'ansari',   'a', 0, 'logical',
		'lepage',   'l', 0, 'logical',
		'pairfile', 'P', 1, 'character',
		'paired'  , 'p', 0, 'logical',
		'group1'  , '1', 0, 'logical',
		'group2'  , '2', 0, 'logical',
		'filter',   'F', 0, 'logical',
		'stat',     'S', 1, 'character',
		'dist',     'D', 1, 'character',
		'linkage',  'L', 1, 'character',
		'thread',   'm', 1, 'integer',
		'width' ,   'W', 1, 'integer',
		'height',   'H', 1, 'integer',
		'pdf'     , 'Y', 0, 'logical'
		)), pars[4:length(pars)])
	else opt = list()

	Format = ifelse(is.null(opt$format) || opt$format == 'auto', 0, ifelse(opt$format == 'sitemethy', 1, ifelse(opt$format == 'sitenread', 2, ifelse(opt$format == 'winmethy', 3, 4))))
# format = site methy
##Key1   Key2    Key3    KO1     KO2     KO3     WT1     WT2     WT3
#chr1    3010894 3010896 0.684   0.588   0.429   0.800   0.714   0.647
# format = site nread
##Key1   Key2    Key3    KO1     KO1     KO2     KO2     KO3     KO3     WT1     WT1     WT2     WT2     WT3     WT3
#chr1    3010894 3010896 6       19      7       17      4       7       2       10      2       7       6       17
# format = win methy
##Key1   Key2    Key3    KO1     KO1     KO1     KO2     KO2     KO2     KO3     KO3     KO3     WT1     WT1     WT1     WT2     WT2     WT2     WT3  WT3      WT3
#chr1    3010701 3010900 1       0.684   0.000   1       0.588   0.000   1       0.429   0.000   1       0.800   0.000   1       0.714   0.000   1    0.647    0.000
# format = win nread
##Key1   Key2    Key3    KO1     KO1     KO1     KO2     KO2     KO2     KO3     KO3     KO3     WT1     WT1     WT1     WT2     WT2     WT2     WT3  WT3      WT3
#chr1    3010701 3010900 1       13      19      1       10      17      1       3       7       1       8       10      1       5       7       1    11       17

	if (!is.null(opt$fdr     )) FdrCutoff    = as.numeric(opt$fdr)
	if (!is.null(opt$pval    )) PvalCutoff   = as.numeric(opt$pval)
	if (!is.null(opt$cpg     )) MinCpgNum    = as.integer(opt$cpg)
	if (!is.null(opt$diff    )) MetDiff      = as.numeric(opt$diff)
	if (!is.null(opt$top     )) Top          = as.numeric(opt$top)
	if ( is.null(opt$ttest   )) opt$ttest    = F
	if ( is.null(opt$wilcoxon)) opt$wilcoxon = F
	if ( is.null(opt$raoscott)) opt$raoscott = F
	if ( is.null(opt$kstest  )) opt$kstest   = F
	if ( is.null(opt$ansari  )) opt$ansari   = F
	if ( is.null(opt$vartest )) opt$vartest  = F
	if ( is.null(opt$lepage  )) opt$lepage   = F
	DistMethod = ifelse(is.null(opt$dist), 'euclidean', opt$dist)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	Paired = ifelse(is.null(opt$paired), F, T)
	Filter = ifelse(is.null(opt$filter), F, T)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}

suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
suppressMessages(library(aod))
if (opt$lepage) suppressMessages(library(NSM3))
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
	data <- read.table(file, sep="\t", header=header, check.names=F, comment.char=comment.char, ...)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

drawHeatmap <- function(x, col=bluered(100), file=NA, ...)
{
	Width = ifelse(!is.null(opt$width), opt$width, 300 + ncol(x)*ifelse(ncol(x) < 100, 8, 7))
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

	hv <- heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, cexCol=1.3, cexRow=1, margin=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ColSideColors=c(rep('green3', length(methy.igrp1)), rep('gold', length(methy.igrp2))), ...)

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
	mq.md    <- sum(mq$m < -MetDiff)
	mq.mu    <- sum(mq$m >  MetDiff)
	mq.md0q  <- sum(mq$m <  0       & mq$q < cutoff)
	mq.mu0q  <- sum(mq$m >  0       & mq$q < cutoff)
	mq.mdq   <- sum(mq$m < -MetDiff & mq$q < cutoff)
	mq.muq   <- sum(mq$m >  MetDiff & mq$q < cutoff)
	mq <- mq[mq$q != 0 & is.finite(mq$q),]
	if (!is.na(path)) openGraphics(path, 600, 500)

	max.m <- min(max(abs(mq$m)), 5)
	min.q <- min(mq$q)
	max.q <- min(-log10(min.q), 10)

	smoothScatter(mq$m, -log10(mq$q), nbin=512, nrpoints=0, xlab="Difference", ylab=paste("-log10(", prefix, "-value)", sep=''), main=title, xlim=c(-max.m,max.m), ylim=c(0,max.q), cex=ifelse(PdfFlag, 0.5, 1.5), cex.lab=1.2, cex.axis=1.2)
	abline(h=-log10(cutoff), col="darkgrey", lt=2)
	abline(v=-MetDiff, col="darkgrey", lt=2)
	abline(v= MetDiff, col="darkgrey", lt=2)

	if (max.m > MetDiff)
	{
		text(-max.m, 0, adj=c(0,0), mq.md, col='darkgray', cex=1.2)
		text( max.m, 0, adj=c(1,0), mq.mu, col='darkgray', cex=1.2)
	}
	if (min.q < cutoff)
	{
		sig <- mq[mq$q < cutoff & abs(mq$m) > MetDiff,]
		points(sig$m, -log10(sig$q), pch=20, cex=0.3, col='red')
		ypos <- max.q
		if (max.m > MetDiff)
		{
			text(-max.m, ypos, adj=c(0,1), mq.mdq, col='red', cex=1.2)
			text( max.m, ypos, adj=c(1,1), mq.muq, col='red', cex=1.2)
		}
		text(0, ypos, adj=c(MetDiff,1), paste(mq.md0q, ':', mq.mu0q), col='blue', cex=1.2)
	}

	if (!is.na(path)) closeGraphics()
}

ttest <- function(data)
{
	x <- data[as.character(grp$sample[Grp1 == grp$group])]
	y <- data[as.character(grp$sample[Grp2 == grp$group])]
	mean1 <- mean(as.matrix(x))
	mean2 <- mean(as.matrix(y))

	if (sd(x) == 0 & sd(y) == 0) x[1] = x[1] + 1e-10
	tval <- t.test(x, y, paired=Paired)

	return(c(mean1, mean2, mean1-mean2, tval$p.value))
}

wilcoxon <- function(tab)
{
	x <- tab[as.character(grp$sample[Grp1 ==  grp$group])]
	y <- tab[as.character(grp$sample[Grp2 ==  grp$group])]
	mean1 <- mean(as.matrix(x))
	mean2 <- mean(as.matrix(y))

	w <- wilcox.test(x, y, paired=Paired)
	vtest <- w$stat
	if (is.null(vtest))
	{
		vtest <- NA
		pvalue <- NA
	}
	else
	{
		pvalue <- w$p.value
	}

	return(c(mean1, mean2, mean1-mean2, pvalue))
}

kolmogorovsmirnov <- function(tab, exact=T)
{
	x <- tab[as.character(grp$sample[Grp1 ==  grp$group])]
	y <- tab[as.character(grp$sample[Grp2 ==  grp$group])]
	mean1 <- mean(as.matrix(x))
	mean2 <- mean(as.matrix(y))

	if (exact) {
		z <- noPairDuplicates(noDuplicates(x), noDuplicates(y))
		x <- z[,1]
		y <- z[,2]
	}

	ks <- ks.test(x, y)
	dtest <- ks$stat
	if (is.null(dtest))
	{
		dtest <- NA
		pvalue <- NA
	}
	else
	{
		pvalue <- ks$p.value
	}

	return(c(mean1, mean2, mean1-mean2, pvalue))
}

vartest <- function(data)
{
	x <- data[as.character(grp$sample[Grp1 == grp$group])]
	y <- data[as.character(grp$sample[Grp2 == grp$group])]
	sd1 <- sd(as.matrix(x))
	sd2 <- sd(as.matrix(y))

	if (sd(x) == 0 & sd(y) == 0) x[1] = x[1] + 1e-10
	vval <- var.test(x, y)

	return(c(sd1, sd2, sd1-sd2, vval$p.value))
}

ansari <- function(data)
{
	x <- data[as.character(grp$sample[Grp1 == grp$group])]
	y <- data[as.character(grp$sample[Grp2 == grp$group])]
	mean1 <- mean(as.matrix(x))
	mean2 <- mean(as.matrix(y))

	aval <- ansari.test(x, y)

	return(c(mean1, mean2, mean1-mean2, aval$p.value))
}

lepage <- function(data)
{
	x <- data[as.character(grp$sample[Grp1 == grp$group])]
	y <- data[as.character(grp$sample[Grp2 == grp$group])]
	mean1 <- mean(as.matrix(x))
	mean2 <- mean(as.matrix(y))

	lval <- pLepage(x, y)

	return(c(mean1, mean2, mean1-mean2, lval$p.val))
}

raotest <- function(var)
{
	i1 <- which(names(var) %in% as.character(grp$sample[Grp1==grp$group]))
	i2 <- which(names(var) %in% as.character(grp$sample[Grp2==grp$group]))
	z  <- data.frame(group=c(rep(Grp1, length(i1)), rep(Grp2, length(i2))))
	if (Format == 4) {
		z$y <- c(var[i1+1], var[i2+1])
		z$n <- c(var[i1+2], var[i2+2])
	} else if (Format == 2) {
		z$y <- c(var[i1+0], var[i2+0])
		z$n <- c(var[i1+1], var[i2+1])
	}
	rval <- raoscott(y/n ~ group, weights=n, data=z)
	pvalue <- pchisq(q=rval@X2, df=nrow(rval@tab)-1, lower.tail=FALSE)
	return(c(rval@tab$p[1], rval@tab$p[2], rval@tab$p[1]-rval@tab$p[2], pvalue))
}

noDuplicates <- function(x)
{
	for (j in 1:10) {
		i = which(duplicated(x))
		if (length(i) == 0) break
		x[i] = x[i]-1e-5
	}

	return(x)
}

noPairDuplicates <- function(x, y)
{
	z <- expand.grid(x, y)
	for (j in 1:10) {
		i = which(z[,1]==z[,2])
		if (length(i) == 0) break
		z.i = z[2]==z[i,2]
		z[z.i,2] = z[z.i,2]-1e-5
	}

	x2 <- z[1:length(x),1]
	y2 <- z[seq(1,nrow(z),length(x)),2]
	stopifnot(length(x) == length(x2))
	stopifnot(length(y) == length(y2))

	return(cbind(x2,y2))
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

	hv <- drawHeatmap(data, Rowv=T, Colv=T, file=paste(OutPre, suffix, 'heatmap', sep='.'))
	write.table(cbind(info[rev(hv$rowInd), ], data[rev(hv$rowInd), hv$colInd]), file=paste(OutPre, suffix, 'heatmap.csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
}

output <-function(pct, stat, ext)
{
	pext = paste(ext, paste(PvalCutoff, MetDiff, sep='-'), 'pval', sep='.')
	qext = paste(ext, paste(FdrCutoff, MetDiff, sep='-'), 'qval', sep='.')
	if (is.null(opt$stat)) {
		colnames(stat) = c(Grp1, Grp2, 'diff', 'pvalue')
		stat$pvalue[is.na(stat$pvalue)] = 1
		stat$diff[is.na(stat$diff)] = 0
		stat$padjust = p.adjust(stat$pvalue, 'fdr')
	}
	out <- cbind(pct[1:(DatCol.new-1)], stat)
	if (is.null(opt$stat)) write.table(out, file=paste(OutPre, ext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)

	# for TAB-seq
	max.diff = max(abs(stat$diff))
	max.diff = if (max.diff > 1) max.diff+0.05 else 1

	if (is.null(opt$stat)) {
		draw(hist, stat$diff   , breaks=seq(-max.diff,max.diff,0.05), xlab='Methylation difference', ylab='Number of features', xlim=c(-1, 1), file=paste(OutPre, ext, 'md.hist', sep='.'))
		draw(hist, stat$pvalue , breaks=seq( 0,1,0.01), xlab='P-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'pval.hist', sep='.'))
		draw(hist, stat$padjust, breaks=seq( 0,1,0.01), xlab='Q-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'qval.hist', sep='.'))
	}
	try(drawVolcano(data.frame(m=out$diff, q=out$pvalue ), pair, PvalCutoff, 'p', path=paste(OutPre, pext, 'vol', sep='.')))
	try(drawVolcano(data.frame(m=out$diff, q=out$padjust), pair,  FdrCutoff, 'q', path=paste(OutPre, qext, 'vol', sep='.')))

	log = out$pvalue<PvalCutoff & abs(out$diff)>MetDiff
	if (sum(log) > 0) write.table(cbind(pct[log,1:(DatCol-1)], stat[log,]), file=paste(OutPre, pext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	try(analysis(cbind(pct[log,1:(DatCol-1)], stat[log,]), pct[log,DatCol.new:length(pct)], pext))

	log = out$padjust<FdrCutoff & abs(out$diff)>MetDiff
	if (sum(log) > 0) write.table(cbind(pct[log,1:(DatCol-1)], stat[log,]), file=paste(OutPre, qext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	try(analysis(cbind(pct[log,1:(DatCol-1)], stat[log,]), pct[log,DatCol.new:length(pct)], qext))

	if (sum(log) > Top*2)
	{
		log.2 = head(order(abs(stat[log, 'diff']), decreasing=T), Top)
#		log.2 = head(order(stat[log, 'padjust']), Top)
		try(analysis(cbind(pct[log,][log.2, 1:(DatCol-1)], stat[log,][log.2,]), pct[log,][log.2, DatCol.new:length(pct)], paste(qext, 'top', Top, sep='')))
#		if (sum(log.2) > 0) write.table(cbind(pct[log,][log.2,1:(DatCol-1)], stat[log,][log.2,]), file=paste(OutPre, ext, paste('q', FdrCutoff, 'diff', MetDiff, 'top', Top, sep=''), 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	}
}

stat.test <- function(x, fun, ext)
{
	if (!is.null(opt$stat)) {
		com = loadFile(paste(opt$stat, ext, 'csv', sep='.'))
		idx = which(colnames(com) == 'diff')
		stopifnot(length(idx) == 1)
		return(com[(idx-2):(idx+2)])
	} else if (is.null(opt$thread) || opt$thread == 1) {
		return(data.frame(t(apply(x, 1, fun))))
	} else {
		return(data.frame(do.call(rbind, mclapply(split(x, 1:nrow(x)), function(y) fun(unlist(y)), mc.cores=opt$thread))))
	}
}


cat("Start: ", date(), "\n")

grp <- loadFile(GrpFile, header=F, comment.char='#')
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

if (!is.null(opt$group1) || !is.null(opt$group2)) {
	stopifnot(!is.null(opt$group1) && !is.null(opt$group2))
	pairs = as.matrix(t(c(opt$group1, opt$group2)))
} else if (is.null(opt$pairfile)) {
	pairs = t(combn(levels(grp$group), 2))
} else {
	pairs <- loadFile(opt$pairfile, header=F, comment.char='#')
}

org <- loadFile(InFile)
colnames(org)[1:3] <- c('chr', 'start', 'end')
if (colnames(org)[DatCol] == paste('Key', DatCol, sep='')) {
	colnames(org)[DatCol] <- 'ID'
	DatCol = DatCol+1
}

org = na.omit(org)
cat('Removed features with NA: ', length(na.action(org)), '\n')
org.row = nrow(org)
if (colnames(org)[DatCol-1] == 'ID') org$ID = factor(org$ID)
col.tab = table(colnames(org)[DatCol:length(org)])

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
#	icpg <- seq(DatCol+0, length(org), 3)
	imet <- seq(DatCol+1, length(org), 3)
#	itot <- seq(DatCol+2, length(org), 3)
}

if (Format == -1) {
	Format = ifelse(max(org[imet[1]], na.rm=T) > 1, 4, 3)
} else if (Format == 2) {
	org = org[apply(org[seq(DatCol+1, length(org), 2)], 1, function(x) all(x!=0)), ]
	cat('Removed features with 0: ', org.row-nrow(org), '\n')
}

for (p in 1:nrow(pairs)) {
	Grp1 = as.character(pairs[p,1])
	Grp2 = as.character(pairs[p,2])
	cat(Grp1, ' vs ', Grp2, '\n')

	igrp1 <- apply(as.matrix(as.character(grp$sample[Grp1 == grp$group])), 1, function(x) which(x == colnames(org))[1])
	stopifnot(length(igrp1) != 0)
	igrp2 <- apply(as.matrix(as.character(grp$sample[Grp2 == grp$group])), 1, function(x) which(x == colnames(org))[1])
	stopifnot(length(igrp2) != 0)

	if (Format == 1) {
		methy = org[c(1:(DatCol-1), igrp1, igrp2)]
		colnames(methy)[DatCol:length(methy)] = colnames(org)[c(igrp1, igrp2)]
		DatCol.new = DatCol
	} else if (Format == 2) {
		count = org[c(1:(DatCol-1), as.vector(rbind(igrp1, igrp1+1)), as.vector(rbind(igrp2, igrp2+1)))]
		colnames(count)[seq(DatCol+1, length(count), 2)] = colnames(org)[c(igrp1, igrp2)]
		methy = cbind(count[1:(DatCol-1)], round(count[seq(DatCol, length(count), 2)]/count[seq(DatCol+1, length(count), 2)], 3))
		DatCol.new = DatCol
	} else {
		count = org[c(1:(DatCol-1), as.vector(rbind(igrp1, igrp1+1, igrp1+2)), as.vector(rbind(igrp2, igrp2+1, igrp2+2)))]
		colnames(count)[seq(DatCol+1, length(count), 3)] = colnames(org)[c(igrp1, igrp2)]
		colnames(count)[seq(DatCol+2, length(count), 3)] = colnames(org)[c(igrp1, igrp2)]
		nsite = apply(count[seq(DatCol, length(count), 3)], 1, function(x) ifelse(any(is.na(x)), 0, min(x)))
		cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(nsite < MinCpgNum), '\n')
		count = count[nsite >= MinCpgNum,]
		nsite = nsite[nsite >= MinCpgNum ]
		if (colnames(org)[DatCol-1] == 'ID') count$ID = factor(count$ID)
		print(dim(count))

		# window.methy or window.nread format
		if (Format == 4) {
			methy = cbind(count[1:(DatCol-1)], nsite, round(count[seq(DatCol+1, length(count), 3)]/count[seq(DatCol+2, length(count), 3)], 3))
		} else {
			methy = cbind(count[1:(DatCol-1)], nsite, count[seq(DatCol+1, length(count), 3)])
		}
		DatCol.new = DatCol + 1
		colnames(methy)[DatCol] = 'MinCpgNum'
	}

	print(head(methy))
	print(tail(methy))
	print(dim(methy))

	methy.igrp1 <- apply(as.matrix(as.character(grp$sample[Grp1 == grp$group])), 1, function(x) which(x == colnames(methy))[1])
	cat('Group 1 [', Grp1, '] =', methy.igrp1, '[', colnames(methy)[methy.igrp1], ']\n')
	stopifnot(length(methy.igrp1) != 0)
	methy.igrp2 <- apply(as.matrix(as.character(grp$sample[Grp2 == grp$group])), 1, function(x) which(x == colnames(methy))[1])
	cat('Group 2 [', Grp2, '] =', methy.igrp2, '[', colnames(methy)[methy.igrp2], ']\n')
	stopifnot(length(methy.igrp2) != 0)

	if (Filter) {
		log = apply(methy[DatCol.new:length(methy)], 1, range)
		log = log[2,] - log[1,] < MetDiff
		methy = methy[!log,]
		count = count[!log,]
		cat('Removed features with methylation difference <', MetDiff, ': ', sum(log), '\n')
		cat('Dimension =', dim(methy), '\n')
		type = paste('filter', MetDiff, '.', sep='')
	} else {
		type =''
	}

	pair = paste(Grp1, '-', Grp2, sep='')
	cat('Dimension =', dim(methy), '\n')
	
	if (opt$ttest) {
		type2 = paste(type, 'ttest', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], ttest, ext)
		output(methy, com, ext)
	}
	
	if (opt$wilcoxon) {
		type2 = paste(type, 'wilcoxon', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], wilcoxon, ext)
		output(methy, com, ext)
	}
	
	if (opt$kstest) {
		type2 = paste(type, 'kstest', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], kolmogorovsmirnov, ext)
		output(methy, com, ext)
	}
	
	if (opt$vartest) {
		type2 = paste(type, 'vartest', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], vartest, ext)
		output(methy, com, ext)
	}
	
	if (opt$ansari) {
		type2 = paste(type, 'ansari', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], ansari, ext)
		output(methy, com, ext)
	}
	
	if (opt$lepage) {
		type2 = paste(type, 'lepage', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(methy[DatCol.new:length(methy)], lepage, ext)
		output(methy, com, ext)
	}
	
	if (opt$raoscott && (Format == 2 || Format == 4)) {
		type2 = paste(type, 'raoscott', sep='')
		print(type2)
		ext = paste(type2, tolower(pair), sep='.')
		com <- stat.test(count[DatCol:length(count)], raotest, ext)
		output(methy, com, ext)
	}
}

cat("End: ", date(), "\n")
