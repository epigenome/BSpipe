FdrCutoff <- 0.05
PvalCutoff <- 0.01
MinCpgNum <- 5
MetDiff <- 0.25
Top <- 500
MaxRowForHeatmap <- 10000
PdfFlag = 0
DatCol = 4


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre') || !exists('GrpFile')) {
	if (length(pars) < 3) {
		message("USAGE: Rscript fisherGroup.r in_file out_prefix group_file [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	GrpFile = pars[3]
	if (length(pars) > 4) opt = getopt(matrix(ncol=4, byrow=T, c(
		'diff',     'd', 1, 'numeric',
		'fdr',      'f', 1, 'numeric',
		'pval',     'v', 1, 'numeric',
		'cpg',      'c', 1, 'integer',
		'top',      'T', 1, 'integer',
		'pairfile', 'P', 1, 'character',
		'fraction', 'R', 1, 'numeric',
		'dist',     'D', 1, 'character',
		'linkage',  'l', 1, 'character',
		'thread',   'm', 1, 'integer',
		'width' ,   'W', 1, 'integer',
		'height',   'H', 1, 'integer',
		'pdf'     , 'Y', 0, 'logical'
		)), pars[4:length(pars)])
	else opt = list()

	if (!is.null(opt$fdr     )) FdrCutoff    = as.numeric(opt$fdr)
	if (!is.null(opt$pval    )) PvalCutoff   = as.numeric(opt$pval)
	if (!is.null(opt$cpg     )) MinCpgNum    = as.integer(opt$cpg)
	if (!is.null(opt$diff    )) MetDiff      = as.numeric(opt$diff)
	if (!is.null(opt$top     )) Top          = as.numeric(opt$top)
	if ( is.null(opt$fraction)) opt$fraction = 1
	DistMethod = ifelse(is.null(opt$dist), 'euclidean', opt$dist)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}


suppressMessages(library(plyr))
suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
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

methy <- function(x)
{
	x[is.nan(x)] = NA
	return(x)
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

	hv <- heatmap.2(as.matrix(na.omit(x)), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, cexCol=1.3, cexRow=1, margin=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ColSideColors=c(rep('green3', length(igrp1.name)), rep('gold', length(igrp2.name))), ...)

	if (!is.na(file)) closeGraphics()
	return(hv)
}

draw <- function(fun, x, file=NA, width=600, height=500, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	fun(x, ...)

	if (!is.na(file)) closeGraphics()
}

drawVolcano <- function(mq, title, cutoff, path=NA)
{
	mq.md   <- sum(mq$m < -MetDiff                , na.rm=T)
	mq.mu   <- sum(mq$m >  MetDiff                , na.rm=T)
	mq.md0q <- sum(mq$m <        0 & mq$q < cutoff, na.rm=T)
	mq.mu0q <- sum(mq$m >        0 & mq$q < cutoff, na.rm=T)
	mq.mdq  <- sum(mq$m < -MetDiff & mq$q < cutoff, na.rm=T)
	mq.muq  <- sum(mq$m >  MetDiff & mq$q < cutoff, na.rm=T)
	mq$m[is.na(mq$m)] <- 0
	mq <- mq[mq$q != 0 & is.finite(mq$q),]
	if (!is.na(path)) openGraphics(path, 600, 500)

	max.m <- min(max(abs(mq$m)), 5)
	min.q <- min(mq$q)
	max.q <- min(-log10(min.q), 10)

	smoothScatter(mq$m, -log10(mq$q), nbin=512, nrpoints=0, xlab="Difference", ylab="-log10(q-value)", main=title, xlim=c(-max.m,max.m), ylim=c(0,max.q), cex=ifelse(PdfFlag, 0.5, 1.5), cex.lab=1.2, cex.axis=1.2)
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

analysis <- function(stat, data, path, nelm=0)
{
	if (nrow(data) <= 1) return(invisible())

	if (nelm == 0) {
		nelm = min(nrow(data), MaxRowForHeatmap)
		write.table(cbind(data, stat), file=paste(path, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	}

	if (nrow(data) < nelm) nelm = nrow(data)
	if (nrow(data) > nelm) {
		is = grep('padjust', colnames(stat))
		if (length(is) > 1) {
			pad.mean = apply(stat[is], 1, mean)
			data = head(data[order(pad.mean),], nelm)
			stat = head(stat[order(pad.mean),], nelm)
		} else {
			data = head(data[order(stat$pad),], nelm)
			stat = head(stat[order(stat$pad),], nelm)
		}
	}
	analysis2(cbind(data[1:(DatCol-1)], stat), data[DatCol.new:length(data)], path)
}

analysis2 <- function(info, data, path)
{
	if (nrow(data) <= 1) return(invisible())

	hv <- drawHeatmap(data, Rowv=T, Colv=T, file=paste(path, 'heatmap', sep='.'))
	write.table(cbind(info[rev(hv$rowInd), ], data[rev(hv$rowInd), hv$colInd]), file=paste(path, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
}

fisher <- function(tab, mean=F)
{
	var <- matrix(as.matrix(tab), nrow=2, ncol=2, byrow=T)
	if (all(!is.na(var[,2]) & var[,2] > 0)) {
		var[,2] = var[,2] - var[,1]
		colnames(var) <- c('methy', 'unmethy')
		rownames(var) <- names(tab)[c(1,3)]
		f <- fisher.test(var)
		ratio <- f$est
		pvalue <- f$p.value
	} else {
		ratio <- NA
		pvalue <- NA
	}

	if (mean) return(c(var[,1]/(var[,1]+var[,2]), ratio, pvalue))
	else      return(c(                           ratio, pvalue))
}

mqf <- function(x, op)
{
	if (length(x) == 2) return(x)
	if (sum(is.na(x)) != 0) {
		flag = T
		x[is.na(x)] = 0
	} else {
		flag = F
	}

	y=matrix(x, ncol=2, byrow=T)
	z=apply(y, 1, function(v) v[1]<FdrCutoff & abs(v[2])>MetDiff)
	z.a = subset(y, z==T)
	z.b = subset(y, z==F)
	if (ifelse(op == 'common', all(z), any(z))) ret = z.a[which.min(z.a[,1]),]
	else                                        ret = z.b[which.max(z.b[,1]),]

	m=apply(y, 1, function(v) v[2] > 0)
	if (flag || length(table(m)) > 1) ret[2] = NA
	return(ret)
}

stat.test <- function(x, fun)
{
	if (is.null(opt$thread) || opt$thread == 1) {
		return(data.frame(t(apply(x, 1, fun))))
	} else {
		return(data.frame(do.call(rbind, mclapply(split(x, 1:nrow(x)), function(y) fun(unlist(y)), mc.cores=opt$thread))))
	}
}

output <- function(var, log, ind1, ind2, subtype=NA)
{
	path1 = paste(Path, paste(tolower(colnames(dat)[ind1]), tolower(colnames(dat)[ind2]), sep='-'), sep='.')
	path2 = ifelse(is.na(subtype), path1, paste(path1, tolower(subtype), sep='.'))
	path3 = paste(path2, paste(FdrCutoff, MetDiff, sep='-'), 'qval', sep='.')
	write.table(cbind(out[1:(DatCol.new-1)], out[colnames(dat)[c(ind1,ind2)]], var), file=paste(path2, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
#	cat('#features with NaN for fisher  test:', sum(is.na(out[,'pvalue'])), '\n')

	mq = data.frame(t(apply(var[2:3], 1, mqf, ifelse(is.na(subtype), 'common', subtype))))
	colnames(mq) = c('q', 'm')
	stopifnot(sum(log, na.rm=T) == sum(mq$q < FdrCutoff & abs(mq$m) > MetDiff, na.rm=T))
	draw(hist, mq$q, breaks=20, xlab='Adjusted p-value', ylab='Number of features', file=paste(path2, 'hist', sep='.'))
	drawVolcano(mq, sub('-', ' - ', pair), FdrCutoff, path=paste(path3, 'vol', sep='.'))

	try(analysis(var[log,], cbind(out[log,1:(DatCol.new-1)], out[log,c(igrp1.name,igrp2.name)]), path3))
	if (sum(log, na.rm=T) > Top*2) try(analysis(var[log,], cbind(out[log,1:(DatCol.new-1)], out[log,c(igrp1.name,igrp2.name)]), paste(path2, paste(FdrCutoff, '-', MetDiff, '.top', Top, sep=''), sep='.'), nelm=Top))
}

unionOutput <- function(tes, all.mq, subtype1, subtype2)
{
	path2 = paste(Path, tolower(subtype1), tolower(subtype2), sep='.')
	path3 = paste(path2, paste(FdrCutoff, MetDiff, sep='-'), 'qval', sep='.')

	mq = data.frame(t(apply(if (is.null(all.mq)) tes else all.mq, 1, mqf, subtype1)))
	colnames(mq) = c('q', 'm')
#	write.table(cbind(out[1:(DatCol.new-1)], combined_padjust=mq$q, combined_diff=mq$m, tes), file=paste(path2, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	log = mq$q < FdrCutoff & abs(mq$m) > MetDiff & !is.na(mq$m)
	if (sum(log, na.rm=T) == 0) return(invisible())

	draw(hist, mq$q, breaks=20, xlab='Adjusted p-value', ylab='Number of features', file=paste(path2, 'hist', sep='.'))
	drawVolcano(mq, sub('-', ' - ', pair), FdrCutoff, path=paste(path3, 'vol', sep='.'))
	try(analysis(tes[log,], cbind(out[log,1:(DatCol.new-1)], out[log,c(igrp1.name,igrp2.name)]), path3))
	if (sum(log, na.rm=T) > Top*2) try(analysis(tes[log,], cbind(out[log,1:(DatCol.new-1)], out[log,c(igrp1.name,igrp2.name)]), paste(path2, paste(FdrCutoff, '-', MetDiff, '.top', Top, sep=''), sep='.'), nelm=Top))
#	if (sum(log, na.rm=T) > 0) write.table(cbind(out[log,1:(DatCol-1)], combined_padjust=mq$q[log], combined_diff=mq$m[log], tes[log,]), file=paste(path3, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
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

#org = na.omit(org)
#cat('Removed features with NA: ', length(na.action(org)), '\n')
if (colnames(org)[DatCol-1] == 'ID') org$ID = factor(org$ID)
col.tab = table(colnames(org)[DatCol:length(org)])

if (all(col.tab == 2)) {
	Format = 1
	itot <- seq(DatCol+1, length(org), 2)
	DatCol.new = DatCol
} else if (all(col.tab == 3)) {
	Format = 2
	itot <- seq(DatCol+2, length(org), 3)
	DatCol.new = DatCol + 1
} else {
	stop('File format is not a nread format')
}

#log = apply(org[itot], 1, function(x) any(x==0))
#org = org[!log, ]
#cat('Removed features with 0 reads: ', sum(log), '\n')

if (Format == 2) { # for windows
	cpg = org[itot-2]
	org = org[-(itot-2)]
	itot <- seq(DatCol+1, length(org), 2)
	colnames(org)[itot] = colnames(org)[itot-1]
}

for (p in 1:nrow(pairs)) {
	Grp1 = as.character(pairs[p,1])
	Grp2 = as.character(pairs[p,2])
	Pair = paste(Grp1, Grp2, sep='-')
	Path = paste(OutPre, 'fisher', tolower(Pair), sep='.')

	org.igrp1 <- apply(as.matrix(as.character(grp$sample[Grp1 == grp$group])), 1, function(x) which(x == colnames(org))[1])
	stopifnot(length(org.igrp1) != 0)
	org.igrp2 <- apply(as.matrix(as.character(grp$sample[Grp2 == grp$group])), 1, function(x) which(x == colnames(org))[1])
	stopifnot(length(org.igrp2) != 0)

	igrp1.name = colnames(org[org.igrp1])
	igrp2.name = colnames(org[org.igrp2])

	idx = as.vector(apply(as.matrix(c(org.igrp1, org.igrp2)), 1, function(x) x:(x+1)))
	log = apply(org[DatCol:length(org)], 1, function(x) sum(!is.na(x[org.igrp1-DatCol+2]) & x[org.igrp1-DatCol+2] > 0)/length(org.igrp1) < opt$fraction | sum(!is.na(x[org.igrp2-DatCol+2]) & x[org.igrp2-DatCol+2] > 0)/length(org.igrp2) < opt$fraction)
	dat = cbind(org[!log, 1:(DatCol-1)], org[!log, idx])
	if (Format == 2) cpg = cpg[!log,]
	colnames(dat)[DatCol:length(dat)] = colnames(org)[idx]
	cat('Removed features with sample fraction <', opt$fraction, ': ', sum(log), '\n')

	if (Format == 1) {     # for sites
		out = cbind(dat[1:(DatCol-1)], data.frame(t(apply(dat[DatCol:length(dat)], 1, function(x) methy(x[seq(1,length(x),2)]/x[seq(2,length(x),2)])))))
		colnames(out)[DatCol:length(out)] = colnames(org)[c(org.igrp1, org.igrp2)]
	} else {               # for windows
		cpg.min = apply(cpg[c(igrp1.name,igrp2.name)], 1, min)
		cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(cpg.min < MinCpgNum), '\n')
		dat     = dat    [cpg.min >= MinCpgNum,]
		cpg.min = cpg.min[cpg.min >= MinCpgNum ]
		if (colnames(dat)[DatCol-1] == 'ID') dat$ID = factor(dat$ID)
		print(dim(dat))
		out = cbind(dat[1:(DatCol-1)], MinCpgNum=cpg.min, data.frame(t(apply(dat[DatCol:length(dat)], 1, function(x) methy(x[seq(1,length(x),2)]/x[seq(2,length(x),2)])))))
		colnames(out)[DatCol:length(out)] = c('MinCpgNum', colnames(org)[c(org.igrp1, org.igrp2)])
#		write.table(out, file=paste(OutPre, 'pct.csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	}

	dat.igrp1 <- apply(as.matrix(as.character(grp$sample[Grp1 == grp$group])), 1, function(x) which(x == colnames(dat))[1])
	stopifnot(length(dat.igrp1) != 0)
	dat.igrp2 <- apply(as.matrix(as.character(grp$sample[Grp2 == grp$group])), 1, function(x) which(x == colnames(dat))[1])
	stopifnot(length(dat.igrp2) != 0)

	all.tes = data.frame()
#	all.com.com = logical(nrow(dat))
#	all.uni.com = logical(nrow(dat))
#	all.uni.uni = logical(nrow(dat))
	all.mq.com = data.frame()

	for (i in dat.igrp1) {
		print(paste('Running', colnames(dat)[i], '...'))
		tes = data.frame()
		for (j in dat.igrp2) {
			pair = paste(colnames(dat)[i], colnames(dat)[j], sep='-')
			tem = stat.test(dat[c(i+0:1,j+0:1)], fisher)
			colnames(tem) = c('ratio', 'pvalue')
			tem$padjust = p.adjust(tem$pvalue, 'fdr')
			tem$diff = out[[colnames(dat)[i]]] - out[[colnames(dat)[j]]]
			if (length(tes) == 0) tes = tem[c('padjust', 'diff')]
			else                  tes = cbind(tes, tem[c('padjust', 'diff')])
			colnames(tes)[(length(tes)-1):length(tes)] = paste(pair, c('padjust', 'diff'), sep='.')
			output(tem[c('pvalue', 'padjust', 'diff')], tem$padjust<FdrCutoff & abs(tem$diff)>MetDiff, i, j)
		}

		lcom = apply(tes, 1, function(x) all(apply(as.matrix(seq(1, length(x), 2)), 1, function(k) x[k]<FdrCutoff & abs(x[k+1])>MetDiff)))
		luni = apply(tes, 1, function(x) any(apply(as.matrix(seq(1, length(x), 2)), 1, function(k) x[k]<FdrCutoff & abs(x[k+1])>MetDiff)))

		if (length(all.tes) == 0) all.tes = tes
		else                      all.tes = cbind(all.tes, tes)
#		all.com.com = all.com.com & lcom
#		all.uni.com = all.uni.com | lcom
#		all.uni.uni = all.uni.uni | luni

		mqcom = data.frame(t(apply(tes, 1, mqf, 'common')))
		if (length(all.mq.com) == 0) all.mq.com = mqcom
		else                         all.mq.com = cbind(all.mq.com, mqcom)

	#	output(tes, lcom, i, 'common')
	#	output(tes, luni, i, 'union')
	}

	if (length(dat.igrp1)*length(dat.igrp2) > 1) {
		unionOutput(all.tes, all.mq.com, 'common', 'common')
		if (length(dat.igrp1) > 1 && length(dat.igrp2) > 1)
			unionOutput(all.tes, all.mq.com, 'union', 'common')
		unionOutput(all.tes, NULL, 'union', 'union')
		write.table(cbind(out[1:(DatCol.new-1)], all.tes), file=paste(Path, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
	}
}	

cat("End: ", date(), "\n")
