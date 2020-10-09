MaxRow <- 5000
MinCpgNum <- 5
#MetDiff <- 0.25
#SdDiff <- 0.2
Top <- 500
PdfFlag = 0
DatCol = 4
GrpColors <- c('green3', 'gold', 'purple', 'orange', 'magenta', 'cyan', 'blue3', 'red', 'pink', 'black')


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre')) {
	if (length(pars) < 2) {
		message("USAGE: Rscript clustering.r in_file out_prefix [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	if (length(pars) > 3) opt = getopt(matrix(ncol=4, byrow=T, c(
		'format'      ,  'b', 1, 'character',
		'hierarchical' , 'h', 0, 'logical',
		'cut'          , 'n', 1, 'integer',
		'kmeans'       , 'k', 1, 'integer',
		'pam'          , 'p', 1, 'integer',
		'diff'         , 'd', 1, 'numeric',
		'sd'           , 's', 1, 'numeric',
		'cpg'          , 'c', 1, 'integer',
		'top'          , 't', 1, 'integer',
		'group'        , 'g', 1, 'character',
		'only'         , 'O', 0, 'logical',
		'dist'         , 'D', 1, 'character',
		'linkage'      , 'L', 1, 'character',
		'width'        , 'W', 1, 'integer',
		'height'       , 'H', 1, 'integer',
		'pdf'          , 'Y', 0, 'logical'
		)), pars[3:length(pars)])
	else opt = list()

	Format = ifelse(is.null(opt$format) || opt$format == 'auto', 0, ifelse(opt$format == 'sitemethy', 1, ifelse(opt$format == 'sitenread', 2, ifelse(opt$format == 'winmethy', 3, 4))))
	if (!is.null(opt$cpg   )) MinCpgNum = as.integer(opt$cpg)
	if (!is.null(opt$diff  )) MetDiff   = as.numeric(opt$diff)
	if (!is.null(opt$sd    )) SdDiff    = as.numeric(opt$sd)
	if (!is.null(opt$top   )) Top       = as.numeric(opt$top)
	if (!is.null(opt$group )) GrpFile   = opt$group
	DistMethod = ifelse(is.null(opt$dist), 'euclidean', opt$dist)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}

suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
suppressMessages(library(cluster))
suppressMessages(library(fpc))


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

loadFile <- function(file, header=T, comment.char='')
{
	print(paste("Loading", file, ".."))
	data <- read.table(file, sep="\t", header=header, check.names=F, comment.char=comment.char)
	print(head(data))
	print(tail(data))
	cat("dim=", dim(data), "\n")
	return(data)
}

drawHist <- function(var, file=NA, ...)
{
	if (!is.na(file)) openGraphics(file, 600, 400)
#	r <- range(as.matrix(var))
#	h <- hist(as.matrix(var), breaks=c((r[1]-0.5):(r[2]+0.5)), ...)
	h <- hist(as.matrix(var), breaks=100, ...)
	print(h$counts)
	if (!is.na(file)) closeGraphics()
}

hclust2 <- function(x, method=Linkage, ...) hclust(x, method=method, ...)
dist2 <- function(x, ...) if (DistMethod == 'pearson' | DistMethod == 'spearman' | DistMethod == 'kendall') as.dist(1-cor(t(x), method=DistMethod), ...) else dist(x, method=DistMethod, ...)

drawHeatmap <- function(x, col=bluered(100), file=NA, ...)
{
	Width = ifelse(!is.null(opt$width), opt$width, 300 + ncol(x)*ifelse(ncol(x) < 100, 8, 7))
	Height = ifelse(!is.null(opt$height), opt$height, 300 + sqrt(nrow(x))*10)
	if (!is.na(file)) openGraphics(file, Width, Height)

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

	if (exists('grp.colors')) hv <- heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, cexCol=1.4, cexRow=1, margins=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ColSideColors=grp.colors, ...)
	else                     hv <- heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, cexCol=1.4, cexRow=1, margins=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ...)

	if (!is.na(file)) closeGraphics()
	return(hv)
}

hierarchicalClust <- function(data, suffix)
{
	print(suffix)
	if (nrow(data) <= 1 || nrow(data) > MaxRow) return(invisible())

	hv <- drawHeatmap(data[DatCol.new:length(data)], Rowv=T, Colv=T, file=paste(OutPre, suffix, 'heatmap', sep='.'))
	write.table(data[rev(hv$rowInd), c(1:(DatCol.new-1), DatCol.new-1+hv$colInd)], file=paste(OutPre, suffix, 'heatmap', 'csv', sep='.'), quote=F, sep="\t", row.names=F, col.names=T)

	if (!is.null(opt$cut) && suffix == 'hclust.all') {
		hclustTree(data, suffix)
	}
}

draw <- function(fun, x, y=NULL, file=NA, width=600, height=500, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	fun(x, y, ...)

	if (!is.na(file)) closeGraphics()
}

plotEachCluster <- function(dat, clu, box=F, xlab=NULL, ylab=NULL, cex.lab=1.3, cex.axis=1.4)
{
	old.par = par(no.readonly = TRUE)
	on.exit(par(old.par))
	n = max(clu) + 1
	par(mfrow=c(n,1))
#	par(cex = 0.6)
	par(mar = c(0, 1, 0, 1), oma = c(6, 5, 1, 3))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	ylim = c(-0.1,1.1)

	for (i in 1:n) {
#		plot(1, axes = FALSE, type = "n")
		if (box) {
			if (i == 1) boxplot(dat           , ylim=ylim, xaxt='n', las=2, col=1:ncol(dat), cex.axis=cex.axis)
			else        boxplot(dat[clu==i-1,], ylim=ylim, xaxt='n', las=2, col=1:ncol(dat), cex.axis=cex.axis)
		} else {
			if (i == 1) matplot(t(dat           ), ylim=ylim, type='l', lty=1, xaxt='n', las=2, col=i, cex.axis=cex.axis)
			else        matplot(t(dat[clu==i-1,]), ylim=ylim, type='l', lty=1, xaxt='n', las=2, col=i, cex.axis=cex.axis)
		}

		if (i == n) axis(1, 1:ncol(dat), colnames(dat), cex.axis=cex.axis, las=2)
#		mtext(letters[i], side = 3, line = -1, adj = 0.1, cex = 0.8, col = "grey40")
#		if (i %in% c(4, 5, 6)) axis(1, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
#		if (i %in% c(1, 4)) axis(2, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
		box(col = "grey60")
		mtext(ifelse(i == 1, 'All', paste('Cluster', i-1)), side=4, line=1, col="black", cex=1)
	}
	mtext(xlab, side=1, outer=T, line=2.2, col="black", cex=cex.lab)
	mtext(ylab, side=2, outer=T, line=2.2, col="black", cex=cex.lab)
}

plotClusterMiddle <- function(dat, clu, fun2, xlab=NULL, ylab=NULL, cex.lab=1.2, cex.axis=1)
{
	# get cluster means
	x = aggregate(dat,by=list(clu),FUN=get(fun2))
	n = 1:(length(x)-1)
	matplot(t(x[n+1]), type='b', lty=1, xaxt='n', col=n+1, ylim=c(0,1), xlab=xlab, ylab=paste(sub('(.)', '\\U\\1\\E', fun2, per=T), 'of', tolower(ylab)), cex.lab=cex.lab)
	axis(1, n, colnames(x)[n+1], las=2, cex.axis=cex.axis)
}

drawCluster <- function(pure, clu, file=NA, main=NULL, xlab=NULL, ylab=NULL)
{
	Ylab = 'Methylation'

	if (is.na(file)) {
		file1 = NA
		file2 = NA
		file3 = NA
		file4 = NA
		file5 = NA
		file6 = NA
	} else {
		file1 = paste(file, 'clplot1', sep='.')
		file2 = paste(file, 'clplot2', sep='.')
		file3 = paste(file, 'line', sep='.')
		file4 = paste(file, 'box', sep='.')
		file5 = paste(file, 'mean', sep='.')
		file6 = paste(file, 'median', sep='.')
	}

	clu.num = length(unique(clu))
	Width = ifelse(!is.null(opt$width), opt$width, 200 + ncol(pure)*10)
	Height = 200+(clu.num+1)*ifelse(clu.num <= 10, 100, ifelse(clu.num <= 50, 50, 25))

	draw(clusplot         , pure, clu, file=file1, color=T, shade=T, labels=4, lines=0, main=main)
	draw(plotcluster      , pure, clu, file=file2)
	draw(plotEachCluster  , pure, clu, file=file3, ylab=Ylab, width=Width, height=Height)
	draw(plotEachCluster  , pure, clu, file=file4, ylab=Ylab, width=Width, height=Height, box=T)
	draw(plotClusterMiddle, pure, clu, file=file5, ylab=Ylab, fun2='mean'  , width=Width, height=600)
	draw(plotClusterMiddle, pure, clu, file=file6, ylab=Ylab, fun2='median', width=Width, height=600)
}

hclustTree <- function(data, suffix)
{
	dis <- dist2(data[DatCol.new:length(data)])
	hcl <- hclust2(dis)
	groups = cutree(hcl, opt$cut)
	drawCluster(data[DatCol.new:length(data)], groups, file=paste(OutPre, '.hcl', Ncluster, sep=''))
	write.table(cbind(data[c(1:(DatCol-1))], cluster=fit$cluster), file=paste(OutPre, suffix, opt$kmeans, 'csv', sep='.'), quote=F, sep="\t", row.names=F, col.names=T)
}

kmeansClust <- function(data, suffix)
{
	print(suffix)
	if (nrow(data) <= 1) return(invisible())

	fit <- kmeans(data[DatCol.new:length(data)], opt$kmeans)
	print(table(fit$cluster))
	drawCluster(data[DatCol.new:length(data)], fit$cluster, file=paste(OutPre, suffix, opt$kmeans, sep='.'))
	write.table(cbind(data[c(1:(DatCol-1))], cluster=fit$cluster), file=paste(OutPre, suffix, opt$kmeans, 'csv', sep='.'), quote=F, sep="\t", row.names=F, col.names=T)
}

pamClust <- function(data, suffix)
{
	print(suffix)
	if (nrow(data) <= 1) return(invisible())

	fit <- pam(data[DatCol.new:length(data)], opt$pam)
	print(table(fit$cluster))
	drawCluster(data[DatCol.new:length(data)], fit$cluster, file=paste(OutPre, suffix, opt$pam, sep='.'))
	write.table(cbind(data[c(1:(DatCol-1))], cluster=fit$cluster), file=paste(OutPre, suffix, opt$pam, 'csv', sep='.'), quote=F, sep="\t", row.names=F, col.names=T)
}


cat("Start: ", date(), "\n")

dat <- loadFile(InFile)
dat.row = nrow(dat)
colnames(dat)[1:3] <- c('chr', 'start', 'end')
if (colnames(dat)[DatCol] == paste('Key', DatCol, sep='')) {
	colnames(dat)[DatCol] <- 'ID'
	DatCol = DatCol+1
}

dat = na.omit(dat)
cat('Removed features with NA: ', length(na.action(dat)), '\n')
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
if (Format == 1) {
	DatCol.new = DatCol
} else if (Format == 2) {
	dat = dat[apply(dat[seq(DatCol+1, length(dat), 2)], 1, function(x) all(x!=0)), ]
	cat('Removed features with 0: ', dat.row-nrow(dat), '\n')
	print(dim(dat))
	dat = cbind(dat[1:(DatCol-1)], round(data.frame(t(apply(dat[DatCol:length(dat)], 1, function(x) x[seq(1,length(x),2)]/x[seq(2,length(x),2)]))), 3))
	DatCol.new = DatCol
} else if (Format > 2) {
	dat.min = apply(dat[icpg], 1, function(x) ifelse(any(is.na(x)), 0, min(x)))
	cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(dat.min < MinCpgNum), '\n')
	dat     = dat    [dat.min >= MinCpgNum,]
	dat.min = dat.min[dat.min >= MinCpgNum ]
	if (colnames(dat)[DatCol-1] == 'ID') dat$ID = factor(dat$ID)
	print(dim(dat))

	# window.methy or window.nread format
	if (Format == 4) {
		dat = cbind(dat[1:(DatCol-1)], dat.min, round(dat[imet]/dat[itot], 3))
	} else {
		dat = cbind(dat[1:(DatCol-1)], dat.min, dat[imet])
	}
	DatCol.new = DatCol + 1
	colnames(dat)[DatCol.new-1] = 'MinCpgNum'
}

print(head(dat))
print(tail(dat))
print(dim(dat))

if (exists('GrpFile')) {
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

	if (!is.null(opt$only)) dat = cbind(dat[1:(DatCol.new-1)], dat[colnames(dat) %in% as.character(grp$sample)])

	grp.colors <- apply(as.matrix(colnames(dat)[DatCol.new:length(dat)]), 1, function (x) ifelse(sum(x==grp$sample) > 0, GrpColors[grp$group[x==grp$sample]], 'white'))
	print(grp.colors)
	stopifnot(length(grp.colors) != 0)
}

#for (i in DatCol.new:length(dat)) {
#	drawHist(na.omit(dat[i]), file=paste(OutPre, 'methy', tolower(colnames(dat)[i]), sep='.'))
#}

if (!is.null(opt$hierarchical)) {
	try(hierarchicalClust(dat, 'hclust.all'))
}

if (!is.null(opt$kmeans)) {
	try(kmeansClust(dat, 'kmeans'))
}

if (!is.null(opt$pam)) {
	try(pamClust(dat, 'pam'))
}

if (!is.null(opt$hierarchical) && exists('MetDiff') && MetDiff > 0) {
	dat$dif = apply(dat[DatCol.new:length(dat)], 1, function(x) {y=range(x); y[2]-y[1]})
	print(range(dat$dif))
	cat('Removed features based on difference between maximum and minimum: ', sum(dat$dif < MetDiff), '\n')
	dat.2 = dat[dat$dif >= MetDiff,]
	cat('Dimension =', dim(dat.2), '\n')

	type = paste('hclust.delta', MetDiff, sep='')
	if (nrow(dat.2) > MaxRow) {
		type = paste('hclust.delta', MetDiff, '.max', MaxRow, sep='')
		dat.2 = dat[order(dat$dif, decreasing=T),][1:min(nrow(dat), MaxRow),]
	}

	if (!is.null(opt$hierarchical)) try(hierarchicalClust(dat.2[1:(length(dat.2)-1)], type))

	if (Top > 0) {
		type = paste('hclust.delta', MetDiff, '.top', Top, sep='')
		dat.2 = dat[order(dat$dif, decreasing=T),][1:min(nrow(dat), Top),]
		if (!is.null(opt$hierarchical)) try(hierarchicalClust(dat.2[1:(length(dat.2)-1)], type))
	}

	dat$dif = NULL
}

if (!is.null(opt$hierarchical) && exists('SdDiff') && SdDiff > 0) {
   dat$sd = apply(dat[DatCol.new:length(dat)], 1, sd)
	print(range(dat$sd))
	cat('Removed features based on standard deviation: ', sum(dat$sd < SdDiff), '\n')
	dat.3 = dat[dat$sd >= SdDiff,]
	cat('Dimension =', dim(dat.3), '\n')

	type = paste('hclust.sd', SdDiff, sep='')
	if (nrow(dat.3) > MaxRow) {
		type = paste('hclust.sd', SdDiff, 'max', MaxRow, sep='')
		dat.3 = dat[order(dat$sd, decreasing=T),][1:min(nrow(dat), MaxRow),]
	}

	if (!is.null(opt$hierarchical)) try(hierarchicalClust(dat.3[1:(length(dat.3)-1)], type))
	if (Top > 0) {
		type = paste('hclust.sd', SdDiff, 'top', Top, sep='')
		dat.3 = dat[order(dat$sd, decreasing=T),][1:min(nrow(dat), Top),]
		if (!is.null(opt$hierarchical)) try(hierarchicalClust(dat.3[1:(length(dat.3)-1)], type))
	}

	dat$sd = NULL
}

cat("End: ", date(), "\n")
