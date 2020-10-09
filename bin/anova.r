FdrCutoff <- 0.05
PvalCutoff <- 0.01
MinCpgNum <- 5
MetDiff <- 0.25
MaxRowForHeatmap <- 10000
Top <- 500
PdfFlag = 0
DatCol = 4
GrpColors <- c('green3', 'gold', 'purple', 'orange', 'magenta', 'cyan', 'blue3', 'red', 'pink', 'black')


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre') || !exists('GrpFile')) {
	if (length(pars) < 3) {
		message("USAGE: Rscript anova.r in_file out_prefix group_file [option]\n")
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
		'filter',   'F', 0, 'logical',
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
	DistMethod = ifelse(is.null(opt$dist), 'euclidean', opt$dist)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	Filter = ifelse(is.null(opt$filter), F, T)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}


suppressMessages(library(geneplotter))
suppressMessages(library(gplots))
suppressMessages(library(aod))
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

	hv <- heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, col=col, scale="none", key=TRUE, symkey=FALSE, density.info="density", trace="none", labRow=NA, cexCol=1.4, cexRow=1, margin=c(7,2), lwid=c(wid,10-wid), lhei=c(hei,10-hei), ...)

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

anova <- function(var)
{
	x <- stack(var)
	z <- merge(grp, x, by.x='sample', by.y='ind')
	ano <- aov(values~group, data=z)
	sum <- summary(ano)[[1]]
	ftest <- sum[1, 'F value']

	if (is.null(ftest)) {
		ftest <- NA
		pvalue <- NA
		hsd <- NA
	} else {
		pvalue <- sum[1, 'Pr(>F)']
		hsd <- TukeyHSD(ano)$group
	}

	ret <- c(ftest, pvalue, t(hsd))
	names(ret) <- c('Ftest', 'pvalue', apply(as.matrix(sub('\\.', '-', rownames(hsd))), 1, paste, gsub(' ', '_', colnames(hsd)), sep='.'))
	return(c(tapply(z$values, z$group, mean), ret))
}

analysis <- function(info, data, suffix, grp.col, save=F)
{
	if (nrow(data) <= 1) return(invisible())

	if (save) write.table(cbind(info[ , c(1:(DatCol-1))], data, info[ , DatCol:length(info)]), file=paste(OutPre, suffix, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)

	if (nrow(data) > MaxRowForHeatmap) {
		ind = grep('pvalue', colnames(info))
		if (length(ind) == 0) ind = grep('.p_adj', colnames(info))
		stopifnot(length(ind) > 0)
		data = head(data[order(info[ind]),], MaxRowForHeatmap)
		info = head(info[order(info[ind]),], MaxRowForHeatmap)
		suffix = paste(suffix, '.max', MaxRowForHeatmap, sep='')
	}

	hv <- drawHeatmap(data, Rowv=T, Colv=T, ColSideColors=grp.col, file=paste(OutPre, suffix, 'heatmap', sep='.'))
	write.table(cbind(info[rev(hv$rowInd), ], data[rev(hv$rowInd), hv$colInd]), file=paste(OutPre, suffix, 'heatmap.csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)
}

suboutput <-function(pct, stat, field, cutoff, ext, suf)
{
	cat('Field: ', field, '\n')
	log = stat[[field]] < cutoff
	if (sum(log) == 0) return(invisible())

	ext3 = paste(ext, paste(cutoff, MetDiff, sep='-'), suf, sep='.')
	idiff = grep('\\.diff', colnames(stat))
	stopifnot(length(idiff) > 0)
	log2 = apply(stat, 1, function(x) any(x[[idiff]] > MetDiff))
	try(analysis(cbind(pct[log & log2, 1:(DatCol-1)], stat[log & log2,]), pct[log & log2, DatCol.new:length(pct)], ext3, grp.colors, save=T))

	pairs = t(combn(levels(grp$group), 2))
	for (p in 1:nrow(pairs)) {
		Grp1 = as.character(pairs[p,1])
		Grp2 = as.character(pairs[p,2])
		cat(Grp1, ' vs ', Grp2, '\n')

		pair = paste(Grp1, Grp2, sep='-')
		idx = grep(pair, colnames(stat))
		if (length(idx) == 0) {
			pair = paste(Grp2, Grp1, sep='-')
			idx = grep(pair, colnames(stat))
			if (length(idx) == 0) {
				message('Undefined pair: ', Grp1, ' ', Grp2)
				q()
			}
		}

		idx.q = which(colnames(stat) == paste(pair, 'p_adj', sep='.'))
		idx.d = which(colnames(stat) == paste(pair, 'diff', sep='.'))
		log2 = stat[[idx.q]] < cutoff & abs(stat[[idx.d]]) > MetDiff
		pair.idx = which(colnames(pct) %in% as.character(grp$sample[Grp1 == grp$group | Grp2 == grp$group]))
		ext4 = paste(ext, tolower(pair), paste(cutoff, MetDiff, sep='-'), suf, sep='.')
		try(analysis(cbind(pct[log & log2, 1:(DatCol-1)], stat[log & log2, idx]), pct[log & log2, pair.idx], ext4, grp.colors[pair.idx-DatCol.new+1], save=T))
		try(drawVolcano(data.frame(m=stat[[idx.d]], q=stat[[idx.q]]), pair, cutoff, path=paste(OutPre, ext4, 'vol', sep='.')))

		if (sum(log & log2) > Top*2)
		{
			log.3 = head(order(abs(stat[log & log2, idx.d]), decreasing=T), Top)
			try(analysis(cbind(pct[log & log2,][log.3, 1:(DatCol-1)], stat[log & log2,][log.3, idx]), pct[log & log2,][log.3, pair.idx], paste(ext4, paste('top', Top, sep=''), sep='.'), grp.colors[pair.idx-DatCol.new+1]))
		}
	}
}

output <-function(pct, stat, ext)
{
	stat$pvalue[is.na(stat$pvalue)] = 1
	idx = grep('pvalue', colnames(stat))
	stopifnot(length(idx) == 1)
	stat = data.frame(stat[1:idx[1]], padjust=p.adjust(stat$pvalue, 'fdr'), stat[(idx[1]+1):length(stat)], check.names=F)
	write.table(cbind(pct[1:(DatCol.new-1)], stat), file=paste(OutPre, ext, 'csv', sep='.'), quote = F, sep = "\t", row.names = F, col.names=T)

	draw(hist, stat$pvalue , breaks=seq( 0,1,0.01), xlab='P-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'pval.hist', sep='.'))
	draw(hist, stat$padjust, breaks=seq( 0,1,0.01), xlab='Q-value', ylab='Number of features', xlim=c( 0, 1), file=paste(OutPre, ext, 'qval.hist', sep='.'))

	suboutput(pct, stat, 'pvalue', PvalCutoff, ext, 'pval')
	suboutput(pct, stat, 'padjust', FdrCutoff, ext, 'qval')
}

stat.test <- function(x, fun)
{
	if (is.null(opt$thread) || opt$thread == 1) {
		return(data.frame(t(apply(x, 1, fun)), check.names=F))
	} else {
		return(data.frame(do.call(rbind, mclapply(split(x, 1:nrow(x)), function(y) fun(unlist(y)), mc.cores=opt$thread)), check.names=F))
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

org <- loadFile(InFile)
colnames(org)[1:3] <- c('chr', 'start', 'end')
if (colnames(org)[DatCol] == paste('Key', DatCol, sep='')) {
	colnames(org)[DatCol] <- 'ID'
	DatCol = DatCol+1
}

org = na.omit(org)
cat('Removed features with NA: ', length(na.action(org)), '\n')
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
	org = cbind(org[1:(DatCol-1)], data.frame(t(apply(org[DatCol:length(org)], 1, function(x) x[seq(1,length(x),2)]/x[seq(2,length(x),2)]))))
	org = org[complete.cases(org), ]
}

type = 'anova'
print(type)

if (Format == 1 || Format == 2) {
	pure = org
	DatCol.new = DatCol
} else {
	org.min = apply(org[imet-1], 1, min)
	cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(org.min < MinCpgNum), '\n')
	org     = org    [org.min >= MinCpgNum,]
	org.min = org.min[org.min >= MinCpgNum ]
	if (colnames(org)[DatCol-1] == 'ID') org$ID = factor(org$ID)
	print(dim(org))

	# window.methy or window.nread format
	if (Format == 4) {
		pure = cbind(org[1:(DatCol-1)], org.min, data.frame(t(apply(org[DatCol:length(org)], 1, function(x) round(x[seq(2,length(x),3)]/x[seq(3,length(x),3)], 3)))))
	} else {
		pure = cbind(org[1:(DatCol-1)], org.min, org[imet])
	}
	DatCol.new = DatCol + 1
	colnames(pure)[DatCol:length(pure)] = c('MinCpgNum', colnames(org)[imet])
}

print(head(pure))
print(tail(pure))
print(dim(pure))

dat.grp <- apply(as.matrix(colnames(pure[DatCol.new:length(pure)])), 1, function(x) grp$group[x==grp$sample])
print(dat.grp)
stopifnot(length(dat.grp) != 0)
grp.colors <- GrpColors[dat.grp]

com <- stat.test(pure[DatCol.new:length(pure)], anova)
output(pure, com, type)

cat("End: ", date(), "\n")
