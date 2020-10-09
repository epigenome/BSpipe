MinCpgNum <- 5
Width <- 700
Height <- 700
PdfFlag = 0
DatCol = 4
GrpColors = c('green3', 'blue3', 'magenta', 'brown', 'black', 'purple', 'pink', 'darkcyan', 'red', 'cyan', 'gold', 'green', 'darkred')


suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile') || !exists('OutPre')) {
	if (length(pars) < 2) {
		message("USAGE: Rscript correlation.r in_file out_prefix [option]\n")
		quit(save='no', status=-1)
	}

	InFile = pars[1]
	OutPre = sub('\\.gz$', '', pars[2])
	if (length(pars) > 3) opt = getopt(matrix(ncol=4, byrow=T, c(
		'format',   'b', 1, 'character',
		'pearson' , 'p', 0, 'logical',
		'spearman', 's', 0, 'logical',
		'group'   ,	'R', 1, 'character',
		'only'    ,	'O', 0, 'logical',
		'note'    , 'n', 0, 'logical',
		'linkage' , 'L', 1, 'character',
		'cpg'     , 'c', 1, 'integer',
		'width'   , 'W', 1, 'integer',
		'height'  , 'H', 1, 'integer',
		'pdf'     , 'Y', 0, 'logical'
		)), pars[3:length(pars)])
	else opt = list()

	Format = ifelse(is.null(opt$format) || opt$format == 'auto', 0, ifelse(opt$format == 'sitemethy', 1, ifelse(opt$format == 'sitenread', 2, ifelse(opt$format == 'winmethy', 3, 4))))
	if (!is.null(opt$cpg   )) MinCpgNum = as.integer(opt$cpg)
	if (is.null(opt$pearson)) {
		if (is.null(opt$spearman)) {
			opt$pearson = T
			opt$spearman = T
		} else {
			opt$pearson = F
		}
	} else {
		if (is.null(opt$spearman)) {
			opt$spearman = F
		}
	}
	CellNote = ifelse(is.null(opt$note), F, T)
	Linkage = ifelse(is.null(opt$linkage), 'complete', opt$linkage)
	GrpFile = ifelse(!is.null(opt$group ), opt$group, NA)
	if (!is.null(opt$width )) Width   = as.numeric(opt$width)
	if (!is.null(opt$height)) Height  = as.numeric(opt$height)
	PdfFlag = !is.null(opt$pdf) && opt$pdf
}


suppressMessages(library(geneplotter))
suppressMessages(library(gplots))


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
	dat <- read.table(file, sep="\t", header=header, check.names=F, comment.char=comment.char, ...)
	print(head(dat))
	print(tail(dat))
	cat("dim=", dim(dat), "\n")
	return(dat)
}

saveFile <- function(dat, file, rname=F, cname=T)
{
	write.table(dat, quote=F, sep="\t", row.names=rname, col.names=cname, file=file)
}

drawHist <- function(var, file=NA, ...)
{
	if (!is.na(file)) openGraphics(file, 600, 400)
#	r <- range(as.matrix(var))
#	h <- hist(var, breaks=c((r[1]-0.5):(r[2]+0.5)), ...)
	h <- hist(var, ...)
	print(h$counts)
	if (!is.na(file)) closeGraphics()
}

drawHeatmap <- function(x, dist.method, file=NA, ...)
{
	if (!is.na(file)) openGraphics(file, Width, Height)
	nmargin <- 2 + max(nchar(colnames(x)))/1.5

	hclust2 <- function(x, method=Linkage, ...) hclust(x, method=method, ...)
	dist2 <- function(x, ...) as.dist(1-cor(t(x), method=dist.method, use='pairwise'))
	if (CellNote) cellnote <- x
	else          cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	if (exists('grp.colors')) heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, scale="none", symm=TRUE, key=TRUE, symkey=FALSE, density.info="density", trace="none", Colv=T, cexRow=1.4, cexCol=1.4, cellnote=cellnote, notecol='black', notecex=1.2, margins=c(nmargin,nmargin), ColSideColors=grp.colors, ...)
	else                     heatmap.2(as.matrix(x), distfun=dist2, hclustfun=hclust2, scale="none", symm=TRUE, key=TRUE, symkey=FALSE, density.info="density", trace="none", Colv=T, cexRow=1.4, cexCol=1.4, cellnote=cellnote, notecol='black', notecex=1.2, margins=c(nmargin,nmargin), ...)

	if (!is.na(file)) closeGraphics()
}

draw <- function(fun, x, file=NA, width=600, height=500, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	fun(x, ...)

	if (!is.na(file)) closeGraphics()
}

myImagePlot <- function(x, wid=5, ...){
	min <- min(x)
	max <- max(x)
	yLabels <- rownames(x)
	xLabels <- colnames(x)
	title <-c()
	revFlag <- 0
	ColorRamp <- NULL
	xlab <- NULL
	ylab <- NULL
	las <- NULL

	# check for additional function arguments
	if( length(list(...)) ){
		Lst <- list(...)
		if( !is.null(Lst$min) ){
			min <- Lst$min
		}
		if( !is.null(Lst$max) ){
			max <- Lst$max
		}
		if( !is.null(Lst$zlim) ){
			min <- Lst$zlim[1]
			max <- Lst$zlim[2]
		}
		if( !is.null(Lst$yLabels) ){
			yLabels <- c(Lst$yLabels)
		}
		if( !is.null(Lst$xLabels) ){
			xLabels <- c(Lst$xLabels)
		}
		if( !is.null(Lst$title) ){
			title <- Lst$title
		}
		if( !is.null(Lst$xlab) ){
			xlab <- Lst$xlab
		} else {
			xlab <- ''
		}
		if( !is.null(Lst$ylab) ){
			ylab <- Lst$ylab
		} else {
			ylab <- ''
		}
		if( !is.null(Lst$rev) ){
			revFlag <- Lst$rev
		}
		if( !is.null(col) ){
			ColorRamp <- Lst$col
		}
		if( !is.null(Lst$colmap) ){
			ColorRamp <- Lst$colmap(max-min+1)
		}
		if( !is.null(Lst$las) ){
			las <- Lst$las
		}
	}
	# check for null values
	if( is.null(xLabels) ){
		xLabels <- c(1:ncol(x))
	}
	if( is.null(yLabels) ){
		yLabels <- c(1:nrow(x))
	}

	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(wid,1), heights=c(1,1))

	if ( is.null(ColorRamp) )
	{
		# Red and green range from 0 to 1 while Blue ranges from 1 to 0
		ColorRamp <- rgb( seq(0,1,length=256),  # Red
                      seq(0,1,length=256),  # Green
                      seq(1,0,length=256))  # Blue
 	}
	ColorLevels <- seq(min, max, length=length(ColorRamp))

 	# Reverse Y axis
 	if ( revFlag ) {
		reverse <- nrow(x) : 1
		yLabels <- yLabels[reverse]
		x <- x[reverse,]
	 }

 	# Data Map
 	xlen = max(nchar(xLabels))
 	ylen = max(nchar(yLabels))
 	if ( is.null(xlab) ) {
		mar = c(4,5+xlen/4,2.5,1+ylen/4)
	}
	else {
		mar = c(6,5+xlen/4,2.5,1+ylen/4)
	}
	#print(mar)
	par(mar=mar)
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab=xlab, ylab=ylab, axes=FALSE, zlim=c(min,max))
 	if( !is.null(title) ){
		title(main=title)
 	}
	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.0, tick=F, las=las)
	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1, tick=F, cex.axis=1.0)

 	# Color Scale
	 if ( is.null(xlab) ) {
		par(mar = c(4,1,2.5,3))
	 }
 	else {
		par(mar = c(6,1,2.5,3))
 	}
	image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="", xaxt="n", yaxt="n")
	axis(4, cex.axis=0.9)

	layout(1)
}

plotPair <- function(logs, title, method, ...)
{
	panel.hist <- function(x, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5) )
		h <- hist(x, plot=F)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$counts; y <- y/max(y)
		rect(breaks[-nB], 0, breaks[-1], y, col="gold", ...)
	}

	panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		c <- cor(x, y, method=method, use='pairwise')
		r <- abs(c)
		txt <- format(c(c, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
		text(0.5, 0.5, txt, cex = cex * r, ...)
	}

	panel.scatter <- function(x, y, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5), new=T)
		smoothScatter(x, y, nrpoints=0, cex=ifelse(PdfFlag, 0.5, 1), ...)
	}

	pairs(logs, main=title, cex.labels=1.5, upper.panel=panel.scatter, lower.panel=panel.cor, labels=colnames(logs), ...)
}

plotPair2 <- function(logs, mat, type, title, method, ...)
{
	panel.scatter <- function(x, y, ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5), new=T)
		smoothScatter(x, y, nrpoints=0, cex=ifelse(PdfFlag, 0.5, 1.2), ...)
	}

#	nf <- layout(matrix(c(1,1,2,2),2,2,byrow=TRUE), widths=c(3.5,3.5), heights=c(2,5), TRUE)
	openGraphics(paste(OutPre, type, method, 'pair1', sep='.'), Width, Height)
	len = max(nchar(colnames(mat)))
	par(mar=c(0,5+len/5,5+len/5,1))
	image(mat[,ncol(mat):1], axes=F)
	if (CellNote) text((row(mat)-1)/(nrow(mat)-1), rev(col(mat)-1)/(ncol(mat)-1), mat)
	axis(3, at=(col(mat)[1,]-1)/(ncol(mat)-1), cex.axis=1.2, lab=colnames(mat), las=2)
	axis(2, at=(col(mat)[1,]-1)/(ncol(mat)-1), cex.axis=1.2, lab=rev(rownames(mat)), las=2)
	closeGraphics()
	openGraphics(paste(OutPre, type, method, 'pair2', sep='.'), Width, Height)
	par(mar=c(2,6,1,1))
	pairs(logs, main=title, cex.labels=1.5, upper.panel=panel.scatter, lower.panel=panel.scatter, labels=colnames(logs), xaxt='n', yaxt='n', gap=0)
	closeGraphics()
}

correlationInSample <- function(x, type, method)
{
	ret <- matrix(1, nrow=length(x), ncol=length(x), dimnames=list(colnames(x), colnames(x)))

	for (i in 1:(length(x)-1)) {
		for (j in (i+1):length(x)) {
			ret[i,j] <- round(cor(x[i], x[j], method=method, use='pairwise'), 2)
			ret[j,i] <- ret[i,j]
		}
	}

	saveFile(ret, paste(OutPre, type, method, 'csv', sep='.'))
	openGraphics(paste(OutPre, type, method, 'mat', sep='.'), Width, Height)
	myImagePlot(ret, min=0, max=1, las=2, col=greenred(100))
	closeGraphics()
	openGraphics(paste(OutPre, type, method, 'pair', sep='.'), Width, Height)
	plotPair(x, NULL, method, gap=0, xaxt='n', yaxt='n')
	closeGraphics()
	plotPair2(x, ret, type, NULL, method)
	drawHeatmap(ret, method, file=paste(OutPre, type, method, 'heatmap', sep='.'))
}


cat("Start: ", date(), "\n")

dat <- loadFile(InFile)
colnames(dat)[c(1,2,3)] <- c('Ref', 'Start', 'End')
if (colnames(dat)[DatCol] == paste('Key', DatCol, sep='')) {
	colnames(dat)[DatCol] <- 'ID'
	DatCol = DatCol+1
}

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
if (Format == 2) {
	dat = cbind(dat[1:(DatCol-1)], round(data.frame(t(apply(dat[DatCol:length(dat)], 1, function(x) x[seq(1,length(x),2)]/x[seq(2,length(x),2)]))), 3))
} else if (Format > 2) {
	dat.min = apply(dat[icpg], 1, function(x) ifelse(any(is.na(x)), 0, min(x)))
	cat('Removed features with #CpGs < ', MinCpgNum, ': ', sum(dat.min < MinCpgNum), '\n')
	dat     = dat    [dat.min >= MinCpgNum,]
	dat.min = dat.min[dat.min >= MinCpgNum ]
	if (colnames(dat)[DatCol-1] == 'ID') dat$ID = factor(dat$ID)
	print(dim(dat))

	# window.methy or window.nread format
	if (Format == 4) {
		dat = cbind(dat[1:(DatCol-1)], round(dat[imet]/dat[itot], 3))
	} else {
		dat = cbind(dat[1:(DatCol-1)], dat[imet])
	}
}

print(head(dat))
print(tail(dat))
print(dim(dat))

if (!is.na(GrpFile)) {
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

	if (!is.null(opt$only)) dat = cbind(dat[1:(DatCol-1)], dat[colnames(dat) %in% as.character(grp$sample)])

	grp.colors <- apply(as.matrix(colnames(dat)[DatCol:length(dat)]), 1, function (x) ifelse(sum(x==grp$sample) > 0, GrpColors[grp$group[x==grp$sample]], 'white'))
	print(grp.colors)
	stopifnot(length(grp.colors) != 0)
}

flag <- apply(dat, 1, function(x) sum(is.na(x)))
naFeat <- apply(dat[DatCol:length(dat)], 1, function(x) sum(is.na(x)))
drawHist(naFeat, freq=F, main=NULL, xlab='Number of samples with no reads', ylab='Number of CpGs', file=paste(OutPre, 'nafeat', sep='.'))

naSam <- apply(dat[DatCol:length(dat)], 2, function(x) sum(is.na(x)))
draw(barplot, naSam, ylab='Number of CpGs with no reads', xlab='Samples', las=2, file=paste(OutPre, 'nasam', sep='.'))

sub <- dat[DatCol:length(dat)]
#union <- sub
#union[is.na(union)] <- 0
common <- na.omit(sub)

#correlationInSample(union, 'union', 'pearson')
#correlationInSample(union, 'union', 'spearman')
if (opt$pearson) correlationInSample(common, 'common', 'pearson')
if (opt$spearman) correlationInSample(common, 'common', 'spearman')

cat("End: ", date(), "\n")
