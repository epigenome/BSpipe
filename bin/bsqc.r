suppressMessages(library(getopt))

if (!exists('pars')) pars = commandArgs(trailingOnly = T)
cat('Arguments:', pars, '\n')

if (!exists('InFile')) {
	if (length(pars) < 2) {
		message("USAGE: Rscript bsqc in_file out_prefix [option]\n")
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
	Height = ifelse(!is.null(opt$height), opt$height, 480)
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

drawLine <- function(var, cov=NULL, file=NA, width=Width, height=Height, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

	par(xpd=T, mar=par()$mar+c(0,0.5,0,4.5))
	ind <- 2:length(var)
	last <- max(var[1])
	colors = if (!is.null(cov)) c('blue', 'green4', 'red') else 1:length(ind)

	if (is.null(cov))
	{
		matplot(var[1], var[ind], type='l', lty=3, lwd=3, col=colors, main=Title, xlab='Base position', ylab='Methylation', xlim=c(1,last), ylim=c(0,1), cex.lab=1.1, cex.axis=1.1)
	}
	else
	{
		barplot(t(cov), col=c('white', 'gray90'), border='lightgray', xaxt='n', yaxt='n')
		par(new=TRUE)
		matplot(var[1], var[ind], type='l', lty=3, lwd=3, col=colors, main=Title, xlab='Base position', ylab='Methylation', xlim=c(1,last), ylim=c(0,1), cex.lab=1.1, cex.axis=1.1)
	}

	legend(last+5, 1, names(var[ind]), cex=1.1, lty=3, lwd=3, col=colors)
	#legend(Key, names(var[ind]), cex=1.1, lty=3, lwd=3, col=colors)

	if (!is.na(file)) closeGraphics()
}

drawBox <- function(var, file=NA, width=Width, height=Height, ...)
{
	if (!is.na(file)) openGraphics(file, width, height)

#	par(xpd=T, mar=par()$mar+c(0,0,0,4))
	ind <- 2:length(var)
	last <- max(var[1])
	colors=1:length(ind)
	barplot(var[[2]], names.arg=var[[1]], col=colors[1], density=0)
	for (i in 2:3) barplot(var[[i]], names.arg=var[[1]], col=colors[i], density=0, add=T)

	if (!is.na(file)) closeGraphics()
}

get <- function(dat, id, drop)
{
	var <- reshape(dat, idvar=id, timevar='lab', drop=c('V1', 'V2', drop), direction='wide')
	colnames(var) = c('Position', sub('[^\\.]+\\.', '', colnames(var)[2:length(var)], perl=T))
	if (length(grep('^CN', colnames(var))) > 0) var <- var[-grep('^CN', colnames(var))]
	print(head(var))
	return(var)
}

cat("Start: ", date(), "\n")

PdfFlag = F

dat <- loadFile(InFile)
dat$lab <- apply(cbind(as.character(dat$V1), as.character(dat$V2)), 1, function(x) paste(x[1], substr(x[2], 1, 1), sep=''))

met  <- get(dat, 'V3', c('V4', 'V5'))
mcov <- get(dat, 'V3', c('V5', 'V6'))
ucov <- get(dat, 'V3', c('V4', 'V6'))

for (i in seq(2, ncol(met), 3)) drawLine(met[c(1,i:(i+2))], cov=cbind(mcov[i]+ucov[i], mcov[i+1]+ucov[i+1]), file=paste(OutPre, sub('[frb]$', '', colnames(met[i])), sep='.'))

ind = grep('C[ACT]f', colnames(mcov))
mcov.ch = data.frame(CHf=rowSums(mcov[ind]), CHr=rowSums(mcov[ind+1]), CHb=rowSums(mcov[ind+2]))
ucov.ch = data.frame(CHf=rowSums(ucov[ind]), CHr=rowSums(ucov[ind+1]), CHb=rowSums(ucov[ind+2]))
met.ch = data.frame(met[1], CHf=mcov.ch[1]/(mcov.ch[1]+ucov.ch[1]), CHr=mcov.ch[2]/(mcov.ch[2]+ucov.ch[2]), CHb=mcov.ch[3]/(mcov.ch[3]+ucov.ch[3]))
drawLine(met.ch, cov=cbind(mcov.ch[1]+ucov.ch[1], mcov.ch[2]+ucov.ch[2]), file=paste(OutPre, 'CH', sep='.'))
head(cbind(met.ch, mcov.ch, ucov.ch))

all = cbind(met[1], met['CGb'], met.ch['CHb'], met['Cb'])
colnames(all) = sub('b$', '', colnames(all))
head(all)
drawLine(all, file=paste(OutPre, 'all', sep='.'))

cat("End: ", date(), "\n")
