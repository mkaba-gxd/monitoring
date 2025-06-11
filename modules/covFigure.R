################################################################################################# 
#                                                                                               #
#  To draw coverage plots                                                                       #  
#                                                                                               #
#  input data                                                                                   #
#                                                                                               #
#  1. target.txt : transcript informaion                                                        # 
#                                                                                               #
#     chr          start       end         gene    strand   exon     transcript_ID              #
#   ----------------------------------------------------------------------------------          #
#     chr13         32890598    32890664    BRCA2     +       2       NM_000059                 #
#     chr13         32893214    32893462    BRCA2     +       3       NM_000059                 #
#     chr13	    32899213    32899321    BRCA2     +       4       NM_000059                 #
#                                                                                               #    
#                                                                                               #
#  2. xxx.bam                                                                                   #
#                                                                                               #
#                                                                                               #
#################################################################################################


##########################################################
#                                                        #
#  get target regions by combining amplicons             # 
#                                                        #
##########################################################


getTargetOri <- function(target.txt.file){

	target.txt <- read.csv(target.txt.file, sep='\t', stringsAsFactors=FALSE)
	
	target.ori <- with(target.txt, GRanges(seqnames=chr, ranges=IRanges(start, end),strand=strand))

	mcols(target.ori) <- target.txt[,c("gene","exon","transcript_ID")]

	target.ori

}



getTargetBed <- function(target.txt.file){

	
	target.ori <- getTargetOri(target.txt.file)

	target <- reduce(target.ori)

	start(target) <- start(target)-1 # convert to 0-based

	chrOrder <- paste("chr",c((1:22),"X","Y","M"),sep="")

	chr <- factor(as.vector(seqnames(target)),chrOrder,ordered=TRUE)

	target <- target[order(chr,start(target))]


       # add gene information

	hit <- findOverlaps(target, target.ori)

	hit <- hit[!duplicated(queryHits(hit))]

	gene <- mcols(target.ori)$gene[subjectHits(hit)]

	mcols(target)$gene <- gene

#	write.table(as.data.frame(target),'target.bed', sep='\t',col.names=FALSE, row.names=FALSE, quote=FALSE) 

	target

}




##########################################################
#                                                        #
#  create target data in mpileup format                  # 
#                                                        #
##########################################################

makeTargetMpileup <- function(target){


	target_size_each <- end(target) - start(target)
	
	target_size_total <- sum(target_size_each)        # total target size
	


	gene_byBase <- rep(mcols(target)$gene, target_size_each)
	
	chromosome_byBase <- rep(as.vector(seqnames(target)),target_size_each)

	# basePosition <- unlist(sapply(target,function(x){seq(start(x)+1,end(x))}))
	
	target.df <- as.data.frame(target)

	basePosition <- unlist(apply(target.df,1,function(x){seq(as.numeric(x[2])+1, x[3])}))
	

	# create target data in mpileup_format with 4 coulums;chromosome, pos, name, cov
	# chr : target chromosome
	# pos : target base
	# name : target name
	# cov : target coverage 
		

	target.mpileup<- data.frame(matrix(vector(),target_size_total,4,
			   dimnames=list(c(),c("chr","pos","gene","cov"))), stringsAsFactors= FALSE)
		
	target.mpileup$chr <- chromosome_byBase
		
	target.mpileup$pos <- basePosition

	target.mpileup$gene <- gene_byBase

	target.mpileup$cov <- rep(0,target_size_total)             # initial value of Cov is '0'

	rownames(target.mpileup) <- paste(target.mpileup$chr, target.mpileup$pos,sep='.')

	target.mpileup


}


#########################################################################
#                                                                       #
#     get coverage                                                      # 
#                                                                       #
#########################################################################


getCoverage <- function(depth.file,target){

	

#	flag0 <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE, isSecondaryAlignment=FALSE, isPaired=NA)

#	param0 <- ScanBamParam(flag=flag0, which=target)

#	p_param <- PileupParam(max_depth=1000000, min_base_quality=1,distinguish_nucleotides=FALSE,distinguish_strands=FALSE)

#	mpileup <- pileup(bam.file, scanBamParam=param0, pileupParam=p_param)


	mpileup <- read.csv(depth.file, sep='\t', header=F, stringsAsFactors=FALSE)
	colnames(mpileup) <- c('seqnames','pos','count') 

	rownames(mpileup) <- paste(mpileup$seqnames, mpileup$pos, sep=".")	

	target.mpileup <- makeTargetMpileup(target)

	idx <- intersect(rownames(mpileup), rownames(target.mpileup))


	target.mpileup[idx,]$cov <- mpileup[idx,]$count

	target.mpileup

}


getCoverageStat <- function(target.mpileup, output){


	covStat <- tapply(target.mpileup$cov, target.mpileup$gene, function(x){summary(x, na.rm=TRUE, digits=3)})

	covStat <- sapply(covStat,unlist)

	covStat <- t(covStat)

	write.table(covStat,paste0(output,"_coverage_stat.txt"), quote=F)

}


######################################################################
#                                                                    #
#  estimate the coverage of 2N                                       # 
#                                                                    #
######################################################################


getDNA2CopyCoverage<- function(target.mpileup){


	covered.index <- which(target.mpileup$cov>0)
		
	TwoCopy.cov <- median(target.mpileup$cov[covered.index], na.rm=TRUE)


	TwoCopy.cov

}



#####################################################################
#                                                                   #
#  reverse X positions                                              # 
#                                                                   #
#####################################################################


reverseXposition <- function(X.pos){


	X.pos <- max(X.pos) - X.pos+1

	X.pos


}


#####################################################################
#                                                                   #
#  get exon postions of a gene                                      # 
#                                                                   #
#####################################################################
	
getGeneExonPos <- function(target, gene){


	gene.exon.pos <- target[which(mcols(target)$gene == gene)] 

	gene.exon.pos

}


#####################################################################
#                                                                   #
#  get X positions                                                  # 
#                                                                   #
#####################################################################

getXpos <- function(target, gene, interval, strand){



	gene.target <- target[which(mcols(target)$gene==gene)]

	gene.len <- end(gene.target) - start(gene.target)
	
	X.pos <- 1:sum(gene.len)
		
	add.num <- vector()	

	for(i in 1:length(gene.len)){


		add.num <- c(add.num, rep(interval*(i-1), gene.len[i]))  # offset at each base


	}
	
	X.pos <- X.pos + add.num

	index <- vector()
	
	for(i in 1:length(gene.target)){

		index <- c(index, paste(as.vector(seqnames(gene.target))[i],seq(start(gene.target)[i]+1,end(gene.target)[i]),sep="."))

	}

	names(X.pos) <- index

	if(strand=='-'){

		X.pos <- reverseXposition(X.pos)
	}
		
	X.pos

}


#####################################################################
#                                                                   #
#  get Y positions                                                  # 
#                                                                   #
#####################################################################

getYpos <- function(target.mpileup,X.pos) {


	Y.pos <- target.mpileup[names(X.pos),]$cov
	
	names(Y.pos) <- names(X.pos)

	Y.pos

}

#####################################################################
#                                                                   #
#  get X positions of exons of a transcript                         # 
#                                                                   #
#####################################################################

getExonXPos <- function(transcript, target.ori, X.pos){


	transcript.target <- target.ori[which(mcols(target.ori)$transcript_ID==transcript)]
	
	index <- vector()

	exon.pos.start.index <- paste(as.vector(seqnames(transcript.target)),start(transcript.target),sep=".")
	
	exon.pos.end.index <- paste(as.vector(seqnames(transcript.target)),end(transcript.target),sep=".")
		
	exon.X.pos.start <- X.pos[exon.pos.start.index]

	exon.X.pos.end <- X.pos[exon.pos.end.index]

	exon.X.pos <- data.frame(matrix(c(exon.X.pos.start, exon.X.pos.end),length(exon.X.pos.start),2,

                                dimnames=list(c(),c('exon.X.pos.start','exon.X.pos.end'))), stringsAsFactors=FALSE)
		

	exon.X.pos

}


#####################################################################
#                                                                   #
#  get Y positions of exons                                         # 
#                                                                   #
#####################################################################

getExonYPos <- function(i,Y.pos){


	offset <- max(median(Y.pos,na.rm=TRUE)*1/30,6)

	width <- offset*3
	
	gap <- offset*1.5 

	exon.Y.pos <- c((-1)*offset-(i-1)*(width+gap), (-1)*offset-i*width-(i-1)*gap, (-1)*offset-i*width-(i-1)*gap, (-1)*offset-(i-1)*(width+gap))

	exon.Y.pos

}


#####################################################################
#                                                                   #
#  get strand information                                           # 
#                                                                   #
#####################################################################

	

getStrand <- function(target.ori, gene){


	strand <- unique(as.vector(strand(target.ori))[which(mcols(target.ori)$gene==gene)])

	strand

}


#####################################################################
#                                                                   #
#  get transcript information                                       # 
#                                                                   #
#####################################################################


getTranscript <- function(target.ori, gene){


	transcript <- unique(mcols(target.ori)$transcript_ID[which(mcols(target.ori)$gene==gene)])

	transcript

}


#####################################################################
#                                                                   #
#  get range of Y value                                             # 
#                                                                   #
#####################################################################
	
	
getYlim <- function(Y.pos, transcript){
 
	
	offset <- max(median(Y.pos,na.rm=TRUE)*1/30,6)

	width <- offset*3

	gap <- offset

	ylim.lo <- (-1)*offset-(width+gap)*length(transcript)

	ylim.hi <- max(300, ceiling(max(Y.pos,2,na.rm=TRUE)))

	ylim <- c(ylim.lo, ylim.hi)

	ylim

} 

#####################################################################
#                                                                   #
#  get range of X value                                             # 
#                                                                   #
#####################################################################
	

getXlim <- function(X.pos){

	offset <- min(max(X.pos)*1/9, 600)

	xlim.lo <- (-1)*offset

	xlim.hi <- max(X.pos)

	xlim <- c(xlim.lo, xlim.hi)

}


#####################################################################
#                                                                   #
#  get label of Y axis                                              # 
#                                                                   #
#####################################################################
	

getYaxisLabel <- function(ylim){

	
	ylim.hi <- ylim[2]
	
	yaxis.label <- as.character(seq(0, (ylim.hi %/% 100 ))*100)

	yaxis.label

}

	
#####################################################################
#                                                                   #
#  draw exons                                                       # 
#                                                                   #
#####################################################################
	
drawPolygon <- function(exon.X.pos, exon.Y.pos, exon.info ){


	exon.X.pos <- exon.X.pos[order(exon.X.pos[,1]),]

	exon.X.pos.start <- exon.X.pos[,1]
	
	exon.X.pos.end <- exon.X.pos[,2]
	

	
	col<- rep('grey',length(exon.info))

	col[grep("i|p",exon.info)] <- 'white'


	line.hight <- mean(exon.Y.pos)


	polygon(c(min(exon.X.pos),min(exon.X.pos),max(exon.X.pos),max(exon.X.pos)),c(line.hight+1,line.hight-1,line.hight-1,line.hight+1),col='black')


	for(i in 1:length(exon.X.pos.start)){
		

		polygon(c(rep(exon.X.pos.start[i],2),rep(exon.X.pos.end[i],2)), exon.Y.pos, col=col[i])

		text(sum(exon.X.pos.start[i],exon.X.pos.end[i])/2,line.hight, exon.info[i],cex=0.8)
		
	
	}

}

#####################################################################
#                                                                   #
#  get line color by base                                           # 
#                                                                   #
#####################################################################

getColor <- function(color,Y.pos){



	col<- rep(color,length(Y.pos))

	names(col)<-names(Y.pos)

	col


}

#####################################################################
#                                                                   #
#  get index read number                                            # 
#                                                                   #
#####################################################################
	
		
getIndexCov <- function(target.mpileup, ylim){


	twoCopy <- getDNA2CopyCoverage(target.mpileup)


	index.cov <- vector(mode="numeric", length= ylim[2]+3)


	
	index.cov[seq(3,ylim[2]+3)] <- round(twoCopy * 2^seq(0,ylim[2]),0)

	index.cov[2] <- round(twoCopy/2, 0)
	
	index.cov <- formatC(index.cov, format='d', big.mark=',')
		
	index.cov

}	


#####################################################################
#                                                                   #
#  draw coverage figures                                            # 
#                                                                   #
#####################################################################

drawCoverageFigure <- function(target.ori, target, target.mpileup, gene, interval, output){



	gene.mpileup <- target.mpileup[which(target.mpileup$gene==gene),] 

	strand <- getStrand(target.ori,gene)
	
	transcript <- getTranscript(target.ori,gene)


	######################################
	# information of selected gene       #
	######################################
	
	chr <- gene.mpileup$chr[1]  # chromosome of selected gene
	
	gene.pos.start <- formatC(gene.mpileup$pos[1], format='d', big.mark=',')  # start position of the selected gene
	
	gene.pos.end <- formatC(gene.mpileup$pos[dim(gene.mpileup)[1]], format='d', big.mark=',') # end position of the selected gene

	mainText <- paste(gene,"\n", chr,": ",gene.pos.start," - ",gene.pos.end,sep="")  # text for main in figure


	X.pos <- getXpos(target, gene, interval, strand)

	Y.pos <- getYpos(target.mpileup, X.pos)
	


	ylim <- getYlim(Y.pos, transcript)
	
	xlim <- getXlim(X.pos)

	
	yaxis.label <- getYaxisLabel(ylim)

	col <- getColor('blue', Y.pos)
  

	######################################
	# draw CNV figure                    # 
	######################################



	pdf(paste(output,"_dnacopy_",gene,".pdf",sep=''), w=16, h=6)

	par(mar=c(2,4,4,2))

	plot(x=X.pos,y=Y.pos,type="h", col=col, ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", main=mainText, ylab='',xlab='')
	
	yaxis.pos <- seq(0,ylim[2] %/% 100)*100	

	axis(side=2, at=yaxis.pos, labels=yaxis.label, cex.axis=1.2)
	

	mtext("coverage (bp)", side=2, line=2.5, at=mean(yaxis.pos), cex=1.4, adj=0.5,srt=90)

	abline(h= yaxis.pos[-1], lty=2, col='grey')

	abline(h=100, lty=2, col="red")


	#####################################
	# draw exons                        #
	#####################################
	

	for(i in 1:length(transcript)){


		exon.X.pos <- getExonXPos(transcript[i], target.ori, X.pos)

		exon.Y.pos <- getExonYPos(i,Y.pos)

		exon.info <- mcols(target.ori)$exon[which(mcols(target.ori)$transcript_ID==transcript[i])]


		if(strand=='-'){

			exon.info <- rev(exon.info)
		}
			
	
		drawPolygon(exon.X.pos, exon.Y.pos, exon.info)

		text(xlim[1], mean(exon.Y.pos), transcript[i], pos=4, offset= -1, cex= 0.9)
	

	}

	dev.off()

}



getBoxplotYMax <- function(copyNumber){

	
	Y.max <- ceiling(max(sapply(copyNumber,function(x){quantile(x,3/4) + 1.5*IQR(x)})))

	Y.max	

}

################################################################
#                                                              #
#  draw CNV boxplot                                            # 
#                                                              #
################################################################


drawDNACopyNumberBoxplot <- function(target.mpileup, geneList, output){

	
	copyNumber <-list()

	twoCopy <- getDNA2CopyCoverage(target.mpileup)


	for(i in 1:length(geneList)){
	
		copyNumber[[geneList[i]]] <- 2*(target.mpileup$cov[which(target.mpileup$gene==geneList[i])]/twoCopy)

	}


	gene.count <- length(geneList)

	gene.median <- sapply(copyNumber,function(x){round(median(unlist(x),na.rm=TRUE),1)})

	box <- boxplot(copyNumber)

	Y.max <- max(ceiling(max(box$stat[5,])),10)


	options(bitmapType="cairo")	

	fig.width <- 20 * gene.count

	png(paste(output,'_CNV_boxplot.png',sep=''),width=fig.width, height=450)

	boxplot(copyNumber, col='white', outline=F, boxlwd=0.7,ylim=c(0,Y.max), yaxt='n',xaxt='n')

	abline(h= seq(0,Y.max, by=2),lty=5,col='grey')

	abline(h=2, lty=2, col='red')


	x.lab <- as.character(seq(0,Y.max,by=2))
	
	
	boxplot(copyNumber, notch=F, col='#FDAE61',borider='#9E0142', outline=F, boxlwd=0.7, ylab='DNA copy number',las=2, ylim=c(0,16),xaxt='n', yaxt="n",add=T, cex.lab=1.3,main=output)

	axis(side=2, at=c(0, labels=x.lab), cex.axis=1.3)


	text(x=(1:length(geneList)+0.25), y=-log(Y.max,2)*0.3, geneList,xpd=TRUE, srt=45, pos=2, cex=1)
	dev.off()


	write.table(gene.median, paste0(output,"_CNV.txt"), sep='\t',quote=FALSE,col.names=FALSE)


}	



#####################################################################
#                                                                   #
#  draw all figures                                                 # 
#                                                                   #
#####################################################################

drawCoverageFigureAll <- function(depth.file, target.txt.file, output){

	target.ori <- getTargetOri(target.txt.file)

	target <- getTargetBed(target.txt.file)
	
	geneList <- unique(mcols(target)$gene)

	target.mpileup <- getCoverage(depth.file,target)

#	getCoverageStat(target.mpileup, output)

	interval <- 20	

	for(i in 1:length(geneList)){


		gene <- geneList[i]

		drawCoverageFigure(target.ori, target,target.mpileup, gene, interval, output)

	}


#	drawDNACopyNumberBoxplot(target.mpileup, geneList, output)

}
###################################################
#  R CMD BATCH  \                                 #
#  --no-save \                                    #
#  --slave \                                      #
#  --no-restore \                                 #
#  '--args bam.file target.txt output' \          #
#  covFigure.R                                    # 
###################################################

#source('covFigure.R')

suppressMessages({
  suppressPackageStartupMessages({
    library(GenomicAlignments)
  })
})

args <- (commandArgs(TRUE))

bam.file <- args[1]

target.txt.file <- args[2]

output <- args[3]

drawCoverageFigureAll(bam.file, target.txt.file, output)


