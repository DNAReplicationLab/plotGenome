##  loadBed - loads any bed file (according to the BED file format defined by UCSC: http://genome.ucsc.edu/FAQ/FAQformat).
##  First three columns are mandatory.
##  If supplied with a 'genome' dataframe (which must contain "chrom","chromStart","chromEnd" columns), the data will be cleansed and reordered accordingly.
##  If genome dataframe lacks chrM, then the function will purge any chrM data read from the bed file.
##  The chromosome order will be specified by the order in the genome dataframe "chrom" column.
##  Optional 'fname' string will replace the "name" column content in the resulting bed dataframe.
#~~~~~~~~~~~~~~ example usage:(relative path, unspecified genome):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- loadBed("file.bed")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Run the following to create sacCer3 genome dataframe
# sacCer3<-data.frame(chrom=as.character(paste0("chr",as.roman(1:16))),chromStart=rep.int(0,16), chromEnd=as.integer(c(230208,813136,316613,1531914,576869,270148,1090944,562639,439885,745446,666445,1078173,924430,784328,1091285,948060)), cen=as.integer(c(151524,238265,114443,449766,152046,148569,496979,105645,355687,436366,440188,150888,268090,628817,326643,556015)))
# sacCer3$chrom <- factor(sacCer3$chrom, levels = unique(sacCer3$chrom))
#~~~~~~~~~~~~~~ example usage: (absolute path, specified genome dataframe):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- loadBed("/path/to/bed/file.bed",sacCer3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loadBed <- function(file,genome=NULL,fname=NULL) {
	bed<-read.table(file,sep="",stringsAsFactors=F)
	bed.col.n<-as.integer(dim(bed)[2])
	if (bed.col.n<3) stop("Bed file must contain at least 3 columns!")
	bed.cols<-c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
	colnames(bed)<-bed.cols[1:bed.col.n]
		## check if loaded bed has column names and drop them
	if (is.na(suppressWarnings(as.numeric(bed$chromStart[1])))) {
		bed <- bed[-1,]
	}
		## replace name column if supplied
	if (!is.null(fname)) {
		bed$name<-paste(fname)
	}
	bed$name<-factor(bed$name)
	#attr(bed,"source") <- paste(fname)
		## explicitly state data types
	bed$chromStart <- as.integer(bed$chromStart)
	bed$chromEnd <- as.integer(bed$chromEnd)
	bed$chrom<-factor(bed$chrom,levels=unique(bed$chrom))
		## if genome dataframe is provided, only use data for chromosomes in it
	if (!is.null(genome)) {
		if (!all(bed.cols[1:3] %in% colnames(genome))) { stop('Genome dataframe must have "chrom", "chromStart", and "chromEnd" columns! Aborting.') }
		junk<-setdiff(bed$chrom,genome$chrom)
		for (bad.chr in junk) {
			trash<-which(bed$chrom == bad.chr)
			if (length(trash) > 0) {
				bed<-bed[-trash,]
				message("\n\tExcluded data for ",bad.chr,".\n")
			}
		}
		bed$chrom<-factor(bed$chrom,levels=unique(genome$chrom))
	}
		## order the resulting bed dataframe
	bed<-bed[order(bed$chrom,bed$chromStart),]
	rownames(bed) <- 1:nrow(bed)
	return(bed)
}





##  gplotBed - plots boxplots of values in the "score" column, per chromosome. It highlights outliers based on 
##  interquartile range (IQR), i.e. datapoints below first quartile - 3xIQR and above third quartile + 3xIQR
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		gplotBed(bedFileContent)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  if plotting is set to FALSE, it returns the plot object which can be manipulated further or plotted by calling it
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		plotObject <- gplotBed(bedFileContent,plotting=F)
#~ 		plotObject  ## (to plot using X11 - requires "-X" argument when connected via ssh)
##  To save as a pdf, run the following 3 lines of code (pdf will be saved in the current working directory):
#~		cairo_pdf("plot.pdf"),family="Arial",width=8,height=11)
#~		plotObject
#~		dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gplotBed <- function(bed,plotting=T) {
	require(ggplot2)
	theTheme<-theme(
		text = element_text(family="Arial"),
		plot.title = element_text(size=rel(2),hjust=0.5,face="bold"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0.75),
		panel.background = element_blank(),
		axis.line = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(size = rel(1.8)),
		axis.text = element_text(size = rel(1.5)),
		axis.text.x = element_text(angle=45,hjust=1),
		legend.position='none',
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	nChr <- length(unique(bed$chrom))
	scoreStats <- summary(bed$score)
	maxScore <- scoreStats[6]
	minScore <- scoreStats[1]
	iqr <- as.numeric(scoreStats[5]-scoreStats[2])
	top3 <- as.numeric(scoreStats[5]+3*iqr)
	low3 <- scoreStats[2]-3*iqr
    boxplot <- ggplot(bed,aes(chrom,score))   
	if (as.numeric(scoreStats[6]) > as.numeric(scoreStats[5]+3*iqr)) {
		boxplot <- boxplot + geom_rect(
			xmin=1-nChr*0.03,
			ymin=top3,
			xmax=nChr+nChr*0.03,
			ymax=maxScore+(maxScore-minScore)*0.02,
			fill="gray95",
			size=0
		)
	}
	if (minScore < low3) {
		boxplot <- boxplot + geom_rect(
			xmin=1-nChr*0.03,
			ymin=minScore-(maxScore-minScore)*0.02,
			xmax=nChr+nChr*0.03,
			ymax=low3,
			fill="gray95",
			size=0
		)
	}
	boxplot <- boxplot + geom_boxplot()
	boxplot <- boxplot + geom_hline(yintercept=as.numeric(scoreStats[3]),size=.75,color="orchid")  ##  Median line
	boxplot <- boxplot + coord_cartesian(xlim=c(1,nChr),ylim=c(minScore,maxScore)) 
	boxplot <- boxplot + ggtitle("Raw reads distribution")
	boxplot <- boxplot + labs(y="Score value")
	boxplot <- boxplot + theTheme
	if (plotting) { boxplot } else { return(boxplot) }
}




##  plotCoverage - plots values in the "score" column for the whole genome
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		plotCoverage(bed,plotting=F)  ## returns plot object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotCoverage <- function(bed,region=F,plotting=T) {
	# sanity check
#	bedName <- deparse(substitute(bed))
#	if (!exists(bedName)) {
#		cat("\n\tCannot find ",bedName," dataframe. Make sure it is spelled correctly\n\n",sep="")
#	stop("")
#	}

	##  load ggplot2 library and set plotting theme
	require(ggplot2)
	fontSize <- 12
	pointSize <- 1 / 0.352777778
	gplotTheme<-theme(
		text = element_text(family="Arial",size=fontSize),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.title = element_text(hjust = 0.5,size=rel(1.5)),
		axis.line = element_line(colour="gray50",size=.5),
		axis.title = element_text(size = rel(1.5),face="bold"),
		axis.text.x = element_text(size=rel(1.4)),
		legend.position="none",
		strip.background = element_blank(),
		strip.text = element_blank(),
		plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
	)

	theMedian <- as.numeric(summary(bed$score)["Median"])
	ylims <- c(min(bed$score),max(bed$score))
	##  Create genome dataframe using a helper function
	if (region!=F) {
		plottingRegion <- T
		region <- as.character(region)
		genome <- makeGenome(bed,region)
		if (length(grep(":",region))==1) {
			region <- unlist(strsplit(region,":"))
			tmp <- unlist(strsplit(region[2],"-"))
			region <- c(region[1],tmp)
			bed<-bed[bed$chrom==region[1] & bed$chromStart>=as.integer(region[2]) & bed$chromEnd<=as.integer(region[3]),]
			#bed<-bed[bed$chromStart>=as.integer(region[2]),]
			#bed<-bed[bed$chromEnd<=as.integer(region[3]),]
		} else {
			bed<-bed[bed$chrom==region]
		}
	} else {
		genome <- makeGenome(bed)
		plottingRegion <- F
	}
	xAxisName <- attributes(genome)$axisName

	##  Make nice x-axis labels using a helper function
	theMax <- max(genome$chromEnd)
	theMin <- min(genome$chromStart)
	labels <- makeLabels(theMin,theMax,"b")
#~ 	minorTickFactor <- as.numeric(attributes(labels)$minorTickFactor)

		## Make a dataframe with the ticks for the chromosome lines
	if (class(region)=="logical") {
		label <- labels$ticks[2]-labels$ticks[1]
		ticks <- numeric()
		chrom <- character()
		y <- numeric() 
		for (chr in levels(genome$chrom)) {
			row <- which(genome$chrom == chr)
			tmpTicks <- seq(genome[row,2],genome[row,3],by=label)
			tmpChr <- rep(chr,length(tmpTicks))
			tmpY <- rep(max(bed[bed$chrom==chr,"score"])/12,length(tmpTicks))
			ticks <- append(ticks,tmpTicks)
			chrom <- append(chrom,tmpChr)
			y <- append(y,tmpY)
		}
		chrTicks <- data.frame(chrom=factor(chrom,levels=unique(genome$chrom)),ticks=as.numeric(ticks),y=as.numeric(y))
	}

	##  plotting
	plot<-ggplot(bed)
	plot<-plot+scale_y_continuous(
		name="Reads per bin",
		limits=c(0,NA)
	)  ##  y scale
	plot<-plot+scale_x_continuous(
		breaks=labels$ticks,
		labels=labels$labels,
		name=xAxisName,
		limits=c(theMin-theMax/200,theMax+theMax/100),expand=c(0.01,0.01)
	)  ##  x scale
	plot<-plot+geom_segment(
		aes(x=chromStart,xend=chromEnd,y=theMedian,yend=theMedian),
		data=genome,
		size=.75,
		color="orchid"
	)  ##  Median line
	if (!plottingRegion) {
		plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=0,yend=0),data=genome,size=.7)  ##  Chromosome length line

		plot <- plot+geom_segment(
			aes(x=ticks,xend=ticks,y=0,yend=y),
			data=chrTicks,size=.7,color="gray25") ## Chromosome line ticks


#~ 		plot<-plot+geom_segment(
#~ 			aes(x=chromStart,xend=chromEnd,y=theMedian,yend=theMedian),
#~ 			data=genome,
#~ 			size=.75,
#~ 			color="orchid"
#~ 	)  ##  Major ticks
		
		plot<-plot+facet_grid(chrom ~ .,scales = "free", space = "free_x")  ##  Facet
		plot<-plot+geom_text(
			aes(x=chromEnd+theMax/200,y=midY,label=chrom),
			data=genome,
			angle=270,
			size=fontSize/pointSize,
			vjust=0
		)  ## Chromosome names
	}
	if ("cen" %in% colnames(genome)) {
		plot<-plot+geom_vline(aes(xintercept=cen),data=genome,color="green3")
	}
	bedName<-levels(bed$name)[1]  ##  Retrieve sample name
	plot<-plot+ggtitle(paste0(bedName," sequencing coverage"))  ##  Plot title
	plot<-plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=score),color="blue",size=0.5)  ##  Raw data
	plot<-plot+gplotTheme
	if (plotting) { plot } else { return(plot) }
}




##  rmMax - removes a row from the bed dataframe that contains maximum value in the "score" column
##  defaults to removing a single row, but n=2 will do it twice, i.e. remove two maximum values
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- rmMax(bedFileContent)  ## removes a single row of data
#~		bedFileContent <- rmMax(bedFileContent,3)  ## removes 3 rows of data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rmMax <- function(bed,n=1) {
	for (times in 1:n) {
		bed<-bed[-which(bed$score==max(bed$score)),]
	}
	rownames(bed) <- 1:nrow(bed)
	return(bed)
}





##  rmChr - removes data for the whole chromosome, for example if you want to get rid of chrM data
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- rmChr(bedFileContent,"chrM")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rmChr <- function(bed,chr) {
	bed<-bed[bed$chrom!=chr,]
	rownames(bed) <- 1:nrow(bed)
	bed$chrom <- factor(bed$chrom, levels = unique(bed$chrom))
	return(bed)
}





##  rmOutliers - removes outliers based on IQR (boxplot's way)
##  defaults to 3xIQR if not supplied with different number
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- rmOutliers(bedFileContent)  # removes datapoints outside of IQR +- 3xIQR range
#~		bedFileContent <- rmOutliers(bedFileContent, 1.5)  # removes datapoints of suspected outliers, IQR +- 1.5xIQR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rmOutliers <- function(bed,lim=3) {
	scoreStats <- summary(bed$score)
	iqr <- scoreStats[5]-scoreStats[2]
	bed<-bed[!(bed$score<scoreStats[2]-lim*iqr | bed$score>scoreStats[5]+lim*iqr),]
	return(bed)
}






##  rmOutliersMan - removes outliers using supplied proportions
##  Bins that contain fewer than Min*Median or more that Max*Median reads per bin will be discarded. 
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		bedFileContent <- rmOutliersMan(bedFileContent,min=0.5,max=1.5)  # removes data outside of the specified range
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rmOutliersMan <- function(bed,loLim=NULL,hiLim=NULL) {
	if (!(is.null(loLim) & is.null(hiLim))) {
		median <- as.numeric(summary(bed$score)[3])
		if (!is.null(loLim)) {
			if (!is.na(as.numeric(loLim))) {
				lo <- as.numeric(as.numeric(loLim)*median)
				bed <- bed[bed$score>lo,]
			}
		}
		if (!is.null(hiLim)) {
			if (!is.na(as.numeric(hiLim))) {
				hi <- as.numeric(as.numeric(hiLim)*median)
				bed <- bed[bed$score<hi,]
			}
		}
	}
	return(bed)
}







##  makeRatio - merges two bed dataframes and calculates ratio of their "score" values,
##  replicating over non-replicating, adjusted by the score sum.
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		mergedBed <- makeRatio(replicatingBed,nonreplicatingBed)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

makeRatio <- function(bedRep,bedNonRep) {
	byVector <- c("chrom","chromStart","chromEnd")
	ratioDF <- merge(bedRep,bedNonRep,by=byVector,all=F,sort=F)
		# rename columns
	mergedNames <- (colnames(ratioDF))
	mergedNames <- gsub(".x",".rep",mergedNames)
	mergedNames <- gsub(".y",".nonRep",mergedNames)
	names(ratioDF) <- mergedNames
		#~ Calculate ratio
	repSum <- sum(as.numeric(ratioDF$score.rep))
	nonRepSum <- sum(as.numeric(ratioDF$score.nonRep))
	corrFactor <- repSum/nonRepSum
	ratioDF$ratio <- ratioDF$score.rep/ratioDF$score.nonRep/corrFactor
		#~ Clean up the dataframe
	ratioDF <- ratioDF[,!(names(ratioDF) %in% c("score.rep","score.nonRep"))]
		## order the resulting ratio dataframe
	ratioDF<-ratioDF[order(ratioDF$chrom,ratioDF$chromStart),]
	rownames(ratioDF) <- 1:nrow(ratioDF)
	ratioDF$ratioFactor <- "1.0"
	return(ratioDF)
}






##  gplotRatio - plots histogram of values in a supplied ~~~~~COLUMN~~~~~ using ggplot2 and
##  highlights interval between 1 and 2 in green
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		gplotRatio(mergedBed$fit)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gplotRatio <- function(ratio,plotting=T) {
	theMax <- max(as.numeric(ratio))
	theMin <-min(as.numeric(ratio))
	
	if (theMin > 0.99) theMin <- 0.99 
	if (theMax < 2.01) theMax <- 2.01

	labels <- makeLabels(theMin,theMax)

	require(ggplot2)
	theTheme<-theme(
		text = element_text(family="Arial"),
		plot.title = element_text(size=rel(2),hjust=0.5,face="bold"),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0.75),
		panel.background = element_blank(),
		axis.line = element_blank(),
		axis.title = element_text(size = rel(1.5)),
		axis.text = element_text(size = rel(1.3)),
		legend.position='none',
		plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
	)
	hist <- ggplot() + aes(ratio)
	hist <- hist + annotate("rect", xmin=1, xmax=2, ymin=0, ymax=Inf, alpha=0.5, fill="lightgreen")
	hist <- hist + geom_histogram(bins=100,fill='gray90',color='gray30')
	hist <- hist + ggtitle("Ratio distribution")
	hist <- hist + scale_x_continuous(
						name="Ratio values",
						breaks=labels$ticks,
						labels=labels$labels,
						limits=c(theMin,theMax),
						expand=c(0.02,0.02))
	hist <- hist + scale_y_continuous(name="Counts within bin")
	hist <- hist + theTheme
	if (plotting) { hist } else { return(hist) }
}





##  trimRatio - function to manually remove outliers from the ratio dataframe.
##	Apply after calculationg the ratio. Has conservative default parameters: 
##  low limit of 0.3 and high limit of 1.7. 
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		ratioDF <- trimRatio(ratioDF,0.5,1.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trimRatio <- function(ratioDF,loLim=0.3,hiLim=1.7) {
	loLim <- as.numeric(loLim)
	hiLim <- as.numeric(hiLim)
	if (hiLim > loLim) {
		if ("tmpRatio" %in% colnames(ratioDF)) {
			ratioDF <- ratioDF[!(ratioDF$tmpRatio<loLim | ratioDF$tmpRatio>hiLim),]
		} else {
			ratioDF <- ratioDF[!(ratioDF$ratio<loLim | ratioDF$ratio>hiLim),]
		}
	} else if (hiLim == loLim) {
		stop("Supplied limits are equal. I quit!")
	} else {
		if ("tmpRatio" %in% colnames(ratioDF)) {
			ratioDF <- ratioDF[!(ratioDF$tmpRatio<hiLim | ratioDF$tmpRatio>loLim),]
		} else {
			ratioDF <- ratioDF[!(ratioDF$ratio<hiLim | ratioDF$ratio>loLim),]
		}
	}
	return(ratioDF)
}






##    normaliseRatio - function to normalise ratio based on either manually supplied factor (based on flow cytometry)
##  or automatically by minimizing difference between the sum of datapoints below 1 and sum of datapoints
##  above 2. Not weighted, so number of datapoints below 1 will be roughly twice of those above 2.
##    Normalised data either replaces the original (default), or added into temporary columns. 
##    Automatically calculated normalisation factor passed as the dataframe comment. To extract it,
##  use 'attributes(mergedBed)$comment' command to extract it
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		mergedBed <- normaliseRatio(mergedBed)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

normaliseRatio <- function(ratio,rFactor=NULL,replace=T) {
	if (is.null(rFactor)) {
		f <- function (factor,values,ceiling=2) abs(sum(factor*values[factor*values>ceiling]-ceiling)-sum(1-factor*values[factor*values<1]))
		xmin <- optimise(f,c(1,2),tol = 0.001,values=na.omit(ratio$ratio),maximum=F)
		ratioFactor <- xmin[[1]]
	} else {
		ratioFactor <- as.numeric(rFactor)
	}
	if (replace) { 
		ratio$ratio <- round(ratio$ratio*ratioFactor,3)
	} else {
		ratio$tmpRatio <- ratio$ratio*ratioFactor
	}
	comment(ratio) <- as.character(round(ratioFactor,3))
	return(ratio)
}






##  smoothRatio - groups ratio values based on the supplied min group and split parameters
##  (defaults to 10 and 5), i.e. it needs 10 adjacent bins to have values in order to create a
##  new group and, if values from 5 successive bins are missing, it closes the group.
##  It then applies cubic spline smoothing within groups.
#~~~~~~~~~~~~~~ example usage:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~		ratioBed <- smoothRatio(ratioBed)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

smoothRatio <- function(ratioDF,groupMin=10,split=5) {
	if (groupMin<4) {
		warning('Minimum group size is four. Setting groupMin parameter to 4.')
		groupMin <- 4
	}
	if (split<1) stop('Split must be a positive integral')
	ratioDF$chrom <- factor(ratioDF$chrom, levels = unique(ratioDF$chrom))
#~ 	if ("name" %in% colnames(ratioDF)) {
#~ 		ratioDF$name <- factor(ratioDF$name, levels = unique(ratioDF$name))
#~ 	}
	if ("name.rep" %in% colnames(ratioDF)) {
		ratioDF$name.rep <- factor(ratioDF$name.rep, levels = unique(ratioDF$name.rep))
	}
	if ("name.nonRep" %in% colnames(ratioDF)) {
		ratioDF$name.nonRep <- factor(ratioDF$name.nonRep, levels = unique(ratioDF$name.nonRep))
	}
		## initiate group number (j)
#~ 	if ("group" %in% colnames(ratioDF)) {
#~ 		if (!is.na(suppressWarnings(max(as.numeric(ratioDF$group))))) {
#~ 			j <- max(ratioDF$group)+1
#~ 		}
#~ 	} else {
#~ 		j <- 1
#~ 	}
	j <- 1
		##  create dummy columns "group" and "splineSmooth"
	if (!("group" %in% colnames(ratioDF))) ratioDF$group <- paste("NA")
	if (!("splineSmooth" %in% colnames(ratioDF))) ratioDF$splineSmooth <- paste("NA")

		##  initiate vectors
	bin <- vector(mode="numeric",length=0)
	group <- vector(mode="numeric",length=0)
	
		##  Actions on every ratio in the supplied dataframe
	for (ratio in unique(paste(as.character(ratioDF$name.rep),as.character(ratioDF$name.nonRep)))) {
		if (length(grep(" ",ratio))==1) {
			ratio <- unlist(strsplit(ratio," "))
			rep <- ratio[1]
			nonRep <- ratio[2]
#~ 			currentRatio <- ratioDF[ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep,]
			rows <- which(ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep)
			chroms <- ratioDF[rows,"chrom"]
			chroms <- factor(chroms,levels=unique(chroms))
			allStarts <- ratioDF[rows,"chromStart"]
			anEnd <- ratioDF[rows,"chromEnd"][1]

			bin <- append(bin,anEnd-allStarts[1])
			if (length(unique(bin))!=1) {
				stop("All samples must have the same genomic bin size")
			} else {
				bin <- unique(bin)
			}

				##  Actions on every chromosome of the current ratio
			for (chr in levels(chroms)) {
					## Grouping lines
 				starts <- allStarts[chroms==chr]
				group <- j
				i <- 1
				for (i in 1:(length(starts)-1)) {
					i <- i + 1
					if ( (starts[i] - starts[i-1])/bin > (split-1) ) { j <- j+1 }
					group <- append(group,j)
				}
					## Remove small groups
				for (aGroup in unique(group)) {
					tmp<-which(group == aGroup)
					if (length(tmp) < groupMin) {
						group[tmp] <- NA
					}
				}
				ratioDF[ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep & ratioDF$chrom==chr,"group"] <- as.integer(group)
				j <- j+1
			}
		}
	}

	#~ Smoothing (spline) in groups
	for (aGroup in as.integer(unique(na.omit(ratioDF$group)))) {
		rows <- which(ratioDF$group==aGroup)
		x <- ratioDF[rows,"chromStart"]
		y <- ratioDF[rows,"ratio"]
		splineObject <- smooth.spline(x,y)
		ratioDF[rows,"splineSmooth"] <- round(as.numeric(splineObject$y),3)
	}
	ratioDF$splineSmooth <- suppressWarnings(as.numeric(ratioDF$splineSmooth))
	ratioDF$group <- suppressWarnings(as.integer(ratioDF$group))
	return(ratioDF)
}





plotGenome <- function(ratioDFs,geom="geom_point",ylims=c(1,2),plotting=T,
	genome=NULL,region=F,guide=NULL,
	lines=NULL,circles=NULL,rectangles=NULL,pointers=NULL,
	colourLines='#00FF00',colourCircles='#FFFFFF',
	colourRectangles='#FF0000',colourPointers='#FF7F00'
) {
#~ if (!is.null(guide)) warning(dfHead(guide))

#~ warning(paste("ratioDFs before:",printColTypes(ratioDFs)))
	ratioDFs$chrom <- factor(ratioDFs$chrom, levels = unique(ratioDFs$chrom))
	ratioDFs$name.rep <- factor(ratioDFs$name.rep, levels=unique(ratioDFs$name.rep))
	ratioDFs$name.nonRep <- factor(ratioDFs$name.nonRep, levels=unique(ratioDFs$name.nonRep))
#~ warning(paste("ratioDFs after:",printColTypes(ratioDFs)))
	if ("name.rep" %in% colnames(ratioDFs)) {
		ratioDFs$name.rep <- factor(ratioDFs$name.rep, levels = unique(ratioDFs$name.rep))
	}
	if ("name.nonRep" %in% colnames(ratioDFs)) {
		ratioDFs$name.nonRep <- factor(ratioDFs$name.nonRep, levels = unique(ratioDFs$name.nonRep))
	}
	ylims <- ylims[order(ylims)]
	samples <- unique(ratioDFs[,c("name.rep","name.nonRep")])
	samples$name <- paste0(as.character(samples[,1])," (",as.character(samples[,2]),")")
		##  First ratio
	firstRatioName <- samples[1,3]

#~ 	if (length(grep(" ",firstRatioName))==1) {
#~ 		ratio <- unlist(strsplit(firstRatioName," "))
		rep <- samples[1,1]
		nonRep <- samples[1,2]
		firstRatio <- ratioDFs[ratioDFs$name.rep==rep & ratioDFs$name.nonRep==nonRep,]
#~  		warning(paste0("First ratio: ",printColTypes(firstRatio)))
#~ 		warning(paste0("First ratio: ",dfHead(firstRatio)))
#~ 	}


		##  If region is supplied, create the wee genome dataframe
	if (class(region)=="character") {
			region <- as.character(region)
			if (length(grep(":",region))==1) {
				region <- unlist(strsplit(unlist(strsplit(region,"-")),":"))
				genome <- data.frame(chrom=as.character(region[1]),chromStart=as.integer(region[2]),chromEnd=as.integer(region[3]))
				xAxisName <- paste(region[1],"coordinates")
			}
	} else {
		if (is.null(genome)) {
				##  Genome is not supplied
					## create genome dataframe using a helper function based on the first sample
				genome <- makeGenome(ratioDFs)
		} else { ## genome dataframe is supplied
				##  remove missing chromosomes from the genome
			goods<-intersect(levels(firstRatio$chrom),levels(genome$chrom))
			if (!all(levels(genome$chrom) %in% goods)) {
				genome<-genome[genome$chrom %in% goods,]
				genome$chrom<-factor(genome$chrom,levels=unique(genome$chrom))
				rownames(genome) <- 1:nrow(genome)
			}
		}
		xAxisName <- "Chromosome coordinates"
	}
	dummy <- genome ## dummy dataframe
	dummy$ratio <- NA
#~ 	warning(paste0("Genome: ",printColTypes(genome)))
#~ 	warning(paste0("Genome: ",dfHead(genome)))

		##  Make nice x-axis labels using a helper function
	theMax <- max(genome$chromEnd)
	theMin <- min(genome$chromStart)
	labels <- makeLabels(theMin,theMax,"b")
	
		## Make a dataframe with the ticks for the chromosome lines
	if (class(region)=="logical") {
		label <- labels$ticks[2]-labels$ticks[1]
		ticks <- numeric()
		chrom <- character()
		for (chr in levels(genome$chrom)) {
			row <- which(genome$chrom == chr)
			tmpTicks <- seq(genome[row,2],genome[row,3],by=label)
			tmpChr <- rep(chr,length(tmpTicks))
			ticks <- append(ticks,tmpTicks)
			chrom <- append(chrom,tmpChr)
		}
		chrTicks <- data.frame(chrom=factor(chrom,levels=unique(genome$chrom)),ticks=as.numeric(ticks))
	}
		##  Get the bin size for the first sample
	bin <- firstRatio$chromEnd[1]-firstRatio$chromStart[1]

	##  load ggplot2 library and set plotting theme
	fontSize <- 10
	pointSize <- 1 / 0.352777778
	require(ggplot2)
	gplotTheme<-theme(
		text = element_text(family="Arial",size=fontSize),
		plot.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		plot.title = element_text(hjust = 0.5,size=rel(1.8)),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.line.x = element_line(colour="gray50",size=.5),
		axis.title = element_text(size = rel(1.6),face="bold"),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size=rel(1.4)),
		legend.position="bottom",
		legend.spacing.x=unit(0.1,"cm"),
		legend.title=element_blank(),
		legend.text = element_text(size=rel(1)),
		legend.key = element_blank(),
		legend.key.size = unit(c(0.4,0.4),"cm"),
		strip.background = element_blank(),
		strip.text = element_blank(),
		plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
	)

	##  start creating plot object - labels, scales names etc.
	plot <- ggplot(genome)
	plot <- plot + scale_y_continuous(
		name="Relative copy number",
		limits = c(ylims[1]-((ylims[2]-ylims[1])/5),ylims[2]+((ylims[2]-ylims[1])/5))
	)  ##  y scale
	plot <- plot + annotate(
		geom="segment",
		x=0,
		xend=0,
		y=ylims[1],
		yend=ylims[2],
		colour="gray50",
		lwd=0.5
	)  ##  y scale line
	plot <- plot + scale_x_continuous(
		breaks=labels$ticks,
		labels=labels$labels,
		name=xAxisName,
		limits=c(theMin-(theMax-theMin)/50,theMax+(theMax-theMin)/50),
		expand=c(0.02,0.02)
	)  ##  x scale
	if (class(region)=="logical") {
		plot <- plot + facet_grid(
			chrom ~ .,
			scales = "free",
			space = "free_x"
		)  ##  Facet
		plot <- plot+geom_text(
			aes(x=chromEnd+(theMax-theMin)/100,y=ylims[1]+((ylims[2]-ylims[1])/2),label=chrom),
			angle=270,
			size=fontSize/pointSize,
			vjust=0
		)  ##  Chromosome labels
		plot <- plot + geom_segment(
			aes(x=0,xend=genome$chromEnd,y=(ylims[1]-((ylims[2]-ylims[1])/10)),yend=(ylims[1]-((ylims[2]-ylims[1])/10))),
			size=.6,color="gray50"
		)  ##  Chromosome length line
		
		plot <- plot+geom_segment(
			aes(x=ticks,xend=ticks,y=ylims[1]-((ylims[2]-ylims[1])/6),yend=ylims[1]-((ylims[2]-ylims[1])/10)),
			data=chrTicks,size=.7,color="gray35") ## Chromosome line ticks		
		
	}
	plot <- plot + geom_segment(
		aes(x=0,xend=chromEnd,y=ylims[1],yend=ylims[1]),
		linetype="dashed",
		size=.25,
		color="gray40"
	)  ##  Lower y line
	plot <- plot + geom_segment(
		aes(x=0,xend=chromEnd,y=ylims[2],yend=ylims[2]),
		linetype="dashed",
		size=.25,
		color="gray40"
	)  ##  Upper y line
	plot <- plot + geom_text(
		aes(x=theMin-(theMax-theMin)/100,y=ylims[1],label=as.character(ylims[1])),
		size=fontSize/pointSize,
		hjust=1
	)  ##  Lower y break
	plot <- plot + geom_text(
		aes(x=theMin-(theMax-theMin)/100,y=ylims[2],label=as.character(ylims[2])),
		size=fontSize/pointSize,
		hjust=1
	)  ##  Upper y break

##  initialising legendary vectors here
	color.names <- vector(mode="character",length=0)
	color.values <- vector(mode="character",length=0)
	color.lines <- vector(mode="character",length=0)
	color.shapes <- vector(mode="numeric",length=0)
	color.sizes <- vector(mode="numeric",length=0)
	color.fills<-  vector(mode="character",length=0)
	fill.names <- vector(mode="character",length=0)
	fill.values <- vector(mode="character",length=0)
	fill.shapes <- vector(mode="numeric",length=0)
	fill.sizes <- vector(mode="numeric",length=0)
	fill.colors <- vector(mode="character",length=0)
	fill.lines<- vector(mode="character",length=0)
	line.names <-vector(mode="character",length=0)
	line.values <-vector(mode="character",length=0)
	line.colors <- vector(mode="character",length=0)
	line.shapes <- vector(mode="numeric",length=0)
	line.sizes <- vector(mode="numeric",length=0)

	if (is.null(guide)) {
##  plot everything (scatter + line)
			##  Colors
		dotColors <- c("gray50","deepskyblue4","orange3","darkolivegreen4","hotpink4","purple4","red3","yellow3")
#~ 		lineColors <- c("gray60","deepskyblue3","orange2","darkolivegreen3","hotpink3","purple3","red2","yellow2")
			##  Actions on every ratio in the supplied dataframe
		for (row in nrow(samples)) {
			rep <- samples[row,1]
			nonRep <- samples[row,2]
			currentRatio <- ratioDFs[ratioDFs$name.rep==rep & ratioDFs$name.nonRep==nonRep,]
			if (class(region)=="character") {
				currentRatio <- currentRatio[currentRatio$chrom==region[1] & currentRatio$chromStart>=as.numeric(region[2]) & currentRatio$chromEnd<=as.numeric(region[3]),]
			}
			geom_str <- paste0("geom_point(aes(x=chromStart+bin/2,y=ratio,color='",rep," (",nonRep,")'),data=currentRatio,size=0.1)")
			plot <- plot + eval(parse(text=geom_str)) ##  Raw data
			if (as.character(geom)=="geom_ribbon") {
				fills <- rbind(fills,data.frame(name=samples[row,3],fill=dotColors[row],color=NA,shape=15,line=NA,size=2,stringsAsFactors=F))
			} else {
					color.names<-c(color.names,amples[row,3])
					color.values<-c(color.values,dotColors[row])
					color.shapes<-c(color.shapes,20)
					color.lines<-c(color.lines,"blank")
					color.sizes<-c(color.sizes,5)
					color.fills<-c(color.fills,NA)
			}
			
			if ("splineSmooth" %in% colnames(currentRatio)) {
				if (length(na.omit(currentRatio$splineSmooth)) > 0) {
					geom_str <- paste0("geom_line(aes(x=chromStart+bin/2,y=splineSmooth,group=group),color='",dotColors[row],"',data=currentRatio)")
					plot <- plot + eval(parse(text=geom_str))
				}
			}
		}
	} else {  ## guided plotting (shiny)
		geoms <- c("geom_point","geom_ribbon","geom_segment")
		geom_pars <- c(
			"geom_point(aes(x=chromStart+bin/2,y=ratio,color='changeMe'),size=0.2,data=currentRatio)",
			"geom_ribbon(aes(x=chromStart+bin/2,ymin=ylims[1]-(ylims[2]-ylims[1])/20,ymax=ratio,fill='changeMe'),data=currentRatio)",
			"geom_segment(aes(x=chromStart+bin/2,xend=chromStart+bin/2,yend=ratio,y=ylims[1]-(ylims[2]-ylims[1])/20),data=currentRatio,color=guide$color[i])"
			)
		geom_str <- geom_pars[match(geom,geoms)]
		guide <- guide[order(guide$order,decreasing=T,na.last = T),]
		rownames(guide) <- 1:nrow(guide)
		layers <- length(na.omit(guide$order))
		sampleNames <- character()
		for (i in 1:layers) {
			if (guide[i,"raw"] | guide[i,"smooth"]) {
				rep <- guide$name.rep[i]
				nonRep <- guide$name.nonRep[i]
				currentRatio <- ratioDFs[ratioDFs$name.rep==rep & ratioDFs$name.nonRep==nonRep,]
				ratioFactor <- currentRatio$ratioFactor[1]
				sampleName <- paste0(rep," (",nonRep,",",ratioFactor,")")
				sampleNames <- append(sampleNames,sampleName)
				if (class(region)=="character") {
					currentRatio <- currentRatio[currentRatio$chrom==region[1] & currentRatio$chromStart>=as.numeric(region[2]) & currentRatio$chromEnd<=as.numeric(region[3]),]
				}
				if (xor(as.logical(guide$raw[i]),as.logical(guide$smooth[i]))) { ##  Either raw only or smooth only
					if (as.logical(guide$raw[i])) { # raw only
						geom_string <- gsub("changeMe",sampleName,geom_str)
						plot <- plot + eval(parse(text=geom_string))
					} else { # smooth only
						geom_string <- gsub("changeMe",sampleName,geom_str)
						geom_string <- gsub("=ratio,","=splineSmooth,group=group,",geom_string)
						geom_string <- gsub("geom_point","geom_line",geom_string)
						geom_string <- gsub("size=0.2,","",geom_string)
						plot <- plot + eval(parse(text=geom_string))
					}
				} else if (as.logical(guide$raw[i]) & as.logical(guide$smooth[i])) { ## Both
					geom_string <- gsub("changeMe",sampleName,geom_str)
					plot <- plot + eval(parse(text=geom_string))
					geom_string <- paste0("geom_line(aes(x=chromStart+bin/2,y=splineSmooth,group=group),data=currentRatio,color='",as.character(guide$color[i]),"')")
					plot <- plot + eval(parse(text=geom_string))
				}					
				if (geom == "geom_point") {
					color.names<-c(color.names,sampleName)
					color.values<-c(color.values,guide$color[i])
					color.shapes<-c(color.shapes,20)
					color.lines<-c(color.lines,"solid")
					color.sizes<-c(color.sizes,5)
					color.fills<-c(color.fills,NA)
				} else if (geom == "geom_ribbon") {
					fill.names<-c(fill.names,sampleName)
					fill.values<-c(fill.values,guide$color[i])
					fill.shapes<-c(fill.shapes,22)
					fill.sizes<-c(fill.sizes,4)
					fill.colors<-c(fill.colors,NA)
					fill.lines<-c(fill.lines,"blank")	
				} else if (geom=="geom_segment") {
					geom_string <- paste0("geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,linetype='",sampleName,"'),data=dummy)")
					plot <- plot + eval(parse(text=geom_string))
					line.names<-c(line.names,sampleName)
					line.values<-c(line.values,"solid")
					line.colors<-c(line.colors,guide$color[i])
					line.shapes<-c(line.shapes,NA)
					line.sizes<-c(line.sizes,5)	
				}
			}
		}
	}

		##  optional plot features	
	if ("p.value" %in% colnames(firstRatio)) {
		temp<-which(firstRatio$p.value < 0.001)
		sign.999<-firstRatio[temp,]
		temp<-which(firstRatio$p.value <= 0.01 & firstRatio$p.value > 0.001)
		sign.99<-firstRatio[temp,]
	##  Plot significant p-values
		if (nrow(sign.999)!=0) {
			if (class(region)=="character") {
				sign.999 <- sign.999[sign.999$chrom==region[1] & sign.999$chromStart>=as.numeric(region[2]) & sign.999$chromEnd<=as.numeric(region[3]),]
			}
			plot<-plot+geom_segment(aes(
				x=chromStart+bin/2,
				xend=chromStart+bin/2,
				y=ylims[2]+((ylims[2]-ylims[1])/10),
				yend=ylims[2]+((ylims[2]-ylims[1])/5)
			),data=sign.999,color="gray25")
			plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,linetype="p-value***"),data=dummy)  ##  dummy
#~ 			colors <- rbind(colors,data.frame(name="p-value***",color="gray25",shape=NA,line="solid",size=1,stringsAsFactors=F))
			line.names<-c(line.names,"p-value***")
			line.values<-c(line.values,"solid")
			line.colors<-c(line.colors,"gray25")
			line.shapes<-c(line.shapes,NA)
			line.sizes<-c(line.sizes,2)
		}
		if (nrow(sign.99)!=0) {
			if (class(region)=="character") {
				sign.99 <- sign.99[sign.99$chrom==region[1] & sign.99$chromStart>=as.numeric(region[2]) & sign.99$chromEnd<=as.numeric(region[3]),]
			}
			plot<-plot+geom_segment(aes(
				x=chromStart+bin/2,
				xend=chromStart+bin/2,
				y=ylims[2]+((ylims[2]-ylims[1])/10),
				yend=ylims[2]+((ylims[2]-ylims[1])/5)
			),data=sign.99,color="gray75")
			plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,linetype="p-value**"),data=dummy)  ##  dummy
#~ 			colors <- rbind(colors,data.frame(name="p-value**",color="gray50",shape=NA,line="solid",size=1,stringsAsFactors=F))
			line.names<-c(line.names,"p-value**")
			line.values<-c(line.values,"solid")
			line.colors<-c(line.colors,"gray75")
			line.shapes<-c(line.shapes,NA)
			line.sizes<-c(line.sizes,2)
		}
	}

	if (!is.null(lines)) {
		lines$chrom <- factor(lines$chrom,levels=levels(genome$chrom))
		linesName <- levels(lines$name)[1]
		if (class(region)=="character") {
			lines <- lines[lines$chrom==region[1] & lines$chromStart>=as.numeric(region[2]) & lines$chromEnd<=as.numeric(region[3]),]
		}
		plot<-plot+geom_segment(aes(x=chromStart,xend=chromEnd,y=ylims[1],yend=ylims[2]),data=lines,color=colourLines)
		plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,linetype=linesName),data=dummy)  ##  dummy
		line.names<-c(line.names,linesName)
		line.values<-c(line.values,"solid")
		line.colors<-c(line.colors,colourLines)
		line.shapes<-c(line.shapes,NA)
		line.sizes<-c(line.sizes,3)
	}

	if (!is.null(circles)) {
		circles$chrom <- factor(circles$chrom,levels=levels(genome$chrom))
		circlesName <- levels(circles$name)[1]
		if (class(region)=="character") {
			circles <- circles[circles$chrom==region[1] & circles$chromStart>=as.numeric(region[2]) & circles$chromEnd<=as.numeric(region[3]),]
		}
		plot <- plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=ylims[1]-((ylims[2]-ylims[1])/10),fill=circlesName),
			data=circles,color="black",size=3,shape=21)
		fill.names<-c(fill.names,circlesName)  ##  Circles legend values
		fill.values<-c(fill.values,coloursCircles)
		fill.shapes<-c(fill.shapes,21)
		fill.sizes<-c(fill.sizes,3)
		fill.colors<-c(fill.colors,"black")
		fill.lines<-c(fill.lines,"blank")
	}

	if (!is.null(rectangles)) {  ##  If a dataframe supplied, plot the rectangles
		rectangles$chrom <- factor(rectangles$chrom,levels=levels(genome$chrom))
		rectanglesName <- levels(rectangles$name)[1]
		if (class(region)=="character") {
			rectangles <- rectangles[rectangles$chrom==region[1] & rectangles$chromStart>=as.numeric(region[2]) & rectangles$chromEnd<=as.numeric(region[3]),]
		}
		plot <- plot+geom_rect(aes(
			xmin=chromStart,xmax=chromEnd,
			ymin=(ylims[1]-((ylims[2]-ylims[1])/10)) - (ylims[2]-ylims[1])/20,
			ymax=(ylims[1]-((ylims[2]-ylims[1])/10)) + (ylims[2]-ylims[1])/20),
			data=rectangles,color=colourRectangles,fill=NA)
		plot <- plot+geom_point(aes(x=theMin-(theMax-theMin)/10,y=ylims[1]-((ylims[2]-ylims[1])/2),fill=rectanglesName),data=dummy,shape=22)  ##  dummy
		fill.names<-c(fill.names,rectanglesName)  ##  Rectangles legend values
		fill.values<-c(fill.values,NA)
		fill.shapes<-c(fill.shapes,22)
		fill.sizes<-c(fill.sizes,5)
		fill.colors<-c(fill.colors,colourRectangles)
		fill.lines<-c(fill.lines,"blank")
	}
	if (!is.null(pointers)) {  ##  If a dataframe supplied, plot the pointers
		pointers$chrom <- factor(pointers$chrom,levels=levels(genome$chrom))
		pointersName <- levels(pointers$name)[1]
		if (class(region)=="character") {
			pointers <- pointers[pointers$chrom==region[1] & pointers$chromStart>=as.numeric(region[2]) & pointers$chromEnd<=as.numeric(region[3]),]
		}
		plot <- plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=ylims[1],fill=pointersName),
			data=pointers,color="black",size=3,shape=25)
		fill.names<-c(fill.names,pointersName)  ##  Pointers legend values
		fill.values<-c(fill.values,colourPointers)
		fill.shapes<-c(fill.shapes,25)
		fill.sizes<-c(fill.sizes,4)
		fill.colors<-c(fill.colors,"black")
		fill.lines<-c(fill.lines,"blank")
	}



		my.colors<-color.values  ##  Preparing named vector for colors
		names(my.colors)<-color.names
		my.fills<-fill.values
		names(my.fills)<-fill.names
		my.lines<-line.values
		names(my.lines)<-line.names

		if (length(my.colors) > 0) {
			plot <- plot + scale_colour_manual(
				name="",values=my.colors,guide = guide_legend(override.aes = list(
					linetype=color.lines,
					shape= color.shapes,
					size=color.sizes,
					fill=color.fills)),breaks=names(my.colors)
			)
		}
		if (length(my.fills) > 0) {
			plot <- plot + scale_fill_manual(
				name="",values=my.fills,guide = guide_legend(override.aes = list(
					#size=fill.sizes,
					shape=fill.shapes,
					linetype=fill.lines,
					color=fill.colors)),breaks = names(my.fills)
			)
		}
		if (length(my.lines) > 0) {
			plot <- plot + scale_linetype_manual(
				name="",values=my.lines,guide = guide_legend(override.aes = list(
					color=line.colors,
					#size=line.sizes,
					shape=line.shapes)),breaks = names(my.lines)
			)
		}
		##  custom theme
	plot <- plot + gplotTheme

		##  and we're done
	if (plotting) {
		if (region==F) {
			X11(width=7,height=10,pointsize=15)
		} else {
			X11(width=7,height=3,pointsize=15)
		}
		plot
	} else {
		return(plot)
	}
}




###  ratioStats
###  merges two ratio dataframes horizontally
###  and calculates differences between ratios
###  based on z-score stats

ratioStats <- function(ratio1,ratio2,names=NULL) {
	if (!is.null(names)) {
		nameRatio1 <- as.character(names[1])
		nameRatio2 <- as.character(names[2])
	} else {
		if (typeof(ratio1)=="character" & typeof(ratio2)=="character") {
			nameRatio1<-ratio1
			nameRatio2<-ratio2
			try(ratio1 <- eval(as.name(nameRatio1)))
			try(ratio2 <- eval(as.name(nameRatio2)))
		} else {
			nameRatio1<-deparse(substitute(ratio1))
			nameRatio2<-deparse(substitute(ratio2))
		}
	}
	if (typeof(ratio1) != "list") { stop("Could not find dataframe ",nameRatio1,"!")
	} else if (typeof(ratio2) != "list") { stop("Could not find dataframe ",nameRatio2,"!")
	} else {
		ratios <- merge(ratio1,ratio2,by=c("chrom","chromStart","chromEnd"),all=T,sort=F)
		ratios <- ratios[order(ratios$chrom,ratios$chromStart),]
		rownames(ratios) <- 1:nrow(ratios)
#~ 		cat("\tComparing replication profiles between ",nameRatio1," and ",nameRatio2,"\n",sep="")
		calc <- data.frame(diff=ratios$ratio.x-ratios$ratio.y,means=(ratios$ratio.x+ratios$ratio.y)/2)
		q=qqnorm(calc$diff, plot.it=F)
		qnt_x <- q$x[! quantile(q$x,0.25,na.rm=T)>=q$x & q$x<=quantile(q$x,0.75,na.rm=T)]
		qnt_y <- q$y[! quantile(q$x,0.25,na.rm=T)>=q$x & q$x<=quantile(q$x,0.75,na.rm=T)]
		my_fit=lm(qnt_y ~ qnt_x)
		est_mean=my_fit$coefficients[1]
		est_sd=my_fit$coefficients[2]
		calc$z.score <- (abs(calc$diff-est_mean))/est_sd
		ratios$p.value <- 1-pnorm(calc$z.score)
		
		longRatios <- subset(ratios,!is.na(name.rep.x),select=names(ratios)[-c(10:15)])
		names(longRatios) <- gsub("\\.x$","",names(longRatios))
		tmp <- subset(ratios,!is.na(name.rep.y),select=names(ratios)[-c(4:9)])
		tmp$p.value <- NA
		names(tmp) <- gsub("\\.y$","",names(tmp))
		longRatios <- rbind(longRatios,tmp)
		rownames(longRatios) <- 1:nrow(longRatios)
		return(longRatios)
		}
}




###  makeGenome helper function
###  creates genome dataframe based on a supplied bed/ratio dataframe
###  understands regions in 'chrX:1000-2000' string format
###  supplies axis name as a comment to the genome dataframe

makeGenome <- function(DF,region=F) {
	if (region==F) {
		##  using the whole supplied dataframe
		DF$chrom <- factor(DF$chrom,levels=unique(DF$chrom))
		chroms <- levels(DF$chrom)
		ends <- vector(mode="numeric",length=0)
		midY <- vector(mode="numeric",length=0)
#~ 	warning(paste0("Making genome for chromosomes: "),paste0(chroms,collapse=","))
		if ("score" %in% colnames(DF)) {
			for (chr in chroms) {
				ends <- append(ends,max(DF[DF$chrom==chr,"chromEnd"]))
				midY <- append(midY,max(DF[DF$chrom==chr,"score"])/2)
			}
		} else if ("ratio" %in% colnames(DF)) {
			for (chr in chroms) {
				ends <- append(ends,max(DF[DF$chrom==chr,"chromEnd"]))
				midY <- append(midY,max(DF[DF$chrom==chr,"ratio"])/2)
			}
		}
		genome <- data.frame(
			"chrom"=as.character(chroms),
			"chromStart"=as.integer(rep(0,length(chroms))),
			"chromEnd"=as.integer(ends),
			"midY"=as.numeric(midY))
		genome$chrom <- factor(genome$chrom,levels=unique(genome$chrom))
		attr(genome,"axisName") <- paste("Chromosome coordinates")
	} else {
		genome <- data.frame(chrom=as.character(NA),chromStart=as.integer(NA),chromEnd=as.integer(NA),stringsAsFactors=FALSE)
		##  if region is supplied, first create the region DF
		region <- as.character(region)
		if (length(grep(":",region))==1) {
			region <- unlist(strsplit(region,":"))
			tmp <- unlist(strsplit(region[2],"-"))
			region <- c(region[1],tmp)
		} else {
			tmp <- c(0,max(subset(DF,chrom==region,select="chromEnd")$chromEnd))
			region <- append(region,tmp)
		}
		if (length(region)==3) {
			if (region[3]>region[2]) {
				genome[1,1] <- as.character(region[1])
				genome[1,2] <- as.integer(region[2])
				genome[1,3] <- as.integer(region[3])
			} else {
				stop("Supplied region end should be bigger than region start")
			}
		} else {
			stop("Couldn't understand the supplied region. Please supply it as a string in 'chrX' or 'chrX:1000-2000' format.") 
		}
		attr(genome,"axisName") <- as.character(paste(region[1],"coordinates"))
	}
	return(genome)
}




###  makeLabels helper function
###  creates named vector of labels based on supplied minimum and maximum values
###  units may be added

makeLabels <- function(theMin,theMax,unit="") {
#~ 	warning(paste0("Making labels for ",theMin,":",theMax," range"))
	if (!is.na(as.numeric(theMin)) & !is.na(as.numeric(theMax))) {
		step <- 10^floor(log10(theMax-theMin))
		minorFactor <- 2
		if ((theMax-theMin)/step < 6) {
			step <- step/2
			minorFactor <- 2
		}
		if ((theMax-theMin)/step < 6) {
			step <- step/2
			minorFactor <- 2.5
		}
		if ((theMax-theMin)/step < 6) {
			step <- step/2.5
			minorFactor <- 2
		}
		minTick <- step*ceiling(theMin/step)
		if ((minTick-theMin)/step > 0.5) minTick <- minTick-step
		maxTick <- step*floor(theMax/step)
		if ((theMax-maxTick)/step > 0.5) maxTick <- maxTick+step
		labels <- seq(from = minTick, to = maxTick, by = step)
		if (log10(step) < 0) {
			prefix <- switch(abs(ceiling(log10(step)/3))+1,"","m","μ","n","p","f","a")
			divisor <- switch(abs(ceiling(log10(step)/3))+1,1,0.001,0.000001,0.000000001,0.000000000001,0.000000000000001,0.000000000000000001)
		} else {
			prefix <- switch(floor(log10(step)/3)+1,"","k","M","G","T","P","E")
			divisor <- switch(floor(log10(step)/3)+1,1,1000,1000000,1000000000,1000000000000,1000000000000000,1000000000000000000)
		}
		names <- paste0(labels/divisor,prefix,unit)
		names <- gsub(paste0("^",0,prefix,unit),"0",names)
		labels <- data.frame(ticks=labels,labels=names)
		attr(labels,"minorTickFactor") <- as.character(minorFactor)
		return(labels)
	} else {
		stop("There is a problem with the supplied minimum and maximum!")
	}
}



# use this function to determine hex colors:
getColorHex <- function(color)
	{
		c <- col2rgb(color)
		sprintf("#%02X%02X%02X", c[1],c[2],c[3])
	}

# a function to help debugging. Prints column types in a dataframe
printColTypes <- function(df) {
	string <- paste("Dim:",dim(df)[1],"rows,",dim(df)[2],"columns.")
	colNames <- names(df)
	for (i in 1:ncol(df)) {
		string <- paste0(string,colNames[i],":",typeof(df[,i]),",",class(df[,i]),"|")
	}
	return(string)
}

# a function to help debugging. Prints head of a dataframe data
dfHead <- function(df) {
	string <- paste("Dataframe head (",dim(df)[1],"rows,",dim(df)[2],"columns). ")
	string <- paste0(string,paste0(names(df),collapse=":")," | ")
	if (nrow(df) > 3 ) {
		j <- 3
	} else {
		j <- nrow(df)
	}
	for (i in 1:j) {
		string <- paste0(string,paste0(df[i,],collapse=":")," | ")
	}
	return(string)
}
