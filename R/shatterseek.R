#' Class to store SV data
#' @param chrom1 (character): chromosome for the first breakpoint
#' @param pos1 (character): position for the first breakpoint
#' @param chrom2 (character): chromosome for the second breakpoint
#' @param pos2 (character): position for the second breakpoint
#' @param SVtype (character): type of SV, encoded as: DEL (deletion-like; +/-), DUP (duplication-like; -/+), h2hINV (head-to-head inversion; +/+), and t2tINV (tail-to-tail inversion; -/-).
#' @param strand1 (e.g. + for DEL)
#' @param strand2 (e.g. - for DEL)
#' @return an instance of the class 'SVs' that contains SV data. Required format by the function shatterseek
#' @export

SVs <- setClass("SVs",
				representation(
							   chrom1="character",
							   pos1 = "numeric",
							   chrom2 = "character",
							   pos2 = "numeric",
							   SVtype = "character",
							   strand1 = "character",
							   strand2 = "character"
							   #numSV = "numeric"
							   ) 
				)

setMethod("initialize", "SVs", function(.Object, ...) {
			  .Object <- callNextMethod()
			  if(length(.Object@chrom1) != length(.Object@pos1)|| length(.Object@chrom1) != length(.Object@pos2))
				  stop("slots lengths are not all equal")
			  if(length(.Object@chrom2)==0) .Object@chrom2 = .Object@chrom1
			  if(length(.Object@SVtype)==0) .Object@SVtype = rep("",length(.Object@chrom1))
			  if(length(.Object@chrom1)!=length(.Object@SVtype) || length(.Object@chrom1) != length(.Object@chrom2) ) stop("slots lengths are not all equal")

			  ind = which(.Object@pos1>.Object@pos2)
			  if(length(ind)>0){
				  tmp1 = .Object@chrom2[ind]
				  .Object@chrom2[ind] = .Object@chrom1[ind]
				  .Object@chrom1[ind] = tmp1
				  tmp1 = .Object@pos2[ind]
				  .Object@pos2[ind] = .Object@pos1[ind]
				  .Object@pos1[ind] = tmp1
			  }
			  #.Object@numSV = length(.Object@chrom1)
			  ind.match.sv1 = match(.Object@chrom1,chromNames)
			  ind.match.sv2 = match(.Object@chrom2,chromNames)	

			  if(sum(is.na(ind.match.sv1))!=0 | sum(is.na(ind.match.sv2))!=0){
				  stop(paste("chromosome name must be in \"", paste(chromNames,collapse=" ")),"\"")
			  }

			  .Object
				})

#setMethod("show","SVs",function(object){
#			  print(paste("SVs with", object@numSV, "structural variations"))
#				})


setAs("SVs","data.frame",function(from, to){
		  to = data.frame(chrom1=from@chrom1,
						  pos1=from@pos1,
						  chrom2=from@chrom2,
						  pos2 = from@pos2,
						  strand1 = from@strand1,
						  strand2 = from@strand2,
						  SVtype = from@SVtype,stringsAsFactors=FALSE
						  )

				})

#' Class to store CNV data
#'
#' @param chrom (character): chromosome (also in Ensembl notation)
#' @param start (numeric): start position for the CN segment
#' @param end (numeric): end position for the CN segment
#' @param CN (numeric): integer total copy number (e.g. 2 for unaltered chromosomal regions)
#' @return an instance of the class 'CNVsegs' that contains CNV data. Required format by the function shatterseek
#' @export
CNVsegs = setClass("CNVsegs",
				   representation(
								  chrom="character",
								  start = "numeric",
								  end="numeric",
								  total_cn = "numeric",
								  numSegs = "numeric"
								  )
				   )

setMethod("initialize", "CNVsegs", function(.Object, ...) {
			  .Object <- callNextMethod()
			  if(length(.Object@chrom) != length(.Object@start)|| length(.Object@chrom) != length(.Object@total_cn) || length(.Object@chrom) != length(.Object@end))
				  stop("slots lengths are not all equal")
			  .Object@numSegs = length(.Object@chrom)
			  ind.match.cnv = match(.Object@chrom,chromNames)
			  if(sum(is.na(ind.match.cnv))!=0){
				  stop(paste("chromosome name must be in \"", paste(chromNames,collapse=" ")),"\"")
			  }

			  .Object
				   })

setMethod("show","CNVsegs",function(object){
			  print(paste("CNVsegs with", object@numSegs, "segments"))
				   })

setAs("CNVsegs","data.frame",function(from,to){
		  to = data.frame(chrom=from@chrom,start=from@start,end=from@end,total_cn=from@total_cn,stringsAsFactors=FALSE)
				   }
)


chromoth = setClass("chromoth",
					representation(
								   chromSummary = "data.frame",
								   detail = "list"
								   )
					)

setMethod("show","chromoth",function(object){
			  print(object@chromSummary)
					})

setAs("chromoth","data.frame",function(from,to){
		  to = from@chromSummary
					})

## if chromNames is not NULL, return the maximum number of cluster size in each chromsome specified in chromNames
cluster.SV = function(SV.sample,min.Size=1,chromNames){
	SV.df = as.data.frame(SV.sample)
	SV.df$chromothEvent = rep(0,nrow(SV.df))

	chromNames = as.character(chromNames)
	numSV.byChrom = data.frame(chrom=chromNames,numbSV=rep(0,length(chromNames)),stringsAsFactors=FALSE)
	maxCluster.byChrom = data.frame(chrom=chromNames,clusterSize=rep(0,length(chromNames)),stringsAsFactors=FALSE)

	ind.chr = which(SV.df$chrom1==SV.df$chrom2)
	if(length(ind.chr)<1){
		rt.value = list(SV=SV.df,graph=NULL,connComp=list(),num.chromth=0,
						#chromothripsis=list(),
						#chromothripsis.chr=c(),
						maxSVs=0,degree=list())
		rt.value$numSVByChrom = numSV.byChrom
		rt.value$maxClusterSize = maxCluster.byChrom
		return(rt.value)
	}
	tmp.numSV.byChrom = table(SV.df$chrom1[ind.chr])
	ind.numSV.byChrom = match(names(tmp.numSV.byChrom),numSV.byChrom$chrom)
	numSV.byChrom$numbSV[ind.numSV.byChrom[!is.na(ind.numSV.byChrom)]] = tmp.numSV.byChrom[!is.na(ind.numSV.byChrom)]


	gr.chr =  GRanges(seqnames = Rle(SV.df$chrom1[ind.chr]),ranges = IRanges(SV.df$pos1[ind.chr]-1,SV.df$pos2[ind.chr]+1),strand=Rle("+",length(ind.chr)))
	ovlp.chr = as.data.frame(findOverlaps(gr.chr,gr.chr))
	#ovlp.chr = data.frame(ovlp.chr[[1]],ovlp.chr[[2]])
	ind.ovlp.chr = which(ovlp.chr[,1] < ovlp.chr[,2])
	# to remove redundancies: e.g. 1 -2 ; 2 -1 both indicate that the SVs 1 and 2 overlap
	ovlp.chr = ovlp.chr[ind.ovlp.chr,]


	ovlp.chr.witin = as.data.frame(findOverlaps(gr.chr,gr.chr,type="within"))
	ind.rm = which(ovlp.chr.witin[,1]==ovlp.chr.witin[,2])
	if(length(ind.rm)>0) ovlp.chr.witin = ovlp.chr.witin[-ind.rm,]
	#ovlp.chr.witin = ovlp.chr.witin[-ind.rm,]
	if(nrow(ovlp.chr.witin)>0 & nrow(ovlp.chr)>0){
		ovlp.chr.witin1 = ovlp.chr.witin
		ind.tmp.within = which(ovlp.chr.witin[,2] < ovlp.chr.witin[,1])
		if(length(ind.tmp.within)>0){
			ovlp.chr.witin[ind.tmp.within,1] = ovlp.chr.witin1[ind.tmp.within,2]
			ovlp.chr.witin[ind.tmp.within,2] = ovlp.chr.witin1[ind.tmp.within,1]
		}
		ind.equal.rm = as.data.frame(findOverlaps(IRanges(ovlp.chr[,1],ovlp.chr[,2]),IRanges(ovlp.chr.witin[,1],ovlp.chr.witin[,2]),type="equal"))
		if(length(ind.equal.rm)>0) ovlp.chr = ovlp.chr[-ind.equal.rm[,1],]
	}
	adjMatrix = matrix(0,nrow=length(ind.chr),ncol=length(ind.chr))
	if(nrow(ovlp.chr)>0){
		adjMatrix[ovlp.chr[,1] + (ovlp.chr[,2]-1)*length(ind.chr)] = 1
		adjMatrix[ovlp.chr[,2] + (ovlp.chr[,1]-1)*length(ind.chr)] = 1
	}
	rownames(adjMatrix) = ind.chr
	colnames(adjMatrix) = ind.chr
	g1 <- graphAM(adjMat=adjMatrix)
	gn = as(g1,"graphNEL")
	cmpnt = connComp(gn) ## finds the clusters

	ind.LargComp = which(sapply(cmpnt,length)>=min.Size)
	chromothripsis.Reg = cmpnt[ind.LargComp]
	#chromothripsis.Reg.chr = sapply(chromothripsis.Reg,FUN=function(v){SV.df$chrom1[as.numeric(v[1])]})

	degree.chromtheripsis = list()
	maxSVs = 0
	if(length(ind.LargComp)>0){
		for(i in 1:length(chromothripsis.Reg)){
			gn.sub = subGraph(chromothripsis.Reg[[i]],gn)
			degree.chromtheripsis[[i]] = degree(gn.sub)
		}
		tmp = as.numeric(unlist(chromothripsis.Reg))
		SV.df$chromothEvent[unique(tmp)] = 1
		maxSVs = max(sapply(chromothripsis.Reg,length))
	}

	rt.value=list(SV=SV.df,graph=gn,connComp=cmpnt,num.chromth=length(ind.LargComp),
				  #chromothripsis = chromothripsis.Reg,
				  #chromothripsis.chr=chromothripsis.Reg.chr,
				  maxSVs = maxSVs,degree = degree.chromtheripsis)
	rt.value$numSVByChrom = numSV.byChrom

	tmp.chr = as.character(sapply(cmpnt,FUN=function(v){SV.df$chrom1[as.numeric(v[1])]}))
	tmp.Size = sapply(cmpnt,length)
	tmp.maxSize = aggregate(tmp.Size,by=list(tmp.chr),max)
	ind.match = match(chromNames,tmp.maxSize[,1])
	maxCluster.byChrom[,2][!is.na(ind.match)] = tmp.maxSize[ind.match[!is.na(ind.match)],2]
	maxCluster.byChrom[maxCluster.byChrom[,2]<min.Size,2] = 0
	rt.value$maxClusterSize = maxCluster.byChrom
	#print(rt.value)
	return(rt.value)
}

#' Main ShatterSeek function
#' Identified cluster of interleaved SVs and calculates statistical metrics for each chromosome (chromosomes 1-22 and X)
#' @param SV.sample an instance of class SVs
#' @param seg.sample an instance of class CNVsegs
#' @param min.Size minimum number of inleaved SVs required to report a cluster. Default is 3
#' @param genome reference genome (hg19 or hg38)
#' @export
shatterseek = function(SV.sample,seg.sample,min.Size=1, genome="hg19"){
	cat("Running..\n\n\n")
	if(!is(SV.sample,"SVs")){stop("SV.sample must be a SVs object")}
	if(!missing(seg.sample)){
		if(!is(seg.sample,"CNVsegs")) {stop("seg.sample must be a CNVsegs object")}
	}
  if ( !(as.character(genome) %in% c("hg19","hg38"))){stop("Reference genome assembly is not supported (Use hg19 or hg38)")}
  
	SV.sample = as(SV.sample,"data.frame")
	chromothSample = cluster.SV(SV.sample[SV.sample$chrom1==SV.sample$chrom2,],min.Size=min.Size,chromNames=chromNames) ## pass only intra
	chromothSample$SV = SV.sample[SV.sample$chrom1==SV.sample$chrom2,]
	chromothSample$SVinter = SV.sample[SV.sample$chrom1!=SV.sample$chrom2,]

	chromSummary = data.frame(chromothSample$maxClusterSize)#,SVpvalue=chromothSample$SVpvalue[,2])
	seg.sample = as(seg.sample,"data.frame")
	chromothSample$CNV = seg.sample

	if(!missing(seg.sample)){
		seg.sample = as(seg.sample,"data.frame")
		chromothSample$CNV = seg.sample
	}
	
	out = chromoth(chromSummary=chromSummary,detail=chromothSample)
	cat("Evaluating the statistical criteria\n")
	out@chromSummary = suppressWarnings(statistical_criteria(out, genome))
	cat("Successfully finished!\n")
	return(out)
}
