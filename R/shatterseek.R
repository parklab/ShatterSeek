SVs <- setClass("SVs",
				representation(
							   chrom1="character",
							   pos1 = "numeric",
							   chrom2 = "character",
							   pos2 = "numeric",
							   SVtype = "character",
							   strand1 = "character",
							   strand2 = "character",
							   numSV = "numeric"
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
			  .Object@numSV = length(.Object@chrom1)
			  ind.match.sv1 = match(.Object@chrom1,chromNames)
			  ind.match.sv2 = match(.Object@chrom2,chromNames)	

			  if(sum(is.na(ind.match.sv1))!=0 | sum(is.na(ind.match.sv2))!=0){
				  stop(paste("chromosome name must be in \"", paste(chromNames,collapse=" ")),"\"")
			  }

			  .Object
				})

#setMethod("as.data.frame","SVs",function(.Object, ...){
#	rslt = data.frame(chrom1=.Object@chrom1,
#			  pos1=.Object@pos1,
#			  chrom2=.Object@chrom2,
#			  pos2 = .Object@pos2,
#			  SVtype = .Object@SVtype
#			)
#	rslt
#	})

setMethod("show","SVs",function(object){
			  print(paste("SVs with", object@numSV, "structural variations"))
				})


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

#setMethod("as.data.frame","SVs",function(.Object, ...){
#	rslt = data.frame(chrom=.Object@chrom,start=.Object@start,end=.Object@end,total_cn=.Object@total_cn)
#	rslt
#	})

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
		rt.value = list(SV=SV.df,graph=NULL,connComp=list(),num.chromth=0,chromothripsis=list(),chromothripsis.chr=c(),maxSVs=0,degree=list())
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
	chromothripsis.Reg.chr = sapply(chromothripsis.Reg,FUN=function(v){SV.df$chrom1[as.numeric(v[1])]})

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

	rt.value=list(SV=SV.df,graph=gn,connComp=cmpnt,num.chromth=length(ind.LargComp),chromothripsis = chromothripsis.Reg,chromothripsis.chr=chromothripsis.Reg.chr,maxSVs = maxSVs,degree = degree.chromtheripsis)
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

###---------------------------------------------------------------------------

#cancer_types = c("BLCA", "BOCA", "BRCA", "BTCA", "CESC", "CLLE", "COAD", "DLBC", "EOPC", "ESAD","CMDI",
#				 "GACA", "GBM" , "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG" , "LICA", "LIHC",
#				 "LINC", "LIRI", "LUAD", "LUSC", "MALY", "MELA", "ORCA", "OV"  , "PACA", "PAEN",
#				 "PBCA", "PRAD", "READ", "RECA", "SARC", "SKCM", "STAD", "THCA", "UCEC")

#sampleSVproportionalChr = function(numSV,probChrs,cancer.type){
#
#	if (!missing(cancer.type)){ 
#		if (!(cancer.type %in% cancer_types)){stop("The input cancer type is not supported. Please set the cancer.type argument")}
#		SV.icgc.pancer.now = SV.icgc.pancer[SV.icgc.pancer$cancer_type == cancer.type,]
#	} else {SV.icgc.pancer.now = SV.icgc.pancer}
#
#	#tol.sv = nrow(SV.tcga.pancer)
#	tol.sv = nrow(SV.icgc.pancer.now)
#	if (numSV < 10){
#		indSV.smple = sample.int(tol.sv,numSV)
#	}else{
#		indSV.smple = sample.int(tol.sv,numSV,prob=probChrs)
#	}
#	#SV.tmp = SV.tcga.pancer[indSV.smple,]
#	SV.tmp = SV.icgc.pancer.now[indSV.smple,]
#	return(SV.tmp)
#}


#sampleSV = function(numSV,cancer.type){
#	#print(cancer.type)
#
#	if (cancer.type != ""){ 
#		if (!(cancer.type %in% cancer_types)){stop("The input cancer type is not supported. Please set the cancer.type argument")}
#		SV.icgc.pancer.now = SV.icgc.pancer[SV.icgc.pancer$cancer_type == cancer.type,]
#		#print(dim(SV.icgc.pancer.now))
#	} else {SV.icgc.pancer.now = SV.icgc.pancer}
#
#	#tol.sv = nrow(SV.tcga.pancer)
#	tol.sv = nrow(SV.icgc.pancer.now)
#	#print(tol.sv)
#	if (nrow(SV.icgc.pancer.now) < numSV){
#		indSV.smple = sample.int(tol.sv,numSV,replace=T)} else{
#			indSV.smple = sample.int(tol.sv,numSV)} 
#	#SV.tmp = SV.tcga.pancer[indSV.smple,]
#	SV.tmp = SV.icgc.pancer.now[indSV.smple,]
#	return(SV.tmp)
#}
#
##--------------------------------------------------------------------------------------------------------------------------------------
#assignPvaluePairedEnd = function(chromothSample,B=5000,cores=10,poisson=TRUE,mean.poisson,min.Size=1,method="Ruibin",cancer.type=""){
#		#print(dim(SV.icgc.pancer))
#        numSV.sample = sum(chromothSample$SV$chrom1 == chromothSample$SV$chrom2)
#        registerDoMC(cores=cores)
#				
#		#--------------------- sample SVs with the same per-chromosome distribution as the sample
#		ratios <- table(chromothSample$SV$chrom1) / numSV.sample #interchromosomal SV probs
#		
#		probabilities <- data.frame(chr=chromNames, probs= rep(0,23))
#		probabilities$probs <- ratios[ match(probabilities$chr, names(ratios))]
#		probabilities$probs[is.na(probabilities$probs)] <- 0
#		
#		####ff = table(SV.icgc.pancer[,1])
#		#####ff = table(SV.tcga.pancer[,1])
#		####names(ff)[23] = "23"
#		####ff = ff[order(as.numeric(names(ff)),decreasing=F)]
#		####probs_all <- (rep(probabilities$probs, ff) / rep(ff,ff))
#		#---------------------
#        rslt = foreach(icount(B),.combine=rbind) %dopar%{
#		if(poisson) numSV=rpois(1,mean.poisson) else numSV=numSV.sample
#	#			print(numSV)
#				if (method == "Ruibin"){
#                tmp=sampleSV(numSV=numSV,cancer.type=cancer.type) ## Ruibin
#				} else{
#				tmp=sampleSVproportionalChr(numSV=numSV,prob=probs_all,cancer.type=cancer.type)
#				}
#                tmp1 = cluster.SV(tmp,min.Size=1,chromNames=chromNames)
#                tmp1$maxClusterSize[,2]
#                }
#	#print(dim(rslt))
#	maxClustersizeSample = sapply(1:nrow(rslt),FUN=function(i){return(max(rslt[i,]))})
#
##	ind.largerThanObs = (maxClustersizeSample>=max(chromothSample$maxClusterSize[,2]))
##        pvalue.tmp = sapply(1:ncol(rslt),FUN=function(j){sum(rslt[,j]>=chromothSample$maxClusterSize[j,2] | ind.largerThanObs | chromothSample$maxClusterSize[j,2]<=min.Size)/B})
#	
#	#pvalue.tmp = pvalue.tmp*length(chromNames)
#	#ind = which(pvalue.tmp>1.0)
#	#pvalue.tmp[ind] = 1.0
#
##	chromothSample$SVpvalue = chromothSample$maxClusterSize
##        chromothSample$SVpvalue[,2] = pvalue.tmp
##	names(chromothSample$SVpvalue)[2] = "pvalue"
#	return(chromothSample)
#        }

##--------------------------------------------------------------------------------------------------------------------------------------

##--------------------------------------------------------------------------------------------------------------------------------------
#assignPvalueCNV = function(chromothSample,segSample,B=5000){
#        chromth.size = unlist(sapply(chromothSample$"chromothripsis",length))
#
#        SVs =  chromothSample$SV
#        cmpnt = chromothSample$"connComp"
#        tmp.chr = as.character(sapply(cmpnt,FUN=function(v){SVs$chrom1[as.numeric(v[1])]}))
#        tmp.Size = sapply(cmpnt,length)
#        maxCluster.byChrom = chromothSample$maxClusterSize
#
#        CNVpavlue.tmp.bychrom = maxCluster.byChrom
#        CNVpavlue.tmp.bychrom[,2] = 1.0
#        chromothSample$CNVpvalue = CNVpavlue.tmp.bychrom
#        names(chromothSample$CNVpvalue) = c("chrom","pvalue")
#        for(j1 in 1:nrow(maxCluster.byChrom)){
#                chr = maxCluster.byChrom[j1,1]
#                ind.chr = which(tmp.chr==chr)
#                cmpnt.chr = cmpnt[ind.chr]
#                tmp.Size.chr = sapply(cmpnt.chr,length)
#
#
#                CNVpavlue.tmp = c()
#                if (length(ind.chr)==0){
#                        CNVpavlue.tmp.bychrom[j1,2] = 1.0
#                        }else{
#                        for(j in ind.chr[tmp.Size.chr==max(tmp.Size.chr)]){
#                                ind.connComp = as.numeric(chromothSample$"chromothripsis"[[j]])
#                                chrom.name = chromothSample$SV$chrom1[ind.connComp[1]]
#                                ind.seg.chrom = which(segSample$chr==chrom.name)
#                                segs.tmp = segSample[ind.seg.chrom,2:4]
#                                SV.tmp = chromothSample$SV[ind.connComp,]
#                                min.Brk = min(SV.tmp$pos1,SV.tmp$pos2)
#                                max.Brk = max(SV.tmp$pos1,SV.tmp$pos2)
#                                ind.seg.bw = which(segs.tmp$end>=min.Brk & segs.tmp$start<=max.Brk)
#                                if(length(ind.seg.bw)>0){
#                                        segs.tmp = segs.tmp[ind.seg.bw,]
#                                        p.tmp = peak.count(segs.tmp[,3])
#                                        pp.tmp = replicate(B,peak.count(sample(segs.tmp[,3],length(segs.tmp[,3]))))
#                                        CNVpavlue.tmp = c(CNVpavlue.tmp,sum(pp.tmp>=p.tmp)/B)
#                                        }else{
#                                        CNVpavlue.tmp = c(CNVpavlue.tmp,1)
#                                        }
#                                }
#                                CNVpavlue.tmp.bychrom[j1,2] = min(CNVpavlue.tmp)
#                        }
#                }
#        chromothSample$CNVpvalue = CNVpavlue.tmp.bychrom
#	chromothSample$CNV = segSample
#        return(chromothSample)
#        }

##--------------------------------------------------------------------------------------------------------------------------------------


#filterSV = function(SV.sample){ ##SV.sample should be a SVs object; return a data.frame object
#	freq.cutoff = 0.2
#	size.cutoff = 1e4
#	SV.sample = as(SV.sample,"data.frame")
#	# keep interchromosomal SVs
#	ind = which(SV.sample$chrom1==SV.sample$chrom2)
#	interchr_SVs = SV.sample[SV.sample$chrom1!=SV.sample$chrom2,]
#	SV.sample = SV.sample[ind,]
#	if(length(ind)==0) return(SV.sample)
#	size = abs(SV.sample$pos1-SV.sample$pos2)
#	ind = which(size>size.cutoff)
#	SV.sample = SV.sample[ind,]
#	if(length(ind)==0) return(SV.sample)
#
#	win.size = 200
#	brkLft.range = GRanges(seqnames=SV.tcga.FP$chrom1,ranges=IRanges(SV.tcga.FP$pos1-win.size,SV.tcga.FP$pos1+win.size),seqinfo=Seqinfo(seqnames=chromNames))
#	brkRgt.range = GRanges(seqnames=SV.tcga.FP$chrom2,ranges=IRanges(SV.tcga.FP$pos2-win.size,SV.tcga.FP$pos2+win.size),seqinfo=Seqinfo(seqnames=chromNames))
#
#	query.left = GRanges(seqnames=SV.sample$chrom1,ranges=IRanges(SV.sample$pos1-win.size,SV.sample$pos1+win.size))
#	query.right = GRanges(seqnames=SV.sample$chrom2,ranges=IRanges(SV.sample$pos2-win.size,SV.sample$pos2+win.size))
#
#	ovlp.left = as.matrix(findOverlaps(query.left,brkLft.range))
#	ovlp.right = as.matrix(findOverlaps(query.right,brkRgt.range))
#
#	tmp.left = paste(ovlp.left[,1],ovlp.left[,2],sep="-")
#	tmp.right = paste(ovlp.right[,1],ovlp.right[,2],sep="-")
#	ind.mtch = match(tmp.left,tmp.right)
#	tmp.match = tmp.left[!is.na(ind.mtch)]
#	if(length(tmp.match)==0) return(SV.sample)
#
#	tmp.split = strsplit(tmp.match,"-",fixed=TRUE)
#	ovlp.tmp = as.numeric(unlist(tmp.split))
#	ovlp = matrix(0,nrow=length(tmp.split),ncol=2)
#	ovlp[,1] = ovlp.tmp[seq(from=1,to=length(ovlp.tmp),by=2)]
#	ovlp[,2] = ovlp.tmp[seq(from=2,to=length(ovlp.tmp),by=2)]
#
#	#tmp = aggregate(ovlp[,2],by=list(ovlp[,1]),FUN=function(v){max(SV.tcga.FP$occFreq[v])})
#	#occFreq = rep(0,nrow(SV.sample))
#	#occFreq[tmp[,1]] = tmp[,2]
#	#SV.sample = SV.sample[occFreq<freq.cutoff,]
#	return(SV.sample)
#}


shatterseek = function(SV.sample,seg.sample,min.Size=1){
	#,cores=20,B=5000,filter=TRUE,poisson=TRUE,mean.poisson,method="Ruibin",cancer.type=""){

	cat("Running..\n\n\n")
	if(!is(SV.sample,"SVs")){stop("SV.sample must be a SVs object")}
	if(!missing(seg.sample)){
		if(!is(seg.sample,"CNVsegs")) {stop("seg.sample must be a CNVsegs object")}
	}
	#if(B<100) B = 100
	#if(filter) SV.sample = filterSV(SV.sample) 
	#	else{
	SV.sample = as(SV.sample,"data.frame")
	#SV.sample = SV.sample[SV.sample$chrom1==SV.sample$chrom2,] ##if filter==TRUE, this was already performed
	#	}

	#if(missing(mean.poisson)) mean.poisson = meanSVNum.tcga.pancer
	#if(mean.poisson<=0) stop("mean.poisson must be positive")
	chromothSample = cluster.SV(SV.sample[SV.sample$chrom1==SV.sample$chrom2,],min.Size=min.Size,chromNames=chromNames) ## pass only intra
	#chromothSample = assignPvaluePairedEnd(chromothSample,B=B,cores=cores,poisson=poisson,method=method,
	#										   mean.poisson=mean.poisson,min.Size=min.Size,cancer.type=cancer.type)
	chromothSample$SV = SV.sample[SV.sample$chrom1==SV.sample$chrom2,]
	chromothSample$SVinter = SV.sample[SV.sample$chrom1!=SV.sample$chrom2,]

	chromSummary = data.frame(chromothSample$maxClusterSize)#,SVpvalue=chromothSample$SVpvalue[,2])
	#chromSummary$SVpvalue[chromSummary$SVpvalue==0] = 1/B
	#chromSummary$CNVpvalue[chromSummary$CNVpvalue==0] = 1/B
	#chromSummary$score = -log10(chromSummary$SVpvalue) -log10(chromSummary$CNVpvalue)
	#chromSummary$score = -log10(chromSummary$SVpvalue)
	seg.sample = as(seg.sample,"data.frame")
	chromothSample$CNV = seg.sample

	if(!missing(seg.sample)){
		seg.sample = as(seg.sample,"data.frame")
		chromothSample$CNV = seg.sample
		#chromothSample = assignPvalueCNV(chromothSample,seg.sample,B=B)
		#chromSummary$CNVpvalue =chromothSample$CNVpvalue[,2]
		#chromSummary$CNVpvalue[chromSummary$CNVpvalue==0] = 1/B
		#chromSummary$score = -log10(chromSummary$SVpvalue) -log10(chromSummary$CNVpvalue)
	}
	
	cat("Evaluating the statistical criteria\n")
	#chromSummary = statistical_criteria(chromSummary)
	cat("Successfully finished!\n")
	return(chromoth(chromSummary=chromSummary,detail=chromothSample))
}

#plot.chromoth = function(chromothripsis,chrom,sample){
#	SV.tmp = chromothripsis@detail$SV
#	SVs.tmp = SV.tmp[SV.tmp$chrom1==chrom[1]&SV.tmp$chrom2==chrom[1],]
#	if(nrow(SVs.tmp)==0) stop("Error: there is no SV to plot")
#	xlim = c(min(SVs.tmp$pos1,SVs.tmp$pos2),max(SVs.tmp$pos1,SVs.tmp$pos2))
#	if("CNV" %in% names(chromothripsis@detail)){
#		CNV.tmp = chromothripsis@detail$CNV
#		segs.tmp = CNV.tmp[CNV.tmp$chrom==chrom[1],2:4]
#		xlim = c(min(segs.tmp[,1],SVs.tmp$pos1,SVs.tmp$pos2),max(segs.tmp[,2],SVs.tmp$pos1,SVs.tmp$pos2))
#	}
#
#	rg = xlim[2]-xlim[1]
#	xlim[1] = max(xlim[1]-0.05*rg,0)
#	xlim[2] = xlim[2]+0.05*rg
#	if(missing(sample)){
#		xlab = paste("chrom ",chrom,sep="")
#	}else{
#		xlab=paste(sample,": ","chrom ",chrom,sep="")
#	}
#	ind = which(chromothripsis@chromSummary$chrom==chrom) 
#	par(mfrow=c(2,1))
#	plot.SV(SVs.tmp,chrom,xlim=xlim,col="red")
#	l.text1 = paste("p-value-SV:",chromothripsis@chromSummary$SVpvalue[ind])
#	l.text2 = paste("Cluster Size:",chromothripsis@chromSummary[ind,2])
#	if("CNV" %in% names(chromothripsis@detail)){
#		segment.plot(segs.tmp,xlab=xlab)
#		l.text3 = paste("p-value-CNV:",chromothripsis@chromSummary$CNVpvalue[ind])
#		legend("topright",legend=c(l.text1,l.text3,l.text2))
#	}else{
#		legend("topright",legend=c(l.text1,l.text2))
#	}
#	#dev.off()
#}
#
#
#plot.SV = function(SV,chr,col="black",xlim=c(0,1e8),ylim=c(-1,4),add=FALSE,xlab=""){
#	ind = which(SV$chrom1==chr&SV$chrom2==chr)
#	plot.mbridge(SV$pos1[ind],SV$pos2[ind],col=col,xlim=xlim,ylim=ylim,add=add,xlab=xlab)
#}
#
#
#plot.bridge  = function(start,end,ystart,yend,height=4,col="red"){
#	x.mid = mean(c(start,end))
#	y.mid = max(ystart+height,yend+height)
#	lines(c(start,x.mid,end),c(ystart,y.mid,yend),col=col)
#}
#
#plot.mbridge = function(start,end,col="red",xlim=c(0,1e8),ylim=c(-1,4),add=FALSE,xlab="",ylab=""){
#	if(length(start)!=length(end)) stop("start and end has different length")
#	#plot.new()
#	#plot.window(xlim=xlim,ylim=ylim)
#	if(!add) {plot(xlim,c(0,0),type="l",ylim=ylim,xlab=xlab,ylab=ylab)}
#	for(i in c(1:length(start))){
#		plot.bridge(start[i],end[i],0,0,col=col)
#	}
#}
#
#
##setMethod("plot.chromoth",c("chromoth","character","character"),function(chromothripsis,chrom,sample){
##        plot.chromoth(chromothripsis,chrom,sample)
##        })
#
#getSVs = function(object){
#	if(!is(object,"chromoth")){stop("object must be a chromoth object")}
#	SV = object@detail$SV
#	SV = SVs(chrom1=SV$chrom1,pos1=SV$pos1,chrom2=SV$chrom2,pos2=SV$pos2,SVtype=SV$SVtype,strand1=SV$strand1,strand1=SV$strand2) #XX
#	return(SV)
#}
#getSegs = function(object){
#	if(!is(object,"chromoth")) {stop("object must be a chromoth object")}
#	seg =  object@detail$CNV
#	seg = CNVsegs(chrom=seg$chrom,start=seg$start,end=seg$end,total_cn=seg$total_cn)
#	return(seg)
#}


