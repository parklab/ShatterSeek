plot_chromothripsis <- function(shatterSeek_output, chr=chr,BAF=NULL,sample_name="",
								DEL_color='darkorange1',DUP_color='blue1',
								t2tINV_color="forestgreen",h2hINV_color="black",
								arc_size=.2){

	if ( !(as.character(chr) %in% chromNames)){stop("Chromosome not valid")}

	common_ggplot2 <- theme_bw() + theme(axis.text.x=element_text(size=7,angle=0),
										 axis.text.y=element_text(size=7),
										 axis.title.y=element_text(size=7),
										 axis.title.x=element_blank(),
										 legend.position="none",
										 legend.text = element_text(size=7),
										 legend.key = element_blank(),
										 plot.margin=unit(c(0.1,0.1,0,0.1),"cm"),
										 plot.title=element_blank(),
										 panel.grid.major = element_blank(),
										 panel.grid.minor = element_blank(), 
										 legend.title=element_blank(),
										 plot.background = element_blank(),
										 axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
										 axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))  
	#-------------------------------------------------------------------------------------------------------

	cand = gsub("chr","",chr)
	chr=paste("chr",cand,sep="")
	summary = shatterSeek_output@chromSummary
	candidate_chrs <- shatterSeek_output@chromSummary$chrom 
	cluster_sizes <- sapply(shatterSeek_output@detail$connComp,length)
	# get the SVs for the cluster of SVs in this chromosome
	cand_clust_size <- shatterSeek_output@chromSummary$clusterSize[ shatterSeek_output@chromSummary$chrom == cand]
	idx = which(cluster_sizes == cand_clust_size) 
	SVsnow <- shatterSeek_output@detail$SV##[ as.numeric(unlist(shatterSeek_output@detail$connComp[idx])) ,]
	CNVsnow <- shatterSeek_output@detail$CNV  
	SVsnow <- unique(SVsnow[SVsnow$chrom1 == cand, ]) # remove if there are more
	CNVsnow <- CNVsnow[CNVsnow$chrom == cand, ] # remove if there are more

	#if(nrow(SVsnow) <1 | nrow(CNVsnow) < 1){stop("NO SVs or CN segments in this chromosome")}


	y1=4; y2=12
	df = SVsnow
	df$y1 = rep(y1,nrow(df)); df$y2 = rep(y2,nrow(df))
	min_coord=0;max_coord=10
	if (nrow(df)!=0){
	  min_coord = min(df$pos1)
	  max_coord = max(df$pos2)
  	df$diff = abs(df$pos1 - df$pos2)
  	df$curv = 1 -(df$diff / max(df$diff))
  	max_diff = max(df$diff)
  	df$curv[(df$diff / max_diff) > 0.2] <- .15
  	df$curv[(df$diff / max_diff) > 0.8] <- .08
  	df$curv[(df$diff / max_diff) < 0.2] <- 1
	}
	d = data.frame(x=c(min_coord),y=c(1),leg=c("DEL","DUP","t2tINV","h2hINV"))
	idx = c()
	
	idx1=which(CNVsnow$start >= min_coord)
	idx2=which(CNVsnow$end <= max_coord)   
	idx=intersect(idx1,idx2)
	CNVsnow = CNVsnow[idx,]

	SV_plot = ggplot(d,aes(x=x,y=y)) +
		geom_point(colour="white") + ylim(0,y2+5) + common_ggplot2  +
		geom_line(data=rbind(d,d),aes(x=x,y=y,colour=leg)) 

	#-------------------------------------------------------------------------------------------------------
	# Interchromosomal SVs
	#-------------------------------------------------------------------------------------------------------
	inter <- shatterSeek_output@detail$SVinter   
	inter <- inter[which(inter$chrom1 == cand | inter$chrom2 == cand), ] 

	if (nrow(inter)>0){
		inter$SVtype = factor(inter$SVtype,levels=c("DEL","DUP","h2hINV","t2tINV"))
		inter$y = rep(0,nrow(inter))
		inter$y[which(inter$SVtype %in% c("DUP","DEL"))] = 4
		inter$y[!(inter$SVtype %in% c("DUP","DEL"))] = 12
		inter$type_SV = rep("",nrow(inter))
		inter$type_SV[which(inter$strand1 == "-" & inter$strand2 == "-")] = "t2tINV"
		inter$type_SV[which(inter$strand1 == "-" & inter$strand2 == "+")] = "DUP"
		inter$type_SV[which(inter$strand1 == "+" & inter$strand2 == "-")] = "DEL"
		inter$type_SV[which(inter$strand1 == "+" & inter$strand2 == "+")] = "h2hINV"
		inter$SVtype = inter$type_SV; inter$type_SV=NULL
		inter$colour = rep("",nrow(inter))
		inter$colour[which(inter$SVtype == "DUP")] = "blue1"
		inter$colour[which(inter$SVtype == "DEL")] = "darkorange1"
		inter$colour[which(inter$SVtype == "h2hINV")] = "black"
		inter$colour[which(inter$SVtype == "t2tINV")] = "forestgreen"

		inter = data.frame(pos = c(inter$pos1,inter$pos2), y=c(inter$y,inter$y), SVtype=c(inter$SVtype,inter$SVtype) )
		SV_plot = SV_plot + geom_point(data=inter,size=1,alpha=1,aes(x=pos,y=as.numeric(y),colour=SVtype))
	}
	#-------------------------------------------------------------------------------------------------------

	SV_plot = SV_plot + geom_hline(yintercept=y1,size=0.5) + geom_hline(yintercept=y2,size=0.5) 
	
	if(nrow(df)>300){options(expressions= 100000)} 
	#to avoid "Error: evaluation nested too deeply: infinite recursion / options(expressions=)?"

	now = df[df$SVtype == "DUP",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = now$curv[i],colour="blue1",ncp=8)
		}
		idx= c(idx,1)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))


	now = df[df$SVtype == "DEL",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = -1*now$curv[i],colour="darkorange1") 
		}
		idx= c(idx,2)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))

	now = df[df$SVtype == "t2tINV",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = now$curv[i],colour="forestgreen") 
		}
		idx= c(idx,3)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y2)) + 
		geom_point(data=now,size=.5,aes(x=pos2,y=y2))

	now = df[df$SVtype == "h2hINV",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=arc_size,data = now[i,] , 
										   aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = -1*now$curv[i],colour="black") 
		}
		idx= c(idx,4)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y2)) + geom_point(data=now,size=.5,aes(x=pos2,y=y2))


	SV_plot = SV_plot + theme(axis.ticks.x=element_blank(),panel.border = element_blank(),
							  axis.title.y=element_text(colour="white"),
							  axis.text.y=element_text(colour="white"),
							  axis.ticks.y=element_line(colour="white")) + 
	scale_x_continuous(expand = c(0.01,0.01))+ coord_cartesian(xlim=c(min_coord,max_coord))

	idx = c(1,2,3,4)
	vals = c(DEL_color='darkorange1',DUP_color='blue1',t2tINV_color="forestgreen",h2hINV_color="black")
	labs = c('DEL','DUP',"t2tINV","h2hINV")

	SV_plot = SV_plot +  scale_colour_manual(name = 'SV type', 
										 values =c('darkorange1','blue1',"forestgreen","black"),
										 labels = labs[idx]) + theme(legend.position="none")


	# CN_plot = ggplot() +
	# 	geom_segment(data=CNVsnow, aes(x = start, y =total_cn , xend = end, yend = total_cn),size=2) + 
	# 	common_ggplot2 + 
	# 	ylab("CN")+ xlim(min_coord,max_coord) +  
	# 	scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})+
	# 	#coord_cartesian(xlim=c(min_coord,max_coord)) +
	# 	scale_y_continuous(minor_breaks = NULL,breaks=c(0,sort(unique(CNVsnow$total_cn))),
	# 					   limits=c(min(CNVsnow$total_cn) -0.15,max(CNVsnow$total_cn)+.25)) + xlab(NULL)
	# 	# first xlim y and then expand


	#-------------------------------------------------------------------------------------------------------
	if (max(CNVsnow$total_cn) <=10){
		CNV_plot = ggplot() +
			geom_segment(data=CNVsnow, aes(x = start, y =total_cn , xend = end, yend = total_cn),size=2) + 
			common_ggplot2 + ylab("CN")+ xlab(NULL)+
			scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))
		CNV_plot = CNV_plot + 
			scale_y_continuous(minor_breaks = NULL,breaks=c(0,sort(unique(CNVsnow$total_cn))),limits=c(min(CNVsnow$total_cn) -0.35,max(CNVsnow$total_cn)+.35)) 
	}

	if (max(CNVsnow$total_cn) <=100 & max(CNVsnow$total_cn) >10){
		idx = which(CNVsnow$total_cn ==0);if(length(idx)>0){CNVsnow$total_cn[which(CNVsnow$total_cn ==0)]=0.99}
		mm2 = max(CNVsnow$total_cn)
		if(mm2 <=20){mm = c(0,1,2,4,10,mm2)}
		if(mm2 <=40 & mm2 >20){mm = c(0,1,2,4,10,20,mm2)}
		if(mm2 <=100 & mm2 >40){mm = c(0,1,2,4,10,20,mm2)}

	CNV_plot = ggplot() + geom_segment(data=CNVsnow, aes(x = start, y =log(total_cn,base=2), xend = end, yend = log(total_cn,base=2)),size=2) + 
					common_ggplot2 + ylab("CN")+ xlab(NULL)+ 
					scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))

	CNV_plot = CNV_plot + scale_y_continuous(minor_breaks = NULL,breaks=log(mm,base=2), 
											   limits=c( -.01 + log(min(CNVsnow$total_cn),base=2), 
														log(max(CNVsnow$total_cn) , base=2)+.01 ), 
											   labels=as.character(mm))
	}

	if (max(CNVsnow$total_cn) >100){
		idx = which(CNVsnow$total_cn ==0);if(length(idx)>0){CNVsnow$total_cn[which(CNVsnow$total_cn ==0)]=0.99}
		CNV_plot = ggplot() +
			geom_segment(data=CNVsnow, aes(x = start, y =log(total_cn,base=10) , xend = end, yend = log(total_cn,base=10)),size=2) + 
			common_ggplot2 + #scale_color_manual(values=c("forestgreen","red1","blue","brown")) +
			ylab("CN")+ xlab(NULL)+ 
			scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})#,limits=c(min_coord,max_coord))
		mm2 = max(CNVsnow$total_cn)
		mm = c(0,1,2,4,20,50,100,mm2)
		CNV_plot = CNV_plot + scale_y_continuous(minor_breaks = NULL,breaks=log(mm,base=10), 
												 limits=c( log(min(CNVsnow$total_cn),base=10), 
														  log(max(CNVsnow$total_cn) , base=10)), 
												 labels=as.character(mm))
	}


	#-------------------------------------------------------------------------------------------------------
	common_ggplot2_chrom =  theme(panel.grid.major = element_blank(),
								  panel.grid.minor = element_blank(),
								  panel.border = element_blank(),
								  plot.background = element_blank(),
								  axis.text=element_blank(),
								  axis.title=element_blank(),
								  plot.title=element_blank(), 
								  axis.ticks=element_blank())


	##chr_info = readRDS("/Users/isidro/Dropbox/Park_lab/paper_chromothripsis/chr_info.rds")
	chr_info$color[chr_info$gieStain == "gneg"] = "white"
	chr_info$color[chr_info$gieStain == "gpos25"] = "grey75"
	chr_info$color[chr_info$gieStain == "gpos50"] = "grey50"
	chr_info$color[chr_info$gieStain == "gpos75"] = "grey25"
	chr_info$color[chr_info$gieStain == "gpos100"] = "grey0"
	chr_info$color[chr_info$gieStain == "acen"] = "red"
	chr_info = chr_info[chr_info$seqnames==chr,]
	chr_info$y = rep(1,nrow(chr_info))

	if(nrow(chr_info)>7){chr_info_annot=chr_info}else{
	chr_info_annot=chr_info[seq(3,(nrow(chr_info)-3),3),]}

	ideogram =ggplot(data=chr_info,aes(x=y,y=y)) +geom_point(colour="white")
	ideogram = ideogram +  
		geom_rect(data=chr_info,mapping = aes(xmin = chr_info$start, xmax = chr_info$end,   ymin = 0, ymax = 1),
				  fill = chr_info$color,color="black",size=.1) 

	ideogram = ideogram +ylim(0,4)
	ideogram = ideogram + annotate(geom = "text",x = (chr_info_annot$start+chr_info_annot$end)/2,y=chr_info_annot$y + 1.5,
								   label=chr_info_annot$name,vjust=.5, angle=90,size=2)
	ideogram = ideogram +theme_bw() + common_ggplot2_chrom + #xlim(min_coord,max_coord) + 
		scale_x_continuous(expand = c(0.01,0.01)) +
		coord_cartesian(xlim=c(min_coord, max_coord))

	#----------------------------------------------------------------------
	# BAF
	#----------------------------------------------------------------------
	#XXX 
	if (!is.null(BAF)){
		BAF$ratio1 = BAF$b_ads1/BAF$b_dp

		BAF_plot = ggplot(BAF[BAF$b_dp > 20,],aes(x=pos,y=ratio1,colour=ratio1)) + 
			geom_point(size=.1,alpha=.1,colour="black")+
			#scale_colour_gradient2(mid="orange", high="black", low="black") +
			common_ggplot2 + #xlim(min_coord,max_coord) + 
			ylab("BAF")+ xlim(min_coord,max_coord) +  
			scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")}) +
			scale_y_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"),limits=c(0,1))  

	}

	summary$pos = paste(summary$chrom,":",summary$start,"-",summary$end,sep="")
	summary$oscillations = paste(summary$max_number_oscillating_CN_segments_2_states,";",
								 summary$max_number_oscillating_CN_segments_2_states_3states,sep="")

	cols_sel = c("pos",
				 #"clusterSize",
				 "clusterSize_including_TRA",
				 "number_SVs_sample",
				 "max_number_oscillating_CN_segments_2_states",
				 "number_CNV_segments",
				 "chr_breakpoint_enrichment",
				 "pval_exp_cluster",
				 "pval_fragment_joins",
				 "inter_other_chroms_coords_all")

	table_now = summary[which(summary$chrom == cand),cols_sel]

	names(table_now) = c("Position",
						 #"Nb. Intrachr. SVs", 
						 "Total nb. SVs (intrachr. + transl.)",
						 "SVs in sample",
						 "Oscillating CN (2 and 3 states)","CN segments",
						 "Pval chr. breakp. enrich.",
						 "Pval exponential dist. breakpoints",
						 "Pval fragment joins",
						 "Links with other chrs")

	mytheme <- gridExtra::ttheme_minimal(padding = unit(c(1.8,1.8),"mm"),
										 core = list(fg_params=list(cex = .5)),
										 colhead = list(fg_params=list(cex = .5)),
										 rowhead = list(fg_params=list(cex = .5)))

	for (cc in 1:ncol(table_now)){if( is.numeric(table_now[,cc])){ table_now[,cc] = round(table_now[,cc],digits=2)  }}

	table_now = t(table_now) 
	colnames(table_now)= sample_name
	ss <- tableGrob(table_now, theme=mytheme)


	########################################
	gp1 <- ggplotGrob(ideogram + theme(plot.margin=unit(c(0.5,0.5,0,0),"cm")))
	gp2 <- ggplotGrob(SV_plot+ theme(legend.position = "none") +
					  theme(plot.margin=unit(c(0.5,0.5,0,0), "cm"),axis.text.x=element_blank())) 
	gp3 <- ggplotGrob(CNV_plot+theme(plot.margin=unit(c(0,0.5,0.2,0), "cm"),legend.position="none"))
	gp4 <- ss

	gp3$widths <- gp2$widths
	gp1$widths <- gp2$widths

	if (!is.null(BAF)){
		gp5 <- ggplotGrob(BAF_plot + theme(plot.margin=unit(c(0.5,0.5,0.5,0),"cm"),legend.position="none"))
		gp5$widths <- gp2$widths
		return(k=list(gp1,gp2,gp3,gp4,gp5))
	}else{
		return(k=list(gp1,gp2,gp3,gp4))
	}
	#return(k=list(ideogram,SV_plot,CN_plot,BAF_plot))
}



