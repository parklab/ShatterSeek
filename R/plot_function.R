# plot chromothripsis

plot_chromothripsis <- function(shatterSeek_output, chr=chr,BAF=NULL,sample_name=""){

	if ( !(as.character(chr) %in% chromNames)){stop("Chromosome not valid")}

	cand = chr
	chr=paste("chr",cand,sep="")

	common_ggplot2 <- theme_bw() + theme(axis.text.x=element_text(size=7,angle=0),
										 axis.text.y=element_text(size=7),
										 axis.title.y=element_text(size=7),
										 axis.title.x=element_blank(),
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
	#-----------------------------------------------------------------------------------------------------------------------------------

	summary = shatterSeek_output@chromSummary
	candidate_chrs <- shatterSeek_output@chromSummary$chrom 
	cluster_sizes <- sapply(shatterSeek_output@detail$connComp,length)
	# get the SVs for the cluster of SVs in this chromosome
	cand_clust_size <- shatterSeek_output@chromSummary$clusterSize[ shatterSeek_output@chromSummary$chrom == cand]
	idx = which(cluster_sizes == cand_clust_size) 
	SVsnow <- shatterSeek_output@detail$SV[ as.numeric(unlist(shatterSeek_output@detail$connComp[idx])) ,]
	CNVsnow <- shatterSeek_output@detail$CNV  
	SVsnow <- unique(SVsnow[SVsnow$chrom1 == cand, ]) # remove if there are more
	CNVsnow <- CNVsnow[CNVsnow$chrom == cand, ] # remove if there are more
	

	#if(nrow(SVsnow <1) | nrow(CNVsnow) < 1){stop("NO SVs or CN segments in this chromosome")}

	df = SVsnow

	min_coord = min(df$pos1)
	max_coord = max(df$pos2)
	idx1=which(CNVsnow$start >= min_coord)
	idx2=which(CNVsnow$end <= max_coord)   
	idx=intersect(idx1,idx2)
	CNVsnow = CNVsnow[idx,]

	#df$SVtype <- gsub("SVCLASS=","",df$SVtype)
	y1=4; y2=12
	df$y1 = rep(y1,nrow(df));    df$y2 = rep(y2,nrow(df))
	df$diff = abs(df$pos1 - df$pos2)
	df$curv = 1 -(df$diff / max(df$diff))
	max_diff = max(df$diff)
	df$curv[(df$diff / max_diff) > 0.2] <- .15
	df$curv[(df$diff / max_diff) > 0.8] <- .08
	df$curv[(df$diff / max_diff) < 0.2] <- 1
	d = data.frame(x=c(min_coord),y=c(1),leg=c("DEL","DUP","t2tINV","h2hINV"))
	idx = c()



	SV_plot = ggplot(d,aes(x=x,y=y)) +
	geom_point(colour="white") + ylim(0,y2+5) + common_ggplot2  +
	geom_line(data=rbind(d,d),aes(x=x,y=y,colour=leg)) 

	#-------------------------------------------------------------------------------------------------------
	# Interchromosomal SVs
	#-------------------------------------------------------------------------------------------------------
	inter2 <- shatterSeek_output@detail$SVinter   
	inter2 <- inter2[which(inter2$chrom1 == cand | inter2$chrom2 == cand), ] 


	if (nrow(inter2)>0){
		inter2$SVtype = factor(inter2$SVtype,levels=c("DEL","DUP","h2hINV","t2tINV"))
		inter2$y = rep(0,nrow(inter2))
		inter2$y[which(inter2$SVtype %in% c("DUP","DEL"))] = 4
		inter2$y[!(inter2$SVtype %in% c("DUP","DEL"))] = 12
		inter2$type_SV = rep("",nrow(inter2))
		inter2$type_SV[which(inter2$strand1 == "-" & inter2$strand2 == "-")] = "t2tINV"
		inter2$type_SV[which(inter2$strand1 == "-" & inter2$strand2 == "+")] = "DUP"
		inter2$type_SV[which(inter2$strand1 == "+" & inter2$strand2 == "-")] = "DEL"
		inter2$type_SV[which(inter2$strand1 == "+" & inter2$strand2 == "+")] = "h2hINV"
		inter2$SVtype = inter2$type_SV; inter2$type_SV=NULL
		inter2$colour = rep("",nrow(inter2))
		inter2$colour[which(inter2$SVtype == "DUP")] = "blue1"
		inter2$colour[which(inter2$SVtype == "DEL")] = "darkorange1"
		inter2$colour[which(inter2$SVtype == "h2hINV")] = "black"
		inter2$colour[which(inter2$SVtype == "t2tINV")] = "forestgreen"
		#
		inter2 = data.frame(pos = c(inter2$pos1,inter2$pos2), y=c(inter2$y,inter2$y), SVtype=c(inter2$SVtype,inter2$SVtype) )
		SV_plot = SV_plot + geom_point(data=inter2,size=1,alpha=1,aes(x=pos,y=as.numeric(y),colour=SVtype)) #+
	}
	#-------------------------------------------------------------------------------------------------------

	SV_plot = SV_plot + geom_hline(yintercept=y1,size=0.5) + geom_hline(yintercept=y2,size=0.5) 

	now = df[df$SVtype == "DUP",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=.2,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = now$curv[i],colour="blue1",ncp=8)#,square=T,squareShape=1) 
		}
		idx= c(idx,1)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))


	now = df[df$SVtype == "DEL",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=.2,data = now[i,] , aes(x = pos1, y = y1, xend = pos2, yend = y1), curvature = -1*now$curv[i],colour="darkorange1") 
		}
		idx= c(idx,2)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y1)) + geom_point(data=now,size=.5,aes(x=pos2,y=y1))

	now = df[df$SVtype == "t2tINV",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=.2,data = now[i,] , aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = now$curv[i],colour="forestgreen") 
		}
		idx= c(idx,3)
	}
	SV_plot = SV_plot + geom_point(data=now,size=.5,aes(x=pos1,y=y2)) + geom_point(data=now,size=.5,aes(x=pos2,y=y2))

	now = df[df$SVtype == "h2hINV",]
	if (nrow(now) > 0){
		for (i in 1:nrow(now)){
			SV_plot = SV_plot + geom_curve( size=.2,data = now[i,] , aes(x = pos1, y = y2, xend = pos2, yend = y2), curvature = -1*now$curv[i],colour="black") 
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
vals = c('DEL'='darkorange1','DUP'='blue1',"t2tINV"="forestgreen","h2hINV"="black")
labs = c('DEL','DUP',"t2tINV","h2hINV")

SV_plot = SV_plot +  scale_colour_manual(name = 'SV type', 
										 values =c('darkorange1','blue1',"forestgreen","black"),
										 labels = labs[idx]) + theme(legend.position="none")

## print gene information 
##if (!is.null(genes)){
##  ## create a data.frame similar to the CNV to plot the genes as bands on top of the CN 
##  for (k in 1:nrow(genes)){
##    idx = which(CNVsnow$start <= genes$start[k] & CNVsnow$end >= genes$end[k])
##    if(length(idx)!=0){
##    genes$total_cn[k] <- CNVsnow$total_cn[idx]} else{
##      genes$total_cn[k] <- NA
##    }
##  }
##
##  genes = genes[!is.na(genes$total_cn),]
##  
##  genes$position= genes$total_cn# +.3
##}
##################################################################################################

CN_plot = ggplot() +
geom_segment(data=CNVsnow, aes(x = start, y =total_cn , xend = end, yend = total_cn),size=2) + 
common_ggplot2 +#xlim(min_coord,max_coord) + 
#    geom_segment(data=genes, aes(x = start, y =total_cn , xend = end, yend = total_cn),size=2) + # colour=type
#    geom_text_repel(min.segment.length = unit(0.05, "lines"),point.padding = unit(1.6, 'lines'),arrow = arrow(length = unit(0.01, 'npc')),
#                    segment.alpha = .5,segment.size = .2,
#                    data=genes[genes$name=="MLH1",],aes(x=(start+end)/2,y=position,label=name),size=2,angle = 0,colour="blue",
#                    fontface="italic") +
ylab("CN")+ xlim(min_coord,max_coord) +  
scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")})+
coord_cartesian(xlim=c(min_coord,max_coord)) +
scale_y_continuous(minor_breaks = NULL,breaks=c(0,sort(unique(CNVsnow$total_cn))),
				   limits=c(min(CNVsnow$total_cn) -0.15,max(CNVsnow$total_cn)+.25)) +
xlab(NULL)
# tiene que estar xlim y luego expand

#########################################


common_ggplot2_chrom =  theme(panel.grid.major = element_blank(),
							  panel.grid.minor = element_blank(),
							  panel.border = element_blank(),
							  plot.background = element_blank(),
							  axis.text=element_blank(),
							  axis.title=element_blank(),
							  plot.title=element_blank(), 
							  axis.ticks=element_blank())


##chr_info = readRDS("/Users/icortes/Dropbox/Park_lab/paper_chromothripsis/chr_info.rds")
chr_info$color[chr_info$gieStain == "gneg"] = "white"
chr_info$color[chr_info$gieStain == "gpos25"] = "grey75"
chr_info$color[chr_info$gieStain == "gpos50"] = "grey50"
chr_info$color[chr_info$gieStain == "gpos75"] = "grey25"
chr_info$color[chr_info$gieStain == "gpos100"] = "grey0"
chr_info$color[chr_info$gieStain == "acen"] = "red"

chr_info = chr_info[chr_info$seqnames==chr,]

chr_info$y = rep(1,nrow(chr_info))

## XXX check
chr_info_annot=chr_info[seq(3,(nrow(chr_info)-3),5),]

ideogram =ggplot(data=chr_info,aes(x=y,y=y)) +geom_point(colour="white")
ideogram = ideogram +  
geom_rect(data=chr_info,mapping = aes(xmin = chr_info$start, xmax = chr_info$end,   ymin = 0, ymax = 1),
		  fill = chr_info$color,color="black",size=.2) 

ideogram = ideogram +ylim(0,4)
ideogram = ideogram + annotate(geom = "text",x = (chr_info_annot$start+chr_info_annot$end)/2,y=chr_info_annot$y + 1,
							   label=chr_info_annot$name,angle=35,size=2)
ideogram = ideogram +theme_bw() + common_ggplot2_chrom + #xlim(min_coord,max_coord) + 
scale_x_continuous(expand = c(0.01,0.01)) +
coord_cartesian(xlim=c(min_coord, max_coord))


#----------------------------------------------------------------------
# BAF
#----------------------------------------------------------------------
#XXX check input data and hence aes 
if (!is.null(BAF)){
	BAF$ratio1 = BAF$b_ads1/BAF$b_dp

	BAF_plot = ggplot(BAF[BAF$b_dp > 20,],aes(x=pos,y=ratio1,colour=ratio1)) + 
	geom_point(size=.1,alpha=.1,colour="black")+
	#scale_colour_gradient2(mid="orange", high="black", low="black") +
	common_ggplot2 + #xlim(min_coord,max_coord) + 
	ylab("BAF")+ xlim(min_coord,max_coord) +  
	scale_x_continuous(expand = c(0.01,0.01),labels = function(x){paste(x/1000000,"MB")},limits=c(min_coord,max_coord)) +
	scale_y_continuous(breaks=c(0,0.5,1),labels=c("0","0.5","1"),limits=c(0,1))  

}

#----------------------------------------------------------------------
# Table 
#----------------------------------------------------------------------


#   [1] "chrom"                                               "clusterSize"                                        
#   [3] "start"                                               "end"                                                
#   [5] "number_DEL"                                          "number_DUP"                                         
#   [7] "number_h2hINV"                                       "number_t2tINV"                                      
#   [9] "number_TRA"                                          "clusterSize_including_TRA"                          
#   [11] "number_SVs_sample"                                   "PGA"                                                
#   [13] "CNV_peaks"                                           "number_CNV_segments"                                
#   [15] "fragment_joints"                                     "chr_breakpoint_enrichment"                          
#   [17] "pval_exp_chr"                                        "pval_exp_cluster"                                   
#   [19] "number_CN_segments"                                  "max_number_oscillating_CN_segments_2_states"        
#   [21] "max_number_oscillating_CN_segments_2_states_3states" "number_CN_segments_chr"                             
#   [23] "max_number_oscillating_CN_segments_2_states_chr"     "max_number_oscillating_CN_segments_3_states_chr"    
#   [25] "inter_main_chrom"                                    "inter_number_DEL"                                   
#   [27] "inter_number_h2hINV"                                 "inter_number_t2tINV"                                
#   [29] "inter_number_TRA"                                    "inter_number_DUP"                                   
#   [31] "inter_fragment_joints"                               "inter_other_chroms"                                 
#   [33] "inter_other_chroms_coords_all"                      



summary$pos = paste(summary$chrom,":",summary$start,"-",summary$end,sep="")
summary$oscillations = paste(summary$max_number_oscillating_CN_segments_2_states,";",summary$max_number_oscillating_CN_segments_2_states_3states,sep="")

cols_sel = c("pos",
		 "clusterSize",
		 "clusterSize_including_TRA",
		 "number_SVs_sample",
		 "oscillations",
		 "number_CNV_segments",
		 "chr_breakpoint_enrichment",
		 "pval_exp_cluster",
		 "fragment_joints",
		 "inter_other_chroms_coords_all")
print(summary)
table_now = summary[which(summary$chrom == cand),cols_sel]

names(table_now) = c("Position",
					 "Nb. Intrachr. SVs", 
					 "Total nb. SVs (intrachr. + transl.)",
					 "SVs in sample",
					 "Oscillating CN (2 and 3 states)","CN segments",
					 "Pval chr. breakp. enrich.",
					 "Pval exponential dist. breakpoints",
					 "Pval fragment joints",
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
gp2 <- ggplotGrob(SV_plot+ theme(legend.position = c(.5,1),legend.direction = "horizontal",legend.background=element_blank()) +
				  theme(plot.margin=unit(c(0.5,0.5,0,0), "cm"),axis.text.x=element_blank())) 
gp3 <- ggplotGrob(CN_plot+theme(plot.margin=unit(c(0,0.5,0.2,0), "cm"),legend.position="none"))
gp4 <- ss #ggplotGrob(ss + theme(plot.margin=unit(c(0.5,0.5,0,0),"cm")))

gp3$widths <- gp2$widths
gp1$widths <- gp2$widths
#gp4$widths <- gp2$widths


if (!is.null(BAF)){
gp5 <- ggplotGrob(BAF_plot + theme(plot.margin=unit(c(0.5,0.5,0.5,0),"cm"),legend.position="none"))
gp5$widths <- gp2$widths
return(k=list(gp1,gp2,gp3,gp4,gp5))
}else{
return(k=list(gp1,gp2,gp3,gp4))
}

#return(k=list(ideogram,SV_plot,CN_plot,BAF_plot))

}
