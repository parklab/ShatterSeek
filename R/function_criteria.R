
statistical_criteria = function(input){

    summary <- input@chromSummary$chrom  
    summary$start <- rep(NA,23)
    summary$end <- rep(NA,23)
    summary$number_DEL <- rep(0,23)
    summary$number_DUP <- rep(0,23)
    summary$number_h2hINV <- rep(0,23)
    summary$number_t2tINV <- rep(0,23)
    summary$number_TRA <- rep(0,23)
    summary$clusterSize_including_TRA <- rep(0,23)
    summary$number_SVs_sample <- rep(0,23)
    summary$PGA <- rep(0,23)
    #summary$PGA_region <- rep(0,23)

    summary$CNV_peaks <- rep(NA,23)
    summary$number_CNV_segments <- rep(NA,23)

    summary$fragment_joints <- rep(NA,23)
    summary$chr_breakpoint_enrichment <- rep(NA,23)
    summary$pval_exp_chr <- rep(NA,23)
    summary$pval_exp_cluster <- rep(NA,23)
    #summary$canonical_score <- rep(NA,23)
    #summary$i_i2 <- rep(NA,23)

    # oscillating
    summary$number_CN_segments <- rep(NA,23)
    summary$max_number_oscillating_CN_segments_2_states <- rep(NA,23)
    summary$max_number_oscillating_CN_segments_2_states_3states <- rep(NA,23)
    #summary$i_i2_pval <- rep(NA,23)
    #summary$i_i2_pval_max <- rep(NA,23)
    #summary$i_i2_cutoff1 <- rep(NA,23)
    #summary$i_i2_sequential <- rep(NA,23)
    #summary$i_i2_sequential_8 <- rep(NA,23)
    #summary$i_i2_sequential_8_3states <- rep(NA,23)

    #summary$i_i2_chr <- rep(NA,23)
    #summary$i_i2_cutoff1_chr <- rep(NA,23)
    summary$number_CN_segments_chr <- rep(NA,23)
    summary$max_number_oscillating_CN_segments_2_states_chr <- rep(NA,23)
    summary$max_number_oscillating_CN_segments_3_states_chr <- rep(NA,23)

    #summary$CN_0 <- rep(0,23)
    #summary$CN_no_0 <- rep(0,23)
    #summary$PGNA <- rep(0,23)
    #summary$mode_CN_1 <- rep(0,23)
    #summary$mode_CN_2 <- rep(0,23)


    ####----------------------------------------------------------------------------------------------------------------------------------------------
    candidate_chrs <- input@chromSummary$chrom  
    # contains the IDs for the SVs that are in a cluster.
    cluster_sizes <- sapply(input@detail$connComp,length)
    ####----------------------------------------------------------------------------------------------------------------------------------------------
    ####----------------------------------------------------------------------------------------------------------------------------------------------
    summary_inter  = data.frame(main_chrom=rep(0,length(candidate_chrs)),
                                #other_chroms=rep("pepe",length(candidate_chrs)),
                                number_DEL=rep(0,length(candidate_chrs)),
                                number_h2hINV=rep(0,length(candidate_chrs)),
                                number_t2tINV=rep(0,length(candidate_chrs)),
                                number_TRA=rep(0,length(candidate_chrs)),
                                number_DUP=rep(0,length(candidate_chrs)),
                                #nb_SVs=rep(0,length(candidate_chrs)),
                                #CNV_peaks_sum=rep(0,length(candidate_chrs)),
                                #CNVpeaks=rep(0,length(candidate_chrs)),
                                #i_i2_sequential_8_3states=rep(0,length(candidate_chrs)),
                                #i_i2_sequential_8=rep(0,length(candidate_chrs)),
                                fragment_joints=rep(0,length(candidate_chrs)))

    other_chroms=rep(" ",length(candidate_chrs))
    other_chroms_coords_all = rep(" ",length(candidate_chrs))



    for (cand in candidate_chrs){
        separator_index=1
        cand_clust_size <- input@chromSummary$clusterSize[ input@chromSummary$chrom == cand]
        idx = which(cluster_sizes == cand_clust_size) #ojo, puede haber clusts del mismo size en chrs diferentes
        SVsnow <- input@detail$SV[ as.numeric(unlist(input@detail$connComp[idx])) ,]
        indice_chromosome = which(summary$chrom == cand)
        SVsnow <- SVsnow[SVsnow$chrom1 == cand, ] # remove if there are more

        # get the copy number data
        CNVsnow <- input@detail$CNV
        # get PGA
        CNVsnow_mapp = CNVsnow[which(CNVsnow$total_cn != 2),]
        summary$PGA[indice_chromosome] = 100 * (sum(as.numeric(CNVsnow_mapp$end - CNVsnow_mapp$start)) / sizes_tot)
        #CNVsnow_mapp = CNVsnow[which(CNVsnow$total_cn != 2 & CNVsnow$chrom == cand),]
        #summary$PGA_region[indice_chromosome] = 100 * (sum(CNVsnow_mapp$end - CNVsnow_mapp$start) / sizes$V2[which(sizes$V1 == cand)])
        #CNVsnow_mapp = CNVsnow[which(CNVsnow$total_cn == 2 & CNVsnow$chrom == cand),]

        #summary$PGNA[indice_chromosome] = 100 * (sum(CNVsnow_mapp$end - CNVsnow_mapp$start) / sizes$V2[which(sizes$V1 == cand)])
        # get CNV data for the current chromosome
        CNVsnow = CNVsnow[which(CNVsnow$chrom == cand),]

        # Oscillations at the chromosome level
        if (nrow(CNVsnow) >= 4){
            k = nrow(CNVsnow)
            v = (CNVsnow$total_cn[1:(k-2)] - CNVsnow$total_cn[3:k])
            v_3states = 1*(abs(v) %in% c(0,1))
            v = 1*(v ==0)
            #i_i2 = fmax(1,v)
            number_CN_segments = length(v)
            #i_i2_cutoff1=fmax(1,v,cutoff=1); 
            max_number_oscillating_CN_segments_2_states = fmaxmax(1,v)
            max_number_oscillating_CN_segments_2_states_3states = fmaxmax(1,v_3states)
            #summary$i_i2_chr[indice_chromosome] = i_i2
            #summary$i_i2_cutoff1_chr[indice_chromosome] = i_i2_cutoff1
            summary$number_CN_segments_chr[indice_chromosome] = number_CN_segments
            summary$max_number_oscillating_CN_segments_2_states_chr[indice_chromosome] = max_number_oscillating_CN_segments_2_states
            summary$max_number_oscillating_CN_segments_3_states_chr[indice_chromosome] = max_number_oscillating_CN_segments_2_states_3states
        }



        if(nrow(SVsnow) !=0){
            min_now = as.numeric(min(c(SVsnow$pos1,SVsnow$pos2)))
            max_now = max(c(SVsnow$pos1,SVsnow$pos2))
            summary$start[indice_chromosome] <- min_now  
            summary$end[indice_chromosome] <- max_now  

            if(nrow(CNVsnow) > 0){
                idxa = which(CNVsnow$start <= summary$start[indice_chromosome])
                idxb = which(CNVsnow$end >= summary$end[indice_chromosome])
                if(length(idxa) ==0){idxa = min(which(CNVsnow$start >= summary$start[indice_chromosome])) }
                if(length(idxb) ==0){idxb = max(which(CNVsnow$end <= summary$end[indice_chromosome]))}

                if(length(idxa) !=0 & length(idxb) !=0  & (max(idxa)-min(idxb))<=0 ){
                    CNVsnow = CNVsnow[seq(max(idxa),min(idxb),1),]

                    summary$CNV_peaks[indice_chromosome] <- peak.count(CNVsnow$total_cn)
                    i_i2_sequential = 0
                    if (nrow(CNVsnow) >= 4){
                        k = nrow(CNVsnow)
                        v = (CNVsnow$total_cn[1:(k-2)] - CNVsnow$total_cn[3:k])
                        v_3states = 1*(abs(v) %in% c(0,1))
                        v = 1*(v ==0)
                        max_number_oscillating_CN_segments_2_states_3states = fmaxmax(1,v_3states)
                        #i_i2 = fmax(1,v)
                        number_CN_segments = length(v)
                        #i_i2_cutoff1=fmax(1,v,cutoff=1); 
                        max_number_oscillating_CN_segments_2_states = fmaxmax(1,v)

                        #########################################################

                        #i_i2_pval = fmaxsample(CNVsnow, cutoff=2,times=5000)
                        #i_i2_pval_max = fmaxsampleb(CNVsnow, times=5000)
                        #summary$i_i2[indice_chromosome] = i_i2
                        summary$max_number_oscillating_CN_segments_2_states_3states[indice_chromosome] = max_number_oscillating_CN_segments_2_states_3states
                        #summary$i_i2_cutoff1[indice_chromosome] = i_i2_cutoff1
                        summary$number_CN_segments[indice_chromosome] = number_CN_segments
                        summary$max_number_oscillating_CN_segments_2_states[indice_chromosome] = max_number_oscillating_CN_segments_2_states
                        #summary$i_i2_pval[indice_chromosome] = i_i2_pval
                        #summary$i_i2_pval_max[indice_chromosome] = i_i2_pval_max
                        #summary$i_i2_sequential[indice_chromosome] = max(i_i2_sequential)
                    }

                    #-----------------------------------------------------
                    # check for candidates with CN=0
                    #-----------------------------------------------------
                    #CN_zero = sum(as.numeric(CNVsnow$end[which(CNVsnow$total_cn == 0)] - CNVsnow$start[which(CNVsnow$total_cn == 0)]))
                    #CN_no_zero = sum(as.numeric(CNVsnow$end[which(CNVsnow$total_cn != 0)] - CNVsnow$start[which(CNVsnow$total_cn != 0)]))
                    #summary$CN_0[indice_chromosome] = CN_zero
                    #summary$CN_no_0[indice_chromosome] = CN_no_zero

                    #-----------------------------------------------------
                    # Assign chromo category (canonical, canonical with 
                    #-----------------------------------------------------
                    #summary$canonical_score[indice_chromosome] = sum(CNVsnow$total_cn%in%c(1,2,3))
                    max_CN = max(CNVsnow$total_cn)
                    min_CN = min(CNVsnow$total_cn)
                    modes_states = as.numeric(names(sort(-table(CNVsnow$total_cn)))[1:2])
                    summary$number_CNV_segments[indice_chromosome] = nrow(CNVsnow)
                    #summary$mode_CN_1[indice_chromosome] = modes_states[1]
                    #summary$mode_CN_2[indice_chromosome] = modes_states[2]
                }

            } 
        }

        #-----------------------------------------------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------------------------------------------
        cat("inter\n\n")
        #-----------------------------------------------------
        # Look for interchromosomal events
        #-----------------------------------------------------
        SVs_all <- rbind(input@detail$SV,input@detail$SVinter)  #SV.sample_backup
        #tryCatch(read.csv(SV_Sample,sep="\t",header=T,stringsAsFactors=F), error=function(e) NULL)
        #XXX should take from summary
        # define candidate interchromosomal regions
        inter = SVs_all[which(SVs_all$chrom1 != SVs_all$chrom2),]
        print(head(inter))

        if(nrow(SVsnow) !=0 & nrow(inter) > 0){
            cat("define windows\n\n")
            # get whether any tranlocation involves the cand chromo and is within the cluster of SVs previously identified
            min_now = as.numeric(min(c(SVsnow$pos1,SVsnow$pos2)))
            max_now = max(c(SVsnow$pos1,SVsnow$pos2))
            idx_inter1 = which(inter$chrom1 == cand & inter$pos1 >= (min_now - 10000) & inter$pos1 <= (max_now + 10000))
            idx_inter2 = which(inter$chrom2 == cand & inter$pos1 >= (min_now - 10000) & inter$pos2 <= (max_now + 10000))
            #idx_inter1 = which(inter$chrom1 == cand & inter$start1 >= (min_now - 10000) & inter$end1 <= (max_now + 10000))
            #idx_inter2 = which(inter$chrom2 == cand & inter$start2 >= (min_now - 10000) & inter$end2 <= (max_now + 10000))
            idx_inter = c(idx_inter1,idx_inter2)
            window = c() # save the SVs in the window comprising multiple chrs
            CNV_window =c()
            selection_chrs = c()
            nb_TRA_tot = c()
            selection_chr_coords = c()

            if (length(idx_inter) >= 1){  #es deecir, hay translocations en la zona del cluster para cand
                inter = inter[idx_inter,]
                cand_inter = unique(c(inter$chrom1,inter$chrom2))
                cand_inter = cand_inter[which(cand_inter != cand)]

                for (chr_inter in sort(cand_inter)){ # para cada chr con translocation, mirar si la tranlocation toca un cluster en el otro chr (del tamano que sea y sin importar que fuera significativo en los test a nivel de chromosoma)
                    # get the SV cluster for the other partner (chr) in the translocation
                    cand_clust_size_window <- input@chromSummary$clusterSize[ input@chromSummary$chrom == chr_inter]
                    idx = which(cluster_sizes == cand_clust_size_window) #ojo, puede haber clusts del mismo size en chrs diferentes
                    SVsnow_window <- input@detail$SV[ as.numeric(unlist(input@detail$connComp[idx])) ,]
                    coords_window = paste(SVsnow_window$chrom1,SVsnow_window$pos1, SVsnow_window$chrom2, SVsnow_window$pos2)#-1)
                    coords_window_from_SVs_all = paste(SVs_all$chrom1,SVs_all$pos1, SVs_all$chrom2, SVs_all$pos2)
                    #                coords_window_from_SVs_all = paste(SVs_all$chrom1,SVs_all$start1, SVs_all$chrom2, SVs_all$start2)

                    idx = which(coords_window_from_SVs_all %in% coords_window)
                    SVsnow_window = SVs_all[idx,]#[,c(1,2,4,5,9,10,11)]
                    #names(SVsnow_window)[1:4] = c("chrom1","pos1","chrom2","pos2")

                    min_now_window = as.numeric(min(c(SVsnow_window$pos1[which(SVsnow_window$chrom1 == chr_inter)],SVsnow_window$pos2[which(SVsnow_window$chrom2 == chr_inter)])))
                    max_now_window = max(c(SVsnow_window$pos1[which(SVsnow_window$chrom1 == chr_inter)],SVsnow_window$pos2[which(SVsnow_window$chrom2 == chr_inter)]))
                    # ver si las translocations estan dentro del cluster of SVs
                    inter_chr_inter1 = inter[which(inter$chrom1 == chr_inter),1:2]#1:3]
                    inter_chr_inter2 = inter[which(inter$chrom2 == chr_inter),3:4]#4:6]
                    names(inter_chr_inter1) = c("chr","pos")
                    names(inter_chr_inter2) = c("chr","pos")
                    inter_chr_inter = rbind(inter_chr_inter1,inter_chr_inter2)
                    print(inter_chr_inter)
                    min_trans = min(c(inter_chr_inter$pos))
                    max_trans = max(c(inter_chr_inter$pos))
                    #min_trans = min(c(inter_chr_inter$pos1,inter_chr_inter$pos2))
                    #max_trans = max(c(inter_chr_inter$pos1,inter_chr_inter$pos2))
                    #min_trans = min(c(inter_chr_inter$start1,inter_chr_inter$start2))
                    #max_trans = max(c(inter_chr_inter$start1,inter_chr_inter$start2))

                    if ( (min_trans >= (min_now_window-10000) | max_trans <= (max_now_window+10000)) & nrow(inter_chr_inter) >=2 ){ #the translocation is within the cluster in the other chr    ;    at least 2 TRA
                        # move down window = rbind(window,SVsnow_window)
                        cat("oeoeoeoeoeoe\n\n")
                        CNV_window_now <- input@detail$CNV
                        idxa = which(CNV_window_now$chrom == chr_inter & CNV_window_now$start <= (min_now_window -10000))
                        if (length(idxa)==0){idxa = which(CNV_window_now$chrom == chr_inter & CNV_window_now$start <= (min_now_window))}
                        idxb = which(CNV_window_now$chrom == chr_inter & CNV_window_now$end >= (max_now_window +10000))
                        if (length(idxb)==0){idxb = which(CNV_window_now$chrom == chr_inter & CNV_window_now$end >= (max_now_window +10000))}
                        if (length(idxb)==0){idxb = max( which(CNV_window_now$chrom == chr_inter & CNV_window_now$end <= (max_now_window )) )}

                        if (length(idxb)!=0 & length(idxa)!=0 ){
                            nb_oscill = fmax2(CNV_window_now[seq(max(idxa),min(idxb),1),],cutoff=1)
                            if (nb_oscill >= 0){  
                                CNV_window = rbind(CNV_window,CNV_window_now[seq(max(idxa),min(idxb),1),])
                                selection_chrs = c(selection_chrs,chr_inter)
                                if (separator_index%%2 == 0){separator ="\n"}else{separator=";"}
                                separator_index = separator_index +1
                                selection_chr_coords= c(selection_chr_coords,paste(chr_inter,":",min_now_window,"-",max_now_window,separator,sep=""))
                                window = rbind(window,SVsnow_window)
                            }
                        }

                    }
                }
                cat("wind")
                print(window)
                if (!is.null(window)){
                    inter = inter#[,c(1,2,4,5,9,10,11)]
                    names(inter)[1:4] = c("chrom1","pos1","chrom2","pos2")
                    window = rbind(window,inter[which(inter$chrom1 %in% selection_chrs | inter$chrom2 %in% selection_chrs),])
                    cat("window\n")
                    print(window)
                    #---------------------------------------------
                    # Randomness of joints
                    #---------------------------------------------
                    # check multinomial distribution for SV types
                    # add SV types for translocations
                    #window$SVtype = window$svclass
                    window$SVtype2 = window$SVtype
                    nb_TRA = length(which(window$SVtype == "TRA"))
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "+" & window$strand2 == "-")] = "DEL"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "+" & window$strand2 == "+")] = "h2hINV"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "-" & window$strand2 == "-")] = "t2tINV"
                    window$SVtype2[which(window$SVtype == "TRA" & window$strand1 == "-" & window$strand2 == "+")] = "DUP"
                    obs = c(0,0,0,0)
                    names(obs) <- c("del","inv","inv2","dup")
                    obs[1] <- sum(window$SVtype2 == "DEL")
                    obs[2] <- sum(window$SVtype2 %in% c("h2hINV"))
                    obs[3] <- sum(window$SVtype2 %in% c("t2tINV"))
                    obs[4] <- sum(window$SVtype2 == "DUP")

                    idxx = which(candidate_chrs == cand)

                    summary_inter$number_DEL[idxx] = obs[1]
                    summary_inter$number_h2hINV[idxx] = obs[2]
                    summary_inter$number_t2tINV[idxx] = obs[3]
                    summary_inter$number_DUP[idxx] = obs[4]
                    summary_inter$number_TRA[idxx] = nb_TRA

                    summary_inter$main_chrom[idxx] = cand
                    other_chroms[idxx] = as.character(paste(as.vector(selection_chrs),collapse="_"))
                    other_chroms_coords_all[idxx] = as.character(paste(as.vector(selection_chr_coords),collapse=""))

                    signif <- chisq.test(obs, p=rep(1/4,4))$p.val
                    summary_inter$fragment_joints[idxx] <- signif 
                    cat('kiol\n\n')

                    #---------------------------------------------
                    # CNV peaks
                    #---------------------------------------------
                    library(plyr)
                    print("peaks")
                    # add data from cand
                    CNV_window_now <- input@detail$CNV
                    idxa = which(CNV_window_now$chrom == cand & CNV_window_now$start <= (min_now- 10000))
                    if(length(idxa)==0){idxa = which(CNV_window_now$chrom == cand & CNV_window_now$start <= (min_now))}
                    if(length(idxa)==0){idxa = min(which(CNV_window_now$chrom == cand & CNV_window_now$start >= (min_now)))}

                    idxb = which(CNV_window_now$chrom == cand & CNV_window_now$end >= (max_now +10000))
                    if(length(idxb)==0){idxb = which(CNV_window_now$chrom == cand & CNV_window_now$end >= (max_now))}
                    if(length(idxb)==0){idxb = max(which(CNV_window_now$chrom == cand & CNV_window_now$end <= (max_now)))}

                    if(length(idxb)!=0 &length(idxa)!=0 & !is.na(idxb) & !is.na(idxa)){
                        CNV_window = rbind(CNV_window,CNV_window_now[seq(max(idxa),min(idxb),1),])
                        peaks = ddply(CNV_window,.(chrom),summarize,p=peak.count(total_cn))
                        #summary_inter$CNVpeaks[idxx] <- paste(peaks$p,collapse="_") #peak.count(CNV_window$total_cn)
                        #summary_inter$CNV_peaks_sum[idxx] <- sum(peaks$p)
                        #summary_inter$nb_SVs[idxx] <- nrow(window)

                        #---------------------------------------------
                        # Sequential segments inter
                        #---------------------------------------------
                        i_i2_sequential = 0
                        if (nrow(CNV_window) >=4){
                            v = (CNV_window$total_cn[1:(nrow(CNV_window)-2)] - CNV_window$total_cn[3:nrow(CNV_window)])
                            if (length(v) >= 4){
                                v = 1*(v ==0)
                                i_i2 = fmax(1,v)
                                i_i2_sequential = fmax(0,v)
                            }else{

                                i_i2 = 0; i_i2_sequential=0
                            }
                        }
                        #summary_inter$i_i2_sequential_8[indice_chromosome] = max(i_i2_sequential)
                        print("seq inter")
                        print(i_i2_sequential)
                        #-----------------------------------------------------
                        # check that i == i+2. Purity of the CNV profile  8  ALLOWING FOR THREE STATES
                        #-----------------------------------------------------
                        i_i2_sequential = 0
                        if (nrow(CNV_window) >= 6){
                            v = (CNV_window$total_cn[1:(nrow(CNV_window)-2)] - CNV_window$total_cn[3:nrow(CNV_window)])
                            if (length(v) >= 6){
                                v = 1*(abs(v) %in%c(0,1))
                                i_i2_sequential = c()
                                for(csum in 1:(length(v)-5)){ i_i2_sequential = c(i_i2_sequential, cumsum(v[csum:(csum+5)])[6])}
                            }
                            i_i2 = sum( (CNV_window$total_cn[1:(nrow(CNV_window)-2)] - CNV_window$total_cn[3:nrow(CNV_window)]) == 0)} else{
                                i_i2 = 0; i_i2_sequential=0
                        }
                        #summary_inter$i_i2_sequential_8_3states[indice_chromosome] = max(i_i2_sequential)


                    }
                }
            }
        }

        summary_inter$other_chroms = other_chroms
        summary_inter$other_chroms_coords_all = other_chroms_coords_all
        #summary_inter$CNVpeaks = CNVpeaks

        cat("end inter\n\n\n")
        #-----------------------------------------------------
        # Exponential distribution (CLUSTER SVS)
        #-----------------------------------------------------
        # the mean of an exp distribution is equal to 1/lambda (lamba ==rate)
        SVsnow_exp <-  rbind(input@detail$SV,input@detail$SVinter)  #SV.sample_backup 
        #SV.sample_backup #tryCatch(read.csv(SV_Sample,sep="\t",header=T,stringsAsFactors=F), error=function(e) NULL)    
        ## ojo, en este es start y en shatter pos1 .
        minnow = summary$start[indice_chromosome]
        maxnow = summary$end[indice_chromosome]

        idx_inter1 = which(SVsnow_exp$chrom1 == cand & SVsnow_exp$pos1 >= (minnow) & SVsnow_exp$pos1 <= (maxnow))
        #    idx_inter1 = which(SVsnow_exp$chrom1 == cand & SVsnow_exp$start1 >= (minnow) & SVsnow_exp$end1 <= (maxnow))
        breaks_inter1= SVsnow_exp$start1[idx_inter1]
        idx_inter2 = which(SVsnow_exp$chrom2 == cand & SVsnow_exp$pos2 >= (minnow) & SVsnow_exp$pos2 <= (maxnow))
        #    idx_inter2 = which(SVsnow_exp$chrom2 == cand & SVsnow_exp$start2 >= (minnow) & SVsnow_exp$end2 <= (maxnow))
        breaks_inter2= SVsnow_exp$pos2[idx_inter2]
        #    breaks_inter2= SVsnow_exp$start2[idx_inter2]
        # add these breakpoints to breaks below

        # not take intra..
        SVsnow_exp <- SVsnow_exp[which(SVsnow_exp$chrom1 == SVsnow_exp$chrom2),]  ##input@detail$SV
        cand_clust_size <- input@chromSummary$clusterSize[ input@chromSummary$chrom == cand]
        idx = which(cluster_sizes == cand_clust_size) #ojo, puede haber clusts del mismo size en chrs diferentes
        SVsnow_exp <- input@detail$SV[ as.numeric(unlist(input@detail$connComp[idx])) ,]
        SVsnow_exp <- SVsnow_exp[SVsnow_exp$chrom1 == cand, ] # remove if there are more

        breaks = sort(c(SVsnow_exp$pos1,SVsnow_exp$pos2, breaks_inter1, breaks_inter2))

        if (length(breaks) >= 6){
            rate = 1 / ( sum(as.numeric(breaks[2:length(breaks)] - breaks[1:(length(breaks)-1)])) / (length(breaks)-1) )
            exponential <- rexp(n=(length(breaks)-1), rate = rate) # generate exponential distribution
            # estimate the parameters
            fit1 <- fitdistr(exponential, "exponential") 
            # goodness of fit test
            pval_exp = ks.test(breaks, "pexp", fit1$estimate)$p.val
            summary$pval_exp_cluster[indice_chromosome] = pval_exp} else{
                summary$pval_exp_cluster[indice_chromosome] = NA}

        #-----------------------------------------------------
        # Exponential distribution to test clustering of breakpoints in a given chromosome 
        #-----------------------------------------------------
        # the mean of an exp distribution is equal to 1/lambda (lamba ==rate)

        SVsnow_exp <- rbind(input@detail$SV,input@detail$SVinter)  #SV.sample_backup
        #input@detail$SV #SV.sample_backup # tryCatch(read.csv(SV_Sample,sep="\t",header=T,stringsAsFactors=F), error=function(e) NULL)    

        idx_inter1 = which(SVsnow_exp$chrom1 == cand) # & SVsnow_exp$start1 >= (minnow) & SVsnow_exp$end1 <= (maxnow))
        breaks1= SVsnow_exp$pos2[idx_inter1]
        #    breaks1= SVsnow_exp$start1[idx_inter1]
        idx_inter2 = which(SVsnow_exp$chrom2 == cand) # & SVsnow_exp$start2 >= (minnow) & SVsnow_exp$end2 <= (maxnow))
        breaks2= SVsnow_exp$pos2[idx_inter2]
        #    breaks2= SVsnow_exp$start2[idx_inter2]
        breaks = sort(unique(c(breaks1,breaks2)))

        if (length(breaks) >= 6){
            rate = 1 / ( sum(as.numeric(breaks[2:length(breaks)] - breaks[1:(length(breaks)-1)])) / (length(breaks)-1) )
            exponential <- rexp(n=(length(breaks)-1), rate = rate) # generate some exponential distribution
            # estimate the parameters
            fit1 <- fitdistr(exponential, "exponential") 
            # goodness of fit test
            pval_exp = ks.test(breaks, "pexp", fit1$estimate)$p.val
            summary$pval_exp_chr[indice_chromosome] = pval_exp} else{
                summary$pval_exp_chr[indice_chromosome] =NA
        }

        #-----------------------------------------------------
        # Are there more SVs in a chrs than expected by chance?
        #-----------------------------------------------------
        # ojo, no considerando interchromo events
        SVsnow_exp <- rbind(input@detail$SV,input@detail$SVinter)  #SV.sample_backup
        #input@detail$SV # SV.sample_backup #tryCatch(read.csv(SV_Sample,sep="\t",header=T,stringsAsFactors=F), error=function(e) NULL)    
        idx_inter1 = which(SVsnow_exp$chrom1 == cand) 
        idx_inter2 = which(SVsnow_exp$chrom2 == cand) 
        nb_SVs_cand = SVsnow_exp[unique(c(idx_inter1,idx_inter2)),]

        if(nrow(nb_SVs_cand) !=0){
            nb_SVs_cand = nrow(nb_SVs_cand)
            nb_SVs_all_sample = nrow(SVsnow_exp) #input@detail$SV)
            prob_cand = info_mappa$tot[which(info_mappa$V1 == cand)] / sum(as.numeric(info_mappa$tot))
            chr_enrich = binom.test(nb_SVs_cand,nb_SVs_all_sample,p=prob_cand)$p.val
            summary$chr_breakpoint_enrichment[indice_chromosome] = chr_enrich} else {
                summary$chr_breakpoint_enrichment[indice_chromosome] = NA
        } 
        #--------------------------------------------------
        # Randomness of DNA fragment joints
        #--------------------------------------------------
        # check multinomial distribution for SV types
        cand_clust_size <- input@chromSummary$clusterSize[ input@chromSummary$chrom == cand]
        idx = which(cluster_sizes == cand_clust_size) #ojo, puede haber clusts del mismo size en chrs diferentes
        SVsnow <- input@detail$SV[ as.numeric(unlist(input@detail$connComp[idx])) ,]
        SVsnow <- SVsnow[SVsnow$chrom1 == cand, ] # remove if there are more

        obs = c(0,0,0,0)
        names(obs) <- c("del","inv","inv2","dup")
        obs[1] <- sum(SVsnow$SVtype == "DEL")
        obs[2] <- sum(SVsnow$SVtype %in% c("h2hINV"))
        obs[3] <- sum(SVsnow$SVtype %in% c("t2tINV"))
        obs[4] <- sum(SVsnow$SVtype == "DUP")
        # get the number of translocations

        SVs_all <-rbind(input@detail$SV,input@detail$SVinter)  #SV.sample_backup
        #input@detail$SV #SV.sample_backup # tryCatch(read.csv(SV_Sample,sep="\t",header=T,stringsAsFactors=F), error=function(e) NULL)
        summary$number_SVs_sample[indice_chromosome] = nrow(SVs_all)

        inter = SVs_all[which(SVs_all$chrom1 != SVs_all$chrom2),]
        obs2 = obs
        if(nrow(inter) > 0){
            # get whether any tranlocation involves the cand chromo and is within the cluster of SVs previously identified
            min_now = as.numeric(min(c(SVsnow$pos1,SVsnow$pos2)))
            max_now = max(c(SVsnow$pos1,SVsnow$pos2))
            idx_inter1 = which(inter$chrom1 == cand & inter$pos1 >= (min_now) & inter$pos1 <= (max_now))
            idx_inter2 = which(inter$chrom2 == cand & inter$pos2 >= (min_now) & inter$pos2 <= (max_now))
            #idx_inter1 = which(inter$chrom1 == cand & inter$start1 >= (min_now) & inter$end1 <= (max_now))
            #idx_inter2 = which(inter$chrom2 == cand & inter$start2 >= (min_now) & inter$end2 <= (max_now))
            idx_inter = unique(c(idx_inter1,idx_inter2))
            inter = inter[idx_inter, ]

            obs2[1] = obs2[1] + length( which(inter$strand1 == "+" & inter$strand2 == "-"))
            obs2[2] = obs2[2] + length( which(inter$strand1 == "+" & inter$strand2 == "+"))
            obs2[3] = obs2[3] + length( which(inter$strand1 == "-" & inter$strand2 == "-"))
            obs2[4] = obs2[4] + length( which(inter$strand1 == "-" & inter$strand2 == "+"))
            summary$number_TRA[indice_chromosome] = length(idx_inter)
        }

        summary$number_DEL[indice_chromosome] = obs[1]
        summary$number_h2hINV[indice_chromosome] = obs[2]
        summary$number_t2tINV[indice_chromosome] = obs[3]
        summary$number_DUP[indice_chromosome] = obs[4]
        summary$clusterSize_including_TRA[indice_chromosome] = sum(obs2)

        if(nrow(SVsnow) != 0){
            signif <- chisq.test(obs2, p=rep(1/4,4))$p.val
            summary$fragment_joints[indice_chromosome] <- signif} else{
                summary$fragment_joints[indice_chromosome] <- NA
        }
    }


    names(summary_inter) = paste0("inter_",names(summary_inter))
    return(cbind(summary,summary_inter))
}
