fmaxmax <- function(x, vec){ #segments of at least 4 alternating 
  v <- rle(vec)
if (sum(vec==0) == length(vec)){return(0)} else{
  max(v$lengths[v$values == x])}
}


fmax <- function(x, vec,cutoff=2){ #segments of at least 4 alternating 
  v <- rle(vec)
  sum(v$lengths[which(v$values == 1 & v$lengths >=cutoff)])
  }
fmax2 <- function(CNVsnow, cutoff=1){ 
if (nrow(CNVsnow) >= 4){
k = nrow(CNVsnow)
v = (CNVsnow$total_cn[1:(k-2)] - CNVsnow$total_cn[3:k])
v = 1*(v ==0)
  v <- rle(v)
  max(v$lengths[which(v$values == 1 & v$lengths >=cutoff)])} else{0}
  #sum(v$lengths[which(v$values == 1 & v$lengths >=cutoff)])} else{NA}

}

fmaxsample <- function(CNVsnow, cutoff=2,times=5000){ 
if (nrow(CNVsnow) >= 4){
o=c()
for (t in 1:times){
idx=sample(1:nrow(CNVsnow))
CN = CNVsnow[idx,]
k = nrow(CN)
v = (CN$total_cn[1:(k-2)] - CN$total_cn[3:k])
v = 1*(v ==0)
  v <- rle(v)
  o = c(o,sum(v$lengths[which(v$values == 1 & v$lengths >=cutoff)]))
}
return(sum(o>fmax2(CNVsnow))/times)
} else{NA}

}
fmax2b <- function(CNVsnow, cutoff=1){ 
if (nrow(CNVsnow) >= 4){
k = nrow(CNVsnow)
v = (CNVsnow$total_cn[1:(k-2)] - CNVsnow$total_cn[3:k])
v = 1*(v ==0)
  v <- rle(v)
  #sum(v$lengths[which(v$values == 1 & v$lengths >=cutoff)])} else{0}
max(v$lengths[v$values == 1])}else{NA}

}

fmaxsampleb <- function(CNVsnow, cutoff=1,times=5000){ 
if (nrow(CNVsnow) >= 4){
o=c()
for (t in 1:times){
idx=sample(1:nrow(CNVsnow))
CN = CNVsnow[idx,]
k = nrow(CN)
v = (CN$total_cn[1:(k-2)] - CN$total_cn[3:k])
v = 1*(v ==0)
  v <- rle(v)
  o = c(o,max(v$lengths[which(v$values == 1)]))
#& v$lengths >=cutoff)]))
  #o = c(o,sum(v$lengths[which(v$values == 1 & v$lengths >=cutoff)]))
}
return(sum(o>fmax2b(CNVsnow))/times)
} else{NA}

}




