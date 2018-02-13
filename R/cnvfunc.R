## identify chromthripsis with array-based segments
## seg should be a n by 3 data frame
## 1st conlumn: start position; 2nd column: end position; 3rd column: segmeans.
## significance: the significance level
## if nrow(seg)<=2 or no significant subsegment indentified, return NULL
## otherwise return the most significant segments (also a n by 3 data frame).

#segment.plot = function(seg,plot=TRUE,scale=1,xlab="",ylab="Log2 Copy Ratio",ylim=NULL,xlim=xlim,...)
#		{ if(nrow(seg)<=0||ncol(seg)!=3) {stop("seg must a n by 3 matrix or data frame (n>0)")}
#		  xlim = c(min(seg[,1]),max(seg[,2]))
#		  rg = xlim[2]-xlim[1]
#		  xlim[1] = max(xlim[1]-0.05*rg,0)
#		  xlim[2] = xlim[2]+0.05*rg
#		  abs.max = max(abs(seg[,3]))
#		  if(is.null(ylim)) {ylim=c(-abs.max,abs.max)*1.2}
#		  
#		  data.show = matrix(0,nrow = 2*nrow(seg)[1],ncol=2)
#   	  	  data.show[,1] = rep(seg[,3],each=2)
#		  tmp.ind = 2*c(1:nrow(seg))-1
#		  data.show[tmp.ind,2] = seg[,1]
#		  tmp.ind = tmp.ind+1
#   	  	  data.show[tmp.ind,2]= seg[,2]
#
#		  if(plot){
#			plot(data.show[,2]/scale,data.show[,1],type="l",xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim)
#			}
#		  else{
#			lines(data.show[,2]/scale,data.show[,1],...)
#			}
#		}



peak.count = function(x){
		if(length(x)<3) return(0);
		x1 = x[1:(length(x)-2)]
		x2 = x[2:(length(x)-1)]
		x3 = x[3:length(x)]
		
		cnt = sum((x1-x2)*(x2-x3)<0)
		return(cnt)
		}




