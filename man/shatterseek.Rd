\name{shatterseek}

\alias{shatterseek}



\title{The main function of the package shatterSeek}

\description{
	Detect chromothripsis with structural variations and copy number variations
	}

\usage{ 
        shatterseek(SV.sample,seg.sample,cores=20,B=5000)
        } 
\arguments{
	\item{SV.sample}{
		A \link{SVs} object that stores all the structural variations of the sample; ShatterSeek only considers structural variations with two breakpoints on the same chromosome. Therefore, inter-chromosomal structural variations will be ignored.
		}
	\item{seg.sample}{
		A \link{CNVsegs} object that contains the copy number segmentation results of the sample.
		}
	\item{cores}{
		A positive integer. The number of cores to use for parallel computing.
		}
	\item{B}{
		A positive integer. The number of permutations for assigning the p-values.
		}
	}

\details{
	
	}


\value{
	A \link{chromoth} object. 
	}

\seealso{ \link{plot.chromoth}
	}

\references{
	The Cancer Genome Atlas Research Network (2015) Comprehensive Molecular Characterization of Cutaneous Melanoma. Cell, 161(7):1681-96.
	}

\examples{
	set.seed(1234)
	data(SKCM) ### load the SKCM data
	s = "TCGA-D3-A2JC-06A-11D-A18Z" ### look at this sample

	SV.sample = SV.SKCM[SV.SKCM$sample.name==s,]
	segSample = seg.SKCM[seg.SKCM$sample.name==s,]
	segSample = segSample[segSample$chr!="Y",]
	SV.sample = SVs(chrom1=SV.sample$chrom1,pos1=SV.sample$pos1,chrom2=SV.sample$chrom2,pos2=SV.sample$pos2) ### structrual variations calls
	segSample = CNVsegs(chrom=segSample$chr,start=segSample$start,end=segSample$end,log2copyratio=segSample$segmean) ### segmentation 
	chromthripsis = shatterseek(SV.sample,segSample,B=1000) ##with both SV and CNV data
	chromthripsis = shatterseek(SV.sample,B=1000) ##with only SV data
	plot.chromoth(chromthripsis,chrom="15")
	}

\author{Ruibin Xi}
