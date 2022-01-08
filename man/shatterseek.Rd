\name{shatterseek}
\alias{shatterseek}

\title{The main function of the package shatterSeek}

\description{
	Detect chromothripsis events from structural variations and copy number variation data
}

\usage{ 
shatterseek(SV.sample,seg.sample, min.Size = 3, genome = "hg19")
} 
\arguments{
	\item{SV.sample}{
		A \link{SVs} object that stores all the structural variations of the sample.
		To detect clusters of interleaved SVs, ShatterSeek only considers structural variations with two breakpoints on the same chromosome.
		Inter-chromosomal structural variations are used in the subsequent steps to detect multichromosomal events and to evaluate the set of statistical criteria described in the tutorial.
	}
	\item{seg.sample}{
		A \link{CNVsegs} object that contains the copy number segmentation results for the sample.
	}
	\item{min.Size}{
		Minimum number of inleaved SVs required to report a cluster. Default is 3.
	}
	\item{genome}{
		Reference human genome assembly. hg19 or hg38 (Default is hg19).
	}
}


\value{
	A \link{chromoth} object. 
}

\seealso{ \link{plot_chromothripsis}
}

\references{
	Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. Cortes-Ciriano et al. 2018
}

\examples{
data(DO17373)
SV_data <- SVs(chrom1=as.character(SV_DO17373$chrom1), 
			   pos1=as.numeric(SV_DO17373$start1),
			   chrom2=as.character(SV_DO17373$chrom2), 
			   pos2=as.numeric(SV_DO17373$end2),
			   SVtype=as.character(SV_DO17373$svclass), 
			   strand1=as.character(SV_DO17373$strand1),
			   strand2=as.character(SV_DO17373$strand2))

CN_data <- CNVsegs(chrom=as.character(SCNA_DO17373$chromosome),
				   start=SCNA_DO17373$start,
				   end=SCNA_DO17373$end,
				   total_cn=SCNA_DO17373$total_cn)

chromothripsis <- shatterseek(SV.sample=SV_data,
                seg.sample=CN_data, 
                genome="hg19")
plots_chr3 <- plot_chromothripsis(ShatterSeek_output = chromothripsis, 
              chr = "chr3",
              sample_name="DO17373",
              genome="hg19")

}

