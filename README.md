# ShatterSeek 

ShatterSeek is an R package that provides utilities to detect chromothripsis events from 
next-generation sequencing (NGS) data.
It takes as input copy number (CN) and structural variation (SV) calls calculated
with the user preferred method. 
Hence, it is compatible with virtually any CN and SV caller.

The main function of the package, shatterseek,
first uses intrachromosomal SVs to detect clusters of interleaved
rearrangements.
Next, it evaluates a set of statistical criteria in each of these regions.
The output consists of a data frame reporting the value for the statistical criteria used 
and additional information for each chromosome.
ShatterSeek also provides functionalities for the visualization of SV and CN profiles,
thus enabling the visual inspection of the candidate chromothripsis regions detected.
ShatterSeek is fast, with an average running time of ~20s per sample.


Further details about the implementation and validation of ShatterSeek
using ~2,600 whole-genome sequencing data sets from The Pan-Cancer Analysis of Whole Genomes (PCAWG) project
can be found in the following publication:

**Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing**

[Cortes-Ciriano et al. Nature Genetics, 2020](https://www.nature.com/articles/s41588-019-0576-7)

All the chromothripsis calls reported in that publication can be visualized at: http://compbio.med.harvard.edu/chromothripsis/

A detailed tutorial illustrating how to use ShatterSeek is provided in ./inst/tutorial/tutorial.pdf

ShatterSeek is free for academic use **only**. For non-academic use, please email Dr. Adrian Li at Harvard University Office of Technology Development (adrian\_li@harvard.edu)

# Prerequisites

ShatterSeek is written entirely in R (>= 3.0.1) and depends on the following packages:
methods, BiocGenerics, graph, S4Vectors, GenomicRanges, IRanges, MASS, ggplot2, grid, and gridExtra

# Installation
`$ git clone https://github.com/parklab/ShatterSeek.git`

`$ R CMD INSTALL ShatterSeek`

Alternatively you can install ShatterSeek directly from R using devtools:
```R
library(devtools)
install_github("parklab/ShatterSeek")
```

# How to use
Please see the package tutorial (especially section "How to use ShatterSeek").

A few considerations that have caused some confusion in active users:

1. Please note that the parameter min.Size of the main function of the package (i.e., shatterseek) is set by default to 1. This parameter defines the minimum number of SVs required to identify a cluster of SVs in each chromosome. The statistical criteria are used to identify which clusters of SVs conform with the features of chromothripsis. For example, we require a minimum number of 6 SVs for a cluster to be considered as chromothripsis. We recommend to use min.Size=1 (default) because in some cases a chromothripsis event in a given chromosome (say chrA) might be linked with a cluster of a few SVs in another chromosomome (say chrB). If the value of min.Size is set to a value higher than 1 (say 6), the connection between these two chromosomes (which might be biologically important) would be missed, as the cluster of SVs in chrB would not be detected.  

2. Please note that ShatterSeek expects adjacent copy number segments to have a different copy number value. 
If two adjacent regions have the same copy number value but are considered as two separate entries in your copy number
data.frame, please merge them. You can use the following code:

```R
## d is a data.frame with colums: chr, start, end, total_cn
dd <- d
dd$total_cn[dd$total_cn == 0] <- 150000
dd$total_cn[is.na(dd$total_cn)] <- 0
library(GenomicRanges)
dd <- as(dd,"GRanges")
cov <- coverage(dd,weight = dd$total_cn)
dd1 <- as(cov,"GRanges")
dd1 <- as.data.frame(dd1)
dd1 <- dd1[dd1$score !=0,]
dd1 = dd1[,c(1,2,3,6)]
names(dd1) <- names(d)[1:4]
dd1$total_cn[dd1$total_cn == 150000] <- 0
d= dd1; rm(dd)
```

# Example 
Once installed you can run ShatterSeek with example data.

```R
library(Shatterkeek)
data(DO17373)
```
Prepare structural variation (SV) data:

```R
SV_data <- SVs(chrom1=as.character(SV_DO17373$chrom1), 
			   pos1=as.numeric(SV_DO17373$start1),
			   chrom2=as.character(SV_DO17373$chrom2), 
			   pos2=as.numeric(SV_DO17373$end2),
			   SVtype=as.character(SV_DO17373$svclass), 
			   strand1=as.character(SV_DO17373$strand1),
			   strand2=as.character(SV_DO17373$strand2))
```

Prepare copy number (CN) data:

```R
CN_data <- CNVsegs(chrom=as.character(SCNA_DO17373$chromosome),
				   start=SCNA_DO17373$start,
				   end=SCNA_DO17373$end,
				   total_cn=SCNA_DO17373$total_cn)

```

Run ShatterSeek:

```R
chromothripsis <- shatterseek(
                SV.sample=SV_data,
                seg.sample=CN_data,
                genome="hg19")
```

Plot chromothripsis per chromosome:

```R
library(gridExtra)
library(cowplot)

plots_chr21 <- plot_chromothripsis(ShatterSeek_output = chromothripsis, 
              chr = "chr21", sample_name="DO17373", genome="hg19")
              
plot_chr21 = arrangeGrob(plots_chr21[[1]],
                         plots_chr21[[2]],
                         plots_chr21[[3]],
                         plots_chr21[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))

plots_chrX <- plot_chromothripsis(ShatterSeek_output = chromothripsis, 
              chr = "chrX", sample_name="DO17373", genome="hg19")
plot_chrX = arrangeGrob(plots_chrX[[1]],
                        plots_chrX[[2]],
                        plots_chrX[[3]],
                        plots_chrX[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))
                         
plot_grid(plot_chr21, plot_chrX)

```

# Contact
If you have any questions or suggestions please contact us:

Isidro Cortes Ciriano: isidrolauscher at gmail.com or icortes@ebi.ac.uk

Peter J Park: peter_park at hms.harvard.edu

