\name{chromoth-class}
\docType{class}

% Class:
\alias{class:chromoth}
\alias{chromoth-class}
\alias{chromoth}
\alias{getSVs}
\alias{getSegs}



\title{chromoth objects}

\description{
	The chromoth object is an object returned by the function \link{shatterseek}.
	}

\details{
	This object stores the results returned by the function \link{shatterseek}. 
	This object may be coerced to a data.frame using \code{as(object,"data.frame")}. 
		The obtained data.frame object contains 5 columns, "chrom","clusterSize","SVpvalue","CNVpvalue", and "score". The "chrom" column is the chromosome name. The "clusterSize" column is the maximum cluster size that the SVs can be clustered on a chromosome. The "SVpvalue" column is the pvalue given by performing statistical tests based on structural variations. The "CNVpvalue" is the pvalue given by performing statistical tests based on copy number variaitons. The "score" equals to \code{-log10(SVpvalue)-log10(CNVpvalue)}}.

\section{Accessors}{
	In the following, x is a chromoth object.
	  \describe{
	    \item{}{
	      \code{getSVs(x)}:
	        Get the structrual variations; Note that since \link{shatterseek} first try to filter possible false structrual variations, the structrual variations returned by this function might be different from the structrual variations provided to \link{shatterseek}.
    		}
	  \item{}{
		\code{getSegs(x)}:
		Get the copy number segmentation.
		}
	}

	}

\seealso{
	\link{plot_chromothripsis}, \link{shatterseek}
	}

\author{Ruibin Xi}
