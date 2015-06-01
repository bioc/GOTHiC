\name{GOTHiChicup}
\alias{GOTHiChicup}
\title{Genome Organisation Through HiC from HiCUP output}
\description{
\code{GOTHiChicup} performs a cumulative binomial test to detect interactions between distal genomic loci that have significantly more reads than expected by chance in Hi-C experiments. It takes mapped and filtered paired NGS reads from HiCUP as input and gives back the list of significant interactions for a given bin size in the genome. }
\usage{
GOTHiChicup(fileName, sampleName, res, restrictionFile, cistrans='all', parallel=FALSE, cores=NULL)
}
\arguments{
  \item{fileName}{
	 A character string with the name of the file containing the mapped, filtered reads from HiCUP, after the default HiCUP output is converted to a table containing only the read ID, chromosome, start and end positions columns. Can be gzipped. (Tab separated text format)
}
  \item{sampleName}{
	A character string that will be used to name the quality control plot. It will be saved in the current directory.
}
  \item{res}{
	An integer that gives the required bin size or resolution of the contact map e.g. 1000000, for fragment level use 1.
}
  \item{restrictionFile}{
	A character string with the name of the digest file from HiCUP. It is used to map reads to restriction fragments. (.txt file name)
}
  \item{cistrans}{
	A character string with three possibilities. "all" runs the binomial test on all interactions, "cis" runs the binomial test only on intrachromosomal/cis interactions, "trans" runs the binomial test only on interchromosomal/trans interactions.
}
  \item{parallel}{
	Logical argument. If TRUE the mapping and the binomial test will be performed faster using multiple cores. The default is FALSE.
}
  \item{cores}{
	An integer specifying the number of cores used in the parallel processing if parellel=TRUE. The default is NULL.
}
}

\value{
	A data.frame containing elements

	\item{chr1 / chr2 }{chromosome(s) containing interacting regions 1 and 2}
 	\item{locus1 / locus2 }{start positions of the interacting regions 1 and 2 in the corresponding chromosome(s)}
	\item{relCoverage1 / relCoverage2 }{relative coverage corresponding to regions 1 and 2}
	\item{probability}{expected frequency}
	\item{expected}{expected number of reads}
	\item{readCount }{observed reads number}
	\item{pvalue}{binomial p-value}
	\item{qvalue}{binomial p-value corrected for multi-testing with Benjamini-Hochberg}
	\item{logObservedOverExpected}{observed/expected read numbers log ratio}
}
\author{
	Borbala Mifsud and Robert Sugar
}
\seealso{\code{binom.test}
}

\examples{
library(GOTHiC)
dirPath <- system.file("extdata", package="HiCDataLymphoblast")
fileName <- list.files(dirPath, full.names=TRUE)[4]
restrictionFile <- list.files(dirPath, full.names=TRUE)[3]
binom=GOTHiChicup(fileName, sampleName='lymphoid_chr20', res=1000000, 
restrictionFile, cistrans='all', parallel=FALSE, cores=NULL)
}
\keyword{manip}
\keyword{package}
