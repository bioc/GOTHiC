\name{GOTHiC}
\alias{GOTHiC}
\title{Genome Organisation Through HiC}
\description{
\code{GOTHiC} performs a cumulative binomial test to detect interactions between distal genomic loci that have significantly more reads than expected by chance in Hi-C experiments. It takes mapped paired NGS reads as input and gives back the list of significant interactions for a given bin size in the genome. }
\usage{
GOTHiC(fileName1, fileName2, sampleName, res,
BSgenomeName='BSgenome.Hsapiens.UCSC.hg19',
genome=BSgenome.Hsapiens.UCSC.hg19, restrictionSite='A^AGCTT',
enzyme='HindIII',cistrans='all',filterdist=10000,
DUPLICATETHRESHOLD=1, fileType='BAM', parallel=FALSE, cores=NULL)
}
\arguments{
  \item{fileName1}{
	File containing the mapped reads of the first fragment ends  (BAM or Bowtie format)
}
  \item{fileName2}{
	File containing the mapped reads of the second fragment ends (BAM or Bowtie format)
}
  \item{sampleName}{
	A character string that will be used to name the exported BedGraph file containing the coverage, R object files with paired and mapped reads, and the final data frame with the results from the binomial test. They will be saved in the current directory.
}
  \item{res}{
	An integer that gives the required bin size or resolution of the contact map e.g. 1000000.
}
  \item{BSgenomeName}{
	A character string of the BSgenome package required to make the restriction fragment file containing information for both the organism the experiment was made in, and the genome version the reads were mapped to. The default is the current human genome build 'BSgenome.Hsapiens.UCSC.hg19'.
}
  \item{genome}{
	The BSgenome package required to make the restriction fragment file containing information for both the organism the experiment was made in, and the genome version the reads were mapped to. The default is the current human genome build BSgenome.Hsapiens.UCSC.hg19.
}
  \item{restrictionSite}{
	A character string that specifies the enzymes recognition site, ^ indicating where the enzyme actually cuts. The default is the HindIII restriction site: 'A^AGCTT'.
}
  \item{enzyme}{
	A character string containing the name of the enzyme used during the Hi-C experiment (i.e. "HindIII", "NcoI"). The default is "HindIII".
}
  \item{cistrans}{
	A character string with three possibilities. "all" runs the binomial test on all interactions, "cis" runs the binomial test only on intrachromosomal/cis interactions, "trans" runs the binomial test only on interchromosomal/trans interactions.
}
  \item{filterdist}{
	An integer specifying the distance between the midpoint of fragments under which interactions are filtered out in order to filter for those read-pairs where the digestion was incomplete. The default is 10000.
}
  \item{DUPLICATETHRESHOLD}{
	An integer specifying the maximum amount of duplicated paired-end reads allowed, over that value it is expected to be PCR bias. The default is 1.
}
  \item{fileType}{
	A character string specifying the format of the aligned reads. The default is 'BAM'. Other accepted format is 'Bowtie'.
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
\seealso{\code{binom.test}, \code{pairReads}, \code{mapReadsToRestrictionSites}
}

\examples{
library(GOTHiC)
dirPath <- system.file("extdata", package="HiCDataLymphoblast")
fileName1 <- list.files(dirPath, full.names=TRUE)[1]
fileName2 <- list.files(dirPath, full.names=TRUE)[2]
binom=GOTHiC(fileName1, fileName2, sampleName='lymphoid_chr20', res=1000000, 
BSgenomeName='BSgenome.Hsapiens.UCSC.hg18', genome=BSgenome.Hsapiens.UCSC.hg18, 
restrictionSite='A^AGCTT', enzyme='HindIII',cistrans='all', filterdist=10000, 
DUPLICATETHRESHOLD=1, fileType='Bowtie', parallel=FALSE, cores=NULL)
}
\keyword{manip}
\keyword{package}
