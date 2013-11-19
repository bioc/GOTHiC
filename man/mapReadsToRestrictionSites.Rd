\name{mapReadsToRestrictionSites}
\alias{mapReadsToRestrictionSites}
\title{Function to map aligned and paired reads to the restriction fragments
}
\description{This function takes mapped paired NGS reads in the format of a GenomicRangesList object where the two end of the reads are in the GenomicRanges paired_reads_1 and paired_reads_2. It prepares the digestion file from the genome supplied to it with the given restriction enzyme and specificity and maps the reads to the fragments.
}
\usage{
mapReadsToRestrictionSites(pairedReadsFile, sampleName,
BSgenomeName, genome, restrictionSite, enzyme, parallel=F, cores=1)
}
\arguments{
  \item{pairedReadsFile}{
	R object of GenomicRangesList containing paired_reads_1 and paired_reads_2 GenomicRanges with the paired mapped reads from a Hi-C experiment.
}
  \item{sampleName}{
	A character string that will be used to name the exported R object file with the mapped reads containing a GenomicRangesList with slots locus1 and locus2. It will be saved in the current directory.
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
  \item{parallel}{
	Logical argument. If TRUE the mapping will be performed faster using multiple cores. The default is FALSE.
}
  \item{cores}{
	An integer specifying the number of cores used in the parallel processing if parellel=TRUE. The default is 1.
}
}

\value{
	A GenomicRangesList
	\item{locus1}{GenomicRanges with the coordinates of the start of the fragment where one end of the read mapped}
	\item{locus2}{GenomicRanges with the coordinates of the start of the fragment where the other end of the read mapped}
}

\author{
	Borbala Mifsud and Robert Sugar
}
\seealso{\code{pairReads}, \code{GOTHiC}
}


\examples{
library(GOTHiC)
data(lymphoid_chr20_paired_filtered)
mapped=mapReadsToRestrictionSites(filtered, sampleName='lymphoid_chr20', 
BSgenomeName='BSgenome.Hsapiens.UCSC.hg18', genome=BSgenome.Hsapiens.UCSC.hg18, 
restrictionSite='A^AGCTT', enzyme='HindIII', parallel=FALSE, cores=1)
}
\keyword{manip}

