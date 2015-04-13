\name{pairReads}
\alias{pairReads}
\title{Function pairs aligned paired NGS reads
}
\description{
	This function takes bowtie output files, pairs the reads, only keeps those where both ends mapped, filters for perfect duplicates to avoid PCR bias, and saves and returns a GenomicRangesList object that contains the paired_reads_1 and paired_reads_2 GenomicRanges with the paired reads 
}
\usage{
pairReads(fileName1, fileName2, sampleName, DUPLICATETHRESHOLD = 1,
fileType='BAM')
}

\arguments{
  \item{fileName1}{
	File containing the mapped reads of the first fragment ends  (BAM or Bowtie format)
}
  \item{fileName2}{
	File containing the mapped reads of the second fragment ends (BAM or Bowtie format)
}
  \item{sampleName}{
	A character string that will be used to name the exported BedGraph file containing the coverage, and the R object file with paired reads. They will be saved in the current directory.
}
  \item{DUPLICATETHRESHOLD}{
	An integer specifying the maximum amount of duplicated paired-end reads allowed, over that value it is expected to be PCR bias. The default is 1.
}
  \item{fileType}{
	A character string specifying the format of the aligned reads. The default is 'BAM'. Other accepted format is 'Bowtie'.
}
}
\value{
	A GenomicRangesList called filtered
	\item{paired_reads_1}{GenomicRanges with the coordinates of where one end of the read mapped}
	\item{paired_reads_2}{GenomicRanges with the coordinates of where the other end of the read mapped}
}

\author{
	Borbala Mifsud and Robert Sugar
}
\seealso{\code{mapReadsToRestrictionSites}, \code{GOTHiC}
}

\examples{
library(GOTHiC)
dirPath <- system.file("extdata", package="HiCDataLymphoblast")
fileName1 <- list.files(dirPath, full.names=TRUE)[1]
fileName2 <- list.files(dirPath, full.names=TRUE)[2]
paired <- pairReads(fileName1, fileName2, sampleName='lymphoid_chr20', 
DUPLICATETHRESHOLD = 1, fileType='Table')
}

\keyword{manip}

