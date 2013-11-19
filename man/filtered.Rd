\name{filtered}
\alias{filtered}
\alias{lymphoid_chr20_paired_filtered}
\docType{data}
\title{
A GenomicRangesList object used as an example in the GOTHiC package}
\description{
\code{filtered} is a GenomicRangesList example object used as an example for the binomialHiC package. This GenomicRangesList contains reads from a human lymphoblastoid cell line HiC experiment (Lieberman-Aiden et al. 2009) for chr20, that were mapped to the genome, paired and PCR duplicate-filtered. 
}
\usage{data(lymphoid_chr20_paired_filtered)}
\format{
  The format is:
GenomicRangesList with 2 slots:
$paired_reads_1 contains the coordinates for one end of the paired reads
$paired_reads_2 contains the coordinates for the other end of the paired reads
}
\author{
	Borbala Gerle and Robert Sugar
}
\seealso{ \code{mapReadsToRestrictionSites}}
	
\examples{
data(lymphoid_chr20_paired_filtered)
}
\keyword{datasets}
