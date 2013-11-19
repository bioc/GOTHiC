
#Copyright: Bori Gerle (gerle@ebi.ac.uk), EMBL-EBI, 2011

.onlyPairing <- function(fileName1, fileName2, sampleName, fileType)
{
#DUPLICATETHRESHOLD: maximum amount of duplicated paired-end reads allowed (over that value it is expected to be PCR bias)
	reads1=readAligned(fileName1,type=fileType) ## returns identifier, strand, chromosome, position read, quality
	reads1_1=as(reads1,"GRanges") ## Grange object containing 4 slots, seqnames, ranges, strand, elementMetadata (data.frame)
	reads2=readAligned(fileName2,type=fileType)
	reads2_1=as(reads2,"GRanges")
#match IDs to find paired reads where both ends were mapped
	id1=as.character(elementMetadata(reads1_1)$id)
    id2=as.character(elementMetadata(reads2_1)$id)
	if(fileType=='BAM'){
		ids1 = sapply(strsplit(id1, '\\-'), '[[', 2)
		ids2 = sapply(strsplit(id2, '\\-'), '[[', 2)
		elementMetadata(reads1_1)$id = ids1
		elementMetadata(reads2_1)$id = ids2
	}else
	{
#old fastq files use /1 and /2 to differentiate read pairs, CASAVA 1.8 separates members of a pair with spaces	
		firsttwo=sapply(strsplit(id1[1:2], '[/ ]'), '[[', 1)
		if(firsttwo[1]==firsttwo[2]){
			ids1 = sapply(strsplit(id1, ' '), '[[', 1)
			ids2 = sapply(strsplit(id2, ' '), '[[', 1)
		}else
		{
			ids1 = sapply(strsplit(id1, '[/ ]'), '[[', 1)
			ids2 = sapply(strsplit(id2, '[/ ]'), '[[', 1)
		}
	}
	rm(id1,id2)
	s1=which(ids1%in%ids2)
	s2=which(ids2%in%ids1)
	paired_reads_1=reads1_1[s1]
	rm(reads1_1)
	paired_reads_1=paired_reads_1[order(elementMetadata(paired_reads_1)$id)]
	paired_reads_2=reads2_1[s2]
	rm(reads2_1)
	paired_reads_2=paired_reads_2[order(elementMetadata(paired_reads_2)$id)]
	paired_df=cbind(as.character(seqnames(paired_reads_1)),as.character(start(ranges(paired_reads_1))),as.character(end(ranges(paired_reads_1))),as.character(seqnames(paired_reads_2)),as.character(start(ranges(paired_reads_2))),as.character(end(ranges(paired_reads_2))))
	colnames(paired_df)=c("chr1","start1","end1", "chr2","start2","end2")
#save paired reads in GRanges objects and data.frame	
	paired_df=data.frame(paired_df,frequencies=rep(1,times=nrow(paired_df)),stringsAsFactors=FALSE)
	return(list(paired_reads_1=paired_reads_1, paired_reads_2=paired_reads_2, paired_df=paired_df,s1=s1,s2=s2))
	
}

.countDuplicatesOfRaw <- function(paired_df, sampleName) 
{	
#count how many times we have the exact same read pairs to be able to remove PCR duplicates
	colus=c("chr1","start1","end1", "chr2","start2","end2")
	vec <- do.call("paste",c(paired_df[colus],sep="_"))
	names(vec) <- c(1:length(vec))
	dup <- duplicated(vec)
	uniq <-vec[!dup]
	uniqfreq <- rep(1, times=length(uniq))
	uniq=cbind(uniq,uniqfreq)
	colnames(uniq)=c("interaction","freqs")
	vec1 <- vec[dup]
	vec1 <- sort(vec1)
	freq <- cbind(vec1,as.numeric(sequence(rle(vec1)$lengths))+1)
	colnames(freq)=c("interaction","freqs")
	vec=rbind(uniq, freq)
	vec <- vec[order(as.numeric(rownames(vec))),]
	paired_df$frequencies <- as.numeric(vec[,2])
	freqs <- as.numeric(vec[,2])
	freqstable <- table(freqs)
	l=c()
	for(i in 1:length(freqstable-1)){
		l[i]=freqstable[i]==freqstable[i+1]
	}
	p=freqstable[l]
	p=p[1]
	ths=as.numeric(names(p))
	return(list(paired_df=paired_df,freqs=freqs,ths=ths))
}	
		
		
.filterOutDuplicates <- function(paired_df, paired_reads_1, paired_reads_2, sampleName, DUPLICATETHRESHOLD) 
{	
#DUPLICATETHRESHOLD: maximum amount of duplicated paired-end reads allowed (over that value it is expected to be PCR bias)
#remove readpairs that are present in more replicates than the DUPLICATETHRESHOLD - keep as many of the reads as the threshold
    dth <- which(paired_df$frequencies<(DUPLICATETHRESHOLD+1))
	paired_reads_1<-paired_reads_1[dth]
	paired_reads_2<-paired_reads_2[dth]	
	return(list(paired_reads_1=paired_reads_1,paired_reads_2=paired_reads_2))
}
		

pairReads <- function(fileName1, fileName2, sampleName, DUPLICATETHRESHOLD=1, fileType='BAM')
{		   
	oP=.onlyPairing(fileName1, fileName2, sampleName, fileType)
	paired_reads_1=oP$paired_reads_1
	paired_reads_2=oP$paired_reads_2
	paired_df=oP$paired_df
	s1=oP$s1
	s2=oP$s2	   
	
	cDoR=.countDuplicatesOfRaw(paired_df, sampleName)
	paired_df=cDoR$paired_df
	freqs=cDoR$freqs
	ths=cDoR$ths
		
	filtered=.filterOutDuplicates(paired_df, paired_reads_1, paired_reads_2, sampleName, DUPLICATETHRESHOLD)
	save(filtered,file=paste(sampleName,"paired_filtered",sep="_"))
	return(filtered)		
}

.getRestrictionSitesFromBSgenome <- function(BSgenomeName, genome, restrictionSite, enzyme)
{
#check if the genome needed is already installed and install from BioConductor if not
      iG=installed.genomes()
    if(!BSgenomeName%in%iG){
    source("http://bioconductor.org/biocLite.R")
    biocLite(BSgenomeName)
    }

   library(BSgenomeName,character.only=TRUE)

    chrEnds <- seqlengths(genome)
#find out how many base pairs from the 5 end of the restriction site the enzyme cuts
    restRemain <- length(sapply(strsplit(restrictionSite,'\\^'),'[[',1))
#restriction site as regular expression
    rSite <- paste(strsplit(restrictionSite,'\\^')[[1]][1],strsplit(restrictionSite,'\\^')[[1]][2],sep='')
#find locations where restriction enzyme cuts for all chromosomes
    params <- new("BSParams",X=genome, FUN=matchPattern, exclude = "random",simplify=TRUE)
    resSite <- bsapply(params, pattern = rSite)
    starts <- lapply(resSite, function(x){c(1,(start(x)+restRemain))})
    lengths <- lapply(resSite, function(x){length(x)+1})
    chrs <- lapply(1:length(names(resSite)),function(i){rep(names(resSite)[i],times=lengths[[i]])})
    ends <- lapply(1:length(names(resSite)), function(i){c(start(resSite[[i]]),chrEnds[i])})
    starts <- unlist(starts)
    chrs <- unlist(chrs)
    ends <- unlist(ends)
#save restriction sites as GRanges    
    resGR <- GRanges(seqnames=chrs, ranges=IRanges(start=starts, end=ends))
    save(resGR,file=paste("Digest",BSgenomeName,enzyme,sep="_"))
    return(resGR)

}


#functions for parallel version of findOverlaps
.putRangesOnFirstCircle <- function(x, circle.length)
{
    x_start0 <- start(x) - 1L  # 0-based start
    x_shift0 <- x_start0 %% circle.length - x_start0
    shift(x, x_shift0)
}

.findOverlaps.circle <- function(circle.length, query, subject,
                                 maxgap, minoverlap, type)
{
    if (is.na(circle.length))
        return(findOverlaps(query, subject,
                            maxgap=maxgap, minoverlap=minoverlap,
                            type=type, select="all"))
    if (type != "any")
        stop("overlap type \"", type, "\" is not yet supported ",
             "for circular sequence ", names(circle.length))
    subject0 <- .putRangesOnFirstCircle(subject, circle.length)
    inttree0 <- IntervalTree(subject0)
    query0 <- .putRangesOnFirstCircle(query, circle.length)
    hits00 <- findOverlaps(query0, inttree0,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
    query1 <- shift(query0, circle.length)
    hits10 <- findOverlaps(query1, inttree0,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
    subject1 <- shift(subject0, circle.length)
    hits01 <- findOverlaps(query0, subject1,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
    ## Merge 'hits00', 'hits10' and 'hits01'.
    union(union(hits00, hits10), hits01)
}

.strandAsSignedNumber <- function(x)
{
    tmp <- as.integer(runValue(x))
    idx <- tmp >= 2L
    tmp[idx] <- tmp[idx] - 3L
    runValue(x) <- tmp
    as.vector(x)
}


.findOverlaps.parallel <- function(query, subject, maxgap=0L, minoverlap=1L,
type=c("any", "start", "end", "within", "equal"),
select=c("all", "first", "last", "arbitrary"),
ignore.strand=TRUE, mc.cores=1, mc.preschedule = TRUE)
{
    if (!isSingleNumber(maxgap) || maxgap < 0L)
    stop("'maxgap' must be a non-negative integer")
    type <- match.arg(type)
    select <- match.arg(select)
    
    ## merge() also checks that 'query' and 'subject' are based on the
    ## same reference genome.
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    
    q_len <- length(query)
    s_len <- length(subject)
    q_seqnames <- seqnames(query)
    s_seqnames <- seqnames(subject)
    q_seqlevels <- levels(q_seqnames)
    s_seqlevels <- levels(s_seqnames)
    q_splitranges <- splitRanges(q_seqnames)
    s_splitranges <- splitRanges(s_seqnames)
    q_ranges <- unname(ranges(query))
    s_ranges <- unname(ranges(subject))
    if (ignore.strand) {
        q_strand <- rep.int(1L, q_len)
        s_strand <- rep.int(1L, s_len)
    } else {
        q_strand <- .strandAsSignedNumber(strand(query))
        s_strand <- .strandAsSignedNumber(strand(subject))
    }
    
    common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
    results <- mclapply(common_seqlevels,
    function(seqlevel)
    {
        if (isCircular(seqinfo)[seqlevel] %in% TRUE) {
            circle.length <- seqlengths(seqinfo)[seqlevel]
        } else {
            circle.length <- NA
        }
        q_idx <- q_splitranges[[seqlevel]]
        s_idx <- s_splitranges[[seqlevel]]
        hits <- .findOverlaps.circle(circle.length,
        seqselect(q_ranges, q_idx),
        seqselect(s_ranges, s_idx),
        maxgap, minoverlap, type)
        q_hits <- queryHits(hits)
        s_hits <- subjectHits(hits)
        compatible_strand <- seqselect(q_strand, q_idx)[q_hits] *
        seqselect(s_strand, s_idx)[s_hits] != -1L
        hits <- hits[compatible_strand]
        remapHits(hits, query.map=as.integer(q_idx),
        new.queryLength=q_len,
        subject.map=as.integer(s_idx),
        new.subjectLength=s_len)
    }, mc.cores=mc.cores, mc.preschedule=mc.preschedule)
    
    ## Combine the results.
    q_hits <- unlist(lapply(results, queryHits))
    if (is.null(q_hits))
    q_hits <- integer(0)
    
    s_hits <- unlist(lapply(results, subjectHits))
    if (is.null(s_hits))
    s_hits <- integer(0)
    
    if (select == "arbitrary") {
        ans <- rep.int(NA_integer_, q_len)
        ans[q_hits] <- s_hits
        return(ans)
    }
    if (select == "first") {
        ans <- rep.int(NA_integer_, q_len)
        oo <- IRanges:::orderIntegerPairs(q_hits, s_hits, decreasing=TRUE)
        ans[q_hits[oo]] <- s_hits[oo]
        return(ans)
    }
    oo <- IRanges:::orderIntegerPairs(q_hits, s_hits)
    q_hits <- q_hits[oo]
    s_hits <- s_hits[oo]
    if (select == "last") {
        ans <- rep.int(NA_integer_, q_len)
        ans[q_hits] <- s_hits
        return(ans)
    }
    new2("Hits", queryHits=q_hits, subjectHits=s_hits,
    queryLength=q_len, subjectLength=s_len,
    check=FALSE)
}


.exportCoverage <- function(reads, sampleName)
{
	coverage_reads = coverage(reads)
	export.bedGraph(as(coverage_reads,"RangedData"), paste("coverage_", sampleName, ".bed", sep =""))
}

mapReadsToRestrictionSites <- function(pairedReadsFile, sampleName,BSgenomeName,genome,restrictionSite,enzyme,parallel=F, cores=1)
{
	paired_reads_1 <- pairedReadsFile$paired_reads_1
	paired_reads_2 <- pairedReadsFile$paired_reads_2

#calculate coverage
	all_reads <- c(paired_reads_1, paired_reads_2)
	.exportCoverage(all_reads, sampleName)
#take restriction sites	
    digestFile <- paste("Digest",BSgenomeName,enzyme,sep="_")
    if(file.exists(digestFile))
    {
    load(paste("Digest",BSgenomeName,enzyme,sep="_"))
    hindIIIRanges <- resGR
    }else
    {

    hindIIIRanges <- .getRestrictionSitesFromBSgenome(BSgenomeName, genome,restrictionSite,enzyme)
    }

#find HindIII sites
	if(parallel)
	{
		library(parallel)
		hindIII_1 <- .findOverlaps.parallel(paired_reads_1, hindIIIRanges, mc.cores=cores, select="first")
		hindIII_2 <- .findOverlaps.parallel(paired_reads_2, hindIIIRanges, mc.cores=cores, select="first")
	}else
	{
	hindIII_1 <- findOverlaps(paired_reads_1, hindIIIRanges, select="first")
	hindIII_2 <- findOverlaps(paired_reads_2, hindIIIRanges, select="first")
	}
	rm(paired_reads_1,paired_reads_2)
	
#take out reads where we cant assign hindIII site
    validInteractions=!is.na(hindIII_1) & !is.na(hindIII_2)
    hindIII_1=hindIII_1[validInteractions]
    hindIII_2=hindIII_2[validInteractions]
	
    hindIII1=pmin(hindIII_1, hindIII_2)
    hindIII2=pmax(hindIII_1, hindIII_2)
    
#assign hindIII restriction sites to reads
	loci1 <- hindIIIRanges[hindIII1]
	loci2 <- hindIIIRanges[hindIII2]
	interactingLoci <- GRangesList(locus1 = loci1, locus2 = loci2)
	outputfilename <- paste("interactingLoci", sampleName, sep = "_") 
	save(interactingLoci, file=outputfilename)
	return(interactingLoci)
}




#Copyright: Robert Sugar (robert.sugar@ebi.ac.uk), EMBL-EBI, 2011

#counts duplicate rows of a data frame
#returns a data frame with unique rows and the frequencies in the new row "frequencies"
#warning: data frame should be ordered!
.countDuplicates <- function(df, frequencyCol = c("frequencies", "frequency"), considerExistingFrequencies = any(frequencyCol %in% names(df)), ordered=FALSE) 
{	
	if(considerExistingFrequencies)
	{
		frequencyCol <- frequencyCol[frequencyCol %in% names(df)][1]		
		if(is.null(df[[frequencyCol]]))
			stop("frequencies column missing")
	} else if(!considerExistingFrequencies)
	{
		frequencyCol <- frequencyCol[1]
		df[, frequencyCol] <- 1 #one for all
	}
	
#order
	chrLocusColumns <- c("chr1", "locus1", "chr2", "locus2")
	coverageColumns <- c("chr", "locus")
	startEndColumns <- c("seqnames1", "start1", "end1", "strand1", "seqnames2", "start2", "end2", "strand2")
	if(all(chrLocusColumns %in% names(df)))
	{
		df <- df[order(df$chr1, df$locus1, df$chr2, df$locus2), ]
	} else
	if(all(coverageColumns %in% names(df)))
	{
		df <- df[order(df$chr, df$locus), ]
	} else
	if(all(startEndColumns %in% names(df)))
	{
		df <- df[order(df$seqnames1, df$start1, df$end1, df$strand1, df$seqnames2, df$start2, df$end2, df$strand2), ]
	} else #use all columns except for frequency
	{
		df <- df[do.call(order, df[, names(df)[-which(names(df) == frequencyCol)]]), ]
	}		
#add up the frequencies for the original and all duplicates of a row
	
	duplicates <- c(duplicated(df[, -which(names(df) == frequencyCol)])) #look for duplicates (excluding frequencies column)
	cumulative <- c(0, cumsum(df[, frequencyCol])) #the 0 is needed to consider not the first unique but the last duplicate
	beforeFirstUnique <- cumulative[!c(duplicates, FALSE)] #the extra FALSE value is needed for the last group
	
	dfUnique <- df[!duplicates, ]
	dfUnique[, frequencyCol] = diff(beforeFirstUnique)
	return(dfUnique)
	
}

#####bin interactions into desired bin size e.g. 1Mb########
.binInteractions <- function(interactions, resolution, frequencyCol = "frequencies", considerExistingFrequencies = frequencyCol %in% names(interactions))
{
	if(!identical(colnames(interactions)[1:4], c("chr1", "locus1", "chr2", "locus2")))
		stop("expecting columns chr1, locus1, chr2, locus2")
	if(considerExistingFrequencies && sum(names(interactions) == "frequencies") == 0)
		stop(paste("expecting column", frequencyCol))
	
	
#if resolution is given, bins will be calculated from interactions using the resolution
	interactions$locus1 <- (interactions$locus1 %/% resolution) * resolution
	interactions$locus2 <- (interactions$locus2 %/% resolution) * resolution	
	
#put smaller locus first to make sure a-b and b-a interactions look the same
	firstSmaller <- subset(interactions, (as.character(chr1) < as.character(chr2)) | (as.character(chr1) == (as.character(chr2)) & locus1 <= locus2))
	secondSmaller <- subset(interactions, !((as.character(chr1) < as.character(chr2)) | (as.character(chr1) == (as.character(chr2)) & locus1 <= locus2)))
#flip first and second if second is smaller
	cnames <- colnames(interactions)
	cnames[1:4] <- c("chr2", "locus2", "chr1", "locus1")
	setnames(secondSmaller,colnames(secondSmaller),cnames)
#stich them back together
	interactions <- rbind(firstSmaller, secondSmaller)
	
# --------- count the number of interactions between bins -------	
	interactions <- interactions[order(interactions$chr1, interactions$locus1, interactions$chr2, interactions$locus2), ]	
#bin it
	interactions <- .countDuplicates(interactions, frequencyCol, considerExistingFrequencies)	
	
	return(interactions)
}

# Author: borbalagerle
###############################################################################


# take GenomicRangesList with interactions in locus1 and locus2
.binomialHiC=function(interactingLoci, res,sampleName,BSgenomeName, genome, restrictionSite, enzyme, parallel=F, cores=NULL,removeDiagonal=TRUE,cistrans='all',filterdist=10000)
{
#take restriction sites
        digestFile <- paste("Digest",BSgenomeName,enzyme,sep="_")
    if(file.exists(digestFile))
    {
    load(paste("Digest",BSgenomeName,enzyme,sep="_"))
    hindIIIRanges <- resGR
    }else
    {
    hindIIIRanges <- .getRestrictionSitesFromBSgenome(BSgenomeName, genome, restrictionSite, enzyme)
    }
    
    elementMetadata(hindIIIRanges)$score=1:length(hindIIIRanges)
    
#remove self-ligations and interactions between adjacent bins by assigning the fragment number to reads and remove those where the difference is smaller or equal to 1
    x=findOverlaps(interactingLoci[[1]],hindIIIRanges,select='first')
    y=findOverlaps(interactingLoci[[2]],hindIIIRanges,select='first')
    interactingLoci1 <- interactingLoci[[1]]
    interactingLoci2 <-	interactingLoci[[2]]
    interactingLoci1$score=hindIIIRanges$score[x]
    interactingLoci2$score=hindIIIRanges$score[y]
    interactingLoci <- GRangesList(locus1=interactingLoci1,locus2=interactingLoci2)
    score=interactingLoci[[2]]$score-interactingLoci[[1]]$score
    interactingLoci[[1]]=interactingLoci[[1]][abs(score)>1]
    interactingLoci[[2]]=interactingLoci[[2]][abs(score)>1]
    
#filter for reads that are closer than 10kb to get rid of incomplete digest products
#	df_int <- data.frame(as.vector(seqnames(interactingLoci[[1]])), start(ranges(interactingLoci[[1]])), as.vector(seqnames(interactingLoci[[2]])), start(ranges(interactingLoci[[2]])))
	df_int <- data.frame(as.character(seqnames(interactingLoci[[1]])), start(ranges(interactingLoci[[1]])), as.character(seqnames(interactingLoci[[2]])), start(ranges(interactingLoci[[2]])))
	colnames(df_int) <- c("chr1", "locus1", "chr2", "locus2")
	df_filtered <- df_int
	df_filtered$dist <-as.vector(ifelse(df_filtered$chr1 == df_filtered$chr2, abs(df_filtered$locus1 - df_filtered$locus2), Inf))
	df_filtered <-df_filtered[df_filtered$dist>filterdist,]
	df_filtered <-df_filtered[,1:4]
#bin interactions according to resolution
	binned_df_filtered <- .binInteractions(df_filtered, res)
	binned_df_filtered$int1 <-paste(binned_df_filtered$chr1,binned_df_filtered$locus1,sep='_')
	binned_df_filtered$int2 <-paste(binned_df_filtered$chr2,binned_df_filtered$locus2,sep='_')
#diagonal removal - interactions within the same bin
	if(removeDiagonal)
	{
		subs <-which(binned_df_filtered$int1!=binned_df_filtered$int2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}
	if(cistrans=='cis'){
		subs <-which(binned_df_filtered$chr1==binned_df_filtered$chr2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}
	if(cistrans=='trans'){
		subs <-which(binned_df_filtered$chr1!=binned_df_filtered$chr2)
		binned_df_filtered <- binned_df_filtered[subs,]	
	}

#all read pairs used in binomial
	numberOfReadPairs <- sum(binned_df_filtered$frequencies)
#calculate coverage 
	all_bins <- unique(c(binned_df_filtered$int1,binned_df_filtered$int2))
	all_bins <- sort(all_bins)
	
	binned_dt=data.table(binned_df_filtered)
	covA <- binned_dt[,sum(c(frequencies)),by=int1]	
	covB <- binned_dt[,sum(c(frequencies)),by=int2]
	covA <- setkey(covA,key='int1')
	setnames(covB,'int2','int1')
	covB <- setkey(covB,key='int1')
	cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
	cov$V1.x[is.na(cov$V1.x)]=0
	cov$V1.y[is.na(cov$V1.y)]=0
	cov$coverage=cov$V1.x+cov$V1.y
	coverage=cov$coverage
	names(coverage)=cov$int1
	
	sumcov <- sum(coverage)
	relative_coverage <- coverage/sumcov
	names(relative_coverage)=names(coverage)
	binned_df_filtered$cov1 <- relative_coverage[binned_df_filtered$int1]
	binned_df_filtered$cov2 <- relative_coverage[binned_df_filtered$int2]
#probability correction assuming on average equal probabilities for all interactions
	numberOfAllInteractions <- length(all_bins)^2
	upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
	chromos <- unique(binned_df_filtered$chr1)
	chrlens <- c()
	for(cr in chromos){ 
		chrlens[cr] <- max(length(unique(binned_df_filtered$locus1[binned_df_filtered$chr1==cr])),length(unique(binned_df_filtered$locus2[binned_df_filtered$chr2==cr])))
	}
	cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
	transBinNumber <- upperhalfBinNumber-cisBinNumber
	
	diagonalProb <- sum(relative_coverage^2)
	if(cistrans=='all'){
		probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
	}
	if(cistrans=='cis'){
		probabilityCorrection <- upperhalfBinNumber/cisBinNumber
	}
	if(cistrans=='trans'){
		probabilityCorrection <- upperhalfBinNumber/transBinNumber
	}
	
#calculate probability of random interaction as a product of the relative coverages of individual bins*2
	binned_df_filtered$probability <- binned_df_filtered$cov1*binned_df_filtered$cov2*2*probabilityCorrection

#number of expected reads by chance
	binned_df_filtered$predicted <- binned_df_filtered$probability * numberOfReadPairs

#calculate cumulative binomial test for observed number of reads given the probability of seeing a random read
	if(parallel)
	{
		binomParams <- as.data.frame(t(cbind(as.numeric(binned_df_filtered[["frequencies"]]), as.numeric(binned_df_filtered[["probability"]]))))

		binned_df_filtered$pvalue <- unlist(mclapply(binomParams, function(x)
			{
					binom.test(x[1], numberOfReadPairs, x[2], alternative = "greater")$p.value
			},
			mc.cores=cores))
	} else
	{
	binned_df_filtered$pvalue <- apply(binned_df_filtered, 1, function(x)
			{
				binom.test(as.numeric(x[["frequencies"]]), numberOfReadPairs, as.numeric(x[["probability"]]), alternative = "greater")$p.value
			}	
	)
	}
	
#observed over expected log ratio
	binned_df_filtered$logFoldChange <- log2(binned_df_filtered$frequencies/binned_df_filtered$predicted)
#multiple testing correction separately for matrices with all interactions/only cis/only transs
	
	if(cistrans=='all'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=upperhalfBinNumber+length(all_bins))}
	}
	if(cistrans=='cis'){
		binned_df_filtered$qvalue <- if(removeDiagonal){p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber)}else{p.adjust(binned_df_filtered$pvalue, method = "BH", n=cisBinNumber+length(all_bins))}	
	}
	if(cistrans=='trans'){
		binned_df_filtered$qvalue <- p.adjust(binned_df_filtered$pvalue, method = "BH", n=transBinNumber)	
	}
   
    test <- data.frame(binned_df_filtered)
    test[,"pvalue"] <- test$pvalue
    pval.plot <- ggplot(test,aes(x=pvalue))
	tryCatch(
			 {
			 dev.new()
			 print(pval.plot + geom_density())
			 },
			 error=function(cond) {
			 message("No interactive plot, try saving image")
			 message(cond)
			 return(tryCatch(
							 {
							 pdf(file=paste(sampleName,"pvalue_distribution.pdf",sep="_"))
							 print(pval.plot + geom_density())
							 dev.off()
							 },
							 error=function(cond2) {
							 message("No pdf output, quality assesment plot is not produced")
							 message(cond2)
							 return(NA)
							 },
							 warning=function(cond2) {
							 message("No pdf output, quality assesment plot is not produced")
							 message(cond2)
							 return(NA)
							 }
							 ))
			 },
			 warning=function(cond) {
			 message("No interactive plot, try saving image")
			 message(cond)
			 return(tryCatch(
							 {
							 pdf(file=paste(sampleName,"pvalue_distribution.pdf",sep="_"))
							 print(pval.plot + geom_density())
							 dev.off()
							 },
							 error=function(cond2) {
							 message("No pdf output, quality assesment plot is not produced")
							 message(cond2)
							 return(NA)
							 },
							 warning=function(cond2) {
							 message("No pdf output, quality assesment plot is not produced")
							 message(cond2)
							 return(NA)
							 }
							 ))
			 })
   
    binned_df_filtered=binned_df_filtered[,c('chr1','locus1','chr2','locus2','cov1','cov2','probability', 'predicted','frequencies', 'pvalue','qvalue','logFoldChange')]
    colnames(binned_df_filtered)=c('chr1','locus1','chr2','locus2','relCoverage1','relCoverage2','probability', 'expected','readCount', 'pvalue','qvalue','logObservedOverExpected')
	return(binned_df_filtered)
}

####################################################################################################################################
####################################################################################################################################
## run GOTHiC

GOTHiC <- function(fileName1, fileName2, sampleName, res, BSgenomeName='BSgenome.Hsapiens.UCSC.hg19', genome=BSgenome.Hsapiens.UCSC.hg19, restrictionSite='A^AGCTT', enzyme='HindIII',cistrans='all',filterdist=10000, DUPLICATETHRESHOLD=1, fileType='BAM', parallel=FALSE, cores=NULL){
	if(file.exists(paste(sampleName,"paired_filtered",sep="_"))){
		message("Loading paired reads file ...")
		load(paste(sampleName,"paired_filtered",sep="_"))
		pairedReadsFile <- filtered

	}else{
		message("Pairing reads")
		
		pairedReadsFile <- pairReads(fileName1, fileName2, sampleName, DUPLICATETHRESHOLD=1, fileType=fileType)
	}
	if (file.exists(paste("interactingLoci",sampleName,sep="_"))){
	   message("Loading mapped reads file ...")
	   load(paste("interactingLoci",sampleName,sep="_"))
	   interactions <- interactingLoci
	}else{
		message("Mapping reads to restriction sites")
		interactions <- mapReadsToRestrictionSites(pairedReadsFile, sampleName, BSgenomeName, genome, restrictionSite, enzyme, parallel, cores)
	}

	message("Computing binomial ...")
	binom <- .binomialHiC(interactions, res,sampleName, BSgenomeName, genome, restrictionSite, enzyme, parallel, cores,cistrans=cistrans,filterdist=filterdist)
	return(binom)
}

