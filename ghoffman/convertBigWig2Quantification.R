#!/bin/env Rscript
#
# Gabriel Hoffman
# August 18, 2017
#
# Given a BED file of genome intervals, and a set of bigWig files
# compute the counts within each interval.  
# Save counts, CPM, RPKM, and expression BED file for downstream processing

library(getopt)

spec = matrix(c(
'help', 	'h', 0, "integer",
'bigWigList', 	'w',1, "character",
'bed', 			'b', 1, "character",
'nthreads', 	't', 1, "integer",
'out', 			'o', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

mssg = "# Given a BED file of genome intervals and a set of bigWig files,
# compute the counts within each interval.
# Save counts, CPM, RPKM, and expression BED file for downstream processing
# Gabriel Hoffman @ Icahn School of Medicine at Mount Sinai

--help\t\t show this message
--bigWigList\t text file with col1: BigWigs, col2: identifier, col3: readlength
--bed\t\t bed file with 6 cols 
--nthreads\t limited by memory usage
--out\t\t prefix of out files\n\n"

if( !is.null(opt$help) ){
	cat(mssg)
	q()
}

if( is.null(opt$bigWigList) || !file.exists(opt$bigWigList)){

	cat(mssg)
	cat("Define --bigWigList\n")
	q()
}

if( is.null(opt$bed) || !file.exists(opt$bed)){
	cat(mssg)
	cat("Define --bed\n")
	q()
}

if( is.null(opt$nthreads) ){
	cat(mssg)
	cat("Define --nthreads\n")
	q()
}

if( is.null(opt$out) ){
	cat(mssg)
	cat("Define --out\n")
	q()
}


suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(edgeR))

registerDoParallel( opt$nthreads )

# read bed file
bed2Granges = function(bed){
	colnames(bed) <- c('chr','start','end','id','score','strand')
	with(bed, GRanges(chr, IRanges(start+1, end), gsub('.', '*', strand), score, id=id))
}

regions = bed2Granges( read.table(opt$bed, stringsAsFactors=FALSE) )

# read BigWig list
bwdf = read.table( opt$bigWigList, stringsAsFactors=FALSE) 

if(ncol(bwdf) != 3){
	c("bigWigList must have 3 cols:\n
		1) name of bigwig file\n
		2) sample name\n
		3) readLength\n")
	q()
}

colnames(bwdf) = c("file", "name", "readLength")

bw <- BigWigFileList( bwdf$file)
names(bw) = bwdf$name

# Compute counts
counts = foreach( i = seq_len(length(bw)), .combine=cbind ) %dopar% {
	data = import(bw[[i]], as = 'RleList')

	counts = foreach( chrom = names(data) ) %do% {
	    sum(Views(data[[chrom]], ranges( regions[regions@seqnames == chrom,])))
	}
	rm(data)
	gc()
	unlist(counts)
}
colnames(counts) = names(bw)
rownames(counts) = regions$id



# Divide by read length and round to integer numbers
counts <- round(counts / bwdf$readLength, 0)

# Normalize by library size
d <- DGEList( counts )
d = calcNormFactors(d)    

quantify = list()

quantify[['cpm']] = cpm(d, log=FALSE)
quantify[['cpm_log2']] = cpm(d, log=TRUE)
quantify[['rpkm']] = rpkm(d, gene.length=regions@ranges@width, log=FALSE)
quantify[['rpkm_log2']] = rpkm(d, gene.length=regions@ranges@width, log=TRUE)

write.table( counts, quote=FALSE, file = paste0(opt$out, "_counts.tsv"), sep='\t')
write.table( quantify[['cpm']], quote=FALSE, file = paste0(opt$out, "_cpm.tsv"), sep='\t')
write.table( quantify[['cpm_log2']], quote=FALSE, file = paste0(opt$out, "_cpm_log2.tsv"), sep='\t')
write.table( quantify[['rpkm']], quote=FALSE, file = paste0(opt$out, "_rpkm.tsv"), sep='\t')
write.table( quantify[['rpkm_log2']], quote=FALSE, file = paste0(opt$out, "_rpkm_log2.tsv"), sep='\t')

# write QTLtools bed file

beddf = data.frame(regions)
annotation = with(beddf, data.frame(chr=seqnames, start, end, gene=id, length=width, strand))
df = data.frame(annotation, quantify[['rpkm_log2']])
colnames(df)[1] = "#chr"

write.table( df, row.names=FALSE, quote=FALSE, file = paste0(opt$out, "_rpkm_log2.bed"), sep='\t')