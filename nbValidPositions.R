# this script scans the mothers' bam file in order to count bases that are covered

# this is required to estimate heterozygosity in SNPs by kb, as we should not consider
# genomic regions that have inadequate sequencing depth (in which there was no chance
# to find heterozygous SNPs)

# the only argument if the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

source("WZ_functions.R")

# we import contig names, useful to call samtools on 
# separate contigs in parallel
contigs = fread("contigLengths.txt")$V1

setwd("0-F0")

# we generate the names of the bams
# the bams must be indexed
bams = stri_c("f", c("WXA","ZM"),"_recalibrated.bam")

# we retrieve the sequencing depth quantiles generated at step 01
# we use them as we must ignore regions of excessive sequencing depth
# (like we did to record heterozygous SNPs)
depthQuantiles = fread("DPquantile.txt")
depthQuantiles = unlist(depthQuantiles[,.(mother)])
names(depthQuantiles) = bams


nbPositions = function(contig) {
  # this function counts the number of positions in the range of sequencing depths, per contig
  
  # we do it separately for each bam
  depthBam = function(bam) {
    
    # path of temporary output file from samtools
    out = paste(contig, bam, sep ="_")
    system(paste("samtools view -F 1024 -F 256 -F 512 -F 4 -h", bam, contig, "| samtools depth -Q 1 -q 10 - >", out))
    
    # if there is no read for the contig, we return 0
    # we have to do that as fread() does not like empty files
    if(file.size(out) == 0) {
      file.remove(out)
      return(0L)
    }
    
    dep = fread(out, col.names = c("contig","pos","depth"))
    file.remove(out)
    
    # we report the number of unique positions in the depth range
    dep[depth >=5L & depth < depthQuantiles[bam], unique(pos)]
  }
  
  # we apply the above for the two bams
  pos = lapply(bams, depthBam)
  
  # and report the number of unique positions combining both bams
  length(unique(unlist((pos), use.names = F)))
}

# we apply the function in parallel for different contigs
res = mclapply(contigs, nbPositions, mc.cores = nCPUs)
res = data.table(CHROM = contigs, nPos = unlist(res))

writeT(res, "nbValidPositionsMothers.txt")