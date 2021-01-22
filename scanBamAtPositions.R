
##%######################################################%##
#                                                          #
####            this script call samtools to            ####
####            scan a bam, and obtains the             ####
####            bases at all reads covering             ####
####       positions specified in a tabular file        ####
#                                                          #
##%######################################################%##


# Note: this script takes a lot of memory

args <- commandArgs(trailingOnly = TRUE)

# the inputs are: 
# path the file containing the position of SNPs (2 columns, CHROM and POS)
vcf <- args[1]

# path to the bam file to scan. Must be indexed or samtools will not run
bam <- args[2]

# number of CPUs (processes) to use as we process batches of SNPs in parallel
nCPUs <- as.integer(args[3])

# the output is a data.table with one row per position AND per read pair
# (so a read pair may appear several times). The columns are:
# - the position in absolute genome coordinate (no contig name, to save space) in column "pos"
# - the read pair identifier (integer) in column "readPair"
# - the mapping quality (integer) "mapQ" It is negative if the reads are not in "proper pair" (to save space)
# - the base the read carries at the position in column "base", coded as integer (1:4), again to speed things

source("WZ_functions.R")

# imports the VCF of informative SNPs
vcf <- fread(vcf)

# we will convert contig coordinates into absolute genome coordinates, which will be useful later on
# for this, we import a tabular file indicating the length of each contig
contigLengths <- fread("contigLengths.txt")

# which we use to obtain the start position of each contig, in absolute genome coordinates
contigStarts <- setNames(c(0L, cumsum(contigLengths$V2[-nrow(contigLengths)])),  contigLengths$V1)

# we will split the snp file in chunks to work in parallel
setorder(vcf, CHROM, POS)

# we count SNPs per contig
nSNPs <- vcf[, .N, by = CHROM]

# and create bins containing this number of SNPs at most (spanning several contigs)
# here a bin contains 80k SNPs, as our system has a lot of RAM
nSNPs[, bin := as.integer(cumsum(N) / 80000L)]

# for each SNP, we add a column indicating its bin
bin <- nSNPs[, rep(bin, N)]


# which allows splitting the snp coordinate file into these bins
bins <- split(vcf[, .(CHROM, POS)], bin, drop = T)
rm(vcf, bin, nSNPs, contigLengths)

# we list the cigar operations that we need to account for in our base counts below
operations <- c("M", "I", "D", "S", "H", "N")

# at some point will encode bases as raw bytes for speed and memory reasons
nums <- as.integer(charToRaw("ACGT"))

# and convert these bytes to integer from 1 to 4, this vector will allow this
baseNum <- integer(max(nums))
baseNum[nums] <- 1:4



# retreives bases and positions for a batch of SNP coordinates
getBasesAtSNPs <- function(snps) {

  # this will be the name of the bed file used by samtools to retreive the data we want
  bed <- snps[1, stri_c(bam, CHROM, POS, sep = "-")]
  writeT(snps, bed, col.names = F)

  # we ingest the relevant positions of the bam file below. This is
  # quite a bit faster than the equivalent command in Rsamtools. We cannot
  # ingest with fread() as the number of columns is variable and fread doesn't
  # like it. Note that the bam must be indexed.
  # note also that we ignore PCR duplicate, secondary aligmnets, unmapped reads
  # (probably not necessary) and reads that failed QC
  sam <- system(paste("samtools view -F 1024 -F 256 -F 512 -F 4 -ML", bed, bam), intern = T)
  
  if (length(sam) == 0L) {
    return(NULL)
  }

  # We split and retains the columns we want. We add rbind in case there is a
  # single SNPs, to force create a matrix
  sam <- data.table(rbind(stri_split(sam, fixed = "\t", simplify = T)[, c(1:6, 10)]))
  setnames(sam, c("qname", "flag", "refname", "pos", "mapq", "cigar", "seq"))

  # We convert read names to integers for speed and memory, and converts flag and
  # mapq to integers as they're still coded as strings
  sam[, c("qname", "flag", "mapq") := .(toInteger(qname), as.integer(flag), as.integer(mapq))]

  # for reads not mapped in proper pairs (sam flag 2), we retain this piece of
  # info by setting the mapping quality negative (avoids making a dedicated
  # column for that, as objects will be quite big)
  # In the paper, we don't use this piece of information
  sam[bitwAnd(flag, 2L) != 2L, mapq := -mapq]
  sam[, flag := NULL]

  # we convert alignment start positions to genomic coordinates
  sam[, pos := as.integer(pos) + contigStarts[refname]]
  snpPos <- snps[, POS + contigStarts[CHROM]]

  # fast way to extract read bases and convert them to integer in a single vector
  bases <- charToRaw(stri_flatten(sam$seq))

  # now converted to more "natural" integer from 1 ot 4
  bases <- baseNum[as.integer(bases)]

  # length of each aligned read part
  nc <- nchar(sam$seq)

  # we use it to create a data table with one row per base across all
  # alignments, "aln" is the alignment number = the orignal sam row number
  aln <- rep(1:nrow(sam), nc)
  bases <- data.table(base = bases, sam[aln, .(qname, mapq)], aln)
  
  # to calculate the genomic position of each base in reads, we need to extract
  # cigar operations. We extract those in two columns, the operation type (a
  # letter) and its "size" (an integer). Types are converted to integer via
  # chmatch(), for speed and memory reasons
  cigar <- data.table(
      operation = chmatch(unlist(stri_extract_all(sam$cigar, regex = "[:upper:]")), operations), 
      size = as.integer(unlist(stri_extract_all(sam$cigar, regex = "[0-9]+")))
      )

  # we ignore hard clipping operations as they correspond to bases not even in the reported read sequence
  cigar <- cigar[operation != 5L]

  # TRUE for deletion in reads
  deletion <- cigar[, operation == 3L | operation == 6L]

  # for these, the operation consumes 0 base on the read
  cigar[, sizeOnRead := ifelse(deletion, 0L, size)]

  # this new column will therefore indicate the start position in the read of
  # each remaining operation (this position is not reset between reads)
  cigar[, start := cumsum(c(1L, sizeOnRead)[1:.N])]

  # We add a logical that is TRUE during insertions and soft clipping operations that
  # "consume" no base in the reference. As we use the "size" of operations, the
  # length of this logical is the same as the "bases" table
  notInRef <- cigar[operation == 2L | operation == 4L, unlist(Map(":", start, start + size - 1L))]

  # this column indicates the size "consumed" in the reference by each read base, defaults to 1
  bases[, toAdd := 1L]

  # bases in soft clipped regions and insertions consume no reference base as they are not on the reference
  bases[notInRef, toAdd := 0L]

  # however, bases at the begining of deletions will consume more than one base
  # (= 1+ the size of the deletion)
  bases[cigar$start[deletion], toAdd := cigar$size[deletion] + 1L]

  # this allows computing the position of each base of the read in reference
  # coordinates, but starting at 1 for each alignment
  bases[, pos := bases[, cumsum(toAdd), by = aln]$V1]

  # which we convert to genomic coordinate by adding the position of the
  # alignment start that we recorded earlier
  bases[, pos := pos + sam[aln, pos] - 1L]

  # we only retain bases at SNP coordinates and remove those corresponding to
  # insertions in reads that have no equivalent in the reference
  # we also ignore bases not correspnding to ACGT (the base would be 0L)
  bases <- bases[toAdd > 0L & pos %in% snpPos & base != 0L, .(readPair = qname, mapq, A = abs(mapq), pos, base)]

  # and we will retain only one alignment per read pair per position (there can
  # be several due to overlap between forward and reverse read in short
  # fragments). We favor alignents of better quality
  # so we sort read pairs by decreasing mapQ
  setorder(bases, pos, -A)
  cat(".")
  file.remove(bed)
  res = bases[!duplicated(data.table(readPair, pos)), -"A"]
}

# run in parallel for differents batches of snps
res <- mclapply(bins, getBasesAtSNPs, mc.cores = nCPUs, mc.preschedule = F)
res = rbindlist(res)

# we save an RDS to save space, as numbers can be quite long
# they take 4 bytes in an RDS, but could take much more in a txt file
# and the RDS is compressed
saveRDS(res, file = stri_c(file_path_sans_ext(bam), ".basesAtSNPs.rds"))
