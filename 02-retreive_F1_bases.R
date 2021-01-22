##%######################################################%##
#                                                          #
####           Stage 2: scanning the F1 bams            ####
#                                                          #
##%######################################################%##



# this stage scans the F1 bams to retreive bases and reads at informative SNPs,
# and adds read counts to the SNP tabular file, wich will be used to determine
# SNP type (1 or 2)
# F1 bams were generated separately (not in a script)


# We scan the F1 bams ------------------------------
# we get the paths of the tabular files listing SNP position of reach family
# these were generated in stage 1. 
posFiles <- list.files(pattern = "informativeSNPs_", full.names = T, recursive = T)

# we get the path of bams file for F1 pools
# here they are stored in a directory "1-F1/"
bamFiles <- list.files("1-F1", pattern = ".bam$", full.names = T)

# note that because list.files() sorts alphabetically, WXa-related files come
# first (given how we named the files) so if one uses this script for its own
# bam files, they must ensure that the family order is preserved

# we launch the script to scan each F1 bam in parallel (using 5 CPUs for each)
m <- mcMap(
  function(posFile, bamFile) {
    system(paste("Rscript scanBamAtPositions.R", posFile, bamFile, 5))
  },
  rep(posFiles, each = 2), # each position file is used for two bam files (two pools)
  bamFiles,
  mc.cores = length(bamFiles)
)


# we retreive result to obtain read counts at SNPs ----------------------------------------
# these counts will be used to infer SNP type 

# we list the resulting files
scans <- list.files("1-F1", pattern = "basesAtSNPs.rds", full.names = T)

# before importing these files, we prepare the tabular file of informative SNPs
vcf <- fread("0-F0/informativeSNPs.txt")

# as we did during the bam scan, we convert SNP positions into absolute genomic positions
# for this we import a table of contig lengths
contigLengths <- fread("contigLengths.txt")

# which we use to obtain the start position of each contig, in absolute genome coordinates
contigStarts <- setNames(c(0L, cumsum(contigLengths$V2[-nrow(contigLengths)])), contigLengths$V1)
vcf[, genomePos := POS + contigStarts[CHROM]]
setorder(vcf, fam, genomePos)

# we import the scans of F1 bams, which are data.tables
F1scans <- mclapply(scans, readRDS, mc.cores = length(scans))
names(F1scans) <- basename(scans)

# we retreive the number of reads covering each base in each pool
# for this we need the genomic position of each SNP in each family
SNPpos <- split(vcf$genomePos, vcf$fam)

# and we count reads
baseCounts <- mcMap(baseMatrix, scannedBam = F1scans, genomePos = rep(SNPpos, each = 2), minMapQ = 20L, mc.cores = length(scans))

# we rbind the counts for each sex (as both families in the vcf table are on
# different rows, ZM coming after WXA) we first extract the sex of the pools,
# which we indicated in the file names between two underscores
sex <- splitToColumns(basename(scans), "_", 2)

sexCounts <- lapply(split(baseCounts, sex), function(x) do.call(rbind, x))
# the rows of these matrices correspond to those of the vcf table

# we can now count the number of reads carrying the maternal allele ("C"),
# and the the number of reads carrying the maternal or paternal alleles ("R")
alleleCounts <- lapply(sexCounts, function(sexCount) data.table(
    C = sexCount[vcf[, cbind(1:.N, A1)]],
    R = sexCount[vcf[, cbind(1:.N, A1)]] + sexCount[vcf[, cbind(1:.N, A2)]]
  )
)

# we add these counts to the vcf table
vcf <- data.table(
  vcf[, 1:13, with = F],
  do.call(cbind, sexCounts),
  do.call(cbind, alleleCounts)
)

# and we rename some columns. We use short names as the table is quite wide
setnames(vcf, 14:ncol(vcf), 
         c("dA", "dC", "dG", "dT", "sA", "sC", "sG", "sT", "Cd", "Rd", "Cs", "Rs"))

writeT(vcf, "tables/vcfAndF1countsLikelihoodQ20.txt")
