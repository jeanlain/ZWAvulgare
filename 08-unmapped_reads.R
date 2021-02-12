## %######################################################%##
#                                                          #
####              searching for W-specific              ####
####             sequences in unmapped reads            ####
#                                                          #
## %######################################################%##


source("WZ_functions.R")

dir.create("unmapped")

# We export fastas of unmapped reads and their mates ------------------------------
bams <- list.files(pattern = ".bam$")
exportFas <- function(bam) {

  # a prefix we use to name files
  prefix <- stri_extract_first(bam, regex = "[^_]+_[^_|^.]+")
  # first, we extract unmapped reads (flag 4)
  system(stri_c("samtools fasta -f4 -1 unmapped/", prefix, "_f4_R1.fas", " -2 unmapped/", prefix, "_f4_R2.fas ", bam))
  # notes that reads 1 and reads 2 are put in separate files

  # then their mates (flag 8)
  # I haven't found how to do both in a single samtools command
  system(stri_c("samtools fasta -f8 -1 unmapped/", prefix, "_f8_R1.fas", " -2 unmapped/", prefix, "_f8_R2.fas ", bam))
  return(NULL)
}

res <- mclapply(bams, exportFas, mc.cores = length(bams))

setwd("unmapped")

# we remove duplicates among reads as the two samtools commands
# have extracted the same read twice if both mates were unmapped
filterFasta <- function(pool) {
  dedup <- function(mate) {
    # "mate" is either 1 or 2
    # we import both fastas from reads 1 (or 2) in a single object
    files <- stri_c(pool, "_f", c(4, 8), "_R", mate, ".fas")
    reads <- readDNAStringSet(files)
    reads <- reads[!duplicated(names(reads))]
    # we add mate number at the end of read names
    names(reads) <- stri_c(names(reads), mate, sep = "_")
    reads
  }

  # we do it for both mates
  reads <- mclapply(1:2, dedup)

  # we combine mates in a single object
  reads <- do.call(c, reads)

  # we sort by read names, which gives us an interleaved fasta
  reads <- reads[order(names(reads))]
  writeXStringSet(reads, stri_c(pool, "_unmapped.fas"))
  return(NULL)
}

pools <- stri_c(rep(c("WXA", "ZM"), each = 2), rep(c("daughters", "sons"), 2), sep = "_")
res <- mclapply(pools, filterFasta, mc.cores = length(pools))



# counting 31-mers and their reverse complement with jellyfish -------------------------------------

out <- list.files(pattern = "_unmapped.fas$")

jelly <- function(fas) {
  prefix <- stri_extract_first(fas, regex = "[^_]+_[^_]+")
  system(stri_c("jellyfish count -m 31 -s 100M -t 10 -C -o ", prefix, ".mer_counts.js", " ", fas))
}

res <- mclapply(out, jelly, mc.cores = length(out))


# dumping k-mers to fasta files -------------------------

out <- list.files(pattern = ".js$")
dump <- function(js) {
  prefix <- stri_extract_first(js, regex = "[^_]+_[^_]+")
  system(stri_c("jellyfish dump ", js, " > ", prefix, ".dumps.fa"))
}

res <- mclapply(out, dump, mc.cores = length(out))


# importing k-mers ------------------------------------


# for daughters, we only import k-mers present at least 5 times
# our trick is to use seqtk since the number of times are encoded as sequence names
# this is faster than importing all k-mers and then filtering within R
# so we create a bed file of the sequence names we want
bed <- "temp.bed"
writeT(data.table(5:10000), bed, col.names = F)

dKmers <- function(fam) {
  # fread(cmd = ) is much faster than system(command =) to ingest stdout
  kmers <- fread(cmd = stri_c("seqtk subseq ", fam, "_daughters.mer.dumps.fa ", bed), header = F, colClasses = "character")$V1

  # we discard sequence names
  kmers[seq(2L, length(kmers), 2L)]
}

sKmers <- function(fam) {
  # for sons, we import all k-mers with fread (it appears to be faster than readDNAStringSet for large files)
  kmers <- fread(stri_c(fam, "_sons.mer.dumps.fa"), header = F, colClasses = "character")$V1
  kmers[seq(2L, length(kmers), 2L)]
}


dWXA <- dKmers("WXA")
sWXA <- sKmers("WXA")
dZM <- dKmers("dZM")
sZM <- sKmers("sZM")

# at this stage, R uses about 150 GB of RAM

# we combine daughter k-mers
K <- unique(c(dWXA, dZM))

# and we compute the various proportions
K_in_dWXA <- K %chin% dWXA
K_in_dZM <- K %chin% dZM
K_n_sWXA <- !K %chin% sWXA
K_n_sZM <- !K %chin% sZM

# the proportion of female specific k-mers in both families
mean(K_in_dZM & K_n_sZM & K_in_dWXA & K_n_sWXA)

# the product of these proportions for both families separately
mean(K_in_dZM & K_n_sZM) * mean(K_in_dWXA & K_n_sWXA)
