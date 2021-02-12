##%######################################################%##
#                                                          #
####           We phase maternal haplotypes             ####
####      and estimate haplotype frequencies using      ####
####               data from the F1 pools               ####
#                                                          #
##%######################################################%##


# NOTE : we denote the "SDR allele to which the maternal allele is linked",
# which is quite a mouthful, simply as "SNP type", and we use integer numbers to
# denote that (for performance reasons). Hence a SNP is of "type 1" if its
# maternal allele is linked to the W and is of "type 2" if this allele is linked
# to the Z
# IMPORTANT : to avoid using different functions for daughters and sons and to simply
# the code, we revert SNP type in sons (a type 1 SNP in daughters becomes a
# type-2 SNP in their brothers and reciprocally). This allows using the same
# equations in both sexes, and is equivalent to the equations in the paper (this
# was actually the solution described in the first version of the manuscript
# before the revision, but this made the biological explanations much more
# complex, so we dropped it in the revision). We refer to the reverted SNP type
# as "SNP group" (in both sexes)

# Because of this, f is the frequency of the Z-linked haplotype in sons, not the
# W-linked! As for the W-linked haplotype daughters, the Z should be inherited by all
# sons for a contig that is totally linked to the SDR (f = 0.5). So we are
# particularly interested in p(f = 0.5) for both sexes, which also simplifies the
# code (we have the same expectations in both sexes)

source("WZ_functions.R")

# our input is the table of informative SNPs generated at step 02
vcf <- fread("tables/vcfAndF1countsLikelihoodQ20.txt")

# and the scans of F1 bases generated at step 03 (imported later)


# we first filter SNPs --------------------------------------------

# assigning types to SNPs = phasing haplotype (see paper)
# 1 means that the maternal allele is linked to the W
# 2 means that it is linked to the Z.
vcf[, type := 2L - ((Cd / Rd) >= (Cs / Rs))]

# when frequencies are the same between the sexes, we attribute types at random
# (actually, based on whether the SNP position is even, to permit perfect
# replication in case this code is ran several times)
vcf[(Cd / Rd) == (Cs / Rs), type := POS %% 2L + 1L]

vcf[, p := pFgivenCountsSNPs(Cd, Rd, group = type) * 
      pFgivenCountsSNPs(Cs, Rs, group = 3L - type)]

# We check if allele frequencies are plausible given parental genotypes. In
# particular, we check if allele 1 is not too frequent in F1s (should be AT MOST
# 50%, considering that a frequency of 0.5 is higly unlikely to being with)
# We do it by computing the p-value of an exact binomial test that f > 0.5
vcf[, plausible := pbinom(Cd - 1L, Rd, 0.5, lower.tail = F) * pbinom(Cs - 1L, Rd, 0.5, lower.tail = F)]


# we flag the SNPs that we retain.
vcf[, retained := Cd + Cs > 0L & Rd > 0L & Rs > 0L # those covered in both pools and where the maternal allele was found
& (plausible > 0.05 | pmin(Cd / Rd, Cs / Rs) < 1 / 20) &
  Rd + Rs <= quantile(Rd + Rs, 0.95) & # those of sequencing depth not higher than the 95% quantile
  father.GQ > 40L & mother.GQ > 10L & # with requirement on genotype quality in parents
  Rd + Rs > (dA + dC + dG + dT + sA + sC + sG + sT) * 0.75] # the maternal + paternal allele must constitute at lest 75% of the F1 reads

# but we retain those that yield not too bad probability of f = 0.5
vcf[p > 0.05, retained := T]

vcf[is.na(retained), retained := F]


# We locate "outlier" SNPs that give a much lower probability that f=0.5
# than the other SNPs of the same contig.
# remember that f should equal 0.5 even in sons if the contig did not recombine 
# with the SDR (see introductory NOTE)
# we reframe the table so as all pools are on different rows and to obtain a
# single column for C and R, but we keep information on the sex with a new
# logical "sex" column which is TRUE for females. Doing so simplifies the script
counts <- vcf[retained == T, data.table(
  CHROM = rep(CHROM, 2), genomePos = rep(genomePos, 2),
  fam = rep(fam, 2), A1 = rep(A1, 2), A2 = rep(A2, 2),
  C = c(Cd, Cs), R = c(Rd, Rs), group = c(type, 3L - type),
  sex = rep(c(T, F), each = .N)
)]


# so we compute the (-log) posterior probabilty that f = 0.5 for each SNP
counts[, p := -log(pFgivenCountsSNPs(C, R, group))]

# and we apply our criterion to locate outlier SNPs
limit <- counts[, pLimit(p), by = .(cr = CHROM, fa = fam, fem = sex)]
counts[, lim := limit[chmatch(stri_c(CHROM, fam, sex), stri_c(cr, fa, fem)), V1]]

# we discard outliers in the vcf table. If a SNPs is an outlier in a pool, we
# discard it from both pools of a family, but not in the other family (if it is
# present)
for (fami in c("WXA", "ZM")) {
  vcf[fam == fami & genomePos %in% counts[fam == fami & p > lim, genomePos], retained := F]
}

setorder(vcf, fam, genomePos)

# we now count c1, r1, c2, r2 and compute P(f|cri) -----------------------------------
# Here, 1 refers to W in daughters and 2 to Z in sons, 
# 2 refers to Z in daughters and to W in sons 
# We do not use the notation cw, cz, rw and rz, as our implementation is not 
# a direct transcription of the method section (see the introductory NOTE)

# we import contig lengths to assign absolute positions to contigs in these tables
contigLengths <- readNamedVector("contigLengths.txt")

# we make tables of the retained SNPs that we will use for reach pool
# at first, one per family
SNPs = split(vcf[retained == T, .(CHROM = chmatch(CHROM, names(contigLengths)), 
                                  genomePos, A1, A2, group = type
                                  )], f = vcf[retained == T, fam])

# we duplicate each table, to have one per pool (2 pools per family)
SNPs = SNPs[c(1, 1, 2, 2)]

# for sons, we revert SNP types (see the NOTE at the begining)
for(i in c(2, 4)) {
  snps = copy(SNPs[[i]])
  snps[,group := 3L-group]
  SNPs[[i]] = snps
}


# we import the scans of F1 bams, which is a list of data.tables
scans <- list.files("1-F1", pattern = "markdup.basesAtSNPs.rds", full.names = T)
F1scans <- mclapply(scans, readRDS, mc.cores = length(scans))

# we name the tables
names(F1scans) <- stri_extract_first(scans, regex = "[A-Z]+_[a-z]")

# we obtain the cri variables (read counts for different SNP groups)
# see the countReads() function in WZ_functions.R
readCounts <- setNames(
  Map(countReads, F1scans, SNPs, 20L, length(contigLengths), T, names(F1scans)), 
  names(F1scans))

# we reclaim some RAM
rm(F1scans, counts, SNPs) ; gc()


# we compute the posterior probabilities of all possible haplotype frequencies
probs <- lapply(readCounts, function(counts) do.call(pFgivenCountsHapl, counts))

# we compute nrec :
# First, the expected value of f (weighted mean of f values given their posterior)
E_f <- vapply(probs, function(mat) colSums(t(mat) * Fset), numeric(nrow(probs[[1]])))

# then nrec
nrec <- rowSums(10 * (1 - 2 * E_f))

# we also record the highest posterior probability for f in each contig
maxP <- lapply(probs, rowMaxs)

# and the probability that f = 0.5 for each contig
p05 <- lapply(probs, function(x) x[, 11])

# we add these to the read count tables
readCounts <- Map(data.table, readCounts, pM = maxP, p05 = p05)

# and we combine the results for all pools
probs <- data.table(
  CHROM = names(contigLengths),
  length = contigLengths,
  do.call(data.table, readCounts),
  nrec,
  p = Reduce("*", p05)
)

writeT(probs, "tables/countsAndPosteriors.txt")
writeT(vcf, "tables/vcfAndF1countsLikelihoodQ20.txt")

