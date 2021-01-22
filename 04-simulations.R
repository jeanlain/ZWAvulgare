##%######################################################%##
#                                                          #
####            here we perform the various             ####
####           simulations used in this study           ####
#                                                          #
##%######################################################%##


source("WZ_functions.R")

setwd("simu")

# we import the data we need for the simulations:
# data.tables listing positions of SNPs, contig, and read pair identifiers at those positions
files <- list.files(pattern = "reads.RDS")
F1reads <- setNames(lapply(files, readRDS), splitToColumns(files, ".", 1))

# for the simulation, we create table listing the read pairs for each contig
# (without duplicates) as a read pair spannig several SNPs appears multiple times
uReads <- lapply(F1reads, function(x) x[!duplicated(readPair), .(readPair, contig)])

# we also create vectors of contig identifiers corresponding to each SNP in each family
orderedContigs <- lapply(F1reads[c(1, 3)], function(x) {

  # we first obtain the last SNP position per contig
  dt <- x[, max(pos), by = contig]

  # we order contigs accoring to ascending SNP positions
  setorder(dt, V1)

  # we compute the number of SNPs per contig
  nSNPs <- dt[, c(V1[1], diff(V1))]

  # and return an integer vector, such that vector[x] returns the contig for SNP x
  # this help simulation running fast
  rep(dt$contig, nSNPs)
})

# we compute the number of SNPs per family
nSNPs <- sapply(orderedContigs, length)
names(nSNPs) <- splitToColumns(names(nSNPs), "_", 1)

# and create of vector for all contigs identifiers, which will be useful
uContigs <- sort(unique(unlist(orderedContigs, use.names = F)))

simulate <- function(d = 0.5) {
  # this is the function that computes nrec given a simulated distance d to the
  # SDR in Morgans it does it for all the contigs and all pools at once (it
  # cannot do it for a single contig or pool)
  
  # d can be a vector, in which case contigs can have different values of d
  d <- rep(d, length.out = max(uContigs))

  # based on d, we sample numbers of chromosomes carrying haplotype Z in pools
  nZ <- setNames(as.data.table(replicate(4, rbinom(max(uContigs), 10, d), simplify = F)), c("WXAd", "WXAs", "ZMd", "ZMs"))

  # we need to correct these for sons (10 sons per pool)
  nZ[, c("WXAs", "ZMs") := .(10L - WXAs, 10L - ZMs)]
  # so for each pool, nZ[x, pool] given the number of
  # chromosomes carrying haplotype Z for the contig in the pool

  # for each family we sample the haplotype carrying each maternale allele (= SNP type)
  # 1 = W, 2 = Z.
  haplo <- lapply(nSNPs, function(x) sample(1:2, x, replace = T))

  # we replicate it for each pool
  haplo <- rep(haplo, each = 2)
  # so haplo[[pool]][x] returns the type of SNP x in pool

  # we determine whether a read carries the maternal allele at each SNP according to nZ
  for (pool in 1:length(F1reads)) {
    
    # we retrive the number of chromosomes carrying haplotype z in the pool
    # it's easyer to turn the nz data table into a list for this
    nZs <- as.list(nZ)[[pool]]
    
    # We first attribute each read to maternal DNA (T) or paternal DNA (F) with same probability
    readHapl <- sample(c(T, F), nrow(uReads[[pool]]), T)
    # so readHapl[x] tells whether read x comes from maternal DNA

    # we assign maternal reads to haplotypes Z and W. To do so, we use a binomial
    # distribution with probabilty = the frequency of chromosomes carrying
    # haplotype Z among the 10 maternally inhertited chromosomes. +1L is used to
    # convert 0,1 output by rbinom() into 1,2 (1 = W, 2 = Z) patternal reads
    # will have number 0
    readHapl[readHapl] <- uReads[[pool]][readHapl, rbinom(.N, 1L, nZs[contig] / 10) + 1L]

    # the read carries the maternal allele if the source haplotype of the
    # maternal allele and of the read are the same
    F1reads[[pool]][, maternal := readHapl[readPair] == haplo[[pool]][pos]]
  }


  # to infer SNP type, we get the number of reads for each allele, per pool
  getCounts <- function(f1reads, snps) {
    counts <- f1reads[, tabulate(pos * 2L - 1L + !maternal, nbins = snps * 2L)]
    matrix(counts, ncol = 2, byrow = T)
  }

  counts <- Map(getCounts, F1reads, rep(nSNPs, each = 2))

  # we split the results per family, as SNP type is relevant to a family
  counts <- split(counts, splitToColumns(names(counts), "_", 1))

  getType <- function(twoPools) {
    # we combine matrices of the two pools in a data.table
    twoPools <- do.call(data.table, twoPools)
    setnames(twoPools, stri_c("V", 1:4))
    type <- twoPools[, 2L - (V1 / (V1 + V2) > V3 / (V3 + V4))]
  }

  types <- lapply(counts, getType)

  # we get the per-contig proportion of SNPs for which SNP type was correctly inferred
  prop <- Map(function(hapl, type, contig) {
    dt <- data.table(test = hapl == type, contig)
    dt <- dt[, .(prop = mean(test)), by = contig]
    prop <- dt[match(uContigs, contig), -1, with = F]
  }, haplo[c(1, 3)], types, orderedContigs)

  # we compute the f posteriors
  computeProb <- function(f1reads, group) {
    # we make the counts of read pairs covering the inferred haplotyes in the same way we did for real data
    readCount <- f1reads[, .(contig = 1:nCont, count = tabulate(contig[!duplicated(readPair)], nbins = nCont)), 
                         by = .(group = group[pos], maternal)]
    readCount <- dcast(readCount, contig ~ ifelse(maternal, "c", "r") + group, value.var = "count", sep = "", fill = 0L)

    # if somes cases, the c2 column may be missing (if nZ = 0 for all contigs)
    if (ncol(readCount) < 5L) readCount[, c2 := 0L]
    
    setcolorder(readCount, c("contig", "c1", "r1", "c2", "r2"))
    readCount[, c("r1", "r2") := .(c1 + r1, c2 + r2)]
    probas <- do.call(pFgivenCountsHapl, readCount[, -"contig"])
  }

  # we compute probabilities with inferred SNP type
  probasInferred <- Map(computeProb, 
                        F1reads, 
                        list(types[[1]], 3L - types[[1]], types[[2]], 3L - types[[2]]))

  # and with simulated SNP type
  probasSimulated <- Map(computeProb, 
                         F1reads, 
                         list(haplo[[1]], 3L - haplo[[2]], haplo[[3]], 3L - haplo[[4]]))

  computeNrec <- function(prob) {
    E_f <- vapply(prob, function(mat) colSums(t(mat) * Fset), numeric(nrow(prob[[1]])))
    nrec <- rowSums(10 * (1 - 2 * E_f))
  }

  nrecs <- lapply(list(nrec = probasInferred, nrecS = probasSimulated), computeNrec)

  res <- data.table(
    contig = uContigs, d = d[uContigs],
    n = nZ[uContigs, WXAd + 10L - WXAs + ZMd + 10L - ZMs], # n is the true number of recombinants
    do.call(data.table, prop), do.call(data.table, nrecs)
  )

  cat(".")
  res
}


# Evaluating the performance of the inferrence of nrec and SNP type ----------------------------

# for this, we simulate contigs at various distances to the SDR.
ds <- 0:20 / 40
res <- rbindlist(lapply(ds, simulate))

# we will only consider contigs with SNPs in both families
twoFams <- res[, !is.na(WXA.prop) & !is.na(ZM.prop)]

# we also retreive the number of informative SNPs per contig for both families
SNPsPerContig <- lapply(orderedContigs, tabulate)
res[, c("nWXA", "nZM") := .(SNPsPerContig[[1]][contig], SNPsPerContig[[2]][contig])]

# we make figure S1
library(ggplot2)
library(gridExtra)

temp <- ggplot(data = res[twoFams & n <= 20L, .(n, WXA.prop, ZM.prop, nrec, nrecS, nWXA, nZM)], mapping = aes(x = factor(n)))

p1 <- temp + geom_abline(intercept = -1, slope = 1, col = "darkgrey") + 
  geom_violin(aes(y = nrec), fill = rgb(1, 0.6, 0.6), col = NA, scale = "width") + xlab("") +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.45), "cm"))

p2 <- temp + geom_abline(intercept = -1, slope = 1, col = "darkgrey") + 
  geom_violin(aes(y = nrecS), fill = rgb(1, 0.6, 0.6), col = NA, scale = "width") + xlab("") +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.45), "cm"))

p3 <- temp + geom_violin(
  aes(y = (WXA.prop * nWXA + ZM.prop * nZM) / (nWXA + nZM), weight = nWXA + nZM), 
  fill = rgb(1, 0.6, 0.6), scale = "width", col = NA
  ) +
  xlab("simulated number of recombinants") + ylab("proportion of SNPs with properly inferred type") +
  theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))

pdf("../../figureS1.pdf", 6, 12)
grid.arrange(p1, p2, p3, nrow = 3)
dev.off()
rm(temp, p1, p2, p3)



# Evaluating the power of the assignments to sex chromosomes ------------------------

# to do so, we simulate autosomal contigs (1000 simulations per contig)
# this is done on a compute server due to the CPU time and RAM required
res <- mclapply(rep(0.5, 1000), simulate, mc.cores = 30)
saveRDS(res, "simuAutosomes_1000repl.RDS")

# we also save nrec only, as we transfer results to a PC with less RAM and over
# a slow internet connextion (thanks Covid) this is also why we use the RDS
# format and convert nrec to integer (32 bits instead of 64 for real numbers)
# these results are also used to assigne "real" contigs to sex chromosomes (at
# step 05)
saveRDS(as.integer(res$nrec * 10^7), stri_c("nrec_1000repl.RDS"))

saveRDS(as.integer(res$nrecS * 10^7), "nrecS_1000repl.RDS")

# we now are on the home PC
# we import all the values of nrec obtained by the simulations assuming autosomes
nrecs <- lapply(list.files(pattern = "nrec[S]*_1000repl"), readRDS)

# we use these results to generate different nrec quantiles for each contig
# these quantiles correspond to the significance levels we investigate

probs <- c(1/1000, 5/1000, 1/100, 5/100)

# we prepare column names corresponding to the above quantiles
qNames <- c("q001", "q005", "q01", "q05")

nrecQuantiles <- function(nrec) {
  # we add a column indicating the contig for each value of nrec
  simu <- data.table(contig = rep(uContigs, 1000), nrec = nrec / 10^7)

  # and we compute the quantiles for each contig
  stats <- simu[, .(q = quantile(nrec, probs)), by = contig]
  stats[, qName := rep(qNames, length.out = .N)]
  
  # we place the different quantiles in different columns
  stats <- dcast(stats, contig ~ qName, value.var = "q")
  stats[match(res$contig, contig), -"contig"]
}

quantiles <- lapply(nrecs, nrecQuantiles)


# We compute the proportion of contigs assigned to sex chromosomes
# depending on d, the significance level an whether SNP type was inferred

# we rbind the data we need for both simulations in a single table
bothRes <- rbind(
  data.table(res[twoFams, .(d, nrec)], quantiles[[1]][twoFams], inferred = T),
  data.table(res[twoFams, .(d, nrec = nrecS)], quantiles[[2]][twoFams], inferred = F)
)

# which facilitates the computation of porportions
props <- bothRes[, .(
  q001 = mean(nrec < q001), q005 = mean(nrec < q005),
  q01 = mean(nrec < q01), q05 = mean(nrec < q05)
), by = .(d, inferred)]

setorder(props, inferred, d)

# we also estimate the porportion of assigned contigs among those at or below a
# certain distance to the SDR
cumu <- lapply(ds, function(x) {
  bothRes[d <= x, .(
    q001 = mean(nrec < q001), q005 = mean(nrec < q005),
    q01 = mean(nrec < q01), q05 = mean(nrec < q05)
  ), by = inferred]
})

cumu <- data.table(d = rep(ds, each = 2), rbindlist(cumu))
setorder(cumu, inferred, d)

pdf("../../figureS3.pdf", 7, 9)

par(mfrow = 2:1, lwd = 1, mai = c(0.8, 0.8, 0.4, 0.4))
ggBackground(range(ds * 100), 0:1, xlab = "distance to the SDR (cM)", ylab = "proportion assigned to sex chromosomes")
with(props[inferred == T, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 1, lwd = 1.5))
with(props[inferred == F, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 3, lwd = 1.5))

legend(
  "topright", lty = c(1, 1, 1, 1, 1, 3), col = c(1:4, 1, 1), 
  legend = c(probs, "inferred SNP type", "true SNP type"), 
  cex = 0.7, bg = grey(0.92), box.col = "white"
  )

text(-7, 1.1, "A", xpd = NA, cex = 1.2, font = 2)

ggBackground(range(ds * 100), 0:1, xlab = "highest distance to the SDR (cM)", 
             ylab = "proportion assigned to sex chromosomes")
with(cumu[inferred == T, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 1, lwd = 1.5))
with(cumu[inferred == F, ], matlines(d * 100, cbind(q001, q005, q01, q05), lty = 3, lwd = 1.5))
text(-7, 1.1, "B", xpd = NA, cex = 1.2, font = 2)

dev.off()

diffProps <- dcast(cumu[, .(d, inferred, q001)], 
                   d ~ ifelse(inferred, "YES", "NO"), value.var = "q001")
diffProps[which.max(abs(YES - NO))]




# investigating the reduction of crossover rates near the SDR -------------------------------
# we simulate contigs assuming uniform recombination rates
# For this, one simulation is enough (each contig is assigned only one d)
res <- simulate(d = runif(max(uContigs), min = 0, max = 0.5))

# we only save nrec, which we use to make figure 4 at stage 05
saveRDS(res$nrec, "nrecContinuousCrossovers.RDS")

