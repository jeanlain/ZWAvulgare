
## %######################################################%##
#                                                          #
####           Investigating sequencing depth           ####
####              ratio males/females (CQ)              ####
#                                                          #
## %######################################################%##


# We do it on a per-window basis

# We compute depth for the 6 F1 pools -------------------------------------

system('Rscript depthPerWindow.R 2000 500 0 10 "-Q 20" dWXA.bam sWXA.bam dZM.bam sZM.bam dBF.bam sBF.bam allContigsWindow2000.Q20.depth.txt')


# we import sequencing depth data on all F1 bams
depth <- fread("1-F1/allContigsWindow2000.Q20.depth.txt")



# Computing and plotting CQ ----------------------------------------------------------

# standardizing depths across bams
# we create a matrix to allow using colMeans() below (not possible on a data.table)
progenyDepths <- depth[, cbind(dWXA, sWXA, dZM, sZM, dBF, sBF)]
meanDepths <- colMeans(progenyDepths)
ratios <- max(meanDepths) / meanDepths
adjusted <- t(t(progenyDepths) * ratios)
depth[, c("dWXA", "sWXA", "dZM", "sZM", "dBF", "sBF") := data.table(adjusted)]

# computing chromosome quotient per F1 family
depth[, c("cqWXA", "cqZM", "cqBF") := .(sWXA / dWXA, sZM / dZM, sBF / dBF)]

# plotting chromosomes quotient and the difference (male - female) / (male + female)
# we also name thise difference with the cq prefix as we use the same code for plotting
CQs <- list(
  depth[, .(cqWXA, cqZM, cqBF)],
  depth[, .(
    cqWXA = (sWXA - dWXA) / (sWXA + dWXA),
    cqZM = (sZM - dZM) / (sZM + dZM),
    cqBF = (sBF - dBF) / (sBF + dBF)
  )]
)

xlim <- list(c(0, 2.2), c(-1, 1))

# the position of ther expected value of CQ using moth metrics
# (we draw a vertica line at this position)
v <- c(1, 0)

xlab <- c("Chromosome Quotient (male / female)", "Sequencing depth difference (male - female)/(male + female)")

pdf("../figureS4.pdf", 6, 8)

par(mfrow = c(2, 1))
for (i in 1:2) {
  dt <- CQs[[i]]


  # ploting CQ distributionq
  dt[depth[, dWXA + sWXA >= 3L] & cqWXA < 3, plot(
    density(cqWXA, from = xlim[[i]][1], to = xlim[[i]][2]),
    main = "", xlab = xlab[i], las = 1, bty = "n", lwd = 2
  )]

  dt[depth[, dWXA + sWXA >= 3L] & cqZM < 3, lines(
    density(cqZM, from = xlim[[i]][1], to = xlim[[i]][2]),
    col = "grey", lwd = 2
  )]

  dt[depth[, dWXA + sWXA >= 3L] & cqBF < 3, lines(
    density(cqBF, from = xlim[[i]][1], to = xlim[[i]][2]),
    col = "salmon", lwd = 2
  )]
  abline(v = v[i], lwd = 0.5, col = "grey")

  legend("topleft", col = c("black", "grey", "salmon"), legend = c("WXA", "ZM", "BF"), bty = "n", lwd = 2)
}

dev.off()


# selecting windows and contigs of low CQ --------------------------------------------------

# we add a logical column to flag windows with low CQ ratio and sufficient female depth
depth[, lowCQ := dWXA > 5L & dZM > 5L & dBF > 5L & cqWXA < 0.3 & cqZM < 0.3 & cqBF < 0.3 & N >= 1500]

# we group successive low-CQ windows into blocks
rows <- 2:nrow(depth)
blocks <- depth[c(F, lowCQ[rows] & !lowCQ[rows - 1L]), .(contig, start = window - 1000L)]
blocks$end <- depth[c(lowCQ[rows - 1L] & !lowCQ[rows], F), window + 1000L]
blocks[start < 0L, start := 0L]


# we retrieve the SNP results for these low-CQ contigs
probs <- fread("tables/countsAndPosteriors.txt")

# the number and length of contigs with low-CQ windows
probs[CHROM %chin% blocks$contig, .(.N, sum(length))]

# the length of low-CQ blocks
CQblockLenghts <- blocks[, .(len = sum(end - start)), by = contig]
sum(CQblockLenghts$len)

# the low-CQ contigs that are linked to the SDR
# actually, there is just one
probs[CHROM %in% blocks$contig & p > 0.5]

# we import gene data for figure 7
gff <- fread(
  "genome.gff",  # this is the A.vulgare gff (from ncbi)
  skip = 1, header = F, drop = c(2, 6, 8),
  col.names = c("contig", "type", "start", "end", "strand", "description")
)

# we extain gene names and compute the mid position of genes
gff[type == "gene", c("gene", "mid") := 
      .(splitToColumns(description, "Name=", 2), (start + end) / 2)]

# we extract gene id
gff[, id := stri_extract_first(description, regex = "Avbf_[0-9]+")]
geneForID <- gff[type == "gene", setNames(gene, id)]

# we add the gene name to each CDS
gff[type == "CDS" | type == "exon", gene := geneForID[id]]


# the genes that intersect with low-CQ blocks
lowCQgenes <- gff[type == "gene" & intersectsRanges(data.table(contig, start, end), blocks)]

# to plot the low-CQ contigs, we also need the data on informative SNPs
SNPs <- fread("tables/SNPcategory.txt")


# we plot data for contigs having genes in low-CQ blocks
pdf("../figure7.pdf", 6, 8)

par(mfrow = c(3, 1), mai = c(0.4, 0.7, 0.3, 0.2))
l <- lapply(lowCQgenes[, unique(contig)], plotContig, depth = depth, gff = gff, SNPs = SNPs)

dev.off()



# comparing results with those from our prevous study based
# on CQ and Y genome scan
# we import the names of previous candidate contigs
CQcandidates <- readLines("CQcandidates.txt")

# assess the assignment of these contigs to sex chromosomes
tab <- probs[, table(CHROM %in% CQcandidates, WZ = nrec < q001)]
fisher.test(tab)

# compare the distance to the SDR non-assigned contigs
# that belong to our preivous selection and those that do not
probs[nrec >= q001, median(nrec / 40 * 100), by = CHROM %in% CQcandidates]
probs[nrec >= q001, wilcox.test(nrec[CHROM %in% CQcandidates], nrec[!CHROM %in% CQcandidates])]

# investigating genes that overlap with the blocks of low CQ
lowCQgenes <- gff[type == "gene" & intersectsRanges(data.table(contig, start, end), blocks), ]
lowCQgeneInfo <- gff[type == "exon" & id %in% lowCQgenes$id, .(nbExons = .N, exonLength = sum(end - start + 1L)), by = id]
lowCQgenes <- merge(lowCQgenes, lowCQgeneInfo)
