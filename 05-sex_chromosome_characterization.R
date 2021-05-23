##%######################################################%##
#                                                          #
####            here we assign contig to sex            ####
####       chromosomes, investigate recombination       ####
####                  and SNP density                   ####
#                                                          #
##%######################################################%##


source("WZ_functions.R")

# We compare nrec to simulated data in order to assign contigs to WZ chromosomes ---------------

# we import the results from stage 03
probs <- fread("tables/countsAndPosteriors.txt")

# we add a logical column indicating whether a contigs has informative SNPs in both families
# we do it indirectly by assessing whether the max posterior is higher than 1/11 (the prior)
probs[, two := pmin(WXA_d.pM, WXA_s.pM, ZM_d.pM, ZM_s.pM) > 1 / 11]

# we do the same for contigs having SNPs at all (in at least one family)
# because only these contigs were used in simulations
probs[, one := pmax(WXA_d.pM, WXA_s.pM, ZM_d.pM, ZM_s.pM) > 1 / 11]


# we import all the values of nrec obtained by simulations
simu <- readRDS("simu/nrec_1000repl.RDS")
simu <- data.table(contig = rep(which(probs$one), 1000), nrec = simu / 10^7)

# we compute the nrec 1/1000 quantile for each contig
quantiles <- simu[, .(q001 = quantile(nrec, 1 / 1000)), by = contig]

# we add this value to the table of actual contigs to assign to sex chromosomes
probs[one == T, q001 := quantiles$q001]

# a contig is assigned to the WZ chromosomes if its nrec is lower than the quantile
probs[, WZ := nrec < q001]

# we plot the distribution of nrec ------------------------------
# the lengths of contig categories to show in the plot
lengths <- probs[two == T, round(c(sum(length), sum(length[WZ])) / 10^6, 1)]

# colors of the polygons (observed, assigned to sex chromosomes and simulated contigs)
colors <- c(grey(0.3), rgb(0.85, 0.85, 1, 1), rgb(0.9, 0.4, 0.4, 0.6))

pdf("../figure3.pdf", 8, 7)

# we prepare the plot area
d <- probs[two == T, density(nrec / 40 * 100, from = 0, to = 55, bw = 0.2)]
probs[two == T, plot(d, type = "n", bty = "n", ylim = c(0, 0.18), las = 1, 
                     main = "", xlab = "Inferred distance to the SDR (cM)")]

# we draw the vertical segments representing integer number of recombinant
abline(v = 0:30 / 40 * 100, col = "grey", lty = 3, lwd = 0.5)

# we plot the density for observed data (all contigs)
closedPolygon(d$x, d$y, border = NA, col = colors[1])

# to plot the distribution for assigned contigs, we need a scale factor for the Y axis
Yscale <- probs[two == T, mean(WZ)]
probs[two == T & WZ == T, densityPolygon(nrec / 40 * 100, from = 0, to = 55, bw = 0.2, 
                                         Yscale = Yscale, border = NA, col = colors[2])]

# we plot the distribution for simulated data
simu[contig %in% which(probs$two), 
     densityPolygon(nrec / 40 * 100, from = 0, to = 55, border = NA, col = colors[3])]

legend(
  "topleft",
  border = NA,
  legend = c("Simulated data", stri_c(c("All analyzed contigs", "Assigned to sex chromosomes"), " (~", lengths, " Mbp)")),
  fill = colors[c(3, 1, 2)], box.col = "white", bg = "white", y.intersp = 1.3, cex = 0.8, text.col = grey(0.3)
)

dev.off()

# we reclaim some RAM
rm(simu)

# we extrapolate chromosome length, given that not all contigs have SNPs in both families.


probs[, .(sum(length[two]), sum(length[which(WZ)]) / sum(length[two]) * sum(length))]



# investigating crossing-over rates ---------------------------------------

nrecSim <- readRDS("simu/nrecContinuousCrossovers.RDS")

# we create a table with the data we need, retrieved from "real" contigs
nrecSim <- data.table(nrecSim, probs[one ==T, .(CHROM, nrec, length, q001, two, WZ)])

# we assign simulated contigs to sex chromosomes and discard non-assigned contigs
wz <- nrecSim[nrecSim <= q001 & two == T]

# we obtain the total length of observed contigs assigned to sex chromosomes
L <- probs[WZ == T & two == T, sum(length)]

# we sample simulated contigs, ensuring that their total length is as close as possible to L
sampleContigs <- function() {

  # randomness is ensured by ordering the selected contigs randomly
  sel <- wz[sample(.N), .(length, nrecSim)]

  # we compute their cumulated length from the start of the table
  sel <- sel[, cl := cumsum(length)]

  # and select a sufficient number of contigs so that this cumulated length is the closest to L
  sel <- sel[1:which.min(abs(L - cl))]

  # we then compute the cumulated length of sampled contigs at every distance to the SDR (converted in cM)
  setorder(sel, nrecSim)
  dt <- sel[, .(cM = nrecSim / 40 * 100, len = cumsum(length) / 10^6)]
}

# we repeat the above 1000 times
samples <- replicate(1000, sampleContigs(), simplify = F)
samples <- rbindlist(samples)

# to compute quantiles for the envelope, we make 100 bins of genetic distance
bins <- quantile(samples$cM, seq(0, 1, 0.01))
samples[, bin := .bincode(cM, bins, include.lowest = T)]

# we compute the quantiles 
limits <- samples[, .(cM  = mean(cM), lower = quantile(len, 0.05), upper = quantile(len, 0.95)), by = bin]
setorder(limits, bin)

# we add a last line so that the envelope reaches the max genetic distance (and not just the mean of the last bin)
limits = limits[c(1:100, 100)]
limits[.N, cM := max(samples$cM)]

# we prepare the curve for observed data (without sampling)
setorder(nrecSim, nrec)
observed <- nrecSim[WZ == T & two == T, .(cM = nrec / 40 * 100, len = cumsum(length) / 10^6)]


pdf("../figure4.pdf", 6, 6)

ggBackground(observed[, .(cM, len)], ylab = "Cumulated contig length (Mbp)",
             xlab = "Inferred distance to the SDR (cM)", las = 1)
with(limits, polygon(x = c(cM, rev(cM)), y = c(lower, rev(upper)), 
                     border = NA, col = rgb(0.9, 0.4, 0.4, 0.6)))
lines(observed[, .(cM, len)])

dev.off()



# investigating the divergence of Z and W chromosomes via female heterozygous SNP density ---------------------------

# we import the number of heterozygous SNPs in mothers, per contig
# This table was generated in stage 1
SNPperContig <- fread("0-F0/heteroSNPs_byContigGQ40.txt")

# to estimate heterozygous SNP density per kb, we need to obtain the number of
# positions in the range of acceptable sequencing depths, for each contig. We use
# the script below
system("Rscript nbValidPositions.R 10")

# we import the results
validPositions <- fread("0-F0/nbValidPositionsMothers.txt")

# we merge both tables
SNPperContig <- merge(SNPperContig, validPositions, by = "CHROM", all = T)
SNPperContig[is.na(mothers), mothers := 0L]

# we add a heterozygous SNP density column to our general result table
probs[, c("het", "nPos") := 
        SNPperContig[match(probs$CHROM, CHROM), .(mothers * 1000 / nPos, nPos)]]

# we select contigs for which the max posterior probability of f is
# higher than 0.5 in all pools
used <- probs[, pmin(WXA_d.pM, WXA_s.pM, ZM_d.pM, ZM_d.pM) > 0.5]


# for the figure, we will compute heterozygous SNP density for different classes
# of nrec (on sex chromosomes) which we delineate by nrec deciles
nrecDeciles <- quantile(probs[WZ == T & used, nrec], probs = seq(0, 1, length.out = 11))

# we add a column to denote the class of nrec
probs[WZ == T & used, nrecClass := .bincode(nrec, nrecDeciles, include.lowest = T)]
meanNrecs <- probs[WZ == T & used, mean(nrec), by = .(class = nrecClass)]
probs[WZ == T & used, meanNrec := meanNrecs[match(nrecClass, class), V1]]

probs[, cM := nrec / 40 * 100]


pdf("../figure6.pdf", 7, 7)

layout(rbind(1:2), widths = 2:1)
par(mai = c(0.8, 0.8, 0.8, 0.4), cex.main = 0.9, oma = c(0, 0, 0, 0))

# we draw the scatter plot for each sex-chromosome contig
with(probs[WZ == T & used == T], ggBackground(
    cM, het, ylab = "Heterozygous SNP density in females (SNPs / kbp)", 
    xlab = "inferred distance to the SDR (cM)", ylim = c(0, max(het)),
    main = stri_c("Assigned to sex chromosomes (n=", sum(WZ), ")")
    ))

with(probs[WZ == T & used == T], points(cM, het, pch = 16, col = "salmon"))

# we add points (diamonds) for the nrec classes
distri <- probs[WZ == T & used == T, weighted.mean(het, nPos), by = meanNrec / 40 * 100]
distri <- probs[WZ == T & used == T, median(het), by = meanNrec / 40 * 100]
points(distri, pch = 23, cex = 1.5, bg = grey(1, 0.5))

# we prepare and draw the violin plot (but we don't use ggplot here)

# to properly "align" with the previous plot,
# we get the max Y value from that plot
maxHet <- probs[WZ == T & used == T, max(het)]

# we compute the density
d <- probs[WZ == F & used == T & nrec > 0, 
           density(het, from = 0, to = maxHet, weights = nPos / sum(nPos))]

# we prepare the XY coordinates of the violon (polygon)
X <- c(d$y, -rev(d$y))
Y <- c(d$x, rev(d$x))
par(mai = c(0.8, 0, 0.8, 0.4))

# we plot the density
ggBackground(
    range(X), range(Y), ylim = c(0, maxHet), axes = F, xlab = "density", 
    ylab = "", xgrid = F, main = stri_c("Other contigs (n=", d$n, ")")
    )

polygon(X, Y, border = NA, col = "salmon")

# we add a boxplot
with(probs[WZ == F & used == T & nrec > 0, ], boxplot(
    het, outline = F, add = T, at = 0, axes = F, col = "white", 
    boxwex = diff(range(X)) / 2, notch = T, lty = 1, staplewex = 0))

dev.off()


# we compute the correlation between heterozygous SNP density and distance to the SDR
probs[WZ == T & used == T, cor(het, cM)]

# we test the correlation
probs[WZ == T & used == T, cor.test(het, cM, method = "spearman", alternative = "less")]


# we report medians of contigs assigned to the SDR and other sex-chromosome contigs
probs[WZ == T & used, median(het), by = p > 0.5]

# we compare these medians with a Mann-Whitney test
probs[WZ == T & used, wilcox.test(het[p > 0.5], het[p <= 0.5], alternative = "greater")]

# and we compare heterozygous SNP density between contigs assigned to sex chromosomes and others
probs[, wilcox.test(het[used & WZ == T], het[used & WZ == F & nrec > 0], alternative = "two.sided")]

writeT(probs, "tables/countsAndPosteriors.txt")
