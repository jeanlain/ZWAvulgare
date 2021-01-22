
##%######################################################%##
#                                                          #
####    We locate genomic regions which are suitable    ####
####   for the SDR, taking advantage of both families   ####
#                                                          #
##%######################################################%##

# we reconstruct parental haplotypes to find SNPS that have recombined
# with the SDR during the divergence of the WXA and ZM lines
# to do this, we retreive useful SNPs in the parents

source("WZ_functions.R")

# we import the parental genotype
# All this could have been performed at stage 1 (to avoid reimporting the
# table), but we did not think of all this at this time
vcf = fread ("0-F0/GATK_filtered_parents_snps.table" , sep = "\t" , header =T , 
             drop = c("ID" , "ALT" , "FILTER" , "INFO"  , "fWXA.PL" , "fZM.PL",  "mWXA.PL" , "mZM.PL" ))

# we select biallelic SNPs (just one coma in the AD column)
vcf = vcf[stri_count(fWXA.AD, fixed = ",") == 1L] 

baseTypes = c("A","C","G","T")

# we split alleles in genotypes and encodes them as integers
bases = lapply(
    vcf[,.(fWXA.GT, mWXA.GT, fZM.GT, mZM.GT)], 
    function(col) cbind(chmatch(stri_sub(col, 1, 1), baseTypes), chmatch(stri_sub(col, 3, 3), baseTypes))
    )

# we put the result in a single matrix
bases = do.call(cbind, bases)

# we name columns
# "m" for "mother" and "x" for first family (WXA). "f" for father and "y" for family ZM
# 1 for first allele, 2 for second allele
colnames(bases) = c("mx1","mx2","fx1", "fx2", "my1","my2","fy1","fy2")

# we only keep SNPs without missing bases
allOk = rowSums(is.na(bases)) == 0L  

# we will retain only SNPs for which each allele is found in several haplotypes
# or else such SNP cannot inform about recombination

# We get one of the two alleles of each SNP
A1 = bases[,1]  

# we compute the number of haplotypes carrying this allele
nHapl = rowSums(bases == A1, na.rm = T)

# and the minimum number of haplotypes carying either allele, for each SNP
nHapl = pmin(nHapl, ncol(bases) - nHapl) 

# so wee keep SNPs for which each allele found in at least 2 haplotypes
# and without missing data
genotypes = data.table(vcf[nHapl >= 2L & allOk, .(
  CHROM, POS, mx.DP = fWXA.DP, fx.DP = mWXA.DP, my.DP = fZM.DP, fy.DP = mZM.GQ, 
  mx.GQ = fWXA.GQ, fx.GQ = mWXA.GQ, my.GQ = fZM.GQ, fy.GQ = mZM.GQ
  )], bases[keep,])

writeT(genotypes, "0-F0/selectedGenotypes.txt")



# we infer if SNPs recombined wiht the WZ locus and find "causal" SNPs ------------------------
# we do it for all SNPs, even those on autosomes !
# we always prefer filtering SNPs at the last moment

genotypes = fread("0-F0/selectedGenotypes.txt")

# for this we need information of SNP type, so we import
# the SNP table we saved at step 03

vcf = fread("tables/vcfAndF1countsLikelihoodQ20.txt")

# to phase maternal haplotypes, we merge the parental genotypes with our table
# of informative SNPs (in particular, we need SNP types). We first put both
# families at the same row. We use merge instead of dcast as we prefer grouping
# columns per family rather than variable type
# we turn the "retained" column into integer to facilitate certain instructions
# r.x = 1 will mean that the SNP was discarded in family x
# r.x = 2 will mean that the SNP was retained in this family
mVCF = merge(vcf[fam == "WXA", .(CHROM, POS, A1, A2, r = as.integer(retained)+1L, type)], 
             vcf[fam == "ZM", .(CHROM, POS, A1, A2, r = as.integer(retained)+1L, type)], 
             by = c("CHROM","POS"), all = T) 
# in the above, WXA-related columns get the suffix ".x" and ZM ".y"


# r.x = 0 will mean that the SNP was not informative in family x
mVCF[is.na(r.x), r.x := 0L] ; mVCF[is.na(r.y), r.y := 0L] 


# we write the SNP positions to disk...
writeT(mVCF[,.(CHROM, POS)], "tables/informativeSNPs.txt")

# ...because we retreive data from F1s of the BF family for these SNPs
# as causal SNPs must be so in all three families (see paper)
# we do this for all SNPs even though only very few
# will be credible candidates. Just in case.

bamFiles = list.files("1-F1", pattern = "BF2875")

# we launch the script to scan each F1 bam in parallel (using 10 CPUs for each)
m <- mcMap(
    function(posFile, bamFile) {
        system(paste("Rscript scanBamAtPositions.R", posFile, bamFile, 10))
    },
    rep(posFiles, each = 2), # each position file is used for two bam files (two pools)
    bamFiles,
    mc.cores = length(bamFiles)
)

# we don't import results just yet


# we merge parental genotypes with the mVCF table
genotypes = merge(genotypes, mVCF, by = c("CHROM","POS"), all.x = T)
genotypes[is.na(r.x), r.x := 0L] ; genotypes[is.na(r.y), r.y := 0L] 

# we exclude informative SNPs that we previously excluded as unreliable
genotypes = genotypes[r.x != 1L & r.y != 1L]

# we will exclude SNPs of excessive sequencing depth in any parent
depthOk = vapply(
    genotypes[,.(mx.DP, fx.DP, my.DP, fy.DP)], 
    function(x) x <= quantile(x, 0.95),logical(nrow(genotypes))  
    )   

# but we don't exclude informative SNPs that we retained earlier
depthOk = rowSums(depthOk) == 4L | genotypes[,r.x + r.y > 0L]

# we now select SNPs according to our criteria on depth and genotype quality.
# Note that we don't exclude informative SNPs that we retained at step 03 (r.x == 2L)
# for non-informative SNPs, mothers must ne homozygous
genotypes = genotypes[depthOk & (r.x == 2L | (mx.GQ >= 40L & fx.GQ >= 40L & mx1 == mx2)) & 
                          (r.y == 2L | (my.GQ >= 40L & fy.GQ >= 40L & my1 == my2))]  

rm(depthOk) ; gc()

# we compute the genomic position of SNPs 
contigLengths = readNamedVector("contigLengths.txt")
contigStarts <- setNames(c(0L, cumsum(contigLengths[-length(contigLengths)])), names(contigLengths))
genotypes[, genomePos := POS + contigStarts[CHROM]]

# we phase maternal haplotypes according to SNP types
# in practice, we place the W-linked allele at the left column
# for type-1 SNPs, allele mx1 (W-linked) becomes allele A1 (maternal allele)
genotypes[type.x == 1L & r.x == 2L, c("mx1", "mx2") := .(A1.x, A2.x)]

# for type-2 SNPs, allele mx1 becomes allele A2 (paternal allele)
genotypes[type.x == 2L & r.x == 2L, c("mx1", "mx2") := .(A2.x, A1.x)]

# we do the same for the other family
genotypes[type.y == 1L & r.y == 2L, c("my1", "my2") := .(A1.y, A2.y)]
genotypes[type.y == 2L & r.y == 2L, c("my1", "my2") := .(A2.y, A1.y)]

# we determine if SNPs have recombined with the SDR. To so wo, we create 2-SNP
# haplotypes where 1 = W and 2 = Z
# so the haplotype is encoded as a 2-digit integer (one digit per allele)
wzHapl = t(c(1, 2, 2, 2, 1, 2, 2, 2)*10 + genotypes[, t(cbind(mx1, mx2, fx1, fx2, my1, my2, fy1, fy2))])

# we add an integer column indicating the category of SNPs
# 1 (TRUE) indicates a recombinant SNP
genotypes[,rec := as.integer(countHapl(wzHapl) == 4L)] 

# 2 when a SNP shows no evidence for recombination and is not informative in both families
genotypes[rec != 1L & r.x + r.y < 4L, rec := 2L]       

# 3 if the SNP is informative in both families
genotypes[rec != 1L & r.x + r.y == 4L, rec := 3L]       

# 4 for potentially causal SNPS
genotypes[rec != 1L & r.x + r.y == 4L & type.x == 1L & type.y == 1L, rec := 4L]    

# if a SNP is not recombinant and is variable in just one family, 
# we consider the SNP as not very relevant with respect to recombination with the SDR
# (the rarer allele may have appeared quite recently in a family)

# this logical is TRUE for SNPs that are variable in both families
in2fams = genotypes[,pmin(mx1, mx2, fx1, fx2) != pmax(mx1, mx2, fx1, fx2) &
                  pmin(my1, my2, fy1, fy2) != pmax(my1, my2, fy1, fy2)]

genotypes[!in2fams & rec != 1L, rec := NA]


# we import data from the BF family for better assessment of causal SNPs
BFscans = lapply(list.files("1-F1/", pattern = "BF2875", full.names = T), readRDS)		

# we retreive the number of reads covering each base in the BF pools
baseCounts = lapply(BFscans, baseMatrix, genotypes$genomePos, minMapQ = 20L)			

# we will count reads carrying the maternal (A1) and paternal (A2) alleles. We
# use the reference alleles of WXA unless a SNP is not informative in this
# family. In the end, it doesn't matter because causal SNPs must have the
# same paternal and maternal alleles in both families
A1 = genotypes[,ifelse(r.x == 2L, A1.x, A1.y)]
A2 = genotypes[,ifelse(r.x == 2L, A2.x, A2.y)]

# if the SNP was not informative, we aribrarily use alleles 1 and 2 (doesn't
# matter since we don't need BF data for these, but it simplifies the code)
A1[is.na(A1)] = 1L
A2[is.na(A2)] = 2L

rows = 1:nrow(genotypes)
alleleCounts = lapply(baseCounts, function(x) {	
  dt = data.table(C = x[cbind(rows, A1)], R = x[cbind(rows, A2)], tot = rowSums(x)) 
  dt[, R := C + R] 
  dt })

# we name data.tables "d" for daugthers and "s" for sons
names(alleleCounts) = c("d","s")   

# so these letters become prefixes by the data.table call below
genotypes = data.table(genotypes, do.call(data.table, alleleCounts))

rm(BFscans, baseCounts, A1, A2) ; gc()

# we compute the product of p05 between pool for the BF SNPs, considering all SNPs as type 1
genotypes[, p := pFgivenCountsSNPs(d.C, d.R, group = 1L) * pFgivenCountsSNPs(s.C, s.R, group = 2L)]	

# if alleles A1 and A2 do not consitute the 80% of reads in BF, we ignore the result
genotypes[d.R < d.tot * 0.8 | s.R < s.tot * 0.8, p := NA]

# potentially causal SNPs are downgraded if the BF data is not supportive
genotypes[rec == 4L & p < 0.01, rec := 3L]

# we add the SNP recombination category to the table of informative SNPs
# by default, we use category 2
mVCF[,rec := genotypes[match(mVCF$genomePos, genomePos), rec]]

# for informative SNPs that were not relevant to investigate recombination
# we use the rec category 2 (not shared, not indicating recombination)
mVCF[is.na(rec), rec := 2L]

# we save the SNP category to disk as we will need itfor figure 7
writeT(mVCF[r.x != 1L & r.y != 1L, .(CHROM, POS, rec)], "tables/SNPcategory.txt")

# we assign SNPs to recombination blocks -----------------------------------------

# we select rows corresponding to contigs that did not recombine with
# the SDR during crosses. We could have done that much earlier, but we always
# prefer applying filters at the last moment

# To select the relevant contigs, we import our main table of results
probs = fread("tables/countsAndPosteriors.txt")

sdr = genotypes[CHROM %in% probs[p > 0.5, CHROM]]

setorder(sdr, CHROM, POS)

# we retreive coordinates of genomic blocks (which are those of SNPs at the edges)
breaks = sdr[, breakPoints(
    data.table(rec, CHROM, genomePos, mx1, mx2, fx1, fx2, my1, my2, fy1, fy2)
    ), by = CHROM]$V1

# we sort the coordinates in ascending order to assign SNPs to blocks
# we do this by creating bins (not forgetting to add the two extreme coordinates)
breaks = sort(unique(c(1L, breaks, sum(contigLengths)+1L)))
sdr[,block := .bincode(genomePos, breaks = breaks, right = F)]

# in case a bin overlap two contis, we make sure it does not correspond to a unique block
contigBlock = sdr[,paste(CHROM, block)]
sdr[,block:= chmatch(contigBlock, unique(contigBlock))]


# we get relevant metrics for blocks
blocks = sdr[,data.table(
    start = min(POS), end = max(POS),              # positions of first and last SNP 
    data.table(rbind(tabulate(rec, nbins = 4))),   # and number of SNPs per category
    used = sum(!is.na(rec))                         # and number of SNPs that inform on recombination
    ), by = .(CHROM, block)]

# we select blocks composed of more than one SNP or that contain causal SNPs
blocks = blocks[start < end | V4 > 0L]


# we retrieve contigs that did not have SNPs that could be used to assess
# recombination with the SDR (hence not in the genotypes table) we will consider
# that there is no proof of recombination with the SDR for these contigs so we
# add them as single blocks constituing whole contigs, with just one SNP of rec
# category 2 (non-shared informative SNP not indicating recombination)
missingContigs = setdiff(probs[p > 0.5, CHROM], blocks$CHROM)
if (length(missingContigs)> 0L)
blocks = rbind(blocks, data.table(
    CHROM = missingContigs, block = 1L, start = 1L, end = 1L, 
    V1 = 0L, V2 = 1L, V3 = 0L, V4 = 0L, used = 1L
    ))


# we compute block start and end positions (as blocks must be contiguous)
n = 2:nrow(blocks)

# block boundaries are the midpoint between the SNP coordinates of adjacent blocks
blocks[,end := c(as.integer((end[n-1] + start[n])/2), end[.N])]
blocks[, start := c(1L, end[n-1]+1L)]

# the start of the first block of a contig is moved to 1 and the end of the last
# block is moved to the length of a contig
blocks[!duplicated(CHROM), start := 1L]
blocks[!duplicated(CHROM, fromLast = T), end := contigLengths[CHROM]]


# we categorize blocks. We attribute number 1 whenever a block harbors 2 or more
# than 50% recombinant SNPs, else:
# 2 if it harbors no shared informative SNPs
# 3 if it harbors shared non-causal SNPs
# 4 if it harbors causal SNPs
blocks[, cat := ifelse(V1 >= 2L | V1 / used >= 0.5, 1L, ifelse(V4 > 0L, 4L, ifelse(V3 > 0L, 3L, 2L)))]
blocks[is.na(cat), cat := 2L]

# the proportion of bocks haboring shared informative SNPs
# among recombinant blocks
blocks[,mean(V1 > 1 & V3 + V4 > 0)]

# as we will plot longer contigs first, we retrieve contig lengths
blocks[,len := contigLengths[CHROM]]
setorder(blocks, -len, CHROM)

# blocks will be plotted as rectangles, their top and bottom
# edges are simply their start and end positions in the contigs 
#(-1 for start to ensure that rectangles are contiguous)
blocks[,c("bottom", "top") := .((start-1L)/1000, end/1000)]

# the X position of the left edges are determined by their 
# host contig (+1 per new contig)
row = 2:nrow(blocks)
blocks[,left := cumsum(c(0L, CHROM[row] != CHROM[row-1]))]

# we prepare the second barplot (the inset) showing total lengths
totLengths = blocks[cat != 0L,.(len = sum(end-start+1L)/1000), by = cat]

# the explanatory labels of this barplot
labels = c("contains SNPs indicating recombination with the SDR",
           "contains no SNPs that are informative in both families",
           "contains SNPs that are informative in both families \nbut no causal SNP",
           "contains potentially causal SNPs")

# this second barplot will appear in the white space of the main plot
# between these Y coordinates
maxY = 200
minY = 100

# so we need to compute the top and bottom of sectors so they fit in there
ratio = (maxY - minY)/sum(totLengths$len)
totLengths[,c("bottom","top") := .(cumsum(c(minY, ratio*len[-.N])), cumsum(ratio*len) + minY)]

# we compute the mid Y position of sectors as we will indicate sizes there
# for the shorter sector, we use the top instead of the middle (as it is too narrow)
totLengths[,mid := ifelse(top -bottom >5, (top + bottom)/2, top)]

# we compute the Y position of the explanatory labels (that are equidistant)
totLengths[,textY := seq(min(mid), max(mid), length.out = .N)]

# the colors corresponding to block categories, from 1 to 4
cols =  c("salmon", "grey", "mediumseagreen","olivedrab2")

pdf("../figure5.pdf", width = 8, height =7)

par(mai = c(0.2, 0.8, 0.4, 0.35))

# the plot showing all contigs with colors for blocks of different categories
# we prepare an empty plot to add rectangles to it
with(blocks, plot(
  x = range(left + 1L), y = range(c(top, bottom)), bty = "n", type = "n",
  ylab = "Contig length (kbp)", las =1,
  xlab = "",  xlim =c(2, max(left) +1),
  xaxt = "n"
))

with(blocks, rect(
  xleft = left, ybottom = bottom, xright = left + 0.9, # we use a 0.1 space between bars
  ytop = top, border = NA, col = cols[cat]
))

# we determine the left X position of the labels of the inset
# so that the longest ends right at the X position of the last contig
textR = max(blocks$left) +1L - max(strwidth(labels, cex = 0.7))

# so we can draw the barplot of the inset 
# (we use rectangles as it is more flexible than adding a regular barplot)
with(totLengths, rect(
  xleft = textR - 30, ybottom = bottom, xright = textR - 10, 
  ytop = top, border =NA, col = cols[cat])
  )

# we draw lines connecting the sectors to the captions
l = with(totLengths, Map(function(middle, Y, col) lines(
  x = textR - 10 + c(0.5, 2, 9.5, 11), 
  y = c(middle, middle, Y, Y), 
  lend = 1, col = col
  ), mid, textY, cols[cat]) 
)

# we write the explanatory captions
totLengths[, text(x = textR, y = textY, labels = labels, pos = 4, col = grey(0.3), cex = 0.7)]

# and the lengths within sectors
totLengths[, text(x = textR -20, y = mid, labels = paste(round(len), "kbp"), col = grey(0.3), cex = 0.7)]

dev.off()


# we compute the length of certain block categories
blocks[, sum(end-start + 1L)]
blocks[cat != 1, sum(end-start + 1L)]



# genes that intersect with regions containing causal SNPs ---------------------------

# importing the gff of A. vulgare 
gff = fread(
    "genome.gff", skip =1, header = F, drop = c(2, 6, 8), 
    col.names = c("contig","type","start", "end", "strand","description"))


causalGenes = gff[type == "gene" & intersectsRanges(
  data.table(contig, start, end), blocks[cat == 4L, .(CHROM, start, end)]
  )]


