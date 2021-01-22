# computes sequencing depth on sliding windows in the F1 (and parents)

source("WZfunctions.R")

args <- commandArgs(trailingOnly = T)

# the first argument is the window size
window = as.integer(args[1])

# the the step size (distance between successive windows)
step = as.integer(args[2])

# then the minimum depth to consider (could be usefil in certain cases)
minDepth = as.integer(args[3])

# the number of processes to use
nCPUs <- as.integer(args[4])

# the arguments passed to samtools depth
arg = args[5]

# the bams to use
bams = args[7:length(args)-1]

# the output file name
out = args[length(args)]

# we retreive contig lengths
header = fread(cmd = paste("samtools view -H", bams[1]), col.names = c("col","contig","len"))
header[,len := as.integer(stri_sub(len, 4, nchar(len)))]
header[,contig := stri_sub(contig, 4, nchar(contig))]

# to selects contigs of sufficient length

header <- header[len >= window*2]


# we generate the samtools command that we will use for each contig. Note the use
# of the -r argument to target a contig via the bam indexes
commands <- header[, paste(
  "/usr/local/bin/samtools depth", arg, "-r",
  paste(contig, ":1-", len, sep = ""),
  paste(bams, collapse = " ")
  
)]



depthPerWindow <- function(command) {
  # processes depth data for a contig

  # we import the results from samtools depth
  depth <- fread(cmd = command, sep = "\t", header = F)
  if(nrow(depth) == 0L) return(data.table())
  # we create sliding windows of width 2000 and sliding by 500 bp. V2 is the position within the contig
  windows <- slidingWindows(depth$V2, window, step)

  # we add windows information to the results. This requires duplicating the rows
  # 4 times since there are four overlapping windows on a given position
  depth <- data.table(depth[rep(1:.N, window/step)], window = as.vector(windows))

  # we compute mean depths per window, and reports also the number of position used (.N)
  
  perWindow = sapply(
    stri_c("V", 4:ncol(depth)-1), 
    function(col) depth[, as.numeric(mean(get(col)[get(col) > minDepth])), by = .(contig = V1, window)]$V1
    )
  
  nPos = depth[,.N, by = .(contig = V1, window)]
  perWindow  = data.table(nPos, perWindow)
  setnames(perWindow, 4:ncol(perWindow), file_path_sans_ext(basename(bams)))
  
  setorder(perWindow, window)
  cat(".")
  perWindow
}


# we run the job in parallel processes
m <- mclapply(commands, depthPerWindow, mc.cores = nCPUs)

ok = sapply(m, is.data.table)

writeT(rbindlist(m[ok]), out)

if(any(!ok)) saveRDS(m[!ok], stri_c(out, ".errors.RDS"))
writeT(header[!ok], stri_c(out, ".errors.txt"))
