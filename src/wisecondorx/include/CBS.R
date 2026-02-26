# -----
# args
# -----

args <- commandArgs(TRUE)
input_file <- paste0(args[which(args == "--infile")+1])
output_file <- paste0(args[which(args == "--outfile")+1])
seed <- as.numeric(paste0(args[which(args == "--seed")+1]))
sex <- paste0(args[which(args == "--sex")+1])
alpha <- as.numeric(paste0(args[which(args == "--alpha")+1]))
binsize <- as.numeric(paste0(args[which(args == "--binsize")+1]))

# -----
# lib
# -----

suppressMessages(library("DNAcopy"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

# load input json
input <- read_json(input_file)
ratios <- as.numeric(unlist(input$ratios))  # Array of ratio per bin, per chromosome
weights <- as.numeric(unlist(input$weights)) # Array of weights per bin, per chromosome

if (sex == "M"){
    chrs = 1:24
} else {
    chrs = 1:23
}

# prepare for CBS

# Count the number of bins per chromosome to determine the length of each chromosome's data
bins_per_chr <- sapply(chrs, FUN = function(x) length(unlist(input$ratios[x])))

# Calculate cumulative sums to define the 0-indexed bounds for each chromosome within the flattened arrays
chr.end.pos <- c(0,cumsum(bins_per_chr))

# Replace ratios of 0 with NA to blacklist these bins, excluding them from the segmentation process
ratios[ratios == 0] = NA # blacklist

# Replace weights of 0 with a near-zero value (1e-99) because DNAcopy cannot handle 0 or NA weights.
# Note: 1^-99 previously evaluated to 1 in R, which incorrectly gave these bins full weight. This is fixed to 1e-99.
weights[weights == 0] = 1e-99 # omit DNAcopy weirdness -- weight cannot be NA or 0

for.cbs <- as.data.frame(ratios)
chr.rep <- c()
chr.rep.2 <- c()
for (chr in chrs){
  chr.rep <- c(chr.rep, rep(chr, chr.end.pos[chr + 1] - chr.end.pos[chr]))
  chr.rep.2 <- c(chr.rep.2, 1:(chr.end.pos[chr + 1] - chr.end.pos[chr]))
}
for.cbs$chromosome <- chr.rep; for.cbs$x <- chr.rep.2
for.cbs <- for.cbs[, c(2,3,1)] ; colnames(for.cbs)[3] <- "y"

# Check for complete NA/chr
cbs.mask <- c()
for (chr in chrs){
    check <- which(for.cbs$chromosome == chr)
    if(!(all(is.na(for.cbs$y[check])))){
        cbs.mask <- c(cbs.mask, check)
    }
}
for.cbs <- for.cbs[cbs.mask,]

# CBS

if (!(is.na(seed) || seed == '')) {
  set.seed(seed)
}
# Setup the CNA object with the data for segmentation
CNA.object <- CNA(for.cbs$y, for.cbs$chromosome, for.cbs$x, data.type = "logratio", sampleid = "X")

# Perform Circular Binary Segmentation. 
# Fixed bug: changed $output_file back to $output as that's what DNAcopy::segment actually returns.
# Added sink() to silence the verbose output of segment() to keep logs clean
f = file()
sink(file=f) ## silence output
CNA.object <- invisible(segment(CNA.object, alpha = as.numeric(alpha), verbose=1, weights=weights[cbs.mask])$output)
sink() ## undo silencing
close(f)

# Remove the sampleID (1) and num.marks (5) columns
CNA.object <- CNA.object[,-c(1,5)]
colnames(CNA.object) <- c("chromosome", "start", "end", "ratio")

# Check if segment covers large NA regions. If so = split

new.CNA.object <- data.frame()

for (row.i in 1:nrow(CNA.object)){
  # Fixed column names: updated `$s`, `$e`, `$chr` to `$start`, `$end`, `$chromosome` 
  # safely matching the explicit `colnames()` assignment above, preventing bugs from partial matching.
  start.i = CNA.object$start[row.i]
  end.i = CNA.object$end[row.i]
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chromosome[row.i], ]
  segment = sub.frame$y[start.i:end.i]
  
  diff.na <- diff(is.na(segment), 1)
  
  start.pos <- which(diff.na == 1) + start.i - 1 # all consecutive NAs (start.pos)
  end.pos <- which(diff.na == -1) + start.i - 1 # all consecutive NAs (end.pos)
  
  selection <- end.pos - start.pos > as.integer((binsize / 2000000) ** -1) # 100 kb -> 20 NA stretch: split
  
  start.pos <- start.pos[selection]
  end.pos <- end.pos[selection]
  
  inverse.start.pos <- c(start.i, end.pos)
  inverse.end.pos <- c(start.pos, end.i)
  
  selection <- inverse.end.pos - inverse.start.pos > 0 # segments should be at least two in length
  if (length(which(selection)) == 0){
      next
  }
  inverse.start.pos <- inverse.start.pos[selection]
  inverse.end.pos <- inverse.end.pos[selection]
  
  # Fixed column names: updated `$chr` and `$r` to explicitly match `$chromosome` and `$ratio` column names.
  sub.frame <- cbind(CNA.object$chromosome[row.i], inverse.start.pos, inverse.end.pos, CNA.object$ratio[row.i])
  new.CNA.object <- rbind(new.CNA.object, sub.frame)
  
}

colnames(new.CNA.object) <- c("chromosome", "start", "end", "ratio")
CNA.object <- new.CNA.object

# Recalculate segmental ratios

for.cbs$w <- weights[cbs.mask]

for (row.i in 1:nrow(CNA.object)){
  # Extract data corresponding to the current segment's chromosome to safely recalculate ratios 
  # Fixed column names: changed `$chr`, `$r`, `$s`, `$e` to full names to avoid partial matching bugs.
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chromosome[row.i], ]
  CNA.object$ratio[row.i] = weighted.mean(sub.frame$y[CNA.object$start[row.i]:CNA.object$end[row.i]],
                                      sub.frame$w[CNA.object$start[row.i]:CNA.object$end[row.i]],
                                      na.rm = T)
}

# Adjust 1-indexed R coordinates to 0-indexed Python coordinates
CNA.object$start <- CNA.object$start - 1 # Make python compliant

# Write output
write_json(CNA.object, output_file)
