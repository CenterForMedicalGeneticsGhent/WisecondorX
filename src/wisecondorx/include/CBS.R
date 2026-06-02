# Parse the command-line arguments passed to the R session
args <- commandArgs(TRUE)

# Helper function to get command-line arguments by flag name
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

# Retrieve input file path from the command-line argument
in.file <- get_arg("--infile")

library("DNAcopy")
library("jsonlite")

# -----
# main: Core segmentation pipeline
# -----

# load input: Read and parse the input JSON file
input <- read_json(in.file)
# Extract and convert the 'ratios' list (containing log ratios per chromosome) into a single numeric vector
ratios <- as.numeric(unlist(input$ratios))
# Extract and convert the 'weights' list (containing weights per bin) into a single numeric vector
weights <- as.numeric(unlist(input$weights))

# Extract configuration values from command-line arguments
id <- get_arg("--id","sample")
seed <- as.numeric(get_arg("--seed", 42))
sex <- get_arg("--ref_sex")
alpha <- as.numeric(get_arg("--alpha", 0.001))
binsize <- as.numeric(get_arg("--binsize", 100000))
out.file <- get_arg("--outfile","cbs_out.json")
nperm_arg <- get_arg("--cbs_nperm")
nperm <- if (is.null(nperm_arg) || is.na(nperm_arg)) 10000 else as.numeric(nperm_arg)
min_width_arg <- get_arg("--cbs_min_width")
min.width <- if (is.null(min_width_arg) || is.na(min_width_arg)) 2 else as.numeric(min_width_arg)
undo_splits_arg <- get_arg("--cbs_undo_splits")
undo.splits <- if (is.null(undo_splits_arg) || is.na(undo_splits_arg)) "none" else as.character(undo_splits_arg)
undo_sd_arg <- get_arg("--cbs_undo_sd")
undo.SD <- if (is.null(undo_sd_arg) || is.na(undo_sd_arg)) 3 else as.numeric(undo_sd_arg)
cbs_smooth_arg <- get_arg("--cbs_smooth")
cbs.smooth <- if (is.null(cbs_smooth_arg) || is.na(cbs_smooth_arg)) FALSE else as.logical(cbs_smooth_arg)

# Determine the number of chromosomes to segment: 24 (1:22, X, Y) for male, or 23 (1:22, X) for female
if (sex == "M"){
    chrs = 1:24
} else {
    chrs = 1:23
}

# prepare for CBS: Set up the chromosome and position metadata for each bin

# Compute the number of bins present in each chromosome list from the input ratios
bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$ratios[x])))
# Compute the cumulative ending index positions of the chromosomes (starting with 0)
chr.end.pos <- c(0,cumsum(bins.per.chr))

# Mask ratio values equal to 0 as NA (to blacklist/ignore them in CBS)
ratios[ratios == 0] = NA # blacklist
# Replace non-finite (e.g. Inf, NaN) or non-positive weights with the minimum positive double value (double.xmin) to prevent division issues
weights[!is.finite(weights) | weights <= 0] = .Machine$double.xmin

# Convert the ratio vector into a data frame
for.cbs <- as.data.frame(ratios)
# Initialize vectors to repeat chromosome numbers and local bin coordinates for all data points
chr.rep <- c()
chr.rep.2 <- c()
for (chr in chrs){
  # Generate repeated chromosome number corresponding to the number of bins in that chromosome
  chr.rep <- c(chr.rep, rep(chr, chr.end.pos[chr + 1] - chr.end.pos[chr]))
  # Generate a 1-based local index position sequence for the bins within the current chromosome
  chr.rep.2 <- c(chr.rep.2, 1:(chr.end.pos[chr + 1] - chr.end.pos[chr]))
}
# Add the chromosome labels and local bin coordinate columns to the data frame
for.cbs$chromosome <- chr.rep; for.cbs$x <- chr.rep.2
# Reorder the columns to (chromosome, x, ratio) and name the ratio column 'y'
for.cbs <- for.cbs[, c(2,3,1)] ; colnames(for.cbs)[3] <- "y"

# Check for complete NA/chr: Omit any chromosomes that have only NA values

# Initialize a mask vector to track rows to keep
cbs.mask <- c()
# Loop over each chromosome to check for non-NA data
for (chr in chrs){
    # Find row indices corresponding to the current chromosome
    check <- which(for.cbs$chromosome == chr)
    # If the chromosome contains at least one non-NA value, add its indices to the mask
    if(!(all(is.na(for.cbs$y[check])))){
        cbs.mask <- c(cbs.mask, check)
    }
}
# Filter the data frame to keep only rows that passed the non-NA mask
for.cbs <- for.cbs[cbs.mask,]

# CBS: Run the Circular Binary Segmentation algorithm

# If a valid random seed was provided, set it to ensure reproducible permutation testing
if (!(is.na(seed) || seed == '')) {
  set.seed(seed)
}
# Construct a Copy Number Array (CNA) object, defining the logratio values, chromosome labels, and local coordinate positions
CNA.object <- CNA(for.cbs$y, for.cbs$chromosome, for.cbs$x, data.type = "logratio", sampleid = id)
# Smooth the CNA data to detect and adjust single-point outliers if cbs.smooth is enabled
if (isTRUE(cbs.smooth)) {
  CNA.object <- smooth.CNA(CNA.object)
}

# Run the DNAcopy segment function silently using the defined parameters and weights
CNA.object <-
  segment(
    CNA.object,
    alpha = as.numeric(alpha),
    nperm = as.numeric(nperm),
    min.width = as.numeric(min.width),
    undo.splits = as.character(undo.splits),
    undo.SD = as.numeric(undo.SD),
    verbose = 1,
    weights = weights[cbs.mask]
  )$output

# Subset the segmentation output columns to exclude sample ID (col 1) and number of markers (col 5)
CNA.object <- CNA.object[,-c(1,5)]
# Rename the columns to standard 'chr', 's' (start bin), 'e' (end bin), and 'r' (ratio)
colnames(CNA.object) <- c("chr", "s", "e", "r")

# Check if segment covers large NA regions. If so = split: Split segments containing long NA stretches

# Initialize a new data frame to hold the updated segments
new.CNA.object <- data.frame()

# Loop through each segment row in the CNA object
for (row.i in 1:nrow(CNA.object)){
  # Extract the start bin index of the current segment
  start.i = CNA.object$s[row.i]
  # Extract the end bin index of the current segment
  end.i = CNA.object$e[row.i]
  # Extract the chromosome data frame corresponding to the current segment's chromosome
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chr[row.i], ]
  # Retrieve the ratio values within the segment boundaries
  segment = sub.frame$y[start.i:end.i]
  
  # Calculate first-order differences of NA presence flags to locate transitions to/from NAs
  diff.na <- diff(is.na(segment), 1)
  
  # Identify positions where NAs start (change from FALSE to TRUE) shifted back to absolute indices
  start.pos <- which(diff.na == 1) + start.i - 1 # all consecutive NAs (start.pos)
  # Identify positions where NAs end (change from TRUE to FALSE) shifted back to absolute indices
  end.pos <- which(diff.na == -1) + start.i - 1 # all consecutive NAs (end.pos)
  
  # Select NA stretches that exceed the gap limit (e.g. 100 kb / binsize ratio)
  selection <- end.pos - start.pos > as.integer((binsize / 2000000) ** -1) # 100 kb -> 20 NA stretch: split
  
  # Filter start and end positions to keep only those exceeding the NA stretch limit
  start.pos <- start.pos[selection]
  end.pos <- end.pos[selection]
  
  # Construct new subsegment starting boundaries (original start and end of NA gaps)
  inverse.start.pos <- c(start.i, end.pos)
  # Construct new subsegment ending boundaries (start of NA gaps and original segment end)
  inverse.end.pos <- c(start.pos, end.i)
  
  # Filter subsegments to ensure they have positive lengths (at least two bins)
  selection <- inverse.end.pos - inverse.start.pos > 0 # segments should be at least two in length
  if (length(which(selection)) == 0){
      next
  }
  inverse.start.pos <- inverse.start.pos[selection]
  inverse.end.pos <- inverse.end.pos[selection]
  
  # Bind the split segments together, retaining chromosome number and segment ratio
  sub.frame <- cbind(CNA.object$chr[row.i], inverse.start.pos, inverse.end.pos, CNA.object$r[row.i])
  # Append the split subsegments to the new data frame
  new.CNA.object <- rbind(new.CNA.object, sub.frame)
  
}

# Rename columns of the updated segments data frame
colnames(new.CNA.object) <- c("chr", "s", "e", "r")
# Replace the original CNA object with the split segments data frame
CNA.object <- new.CNA.object

# Recalculate segmental ratios: Compute final segment means using weights and excluding NAs

# Attach the weights corresponding to the non-NA filtered data
for.cbs$w <- weights[cbs.mask]

# Loop through each row of the updated segments
for (row.i in 1:nrow(CNA.object)){
  # Extract ratios and weights for the segment chromosome
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chr[row.i], ]
  # Compute and update the segmental ratio as the weighted mean ratio of non-NA bins in that segment
  CNA.object$r[row.i] = weighted.mean(sub.frame$y[CNA.object$s[row.i]:CNA.object$e[row.i]],
                                      sub.frame$w[CNA.object$s[row.i]:CNA.object$e[row.i]],
                                      na.rm = T)
}

# Subtract 1 from the start bin index to convert from 1-based R index to 0-based Python index format
CNA.object$s <- CNA.object$s - 1 # Make python compliant
write_json(CNA.object, out.file)
