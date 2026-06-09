#!/usr/bin/env Rscript

# Libraries
library("DNAcopy")
library("jsonlite")

# Helper function to get command-line arguments by flag name
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

# Parse command-line arguments
args <- commandArgs(TRUE)

infile <- get_arg("--infile")
outfile <- get_arg("--outfile", "cbs_out.json")
# CNA parameters
seed <- as.numeric(get_arg("--seed", 42))
id <- get_arg("--id", "sample")
# CBS parameters
alpha <- as.numeric(get_arg("--alpha", 0.001))
nperm_arg <- get_arg("--nperm")
nperm <- if (is.null(nperm_arg) || is.na(nperm_arg)) 10000 else as.numeric(nperm_arg)
pmethod <- get_arg("--p_method", "hybrid")
min_width_arg <- get_arg("--min_width")
min_width <- if (is.null(min_width_arg) || is.na(min_width_arg)) 2 else as.numeric(min_width_arg)
kmax_arg <- get_arg("--kmax")
kmax <- if (is.null(kmax_arg) || is.na(kmax_arg)) 25 else as.numeric(kmax_arg)
nmin_arg <- get_arg("--nmin")
nmin <- if (is.null(nmin_arg) || is.na(nmin_arg)) 200 else as.numeric(nmin_arg)
eta_arg <- get_arg("--eta")
eta <- if (is.null(eta_arg) || is.na(eta_arg)) 0.05 else as.numeric(eta_arg)
trim_arg <- get_arg("--trim")
trim <- if (is.null(trim_arg) || is.na(trim_arg)) 0.025 else as.numeric(trim_arg)
verbose_arg <- get_arg("--verbose")
verbosity <- if (is.null(verbose_arg) || is.na(verbose_arg)) 0 else as.numeric(verbose_arg)

# Load input as dataframe
input <- jsonlite::fromJSON(infile, simplifyDataFrame = TRUE)

# CBS
set.seed(seed)
cna_object <- CNA(
  # a vector or matrix of data from array-CGH, ROMA, or other copy number experiments.
  # If it is a matrix the rows correspond to the markers and the columns to the samples.
  input$ratios,
  # the chromosomes (or other group identifier) from which the markers came. 
  # Vector of length same as the number of rows of genomdat. If one wants the chromosomes to be ordered in the natural order,
  # this variable should be numeric orordered category.
  input$chrom,
  # the locations of marker on the genome. Vector of length same as the number of
  # rows of genomdat. This has to be numeric.
  input$maploc,
  data.type = "logratio",
  sampleid = id,
  presorted = TRUE
)
segments <- segment(
  # an object of class CNA
  cna_object,

  # a vector of weights for the probes. The weights should be inversely proportional
  # to their variances. Currently all weights should be positive i.e. remove probes
  # with zero weight prior to segmentation.
  weights = input$weights,

  # significance levels for the test to accept change-points.
  alpha = alpha,

  # number of permutations used for p-value computation.
  nperm = nperm,

  # method used for p-value computation. For the "perm" method the p-value is
  # based on full permutation. For the "hybrid" method the maximum over the entire
  # region is split into maximum of max over small segments and max over the rest.
  # Approximation is used for the larger segment max. Default is hybrid.
  p.method = pmethod,

  # the minimum number of markers for a changed segment. The default is 2 but
  # can be made larger. Maximum possible value is set at 5 since arbitrary widths
  # can have the undesirable effect of incorrect change-points when a true signal of
  # narrow widths exists.
  min.width = min_width,

  # the maximum width of smaller segment for permutation in the hybrid method.
  kmax=kmax,

  # the minimum length of data for which the approximation of maximum statistic
  # is used under the hybrid method. should be larger than 4*kmax
  nmin=nmin,

  # the probability to declare a change conditioned on the permuted statistic exceeding
  # the observed statistic exactly j (= 1,...,nperm*alpha) times.
  eta=eta,
  
  # the sequential boundary used to stop and declare a change. This boundary is a
  # function of nperm, alpha and eta. It can be obtained using the function "getbdry"
  # and used instead of having the "segment" function compute it every time it is
  # called.
  sbdry=NULL,

  # proportion of data to be trimmed for variance calculation for smoothing outliers
  # and undoing splits based on SD.
  trim = trim,

  # A character string specifying how change-points are to be undone, if at all. De-
  # fault is "none". Other choices are "prune", which uses a sum of squares criterion,
  # and "sdundo", which undoes splits that are not at least this many SDs apart.
  undo.splits = "none",

  # the proportional increase in sum of squares allowed when eliminating splits if
  # undo.splits="prune".
  #undo.prune = undo_prune,

  # the number of SDs between means to keep a split if undo.splits="sdundo".
  #undo.SD = undo_sd,

  # verbosity level. 0 for no output, 1 for current sample, 2 for current chromosome, 3 for current segment.
  verbose = verbosity,
)$out
# -> sample, chrom, loc.start, loc.end, num.mark, seg.mean

# Export to json
write_json(
  segments,
  path = outfile,
  dataframe = "rows",
  pretty = TRUE,
  auto_unbox = TRUE
)
