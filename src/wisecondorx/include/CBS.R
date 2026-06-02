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
id <- get_arg("--id", "sample")
seed <- as.numeric(get_arg("--seed", 42))
alpha <- as.numeric(get_arg("--alpha", 0.001))
outfile <- get_arg("--outfile", "cbs_out.json")
nperm_arg <- get_arg("--cbs_nperm")
nperm <- if (is.null(nperm_arg) || is.na(nperm_arg)) 10000 else as.numeric(nperm_arg)
min_width_arg <- get_arg("--cbs_min_width")
min_width <- if (is.null(min_width_arg) || is.na(min_width_arg)) 2 else as.numeric(min_width_arg)
undo_splits_arg <- get_arg("--cbs_undo_splits")
undo_splits <- if (is.null(undo_splits_arg) || is.na(undo_splits_arg)) "none" else as.character(undo_splits_arg)
undo_sd_arg <- get_arg("--cbs_undo_sd")
undo_sd <- if (is.null(undo_sd_arg) || is.na(undo_sd_arg)) 3 else as.numeric(undo_sd_arg)
cbs_smooth_arg <- get_arg("--cbs_smooth")
cbs_smooth <- if (is.null(cbs_smooth_arg) || is.na(cbs_smooth_arg)) FALSE else as.logical(cbs_smooth_arg)

# Load and pre-process input
input <- read_json(infile)
ratios <- input$ratios
weights <- input$weights

# Keep chromosomes with at least one finite non-zero ratio so ratios and
# weights remain aligned after filtering.
keep_chr <- vapply(
  ratios,
  FUN.VALUE = logical(1),
  FUN = function(chr) {
    chr_vals <- as.numeric(unlist(chr))
    any(is.finite(chr_vals) & chr_vals != 0)
  }
)

ratios <- ratios[keep_chr]
weights <- weights[keep_chr]

# Build flat vectors expected by DNAcopy
bins_per_chrom <- lengths(ratios)
chrom_vec <- rep(1:length(ratios), times = bins_per_chrom)
maploc_vec <- unlist(lapply(bins_per_chrom, seq_len), use.names = FALSE)
genomdat_vec <- as.numeric(unlist(ratios))
genomdat_vec[genomdat_vec == 0] <- NA
weights_vec <- as.numeric(unlist(weights))
weights_vec[!is.finite(weights_vec) | weights_vec <= 0] <-
  .Machine$double.xmin

# CBS: Run the Circular Binary Segmentation algorithm

# Construct a Copy Number Array (CNA) object, defining the logratio values, chromosome labels, and local coordinate positions
cna_object <- CNA(
  genomdat_vec,
  chrom_vec,
  maploc_vec,
  data.type = "logratio",
  sampleid = id
)
# Smooth the CNA data to detect and adjust single-point outliers if cbs.smooth is enabled
if (isTRUE(cbs_smooth)) {
  cna_object <- smooth.CNA(cna_object)
}

# Run the DNAcopy segment function
cna_object <-
  segment(
    cna_object,
    alpha = alpha,
    nperm = nperm,
    min.width = min_width,
    undo.splits = undo_splits,
    undo.SD = undo_sd,
    weights = weights_vec
  )$output

# Export to json
write_json(
  cna_object,
  path = outfile,
  dataframe = "rows",
  pretty = TRUE,
  auto_unbox = TRUE
)