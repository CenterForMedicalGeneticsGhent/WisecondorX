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
undo_prune_arg <- get_arg("--undo_prune")
undo_prune <- if (is.null(undo_prune_arg) || is.na(undo_prune_arg)) 0.05 else as.numeric(undo_prune_arg)
undo_splits_arg <- get_arg("--undo_splits")
undo_splits <- if (is.null(undo_splits_arg) || is.na(undo_splits_arg)) "none" else as.character(undo_splits_arg)
undo_sd_arg <- get_arg("--undo_sd")
undo_sd <- if (is.null(undo_sd_arg) || is.na(undo_sd_arg)) 3 else as.numeric(undo_sd_arg)
verbose_arg <- get_arg("--verbose")
verbosity <- if (is.null(verbose_arg) || is.na(verbose_arg)) 0 else as.numeric(verbose_arg)

# Load and pre-process input
input <- read_json(infile)
genomedat <- input$genomedat
chroms <- input$chroms
maploc <- input$maploc
weights <- input$weights

# CBS
set.seed(seed)
cna_object <- CNA(
  genomedat,
  chroms,
  maploc,
  data.type = "logratio",
  sampleid = id,
  presorted = FALSE
)
segments <- segment(
  cna_object, 
  weights = weights,
  alpha = as.numeric(alpha),
  nperm = nperm,
  p.method = pmethod,
  min.width = min_width,
  kmax=kmax,
  nmin=nmin,
  eta=eta,
  sbdry=NULL,
  trim = trim,
  undo.splits = undo_splits,
  undo.prune = undo_prune,
  undo.SD = undo_sd,
  verbose = verbosity,
)$output

# Export to json
write_json(
  segments,
  path = out.file,
  dataframe = "rows",
  pretty = TRUE,
  auto_unbox = TRUE
)
