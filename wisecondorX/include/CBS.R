# -----
# arg
# -----

arguments = "

--infile          - (mandatory argument) input .json file
--sexchroms       - (mandatory argument) which sex chromosomes will be analyzed? (X or XY)
--outfile         - (mandatory argument) output .json file
--alpha           - (mandatory argument) p-value for segmentation

"

args <- commandArgs(TRUE)


if (all(c(args[1], args[3], args[5], args[7]) %in% c("--infile", "--sexchroms", "--outfile", "--alpha"))){
  in.file <- paste0(args[which(args == "--infile")+1])
  sex.chrom <- paste0(args[which(args == "--sexchroms")+1])
  out.file <- paste0(args[which(args == "--outfile")+1])
  alpha <- paste0(args[which(args == "--alpha")+1])
}

# -----
# param
# -----

if (sex.chrom == "X"){
    exclude.chr = c(24)
} else {
    exclude.chr = c()
}

# -----
# lib
# -----

suppressMessages(library("DNAcopy"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

# process data

input <- read_json(in.file)
ratio <- as.numeric(unlist(input$results_r))

chrs = c(1:24)
chrs = chrs[which(!(chrs  %in% exclude.chr))]
bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$results_r[x])))
chr.end.pos <- c(0)
for (chr in chrs){
  l = bins.per.chr[chr]
  chr.end.pos <- c(chr.end.pos, l + chr.end.pos[length(chr.end.pos)])
}

# prepare for CBS

ratio[ratio == 0] = NA
# Only segments will be further processed with this data.
# Software will not account for 0's (these bins had no reference)

for.cbs <- as.data.frame(ratio)
chr.rep <- c()
chr.rep.2 <- c()
for (chr in chrs){
  chr.rep <- c(chr.rep, rep(chr, chr.end.pos[chr + 1] - chr.end.pos[chr]))
  chr.rep.2 <- c(chr.rep.2, 1:(chr.end.pos[chr + 1] - chr.end.pos[chr]))
}
for.cbs$chromosome <- chr.rep; for.cbs$x <- chr.rep.2
for.cbs <- for.cbs[, c(2,3,1)] ; colnames(for.cbs)[3] <- "y"

# Check for complete NA/chr

remove.this <- c()
for (chr in chrs){
    check.this <- which(for.cbs$chromosome == chr)
    if(all(is.na(for.cbs$y[check.this]))){
        remove.this <- c(remove.this, check.this)
    }
}
if (length(remove.this != 0)){
    for.cbs <- for.cbs[-remove.this,]
}

# CBS

CNA.object <- CNA(for.cbs$y, for.cbs$chromosome, for.cbs$x, data.type = "logratio", sampleid = "X")
f = file()
sink(file=f) ## silence output
CNA.object <- invisible(segment(CNA.object, alpha = as.numeric(alpha), verbose=1)$output)
sink() ## undo silencing
close(f)

# Write output

write_json(t(CNA.object), out.file)