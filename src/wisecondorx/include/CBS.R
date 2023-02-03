# -----
# arg
# -----

args <- commandArgs(TRUE)
in.file <- paste0(args[which(args == "--infile")+1])

# -----
# lib
# -----

suppressMessages(library("DNAcopy"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

# load input

input <- read_json(in.file)
ratio <- as.numeric(unlist(input$results_r))
weights <- as.numeric(unlist(input$results_w))
seed <- as.numeric(input$seed)
gender <- input$ref_gender
alpha <- as.numeric(input$alpha)
binsize <- as.numeric(input$binsize)
out.file <- as.character(input$outfile)

if (gender == "M"){
    chrs = 1:24
} else {
    chrs = 1:23
}

# prepare for CBS

bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$results_r[x])))
chr.end.pos <- c(0,cumsum(bins.per.chr))

ratio[ratio == 0] = NA # blacklist
weights[weights == 0] = 1^-99 # omit DNAcopy weirdness -- weight cannot be NA or 0

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

cbs.mask <- c()
for (chr in chrs){
    check <- which(for.cbs$chromosome == chr)
    if(!(all(is.na(for.cbs$y[check])))){
        cbs.mask <- c(cbs.mask, check)
    }
}
for.cbs <- for.cbs[cbs.mask,]

# CBS

if (is.na(seed) || seed == '') {
  set.seed(seed)
}
CNA.object <- CNA(for.cbs$y, for.cbs$chromosome, for.cbs$x, data.type = "logratio", sampleid = "X")
f = file()
sink(file=f) ## silence output
CNA.object <- invisible(segment(CNA.object, alpha = as.numeric(alpha), verbose=1, weights=weights[cbs.mask])$output)
sink() ## undo silencing
close(f)

CNA.object <- CNA.object[,-c(1,5)]
colnames(CNA.object) <- c("chr", "s", "e", "r")

# Check if segment covers large NA regions. If so = split

new.CNA.object <- data.frame()

for (row.i in 1:nrow(CNA.object)){
  start.i = CNA.object$s[row.i]
  end.i = CNA.object$e[row.i]
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chr[row.i], ]
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
  
  sub.frame <- cbind(CNA.object$chr[row.i], inverse.start.pos, inverse.end.pos, CNA.object$r[row.i])
  new.CNA.object <- rbind(new.CNA.object, sub.frame)
  
}

colnames(new.CNA.object) <- c("chr", "s", "e", "r")
CNA.object <- new.CNA.object

# Recalculate segmental ratios

for.cbs$w <- weights[cbs.mask]

for (row.i in 1:nrow(CNA.object)){
  sub.frame = for.cbs[for.cbs$chromosome == CNA.object$chr[row.i], ]
  CNA.object$r[row.i] = weighted.mean(sub.frame$y[CNA.object$s[row.i]:CNA.object$e[row.i]],
                                      sub.frame$w[CNA.object$s[row.i]:CNA.object$e[row.i]],
                                      na.rm = T)
}

CNA.object$s <- CNA.object$s - 1 # Make python compliant

# Write output
write_json(CNA.object, out.file)