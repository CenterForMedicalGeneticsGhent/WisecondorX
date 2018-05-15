options(warn=-1)

# -----
# arg
# -----

arguments = "

--infile          - (mandatory argument) input .json
--sexchroms       - (mandatory argument) which sex chromosomes will be analyzed? (X or XY)
--beta            - (mandatory argument) which sex chromosomes will be analyzed? (X or XY)
--outdir          - (mandatory argument) output folder

"

args <- commandArgs(TRUE)
len.args <- length(args)

parseArgs <- function(x) strsplit(sub("^--", "", x), " ")
args <- as.matrix(do.call("rbind", parseArgs(args)))

in.file <- paste0(args[,1][which(args[,1] == "infile")+1])
sex.chrom <- paste0(args[,1][which(args[,1] == "sexchroms")+1])
beta <- paste0(args[,1][which(args[,1] == "beta")+1])
out.dir <- paste0(args[,1][which(args[,1] == "outdir")+1])

# -----
# param
# -----

if (sex.chrom == "X"){
  exclude.chr = c(24)
} else {
  exclude.chr = c()
}

gain.cut <- log2(1 + as.numeric(beta) / 4)
del.cut <- log2(1 - as.numeric(beta) / 4)

# -----
# lib
# -----

suppressMessages(library("png"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

dir.create(out.dir, showWarnings = FALSE)

input <- read_json(in.file, na = "string")
binsize <- input$binsize

# get nreads (readable)

nreads <- input$nreads
first.part <- substr(nreads, 1, nchar(nreads) %% 3)
second.part <- substr(nreads, nchar(nreads) %% 3 + 1, nchar(nreads))
nreads <- c(first.part,  regmatches(second.part, gregexpr(".{3}", second.part))[[1]])
nreads <- nreads[nreads != ""]
nreads <- paste0(nreads, collapse = ".")

# get ratios

ratio <- unlist(input$results_r)
ratio[which(ratio == 0)] <- NA #0 were introduced if infinite/na values were seen in python. They're all the same: of no interest

chrs = c(1:24)
chrs = chrs[which(!(chrs  %in% exclude.chr))]
bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$results_r[x])))

labels = as.vector(sapply(chrs, FUN = function(x) paste0("chr", x)))
labels = replace(labels, labels == "chr23", "chrX")
labels = replace(labels, labels == "chr24", "chrY")

# find chromosome positions

chr.end.pos <- c(1)
for (chr in chrs){
  l = bins.per.chr[chr]
  chr.end.pos <- c(chr.end.pos, l + chr.end.pos[length(chr.end.pos)])
}

mid.chr <- c()
for (i in 1:(length(chr.end.pos)-1)){
  mid.chr <- c(mid.chr, mean(c(chr.end.pos[i], chr.end.pos[i+1])))
}

ratio <- ratio[1:chr.end.pos[length(chrs) + 1]]

# get margins

box.list <- list()
l.whis.per.chr <- c()
h.whis.per.chr <- c()

for (chr in chrs){
  box.list[[chr]] <- ratio[chr.end.pos[chr]:chr.end.pos[chr + 1]]
  whis = boxplot(box.list[[chr]], plot = F)$stats[c(1,5),]
  l.whis.per.chr = c(l.whis.per.chr, whis[1])
  h.whis.per.chr = c(h.whis.per.chr, whis[2])
}

chr.wide.upper.limit <- max(0.65, max(h.whis.per.chr), na.rm = T) * 1.25
chr.wide.lower.limit <- min(-0.95, min(l.whis.per.chr), na.rm = T) * 1.25

# plot chromosome wide plot

black = "#3f3f3f"
lighter.grey = "#e0e0e0"
darker.grey = "#9e9e9e"
green = "#00a024"
orange = "#cc8c0c"
red = "#d10606"
blue = "#3f60a0"

png(paste0(out.dir, "/chromosomeWide.png"), width=14,height=10,units="in",res=512)

l.matrix <- matrix(rep(1, 100), 10, 25, byrow = TRUE)
for (i in 1:7){
  l.matrix <- rbind(l.matrix, c(rep(2, 22),rep(3, 3)))
}

layout(l.matrix)

par(mar = c(4,4,4,0), mgp=c(2.2,-0.5,2))

plot(1, main = "", axes=F, # plots nothing -- enables segments function
     xlab="", ylab="", col = "white", xlim = c(chr.end.pos[1], chr.end.pos[length(chr.end.pos)]),
     cex = 0.0001, ylim=c(chr.wide.lower.limit,chr.wide.upper.limit))

for (undetectable.index in which(is.na(ratio))){
  segments(undetectable.index, chr.wide.lower.limit, undetectable.index, chr.wide.upper.limit, col=lighter.grey, lwd = 0.1, lty = 1)
}
par(new = T)

dot.cols = rep(black, length(ratio))
for (ab in input$cbs_calls){
  info = unlist(ab)
  chr = as.integer(info[1])
  if (chr %in% exclude.chr){
    next
  }
  start = as.integer(info[2]) + chr.end.pos[chr]
  end = as.integer(info[3]) + chr.end.pos[chr]
  height = as.double(info[5])
  if (height > gain.cut){
    dot.cols[start:end] = blue
  }
  if (height < del.cut){
    dot.cols[start:end] = red
  }
}
plot(ratio, main = "", axes=F,
     xlab="", ylab=expression('log'[2]*'(ratio)'), col = dot.cols, pch = 21, 
     cex = 0.4, ylim=c(chr.wide.lower.limit,chr.wide.upper.limit), bg = dot.cols)

axis(1, at=mid.chr, labels=labels, tick = F, cex.lab = 3)
axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)

plot.constitutionals <- function(ploidy, start, end){
  segments(start, log2(1/ploidy), end, log2(1/ploidy), col=red, lwd = 2, lty = 3)
  segments(start, log2(2/ploidy), end, log2(2/ploidy), col=orange, lwd = 2, lty = 3)
  segments(start, log2(3/ploidy), end, log2(3/ploidy), col=blue, lwd = 2, lty = 3)
}

genome.len <- chr.end.pos[length(chr.end.pos)]
autosome.len <- chr.end.pos[23]
if ((sex.chrom == "X")){
  plot.constitutionals(2, -genome.len * 0.025, genome.len * 1.025)
} else {
  plot.constitutionals(2, -genome.len * 0.025, autosome.len)
  plot.constitutionals(1, autosome.len, genome.len * 1.025)
}

for (x in chr.end.pos){
  segments(x, chr.wide.lower.limit * 1.03, x, chr.wide.upper.limit * 1.03, col=green, lwd = 1.2, lty = 3)
}

par(xpd=TRUE)
# Legends
legend(x=chr.end.pos[length(chr.end.pos)] * 0.2, 
       y = chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.15, 
       legend = c("Constitutional triploid", "Constitutional diploid", "Constitutional monoploid"),
       text.col = c(blue, orange, red), cex = 1.3, bty="n", text.font = 1.8, lty = c(3,3,3), lwd = 1.5,
       col = c(blue, orange, red))

legend(x=0,
       y = chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.15,
       legend = c("Gain", "Loss", paste0("Number of reads: ", nreads)), text.col = c(blue, red, black),
       cex = 1.3, bty="n", text.font = 1.8, pch = c(16,16), col = c(blue, red, "white"))
par(xpd=FALSE)

# plot segmentation

for (ab in input$cbs_calls){
  info = unlist(ab)
  chr = as.integer(info[1])
  if (chr %in% exclude.chr){
    next
  }
  start = as.integer(info[2]) + chr.end.pos[chr]
  end = as.integer(info[3]) + chr.end.pos[chr]
  bm.score = abs(as.double(info[4]))
  height = as.double(info[5])
  segments(start, height, end, height, col=lighter.grey, lwd = 2, lty = 1)
}

box("figure", lwd = 1)

# boxplots

par(mar = c(4,4,4,0), mgp=c(2.2,-0.5,2))

boxplot(box.list[1:22], ylim=c(min(l.whis.per.chr[1:22]), max(h.whis.per.chr[1:22])), bg=black, 
        axes=F, outpch = 16, ylab = expression('log'[2]*'(ratio)'))
axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)
par(mar = c(4,4,4,0), mgp=c(1,0.5,2))
axis(1, at=1:22, labels=labels[1:22], tick = F, cex.lab = 3)

plot.constitutionals(2, 0, 23)

par(mar = c(4,4,4,0), mgp=c(2.2,-0.5,2))

y.sex.down = min(l.whis.per.chr[23:length(chrs)], na.rm = T)
y.sex.up = max(h.whis.per.chr[23:length(chrs)], na.rm = T)

if(any(is.infinite(c(y.sex.down, y.sex.up)))){
  y.sex.down = 0
  y.sex.up = 0
}

boxplot(box.list[23:length(chrs)], ylim=c(y.sex.down, y.sex.up), 
        bg=black, axes=F, outpch = 16, ylab = expression('log'[2]*'(ratio)'))
axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)
par(mar = c(4,4,4,0), mgp=c(1,0.5,2))
axis(1, at=1:(length(chrs) - 22), labels=labels[23:length(chrs)], tick = F, cex.lab = 3)

if ((sex.chrom == "X")){
  plot.constitutionals(2, 0.6, length(chrs[23:length(chrs)]) + 1)
} else {
  plot.constitutionals(1, 0.6, length(chrs[23:length(chrs)]) + 1)
}


box("outer", lwd = 4)

# write image

invisible(dev.off())

# create chr specific plots

for (c in chrs){
  
  margins <- c(chr.end.pos[c], chr.end.pos[c+1])
  len <- chr.end.pos[c+1] - chr.end.pos[c]
  x.labels <- seq(0, bins.per.chr[c] * binsize, bins.per.chr[c] * binsize / 10)
  x.labels.at <- seq(0, bins.per.chr[c], bins.per.chr[c] / 10) + chr.end.pos[c]
  x.labels <- x.labels[2:(length(x.labels) - 1)]
  x.labels.at <- x.labels.at[2:(length(x.labels.at) - 1)]
  
  mean = mean(box.list[[c]], na.rm = T)
  whis = boxplot(box.list[[c]], plot = F)$stats[c(1,5),]
  
  if (any(is.na(whis))){
    next
  }
  
  png(paste0(out.dir, "/", labels[c],".png"), width=14,height=10,units="in",res=256)

  upper.limit <- 0.6 + whis[2]
  lower.limit <- -1.05 + whis[1]
  
  par(mar = c(4,4,4,0), mgp=c(2.2,-0.2,2))
  
  plot(1, main = "", axes=F, # plots nothing -- enables segments function
       xlab="", ylab="", col = "white", 
       cex = 0.0001, ylim=c(lower.limit,upper.limit), xlim = margins)
  
  for (undetectable.index in which(is.na(ratio))){
    segments(undetectable.index, lower.limit, undetectable.index, upper.limit,
             col=darker.grey, lwd = 1/len * 200, lty = 1)
  }
  par(new = T)
  
  plot(ratio, main = labels[c], axes=F,
       xlab="", ylab=expression('log'[2]*'(ratio)'), col = dot.cols, pch = 21, 
       cex = 0.4, ylim=c(lower.limit,upper.limit),
       xlim = margins, bg=dot.cols)

  for (ab in input$cbs_calls){
    info = unlist(ab)
    chr = as.integer(info[1])
    if (chr %in% exclude.chr){
      next
    }
    start = as.integer(info[2]) + chr.end.pos[chr]
    end = as.integer(info[3]) + chr.end.pos[chr]
    bm.score = abs(as.double(info[4]))
    height = as.double(info[5])
    segments(start, height, end, height, col=lighter.grey, lwd = 2, lty = 1)
  }

  rect(0, lower.limit - 10, chr.end.pos[c], upper.limit + 10, col="white", border=NA)
  rect(chr.end.pos[c+1], lower.limit - 10, chr.end.pos[length(chr.end.pos)], upper.limit + 10, col="white", border=NA)
  
  axis(1, at=x.labels.at, labels=x.labels, tick = F, cex.axis=0.8)
  axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)
  
  if ((sex.chrom == "X")){
    plot.constitutionals(2, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
  } else {
    if (c == 23 | c == 24){
      plot.constitutionals(1, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
    } else {
      plot.constitutionals(2, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
    }
  }
  
  for (x in chr.end.pos){
    segments(x, lower.limit * 1.03, x, upper.limit * 1.03, col=green, lwd = 2, lty = 3)
  }
  for (x in x.labels.at){
    segments(x, lower.limit * 1.02, x, upper.limit * 1.02, col=green, lwd = 1, lty = 3)
  }
  invisible(dev.off())
}

q(save="no")