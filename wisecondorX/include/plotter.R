options(warn=-1)

# -----
# arg
# -----

args <- commandArgs(TRUE)
in.file <- paste0(args[which(args == "--infile")+1])

# -----
# lib
# -----

#suppressMessages(library("png"))
suppressMessages(library("jsonlite"))

# -----
# main
# -----

input <- read_json(in.file, na = "string")
binsize <- as.integer(input$binsize)
out.dir <- input$out_dir

dir.create(out.dir, showWarnings = FALSE)

# param

gender = input$ref_gender
beta = as.numeric(input$beta)

# aberration_cutoff

get.aberration.cutoff <- function(beta, ploidy){
    loss.cutoff = log2((ploidy - (beta / 2)) / ploidy)
    gain.cutoff = log2((ploidy + (beta / 2)) / ploidy)
    return(c(loss.cutoff, gain.cutoff))
}

# get n.reads (readable)

n.reads <- input$n_reads
first.part <- substr(n.reads, 1, nchar(n.reads) %% 3)
second.part <- substr(n.reads, nchar(n.reads) %% 3 + 1, nchar(n.reads))
n.reads <- c(first.part,  regmatches(second.part, gregexpr(".{3}", second.part))[[1]])
n.reads <- n.reads[n.reads != ""]
n.reads <- paste0(n.reads, collapse = ".")

# get ratios

ratio <- unlist(input$results_r)
ratio[which(ratio == 0)] <- NA
weights <- unlist(input$results_w)
weights[which(weights == 0)] <- NA

if (gender == "M"){
  chrs = 1:24
} else {
  chrs = 1:23
}
bins.per.chr <- sapply(chrs, FUN = function(x) length(unlist(input$results_r[x])))

labels = as.vector(sapply(chrs, FUN = function(x) paste0("chr", x)))
labels = replace(labels, labels == "chr23", "chrX")
labels = replace(labels, labels == "chr24", "chrY")

# find chromosome positions

chr.end.pos <- c(0)
for (chr in chrs){
  l = bins.per.chr[chr]
  chr.end.pos <- c(chr.end.pos, l + chr.end.pos[length(chr.end.pos)])
}

mid.chr <- c()
for (i in 1:(length(chr.end.pos)-1)){
  mid.chr <- c(mid.chr, mean(c(chr.end.pos[i], chr.end.pos[i+1])))
}

ratio <- ratio[1:chr.end.pos[length(chrs) + 1]]
weights <- weights[1:chr.end.pos[length(chrs) + 1]]

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
darker.grey = "#545454"

color.A = darker.grey
color.B = rgb(227, 200, 138, maxColorValue = 255)
color.C = rgb(141, 209, 198, maxColorValue = 255)

png(paste0(out.dir, "/genome_wide.png"), width=14,height=10,units="in",res=512)

l.matrix <- matrix(rep(1, 100), 10, 25, byrow = TRUE)
for (i in 1:7){
  l.matrix <- rbind(l.matrix, c(rep(2, 22),rep(3, 3)))
}

layout(l.matrix)

par(mar = c(4,4,4,0), mgp=c(2.2,-0.5,2))

plot(1, main = "", axes=F, # plots nothing -- enables segments function
     xlab="", ylab="", col = "white", xlim = c(chr.end.pos[1], chr.end.pos[length(chr.end.pos)]),
     cex = 0.0001, ylim=c(chr.wide.lower.limit,chr.wide.upper.limit))

plot.constitutionals <- function(ploidy, start, end){
  segments(start, log2(1/ploidy), end, log2(1/ploidy), col=color.B, lwd = 2, lty = 3)
  segments(start, log2(2/ploidy), end, log2(2/ploidy), col=color.A, lwd = 2, lty = 3)
  segments(start, log2(3/ploidy), end, log2(3/ploidy), col=color.C, lwd = 2, lty = 3)
}

genome.len <- chr.end.pos[length(chr.end.pos)]
autosome.len <- chr.end.pos[23]
if (gender == "F"){
  plot.constitutionals(2, -genome.len * 0.025, genome.len * 1.025)
} else {
  plot.constitutionals(2, -genome.len * 0.025, autosome.len)
  plot.constitutionals(1, autosome.len, genome.len * 1.025)
}

for (undetectable.index in which(is.na(ratio))){
  segments(undetectable.index, chr.wide.lower.limit, undetectable.index, chr.wide.upper.limit, col=lighter.grey, lwd = 0.1, lty = 1)
}

par(new = T)

dot.cex = (weights / pi)**0.5
dot.cols = rep(darker.grey, length(ratio))
for (ab in input$results_c){
  info = unlist(ab)
  chr = as.integer(info[1]) + 1
  start = as.integer(info[2]) + chr.end.pos[chr] + 1
  end = as.integer(info[3]) + chr.end.pos[chr]
  height = as.double(info[5])
  ploidy = 2
  if ((chr == 23 | chr == 24) & gender == "M"){
    ploidy = 1
  }

  if (height < get.aberration.cutoff(beta, ploidy)[1]){
    dot.cols[start:end] = color.B
  }
  if (height > get.aberration.cutoff(beta, ploidy)[2]){
    dot.cols[start:end] = color.C
  }
}
plot(ratio, main = "", axes=F,
     xlab="", ylab=expression('log'[2]*'(ratio)'), col = dot.cols, pch = 16,
     ylim=c(chr.wide.lower.limit,chr.wide.upper.limit), cex = dot.cex)

axis(1, at=mid.chr, labels=labels, tick = F, cex.lab = 3)
axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)

for (x in chr.end.pos){
  segments(x, chr.wide.lower.limit * 1.03, x, chr.wide.upper.limit * 1.03, col=black, lwd = 1.2, lty = 3)
}

par(xpd=TRUE)
# Legends
legend(x=chr.end.pos[length(chr.end.pos)] * 0.2, 
       y = chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.15, 
       legend = c("Constitutional triploid", "Constitutional diploid", "Constitutional monoploid"),
       text.col = c(color.C, color.A, color.B), cex = 1.3, bty="n", text.font = 1.8, lty = c(3,3,3), lwd = 1.5,
       col = c(color.C, color.A, color.B))

legend(x=0,
       y = chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.15,
       legend = c("Gain", "Loss", paste0("Number of reads: ", n.reads)), text.col = c(color.C, color.B, black),
       cex = 1.3, bty="n", text.font = 1.8, pch = c(16,16), col = c(color.C, color.B, "white"))
par(xpd=FALSE)

# plot segmentation

for (ab in input$results_c){
  info = unlist(ab)
  chr = as.integer(info[1]) + 1
  start = as.integer(info[2]) + chr.end.pos[chr] + 1
  end = as.integer(info[3]) + chr.end.pos[chr]
  height = as.double(info[5])
  segments(start, height, end, height, col=lighter.grey, lwd = 5 * mean(dot.cex[start:end], na.rm = T), lty = 1)
}

box("figure", lwd = 1)

# boxplots

par(mar = c(4,4,4,0), mgp=c(2.2,-0.5,2))

boxplot(box.list[1:22], ylim=c(min(l.whis.per.chr[1:22], na.rm = T),
                               max(h.whis.per.chr[1:22], na.rm = T)), bg=black, 
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

if (gender == "F"){
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

  if (gender == "F"){
    plot.constitutionals(2, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
  } else {
    if (c == 23 | c == 24){
      plot.constitutionals(1, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
    } else {
      plot.constitutionals(2, chr.end.pos[c] - bins.per.chr[c] * 0.02, chr.end.pos[c+1] + bins.per.chr[c] * 0.02)
    }
  }

  for (undetectable.index in which(is.na(ratio))){
    segments(undetectable.index, lower.limit, undetectable.index, upper.limit,
             col=darker.grey, lwd = 1/len * 200, lty = 1)
  }
  par(new = T)
  
  plot(ratio, main = labels[c], axes=F,
       xlab="", ylab=expression('log'[2]*'(ratio)'), col = dot.cols, pch = 16,
       cex = dot.cex, ylim=c(lower.limit,upper.limit),
       xlim = margins)

  for (ab in input$results_c){
    info = unlist(ab)
    chr = as.integer(info[1]) + 1
    start = as.integer(info[2]) + chr.end.pos[chr] + 1
    end = as.integer(info[3]) + chr.end.pos[chr]
    height = as.double(info[5])
    segments(start, height, end, height, col=lighter.grey, lwd = 6 * mean(dot.cex[start:end], na.rm = T), lty = 1)
  }

  rect(0, lower.limit - 10, chr.end.pos[c], upper.limit + 10, col="white", border=NA)
  rect(chr.end.pos[c+1], lower.limit - 10, chr.end.pos[length(chr.end.pos)], upper.limit + 10, col="white", border=NA)
  
  axis(1, at=x.labels.at, labels=x.labels, tick = F, cex.axis=0.8)
  axis(2, tick = T, cex.lab = 2, col = black, las = 1, tcl=0.5)
  
  for (x in chr.end.pos){
    segments(x, lower.limit * 1.03, x, upper.limit * 1.03, col=black, lwd = 2, lty = 3)
  }
  for (x in x.labels.at){
    segments(x, lower.limit * 1.02, x, upper.limit * 1.02, col=black, lwd = 1, lty = 3)
  }
  invisible(dev.off())
}

q(save="no")