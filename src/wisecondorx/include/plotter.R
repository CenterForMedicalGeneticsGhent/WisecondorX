options(warn=-1)

# -----
# arg
# -----

args <- commandArgs(T)
in.file <- paste0(args[which(args == "--infile")+1])


# -----
# lib
# -----

suppressMessages(library("jsonlite"))

# -----
# main
# -----

input <- read_json(in.file, na="string")
binsize <- as.integer(input$binsize)
out.dir <- input$out_dir

dir.create(out.dir, showWarnings=F)

# param

gender = input$ref_gender
beta = as.numeric(input$beta)
zcutoff = as.numeric(input$zscore)
ylim = input$ylim
plot.title = input$plot_title

if (input$cairo) options(bitmaptype='cairo')

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
n.reads <- paste0(n.reads, collapse=".")

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
bins.per.chr <- sapply(chrs, FUN=function(x) length(unlist(input$results_r[x])))

labels = as.vector(sapply(chrs, FUN=function(x) paste0("chr", x)))
labels = replace(labels, labels == "chr23", "chrX")
labels = replace(labels, labels == "chr24", "chrY")

# define chromosome positions

chr.ends <- c(1, cumsum(bins.per.chr))
chr.mids <- chr.ends[2:length(chr.ends)] - bins.per.chr/2

ratio <- ratio[1:chr.ends[length(chrs) + 1]]
weights <- weights[1:chr.ends[length(chrs) + 1]]

# get margins

box.list <- list()
l.whis.per.chr <- c()
h.whis.per.chr <- c()

for (chr in chrs){
  box.list[[chr]] <- ratio[chr.ends[chr]:chr.ends[chr + 1]]
  whis = boxplot(box.list[[chr]], plot=F)$stats[c(1,5),]
  l.whis.per.chr = c(l.whis.per.chr, whis[1])
  h.whis.per.chr = c(h.whis.per.chr, whis[2])
}
chr.wide.upper.limit <- max(0.65, max(h.whis.per.chr), na.rm=T) * 1.25
chr.wide.lower.limit <- min(-0.95, min(l.whis.per.chr), na.rm=T) * 1.25

if (ylim != 'def'){
  ylim=gsub('[', '', ylim, fixed=T) ; ylim=gsub(']', '', ylim, fixed=T)
  ylim=as.numeric(strsplit(ylim, ',', fixed=T)[[1]])
  chr.wide.lower.limit = ylim[1]
  chr.wide.upper.limit = ylim[2]
}

# plot chromosome wide plot

black = "#3f3f3f"
lighter.grey = "#e0e0e0"

color.A = rgb(84, 84, 84, maxColorValue=255)
color.B = rgb(227, 200, 138, maxColorValue=255)
color.C = rgb(141, 209, 198, maxColorValue=255)
color.X <- c(color.C, color.A, color.B)

color.AA = rgb(84, 84, 84, 80, maxColorValue=255)
color.BB = rgb(227, 200, 138, 80, maxColorValue=255)
color.CC = rgb(141, 209, 198, 80, maxColorValue=255)
color.XX = c(color.CC, color.AA, color.BB)

png(paste0(out.dir, "/genome_wide.png"), width=14,height=10,units="in",res=512,pointsize=18)

l.matrix <- matrix(rep(1, 100), 10, 25, byrow=T)
for (i in 1:7){
  l.matrix <- rbind(l.matrix, c(rep(2, 22),rep(3, 3)))
}

layout(l.matrix)

par(mar=c(2.5,4,1,0), mgp=c(2.4,0.5,0), oma=c(0,0,3,0))

plot(1, main="", axes=F, # plots nothing -- enables segments function
     xlab="", ylab="", xlim=c(chr.ends[1], chr.ends[length(chr.ends)]),
     cex=0, ylim=c(chr.wide.lower.limit,chr.wide.upper.limit))

plot.constitutionals <- function(ploidy, start, end){
  segments(start, log2(1/ploidy), end, log2(1/ploidy), col=color.B, lwd=2, lty=3)
  segments(start, log2(2/ploidy), end, log2(2/ploidy), col=color.A, lwd=2, lty=3)
  segments(start, log2(3/ploidy), end, log2(3/ploidy), col=color.C, lwd=2, lty=3)
}

genome.len <- chr.ends[length(chr.ends)]
autosome.len <- chr.ends[23]
if (gender == "F"){
  plot.constitutionals(2, -genome.len * 0.025, genome.len * 1.025)
} else {
  plot.constitutionals(2, -genome.len * 0.025, autosome.len)
  plot.constitutionals(1, autosome.len, genome.len * 1.025)
}

for (undetectable.index in which(is.na(ratio))){
  segments(undetectable.index, chr.wide.lower.limit, undetectable.index, chr.wide.upper.limit, col=lighter.grey, lwd=0.1, lty=1)
}

dot.cex = (weights / pi)**0.5 * .8
dot.cols = rep(color.A, length(ratio))
for (ab in input$results_c){
  info = unlist(ab)
  chr = as.integer(info[1]) + 1
  start = as.integer(info[2]) + chr.ends[chr] + 1
  end = as.integer(info[3]) + chr.ends[chr]
  z = as.double(info[4])
  height = as.double(info[5])
  ploidy = 2
  if ((chr == 23 | chr == 24) & gender == "M"){
    ploidy = 1
  }

  if (!is.na(beta)){
    if (height < get.aberration.cutoff(beta, ploidy)[1]){
    dot.cols[start:end] = color.B
    }
    if (height > get.aberration.cutoff(beta, ploidy)[2]){
      dot.cols[start:end] = color.C
    }
  } else {
    if (is.na(z)){
      dot.cols[start:end] = 'grey'
      next
    }
    if (z < -zcutoff) dot.cols[start:end] = color.B
    if (z > zcutoff) dot.cols[start:end] = color.C
  }
}

par(new=T)
plot(ratio, main="", axes=F,
     xlab="", ylab=expression('log'[2]*'(ratio)'), col=dot.cols, pch=16,
     ylim=c(chr.wide.lower.limit,chr.wide.upper.limit), cex=dot.cex)

par(xpd=NA)
text(chr.mids, par("usr")[3], labels=labels, srt=45, pos=1)
axis(2, tick=T, cex.lab=2, col=black, las=1, tcl=0.5)
par(xpd=F)

for (x in chr.ends){
  segments(x, chr.wide.lower.limit * 1.03, x, chr.wide.upper.limit * 1.03, col=black, lwd=1.2, lty=3)
}

par(xpd=NA)
segments(sum(bins.per.chr) * -.08, par("usr")[3] - abs(chr.wide.upper.limit - chr.wide.lower.limit) * .1,
         sum(bins.per.chr) * 1.03, par("usr")[3] - abs(chr.wide.upper.limit - chr.wide.lower.limit) * .1, lwd=3)
par(xpd=F)

# Legends

par(xpd=NA)
legend(x=chr.ends[length(chr.ends)] * 0.3,
       y=chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.23,
       legend=c("Constitutional 3n", "Constitutional 2n", "Constitutional 1n"),
       text.col=color.X, cex=1.3, bty="n", lty=c(3,3,3), lwd=1.5, col=color.X)
legend(x=0,
       y=chr.wide.upper.limit + (abs(chr.wide.upper.limit) + abs(chr.wide.lower.limit)) * 0.23,
       legend=c("Gain", "Loss", paste0("Number of reads: ", n.reads)), text.col=c(color.C, color.B, black),
       cex=1.3, bty="n", pch=c(16,16), col=c(color.C, color.B, "white"))
par(xpd=F)

if (!is.null(plot.title)) {
    par(xpd=NA)
    title(plot.title, line=1.3, adj=0.8, col.main=color.A)
    par(xpd=F)
}

# plot segmentation

for (ab in input$results_c){
  info = unlist(ab)
  chr = as.integer(info[1]) + 1
  start = as.integer(info[2]) + chr.ends[chr] + 1
  end = as.integer(info[3]) + chr.ends[chr]
  height = as.double(info[5])
  rect(start, height, end, 0, col=color.XX[dot.cols[start] == color.X], border=color.XX[dot.cols[start] == color.X], lwd=0.1)
  segments(start, height, end, height, col=lighter.grey, lwd=5 * mean(dot.cex[start:end], na.rm=T), lty=1)
}

# boxplots

par(mgp=c(2.4,0.5,0))
boxplot(box.list[1:22], ylim=c(min(l.whis.per.chr[1:22], na.rm=T),
                               max(h.whis.per.chr[1:22], na.rm=T)), bg=black,
        axes=F, outpch=16, ylab=expression('log'[2]*'(ratio)'))

par(xpd=NA)
text(1:22, par("usr")[3], labels=labels[1:22], srt=45, pos=1)
axis(2, tick=T, cex.lab=2, col=black, las=1, tcl=0.5)
par(xpd=F)

plot.constitutionals(2, 0, 23)

y.sex.down = min(l.whis.per.chr[23:length(chrs)], na.rm=T)
y.sex.up = max(h.whis.per.chr[23:length(chrs)], na.rm=T)

if(any(is.infinite(c(y.sex.down, y.sex.up)))){
  y.sex.down = 0
  y.sex.up = 0
}

par(mar=c(2.5,3,1,1))
boxplot(box.list[23:length(chrs)], ylim=c(y.sex.down, y.sex.up),
        bg=black, axes=F, outpch=16, ylab='')

par(xpd=NA)
text(1:(length(chrs) - 22), par("usr")[3], labels=labels[23:length(chrs)], srt=45, pos=1)
axis(2, tick=T, cex.lab=2, col=black, las=1, tcl=0.5)
par(xpd=F)

if (gender == "F"){
  plot.constitutionals(2, 0.6, 1.5)
} else {
  plot.constitutionals(1, 0.5, 2.6)
}

# write image

invisible(dev.off())

# create chr specific plots

for (c in chrs){

  margins <- c(chr.ends[c], chr.ends[c+1])
  len <- chr.ends[c+1] - chr.ends[c]
  x.labels <- seq(0, bins.per.chr[c] * binsize, bins.per.chr[c] * binsize / 10)
  x.labels.at <- seq(0, bins.per.chr[c], bins.per.chr[c] / 10) + chr.ends[c]
  x.labels <- x.labels[2:(length(x.labels) - 1)]
  x.labels.at <- x.labels.at[2:(length(x.labels.at) - 1)]

  whis = boxplot(box.list[[c]], plot=F)$stats[c(1,5),]

  if (any(is.na(whis))){
    next
  }

  png(paste0(out.dir, "/", labels[c],".png"), width=14,height=10,units="in",res=256,pointsize=18)

  upper.limit <- 0.6 + whis[2]
  lower.limit <- -1.05 + whis[1]
  upper.limit <- max(upper.limit, max(ratio[margins[1]:margins[2]], na.rm = T))
  lower.limit <- min(lower.limit, min(ratio[margins[1]:margins[2]], na.rm = T))
  if (ylim != 'def'){
    lower.limit = ylim[1] ; upper.limit = ylim[2]
  }
  par(mar=c(2.5,4,1,0), mgp=c(2.4,0.5,0))

  plot(1, main="", axes=F, # plots nothing -- enables segments function
       xlab="", ylab="", cex=0, ylim=c(lower.limit,upper.limit), xlim=margins)

  if (gender == "F"){
    plot.constitutionals(2, chr.ends[c] - bins.per.chr[c] * 0.02, chr.ends[c+1] + bins.per.chr[c] * 0.02)
  } else {
    if (c == 23 | c == 24){
      plot.constitutionals(1, chr.ends[c] - bins.per.chr[c] * 0.02, chr.ends[c+1] + bins.per.chr[c] * 0.02)
    } else {
      plot.constitutionals(2, chr.ends[c] - bins.per.chr[c] * 0.02, chr.ends[c+1] + bins.per.chr[c] * 0.02)
    }
  }

  for (undetectable.index in which(is.na(ratio))){
    segments(undetectable.index, lower.limit, undetectable.index, upper.limit,
             col=color.A, lwd=1/len * 200, lty=1)
  }

  par(new=T)
  plot(ratio, main=labels[c], axes=F,
       xlab="", ylab=expression('log'[2]*'(ratio)'), col=dot.cols, pch=16,
       cex=dot.cex, ylim=c(lower.limit,upper.limit),
       xlim=margins)

  for (ab in input$results_c){
    info = unlist(ab)
    chr = as.integer(info[1]) + 1
    start = as.integer(info[2]) + chr.ends[chr] + 1
    end = as.integer(info[3]) + chr.ends[chr]
    height = as.double(info[5])
    rect(start, height, end, 0, col=color.XX[dot.cols[start] == color.X],border=color.XX[dot.cols[start] == color.X], lwd=0.1)
    segments(start, height, end, height, col=lighter.grey, lwd=6 * mean(dot.cex[start:end], na.rm=T), lty=1)
  }

  rect(0, lower.limit - 10, chr.ends[c], upper.limit + 10, col="white", border=NA)
  rect(chr.ends[c+1], lower.limit - 10, chr.ends[length(chr.ends)], upper.limit + 10, col="white", border=NA)

  par(xpd=NA)
  text(x.labels.at, par("usr")[3], labels=x.labels, srt=45, pos=1)
  axis(2, tick=T, cex.lab=2, col=black, las=1, tcl=0.5)
  par(xpd=F)

  for (x in chr.ends){
    segments(x, lower.limit * 1.03, x, upper.limit * 1.03, col=black, lwd=2, lty=3)
  }
  for (x in x.labels.at){
    segments(x, lower.limit * 1.02, x, upper.limit * 1.02, col=black, lwd=1, lty=3)
  }
  invisible(dev.off())
}

q(save="no")