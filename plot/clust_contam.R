#! /usr/bin/env Rscript
suppressPackageStartupMessages(library(dendextend))

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    sys.stderr.write("USAGE: clust_contam.R plotname file\n")
    sys.stderr.write("file is a file listing output base names from demux_vcf, one per line\n")
    q()
}

dat <- data.frame()

plotname <- args[1]
plotname <- strsplit(plotname, '.png')[[1]][1]
plotname <- strsplit(plotname, '.pdf')[[1]][1]

outname <- paste(plotname, '.png', sep='')

basenames <- read.table(args[2])

meta <- data.frame(lib=unique(basenames$V1), ncells=0, meancr=0)

idx <- 1
for (basename in unique(basenames$V1)){
    df <- read.table(paste(basename, '.contam.dat', sep=''))
    
    assn <- read.table(paste(basename, '.assignments', sep=''))
    ncells <- length(rownames(assn))
    meta[which(meta$lib==basename),]$ncells <- ncells
    cr <- read.table(paste(basename, '.contam_rate', sep=''))
    meancr <- mean(cr$V2)
    meta[which(meta$lib==basename),]$meancr <- meancr

    colnames(df) <- c("ind1", "type1", "ind2", "type2", basename)
    if (idx == 1){
        dat <- df
    } else{
        dat <- merge(dat, df)
    }
    
    idx <- idx + 1
}

dat <- dat[,-c(1,2,3,4)]
m <- t(as.matrix(dat))

d <- dist(m)
hc <- hclust(d)
dend <- as.dendrogram(hc)

cr <- colorRamp(c('#0267C1', '#EFA00B'))
colFun <- function(x){ return(rgb(cr(x), maxColorValue=255)) }

#cols <- colFun(meta[order.dendrogram(dend),]$ncell / max(meta$ncell))
cols <- colFun(meta[order.dendrogram(dend),]$meancr)

png(outname, width=7, height=7, units='in', res=150)

layout(matrix(1:2,ncol=2), width = c(0.75,0.25),height = c(1,1))
labels_colors(dend) <- cols
par(mai=c(0.5,0.5,0.5,2))
plot(dend, type='rectangle', horiz=TRUE, edgePar=list(lwd=2),)
title("Ambient RNA similarity")

steps <- rev(seq(0,1,0.1))
legend_image <- as.raster(matrix(colFun(steps), ncol=1))
par(mai=c(1,0.5,1,0.25))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Contam rate')
text(x=1.75, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)

dev.off()

