#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3){
    write("USAGE: mt_consensus_haps.R inbase bcmap outbase\n", stderr())
    write("inbase should be base name of a demux_mt run (will load \n", stderr())
    write("    .cellhaps and .vars)\n", stderr())
    write("bcmap should be a tab separated file of barcode<tab>cluster \n", stderr())
    write("    name\n", stderr())
    write("outbase is output prefix (will create .haps, .vars, and .ids)\n", stderr())
    q()
}

library(reshape2)

inbase <- args[1]
bcmap <- args[2]
outbase <- args[3]

tab <- read.table(paste(inbase, '.cellhaps', sep=''), header=T)
bcs <- read.table(bcmap, sep='\t')
colnames(bcs) <- c("bc", "clust")
tab <- merge(tab, bcs, by='bc')

melted <- melt(tab, id.vars=c("bc", "clust"))
melted$count <- 1

cells <- aggregate(melted$count, by=list(clust=melted$clust, site=melted$variable), FUN=sum)
melted <- melted[which(!is.na(melted$value)),]
sums <- aggregate(melted$value, by=list(clust=melted$clust, site=melted$variable), FUN=sum)
tots <- aggregate(melted$count, by=list(clust=melted$clust, site=melted$variable), FUN=sum)

colnames(cells)[3] <- "cells"
colnames(sums)[3] <- "sum"
colnames(tots)[3] <- "tot"

merged <- merge(sums, tots)
merged <- merge(merged, cells)

# Drop sites missing too much from one or more haps
dscutoff <- 0.5
dropsites <- unique(merged[which(merged$tot / merged$cells < dscutoff),]$site)

merged$frac <- merged$sum / merged$tot
merged$val <- 0

ll1 <- dbinom(merged$sum, merged$tot, 0.99, log=TRUE)
ll2 <- dbinom(merged$sum, merged$tot, 0.01, log=TRUE)
ll3 <- dbinom(merged$tot, merged$cells, 0.01, log=TRUE)

merged[which(ll1 > ll2 & ll1 > ll3),]$val <- 1
if (sum(ll3 > ll1 & ll3 > ll2) > 0){
    merged[which(ll3 > ll1 & ll3 > ll2),]$val <- NA
}
#merged[which(merged$frac > 0.9),]$val <- 1

casted <- dcast(merged, clust ~ site, value.var='val')
clusts <- data.frame(clust=casted$clust)
casted <- casted[,-c(1)]
vars <- read.table(paste(inbase, '.vars', sep=''))
vars$idx <- seq(0, length(rownames(vars))-1)

cols <- data.frame(col=colnames(casted))

cols$clean <- as.numeric(gsub("X", "", cols$col))
# If more than one cluster, only keep sites that are variable across clusters.
# If only one cluster, keep all sites.
if (length(unique(bcs$clust)) > 1){
    keeps <- cols[which(cols$col %in% colnames(casted)[which(colSums(casted, na.rm=TRUE) > 0)]),]
    casted <- casted[,which(colSums(casted, na.rm=TRUE) > 0)]
} else{
    keeps <- cols
}

# To do:
# Check for and remove identical haplotypes, ignoring NAs?

# Convert NA to - and 0/1 to string
for (col in colnames(casted)){
    casted[[col]] <- as.character(casted[[col]])
    if (sum(is.na(casted[[col]])) > 0){
        casted[which(is.na(casted[[col]])),][[col]] <- "-"
    }
}

vars <- vars[which(vars$idx %in% keeps$clean),]
vars <- vars[,which(colnames(vars) != "idx")]

write.table(vars, file=paste(outbase, ".vars", sep=""), sep="\t", 
    quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(clusts, file=paste(outbase, ".ids", sep=""), sep="\t", 
    quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(casted, file=paste(outbase, ".haps", sep=""), sep="", 
    quote=FALSE, row.names=FALSE, col.names=FALSE)

colnames(clusts)[1] <- "bc"
casted <- cbind(clusts, casted)

write.table(casted, file=paste(outbase, ".cellhaps", sep=""), 
    sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


