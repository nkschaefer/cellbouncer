#! /usr/bin/env Rscript
library(ggplot2)
library(viridis)
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    write("USAGE:", stderr())
    write("compare_assignments.R <out1.assignments> <out2.assignments> <plotname> (S)", stderr())
    write("Given two .assignments files from cellbouncer programs run on", stderr())
    write("    the same cells, compares them to each other and plots a heatmap", stderr())
    write("    that should show whether the assignments tend to agree well,", stderr())
    write("    and if so, which labels from the first file match which labels", stderr())
    write("    from the second file. An example use case would be running demux_vcf", stderr())
    write("    and demux_mt on the same data, and checking which mitochondrial", stderr())
    write("    haplotypes match which individuals in the VCF.", stderr())
    write("Two plots will be created, with <plotname>.png and <plotname>.pdf.", stderr())
    write("If you provide S as an optional 4th argument, input data will be", stderr())
    write("    filtered to only singlet assignments from both input files.", stderr())
    q()
}

assn1 <- read.table(args[1])
assn2 <- read.table(args[2])
outname <- args[3]

if (length(args) > 3 && (args[4] == "s" | args[4] == "S")){
    assn1 <- assn1[which(assn1$V3=="S"),]
    assn2 <- assn2[which(assn2$V3=="S"),]
}

assn1 <- assn1[,c(1,2,3)]
assn2 <- assn2[,c(1,2,3)]
colnames(assn1) <- c("bc", "assn1", "type1")
colnames(assn2) <- c("bc", "assn2", "type2")

merged <- merge(assn1, assn2)

tab <- as.data.frame(table(merged[,c(2,4)]))
tot1 <- as.data.frame(table(merged[,c(2), drop=FALSE]))
tot2 <- as.data.frame(table(merged[,c(4),drop=FALSE]))
colnames(tot1)[2] <- "tot1"
colnames(tot2)[2] <- "tot2"
tab <- merge(tab, tot1)
tab <- merge(tab, tot2)
tab$jaccard <- tab$Freq / (tab$tot1 + tab$tot2 - tab$Freq)

# Attempt to get these in the same order
# Order both by decreasing counts of singlets
tot1S <- as.data.frame(table(merged[which(merged$type1=="S"),c(2),drop=FALSE]))
tot2S <- as.data.frame(table(merged[which(merged$type2=="S"),c(4),drop=FALSE]))

nS1 <- length(rownames(tot1S))
nS2 <- length(rownames(tot2S))

tot1S <- tot1S[order(tot1S$Freq, decreasing=T),]
tot2S <- tot2S[order(tot2S$Freq, decreasing=T),]

tot1S$rank <- seq(1, length(rownames(tot1S)))
tot2S$rank <- seq(1, length(rownames(tot2S)))

tot1S <- tot1S[,c(1,3)]
tot2S <- tot2S[,c(1,3)]

assn1U <- unique(assn1[,c(2,3)])
assn2U <- unique(assn2[,c(2,3)])

assn1U$name1 <- apply(assn1U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][1] })
assn1U$name2 <- apply(assn1U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][2] })
assn2U$name1 <- apply(assn2U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][1] })
assn2U$name2 <- apply(assn2U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][2] })

colnames(tot1S) <- c("name1", "rank1")
assn1U <- merge(assn1U, tot1S)
colnames(tot1S) <- c("name2", "rank2")
assn1U <- merge(assn1U, tot1S, all.x=TRUE)

colnames(tot2S) <- c("name1", "rank1")
assn2U <- merge(assn2U, tot2S)
colnames(tot2S) <- c("name2", "rank2")
assn2U <- merge(assn2U, tot2S, all.x=TRUE)

assn1U[which(is.na(assn1U$rank2)),]$rank2 <- 0
assn2U[which(is.na(assn2U$rank2)),]$rank2 <- 0

assn1U$rank1b <- assn1U$rank1
assn1U$rank2b <- assn1U$rank2
assn1U[which(assn1U$rank2 < assn1U$rank1),]$rank1b <- 
    assn1U[which(assn1U$rank2 < assn1U$rank1),]$rank2
assn1U[which(assn1U$rank2 < assn1U$rank1),]$rank2b <- 
    assn1U[which(assn1U$rank2 < assn1U$rank1),]$rank1

assn2U$rank1b <- assn2U$rank1
assn2U$rank2b <- assn2U$rank2
assn2U[which(assn2U$rank2 < assn2U$rank1),]$rank1b <- 
    assn2U[which(assn2U$rank2 < assn2U$rank1),]$rank2
assn2U[which(assn2U$rank2 < assn2U$rank1),]$rank2b <- 
    assn2U[which(assn2U$rank2 < assn2U$rank1),]$rank1

assn1U <- assn1U[order(assn1U$rank2b),]
assn1U <- assn1U[order(assn1U$rank1b),]
assn2U <- assn2U[order(assn2U$rank2b),]
assn2U <- assn2U[order(assn2U$rank1b),]

tab$assn1 <- factor(tab$assn1, labels=assn1U$assn1, levels=assn1U$assn1)
tab$assn2 <- factor(tab$assn2, labels=assn2U$assn2, levels=assn2U$assn2)

plt <- ggplot(tab) + 
    geom_tile(aes(x=assn1, y=assn2, fill=jaccard)) + 
    theme_bw() + 
    scale_fill_viridis("Jaccard index", option="inferno") + 
    scale_x_discrete("Assignment 1") + 
    scale_y_discrete("Assignment 2") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14, vjust=2), 
        legend.position="top")

ggsave(plt, file=paste(outname, '.pdf', sep=''), width=8, height=7)
ggsave(plt, file=paste(outname, '.png', sep=''), width=8, height=7, units='in')

