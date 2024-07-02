#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    write("USAGE:", stderr())
    write("compare_assignments.R <out1.assignments> <out2.assignments> <plotname> (S/S1/S2/D/D1/D2)", stderr())
    write("Given two .assignments files from cellbouncer programs run on", stderr())
    write("    the same cells, compares them to each other and plots a heatmap", stderr())
    write("    that should show whether the assignments tend to agree well,", stderr())
    write("    and if so, which labels from the first file match which labels", stderr())
    write("    from the second file. An example use case would be running demux_vcf", stderr())
    write("    and demux_mt on the same data, and checking which mitochondrial", stderr())
    write("    haplotypes match which individuals in the VCF.", stderr())
    write("Two plots will be created, with <plotname>.png and <plotname>.pdf.", stderr())
    write("If you provide S as an optional 4th argument, it can affect how ", stderr())
    write("    assignments are filtered or sorted.", stderr())
    write("    S means limit to only singlets for both input files.", stderr())
    write("    S1 means limit to only singlets in the first input file.", stderr())
    write("    S2 means limit to only singlets in the second input file.", stderr())
    write("    D means treat doublets as unique categories - do not sort based", stderr())
    write("      on frequency of component singlets", stderr())
    write("    D1 means the same behavior as above, but only for the first file", stderr())
    write("    D2 means the same behavior as above, but only for the second file", stderr())
    write("    You can also combine S and D options together for the different files", stderr())
    write("      For example: S1D2 or S2D1. You cannot pass a contradictory combination.", stderr())
    q()
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))

if (nchar(args[1]) < 12 | substr(args[1], nchar(args[1])-11, nchar(args[1])) != ".assignments"){
    args[1] <- paste(args[1], '.assignments', sep='')
}
if (nchar(args[2]) < 12 | substr(args[2], nchar(args[2])-11, nchar(args[2])) != ".assignments"){
    args[2] <- paste(args[2], '.assignments', sep='')
}
assn1 <- read.table(args[1])
assn2 <- read.table(args[2])
outname <- args[3]

if (length(args) > 3){
    if (args[4] == "S" | args[4] == "S1" | args[4] == "S1D2"){
        assn1 <- assn1[which(assn1$V3=="S"),]
    }
    if (args[4] == "S" | args[4] == "S2" | args[4] == "S2D1"){
        assn2 <- assn2[which(assn2$V3=="S"),]
    }
}

assn1 <- assn1[,c(1,2,3)]
assn2 <- assn2[,c(1,2,3)]
colnames(assn1) <- c("bc", "assn1", "type1")
colnames(assn2) <- c("bc", "assn2", "type2")

merged <- merge(assn1, assn2)
tab <- as.data.frame(table(merged[,c(2,4)]))
mfilt <- NA
if (length(args) > 3){
    if (args[4] == "D"){
        mfilt <- merged
    } else if (args[4] == "D1" | args[4] == "S2D1"){
        mfilt <- merged[which(merged$type2=="S"),]
    } else if (args[4] == "D2" | args[4] == "S1D2"){
        mfilt <- merged[which(merged$type1=="S"),]
    } else{
        mfilt <- merged[which(merged$type1=="S" & merged$type2=="S"),]
    }
} else{
    mfilt <- merged[which(merged$type1=="S" & merged$type2=="S"),]
}
tabS <- as.data.frame(table(mfilt[,c(2,4)]))
tot1 <- as.data.frame(table(merged[,c(2), drop=FALSE]))
tot2 <- as.data.frame(table(merged[,c(4),drop=FALSE]))
colnames(tot1)[2] <- "tot1"
colnames(tot2)[2] <- "tot2"

if (length(args) > 3){
    if (args[4] == "S"){
        print(table(merged[which(merged$type1=="S" & merged$type2=="S"),c(2,4)]))
    } else if (args[4] == "S1" | args[4] == "S1D2" | args[4] == "D2"){
        # File 2 has doublets, so will take up more screen space.
        # Switch rows/cols so file 2 prints vertically
        print(t(table(merged[which(merged$type1=="S"),c(2,4)])))
    } else if (args[4] == "S2" | args[4] == "S2D1" | args[4] == "D1"){
        print(table(merged[which(merged$type2=="S"),c(2,4)]))
    } else if (args[4] == "D"){
        print(table(merged[,c(2,4)]))  
    } else{
        print(table(merged[which(merged$type1=="S" & merged$type2=="S"),c(2,4)]))
    }
} else{
    print(table(merged[which(merged$type1=="S" & merged$type2=="S"),c(2,4)]))
}

tab <- merge(tab, tot1)
tab <- merge(tab, tot2)
tab$jaccard <- tab$Freq / (tab$tot1 + tab$tot2 - tab$Freq)

tabS <- merge(tabS, tot1)
tabS <- merge(tabS, tot2)
tabS$jaccard <- tabS$Freq / (tabS$tot1 + tabS$tot2 - tabS$Freq)

# Attempt to get these in the same order
assn1U <- unique(assn1[,c(2,3)])
assn2U <- unique(assn2[,c(2,3)])

if (length(args) > 3 & (args[4] == "D" | args[4] == "D1" | args[4] == "S2D1")){
    assn1U$name1 <- assn1U$assn1
    assn1U$name2 <- NA
} else{
    assn1U$name1 <- apply(assn1U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][1] })
    assn1U$name2 <- apply(assn1U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][2] })
}
if (length(args) > 3 & (args[4] == "D" | args[4] == "D2" | args[4] == "S1D2")){
    assn2U$name1 <- assn2U$assn2
    assn2U$name2 <- NA
} else{
    assn2U$name1 <- apply(assn2U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][1] })
    assn2U$name2 <- apply(assn2U, 1, function(x){ strsplit(x[1], "+", fixed=T)[[1]][2] })
}

tabS <- tabS[order(tabS$jaccard, decreasing=T),]
tabS1 <- data.frame("name1"=unique(tabS$assn1))
tabS1$rank1 <- seq(1, length(rownames(tabS1)))
tabS2 <- data.frame("name1"=unique(tabS$assn2))
tabS2$rank1 <- seq(1, length(rownames(tabS2)))

assn1U <- merge(assn1U, tabS1)
colnames(tabS1) <- c("name2", "rank2")
assn1U <- merge(assn1U, tabS1, all.x=TRUE)

assn2U <- merge(assn2U, tabS2)
colnames(tabS2) <- c("name2", "rank2")
assn2U <- merge(assn2U, tabS2, all.x=TRUE)

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
    scale_fill_viridis("Jaccard index", option="inferno", limits=c(0,1)) + 
    scale_x_discrete("Assignment 1") + 
    scale_y_discrete("Assignment 2") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10), 
        axis.text.y=element_text(size=10), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14, vjust=2), 
        legend.position="top")

ggsave(plt, file=paste(outname, '.pdf', sep=''), width=8, height=7)
ggsave(plt, file=paste(outname, '.png', sep=''), width=8, height=7, units='in')

tot_match <- sum(tab[which(as.character(tab$assn1) == as.character(tab$assn2)),]$Freq)
tot_all <- sum(tab$Freq)
# Only print matching if some labels agree
if (tot_match > 0){
    write(paste("Overall agreement = ", tot_match/tot_all, sep=""), stderr())
}
