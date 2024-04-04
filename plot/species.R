#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Please provide a demux_species output directory.", stderr())
    q()
}

outdir <- args[1]
countsfile <- paste(outdir, 'species_counts.txt', sep='/')
namesfile <- paste(outdir, 'species_names.txt', sep='/')
assnfile <- paste(outdir, 'bc_assignments.txt', sep='/')

counts <- read.table(countsfile)
names <- read.table(namesfile)
assn <- read.table(assnfile)

if (length(colnames(assn)) > 4){
    assn$V4 <- assn$V5
    assn <- assn[,c(1,2,3,4)]
}

colnames(assn) <- c("bc", "species", "type", "llr")

counts <- counts[,seq(1, length(rownames(names))+1)]
colnames(counts) <- c("bc", names$V2)
colsort <- colnames(counts[,-c(1)])
colsort <- c(1, order(colsort)+1)
counts <- counts[,colsort]

assn$bc <- gsub("-1", "", assn$bc)
counts$bc <- gsub("-1", "", counts$bc)

counts$tot <- rowSums(counts[2:length(colnames(counts))])

counts <- merge(counts, assn)
counts <- counts[order(counts$llr),]
counts <- counts[order(counts$species),]
rownames(counts) <- counts$bc
assn[which(assn$type=="D"),]$species <- "Doublet"
spec <- unique(assn[which(assn$type=="S"),]$species)
spec <- spec[order(spec)]
spec <- c(spec, "Doublet")
assn$species <- factor(assn$species, labels=spec, levels=spec)
rownames(assn) <- assn$bc
assn$count <- 1
assn <- merge(assn, counts[,which(colnames(counts) %in% c("bc", "tot"))])
assncounts <- aggregate(assn$count, by=list(species=assn$species), FUN=sum)
assncounts$species <- factor(assncounts$species, labels=spec, levels=spec)
assncounts$x <- assncounts$x / sum(assncounts$x)
assncounts_weight <- aggregate(assn$tot, by=list(species=assn$species), FUN=sum)
assncounts_weight$species <- factor(assncounts_weight$species, labels=spec, levels=spec)
assncounts_weight$x <- assncounts_weight$x / sum(assncounts_weight$x)
colnames(assncounts)[2] <- "frac"
colnames(assncounts_weight)[2] <- "frac_weighted"
assncounts <- merge(assncounts, assncounts_weight, by='species')

print(assncounts[1:10,])
assn <- assn[,c(2), drop=FALSE]

for (n in names$V2){
    counts[[n]] <- counts[[n]] / counts$tot
}
counts <- counts[,-c(1,length(colnames(counts))-3, length(colnames(counts))-2, length(colnames(counts))-1, length(colnames(counts)))]

# Load distribution data
dists <- read.table(paste(outdir, 'dists.txt', sep='/'), header=T, row.names=1)
distann <- dists[,c(1),drop=FALSE]
distann$species <- rownames(distann)
distann[grep(".", distann$species, fixed=T),]$species <- "Doublet"
distann$species <- factor(distann$species, labels=spec, levels=spec)
distann <- distann[,-c(1),drop=FALSE]
dists <- dists[,-c(1)]
dists <- dists[,order(colnames(dists))]
dists <- dists[order(rownames(dists)),]

library(pheatmap)
library(cowplot)
library(ggplot2)

ncolor <- 50
cols=colorRampPalette(c("#255957", "#7C9885", "#E2F1AF"))(ncolor)

cells <- pheatmap(as.matrix(counts),
    cluster_cols=FALSE,
    cluster_rows=FALSE,
    scale='none',
    color=cols,
    breaks=seq(0,1,1/ncolor),
    show_rownames=FALSE,
    show_colnames=FALSE,
    annotation_row=assn)

dists <- pheatmap(as.matrix(dists),
    cluster_rows=FALSE,
    cluster_cols=FALSE,
    scale='none',
    color=cols,
    legend=FALSE,
    annotation_legend=FALSE,
    annotation_names_row=FALSE,
    angle_col=90,
    width=5,
    breaks=seq(0, 1, 1/ncolor),
    show_rownames=FALSE,
    annotation_row=distann)

percs <- ggplot(assncounts) + 
    geom_bar(aes(x="Prop.", fill=species, y=frac), stat='identity') + 
    geom_bar(aes(x="Prop. weighted",, fill=species, y=frac_weighted), stat='identity') + 
    scale_fill_discrete("species") + 
    theme(axis.ticks.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
    theme_nothing()

row2 <- cowplot::plot_grid(dists$gtable, NULL, nrow=1, rel_widths=c(0.8, 0.2))

col1 <- plot_grid(cells$gtable, row2, ncol=1, nrow=2, rel_heights=c(2/3,1/3))
col2 <- plot_grid(percs)

pdf(paste(outdir, '/', 'counts.pdf', sep=''), bg='white', width=8.5, height=7)

#cowplot::plot_grid(cells$gtable, row2, ncol=1, nrow=2, rel_heights=c(2/3,1/3))
cowplot::plot_grid(col1, col2, ncol=2, rel_widths=c(0.8,0.2))

dev.off()


