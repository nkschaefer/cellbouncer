#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(cowplot)
library(viridis)
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Please provide a demux_species output directory.", stderr())
    q()
}

outdir <- args[1]
countsfile <- paste(outdir, 'species_counts.txt', sep='/')
namesfile <- paste(outdir, 'species_names.txt', sep='/')
assnfile <- paste(outdir, 'species.filt.assignments', sep='/')
assnfile_full <- paste(outdir, "species.assignments", sep='/')

counts <- read.table(countsfile)
names <- read.table(namesfile)
assn <- read.table(assnfile)
assnfull <- read.table(assnfile_full)

# Strip extra stuff off of cell barcodes to ensure counts & assignments
# files can match barcodes
assn$V1 <- gsub("[^ACGT]+", "", assn$V1)
assnfull$V1 <- gsub("[^ACGT]+", "", assnfull$V1)

if (length(colnames(assn)) > 4){
    assn$V4 <- assn$V5
    assn <- assn[,c(1,2,3,4)]
}

colnames(assn) <- c("bc", "species", "type", "llr")
colnames(assnfull) <- c("bc", "species", "type", "llr")

snames <- unique(assnfull$species)
snames <- factor(snames, labels=snames, levels=snames)
ncol <- length(snames)
pal <- colorRampPalette(pal_d3("category20")(20))(ncol)
if (ncol < 20){
    pal <- pal_d3("category20")(ncol)
}
name2col <- setNames(pal, snames)

counts <- counts[,seq(1, length(rownames(names))+1)]
colnames(counts) <- c("bc", names$V2)
colsort <- colnames(counts[,-c(1)])
colsort <- c(1, order(colsort)+1)
counts <- counts[,colsort]

assn$bc <- gsub("-1", "", assn$bc)
counts$bc <- gsub("-1", "", counts$bc)

counts$tot <- rowSums(counts[2:length(colnames(counts))])

counts <- merge(counts, assn)

for (cn in colnames(counts)[seq(2, length(colnames(counts))-3)]){
    counts[[cn]] <- counts[[cn]] / counts$tot
}

metacols <- seq(length(colnames(counts))-3, length(colnames(counts)))
meta <- counts[,metacols]
meta <- meta[,c(2), drop=FALSE]
rownames(meta) <- counts$bc

countsm <- as.matrix(counts[,-c(c(1), metacols)])
rownames(countsm) <- counts$bc

basename <- gsub("/$", "", args[1])
basename_split <- strsplit(basename, '/')[[1]]
basename_title <- basename_split[length(basename_split)]

title = ggdraw() + 
    draw_label(basename_title, size=12, fontface='bold', lineheight=0.9, hjust=0.6) +
    theme(plot.margin=margin(0,0,0,0))

name2col2 <- list(species=name2col)
heatm <- pheatmap(countsm, 
    scale='none', 
    cluster_cols=FALSE, 
    show_rownames=FALSE, 
    annotation_row=meta,
    annotation_colors=name2col2,
    annotation_legend=FALSE,
    breaks=seq(0,1,0.01),
    angle_col=90,
    color=mako(100),
    main='Normalized k-mer counts per cell')

assn$count <- 1
assnfull$count <- 1
assnagg <- aggregate(assn$count, by=list(species=assn$species), FUN=sum)
assnagg$x <- assnagg$x / sum(assnagg$x)
assnfullagg <- aggregate(assnfull$count, by=list(species=assnfull$species), FUN=sum)
assnfullagg$x <- assnfullagg$x / sum(assnfullagg$x)

colnames(assnagg)[2] <- "Frac"
colnames(assnfullagg)[2] <- "Frac"
assnagg$type <- "Filtered"
assnfullagg$type <- "All"

assnbothagg <- rbind(assnagg, assnfullagg)

bars <- ggplot(assnbothagg) + geom_bar(aes(x=type, y=Frac, fill=species), stat='identity') + 
    theme_bw() +
    scale_y_continuous("Fraction of barcodes") + 
    scale_x_discrete("Cell group") +
    scale_fill_manual(values=name2col) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

row2 <- plot_grid(heatm$gtable, bars, ncol=2, rel_widths=c(0.6, 0.4))

name_png <- paste(basename, '/species.png', sep='')
name_pdf <- paste(basename, '/species.pdf', sep='')

if (args[1] == '.'){
    name_png <- 'species.png'
    name_pdf <- 'species.pdf'
}

png(name_png, width=8, height=6, bg='white', units='in', res=150)
plot_grid(title, row2, nrow=2, rel_heights=c(0.1,0.9))
dev.off()

pdf(name_pdf, width=8, height=6, bg='white')
plot_grid(title, row2, nrow=2, rel_heights=c(0.1,0.9))
dev.off()

