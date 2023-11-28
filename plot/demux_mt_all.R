#! /usr/bin/env Rscript
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("USAGE:\n", stderr())
    write("demux_mt_all.R <output_prefix>\n", stderr())
    write("Creates PDF and PNG plots named <output_prefix>.unclust.pdf and\n", stderr())
    write("    <output_prefix>.unclust.png.\n", stderr())
    write("In contrast to demux_mt_reclust.R, does not hierarchically cluster\n", stderr())
    write("    cells. Plots cells at all available sites, sorted by assigned\n", stderr())
    write("    haplotype. Useful for visualizing results in cases where\n", stderr())
    write("    hierarchical clustering seems to disagree with haplotype\n", stderr())
    write("    assignments.\n", stderr())
    q()
}

pdf_out <- paste(args[1], '.unclust.pdf', sep="")
png_out <- paste(args[1], '.unclust.png', sep="")

fn <- paste(args[1], ".cellhaps", sep="")
tab <- read.table(fn, header=T)

# Drop sites that are not variable within this subset
#keepsites <- which(apply(tab, 2, function(x){ sum( 
#    as.numeric( x[!is.na(x)] )) }) > length(rownames(tab))/1000)
#keepsites <- c(keepsites, which(colnames(tab) == "bc"))
#tab <- tab[,keepsites]

ann_row <- NA
assnfile <- paste(args[1], '.assignments', sep='')
if (file.exists(assnfile)){
    assn <- read.table(assnfile)
    colnames(assn) <- c("bc", "cluster", "type", "llr")
     
    tab <- merge(tab, assn, by="bc")
    tab <- tab[order(tab$cluster),]
    rownames(tab) <- tab$bc
    tab <- tab[which(! colnames(tab) %in% c("bc", "cluster", "type", "llr"))]

    # Replace all unique doublet combinations with 
    # a single annotation to avoid creating too many
    # colors
    if (length(rownames(assn[which(assn$type=="D"),])) > 0){
        assn[which(assn$type=="D"),]$cluster <- "Doublet"
    }
    assn$cluster <- factor(assn$cluster)
    assn <- assn[,which(colnames(assn) %in% c("bc", "cluster"))]
    rownames(assn) <- assn$bc
    
    assn <- assn[,c(2),drop=FALSE]
    ann_row <- assn
}

clusts <- NA

pdf(pdf_out, bg='white')
res <- pheatmap(as.matrix(tab), 
    cluster_cols=FALSE, 
    cluster_rows=FALSE,
    annotation_row=ann_row,
    legend=FALSE, 
    show_rownames=FALSE, 
    show_colnames=FALSE, 
    fontsize_col=3, 
    fontsize_row=0.5,
    angle_col=90, 
    scale='none', 
    cutree_rows=clusts,
    color=c('#0267C1', '#2AAE24', '#EFA00B'))
dev.off()

png(png_out, bg='white')
res <- pheatmap(as.matrix(tab), 
    cluster_cols=FALSE,
    cluster_rows=FALSE, 
    annotation_row=ann_row,
    legend=FALSE, 
    show_rownames=FALSE, 
    show_colnames=FALSE, 
    fontsize_col=3, 
    fontsize_row=0.5,
    angle_col=90, 
    scale='none', 
    cutree_rows=clusts,
    color=c('#0267C1', '#2AAE24', '#EFA00B'))
dev.off()

