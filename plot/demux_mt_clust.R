#! /usr/bin/env Rscript
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("USAGE:", stderr())
    write("demux_mt_clust.R <output_prefix>", stderr())
    write("OR", stderr())
    write("demux_mt_clust.R <output_prefix> frac", stderr())
    write("    Where frac is a decimal between 0 and 1: the upper limit", stderr())
    write("    on what fraction of sites can be missing in a cell for", stderr())
    write("    it to be included (default 0.5). Higher means more cells", stderr())
    write("    included, but greater risk of messing up clustering.", stderr())
    write("    If you get the \"NA in foreign function call\" error,", stderr())
    write("    try decreasing this number.", stderr())
    write("Given output from a demux_mt run, helps visualize and validate", stderr())
    write("    results by hierarchically clustering cells that pass a", stderr())
    write("    coverage threshold at the sites used to build the", stderr())
    write("    mitochondrial haplotypes. Shows cluster assignments next to", stderr())
    write("    hierarchically clustered cells, so users can evaluate whether", stderr())
    write("    the cluster assignments look reasonable.", stderr())
    write("Creates PDF and PNG plots named <output_prefix>.clust.pdf and", stderr())
    write("    <output_prefix>.clust.png.", stderr())
    write("If you provide D as an additional, final argument, it will keep", stderr())
    write("    from using a single \"Doublet\" label to show all doublet IDs", stderr())
    write("    and will instead show specific doublet types on the plot.", stderr())
    q()
}

frac <- 0.5
if (length(args) > 1){
    frac <- as.numeric(args[2])
}

pdf_out <- paste(args[1], '.clust.pdf', sep="")
png_out <- paste(args[1], '.clust.png', sep="")

fn <- paste(args[1], ".cellhaps", sep="")
tab <- read.table(fn, header=T)

if ("bc" %in% colnames(tab)){
    rownames(tab) <- tab$bc
    tab <- tab[,-c(1)]
}

#if (length(args) > 2){
#    nsites <- as.numeric(args[3])
#    if (nsites > 0 && nsites < length(colnames(tab))){
#        csdf <- data.frame(cs=colSums(tab, na.rm=TRUE))
#        csdf$idx <- seq(1, length(rownames(csdf)))
#        csdf <- csdf[order(csdf$cs, decreasing=T),]
#        idx_keep <- csdf[seq(1,nsites),]$idx
#        tab <- tab[,idx_keep]    
#    }
#}


tab$fracna <- apply(tab, 1, function(x){ sum(is.na(x)) })
tab$fracna <- tab$fracna / (length(colnames(tab))-1)

# Drop cells that don't pass criterion
tab <- tab[which(tab$fracna < frac),-c(length(colnames(tab)))]

# Drop sites that are not variable within this subset
keepsites <- which(apply(tab, 2, function(x){ sum( 
    as.numeric( x[!is.na(x)] )) }) > length(rownames(tab))/1000)
tab <- tab[,keepsites]

ann_row <- NA
assnfile <- paste(args[1], '.assignments', sep='')
if (file.exists(assnfile)){
    assn <- read.table(assnfile)
    # Replace all unique doublet combinations with 
    # a single annotation to avoid creating too many
    # colors
    if (length(args) < 3 | args[3] != "D"){
        if (length(rownames(assn[which(assn$V3=="D"),])) > 0){
            assn[which(assn$V3=="D"),]$V2 <- "Doublet"
        }
    }
    assn$V2 <- factor(assn$V2)
    assn <- assn[,c(1,2)]
    colnames(assn) <- c("bc", "cluster")
    rownames(assn) <- assn$bc
    assn <- assn[,c(2),drop=FALSE]
    ann_row <- assn
}

clusts <- NA

pdf(pdf_out, bg='white')
res <- pheatmap(as.matrix(tab), 
    cluster_cols=FALSE, 
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

