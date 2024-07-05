#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
checks <- c(length(args) < 1, !file.exists(paste(args[1], '.contam_rate', sep='')),
    !file.exists(paste(args[1], '.contam_prof', sep='')))

if (sum(checks) > 0){
    write("USAGE:", stderr())
    write("contam.R <output_prefix>", stderr())
    write("After you have run quant_contam on a data set and have the files", stderr())
    write("    [output_prefix].contam_rate and [output_prefix].contam_prof,", stderr())
    write("    run this program to plot summary information about ambient RNA.", stderr())
    write("It will create a plot in the file [output_prefix].contam.png.", stderr())
    q()
}
basename <- args[1]

cr <- read.table(paste(basename, '.contam_rate', sep=''))
assn <- read.table(paste(basename, '.assignments', sep=''))

# Strip library names, etc. from cell barcodes in case they don't match
assn$V1 <- gsub("[^ACGT]+", "", assn$V1)
cr$V1 <- gsub("[^ACGT]+", "", cr$V1)

names <- unique(assn$V2)

assnorig <- assn

doubopt <- FALSE
if (length(args) > 1 & (args[2] == 'D' | args[2] == 'd')){
    # Include specific doublets.
    names1 <- apply(assn, 1, function(x){ strsplit(x[2], '+', fixed=T)[[1]][1] })
    names2 <- apply(assn, 1, function(x){ strsplit(x[2], '+', fixed=T)[[1]][2] })
    names <- unique(c(names1, names2))
    names <- names[which(!is.na(names))]
    names <- unique(c(assn$V2, names))
    names <- names[order(names)]
    doubopt <- TRUE
} else{
    if (length(rownames(assn[which(assn$V3=="D"),])) > 0){
        names <- c(unique(assn[which(assn$V3=="S"),]$V2, "Doublet"))
        assn[which(assn$V3=="D"),]$V2 <- "Doublet"
    }
}

colnames(cr) <- c("bc", "cr", "cr.se")

p_contam <- pnorm(0, mean(cr[which(cr$cr > 0 & cr$cr < 1),]$cr), 
    sd(cr[which(cr$cr > 0 & cr$cr < 1),]$cr), lower.tail=FALSE)
mu_contam <- mean(cr$cr)
sd_contam <- sd(cr$cr)

colnames(assn) <- c("bc", "id", "type", "llr")
cr <- merge(cr, assn)

names <- c(names, "other_species")
names <- factor(names, labels=names, levels=names)
ncol <- length(names)

pal_default <- c("#f94144", "#f3722c", "#f8961e", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1")
pal_default <- c("#f94144", "#277da1", "#f3722c", "#577590", "#f8961e", "#4d908e", "#f9844a", "#43aa8b", "#f9c74f", "#90be6d")
if (ncol <= 20){
    #pal <- pal_locuszoom("default")(ncol)
    #pal <- pal_default[1:ncol]
    pal <- pal_d3("category20")(ncol)
} else{
    #pal <- colorRampPalette(pal_locuszoom("default")(7))(ncol)
    #pal <- colorRampPalette(pal_default)(ncol)
    pal <- colorRampPalette(pal_d3("category20")(20))(ncol)
}
name2col <- setNames(pal, names)
out_png <- paste(basename, '.contam.png', sep='')
out_pdf <- paste(basename, '.contam.pdf', sep='')

cr <- cr[order(cr$cr, decreasing=T),]
cr <- cr[order(cr$id),]
cr$bc <- factor(cr$bc, labels=cr$bc, levels=cr$bc)
bars <- ggplot(cr) + 
    geom_bar(aes(x=bc, y=cr, fill=id), show.legend=FALSE, stat='identity') + 
    geom_errorbar(aes(x=bc, ymin=cr-cr.se, ymax=cr+cr.se)) +
    coord_flip() + 
    scale_y_continuous("Contamination rate", limits=c(0,1)) +
    scale_x_discrete("Cells") +
    scale_fill_manual(values=name2col) +
    theme_bw() + 
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())

boxes <- ggplot(cr) + 
    geom_boxplot(aes(x=id, y=cr, colour=id), show.legend=FALSE) + 
    coord_flip() + 
    scale_y_continuous("Contamination rate") +
    theme_bw() +
    scale_colour_manual(values=name2col) +
    theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor=element_blank())

dens <- ggplot(cr) + 
    geom_density(aes(x=cr)) +
    theme_bw() + 
    scale_x_continuous("Contam rate", limits=c(0,1)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          panel.background=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_blank())

byllr <- ggplot(cr) + 
    geom_point(aes(x=llr, y=cr, colour=id), show.legend=FALSE) + 
    geom_errorbar(aes(x=llr, ymin=ifelse(cr-cr.se < 0, 0, cr-cr.se),
        ymax=ifelse(cr+cr.se > 1, 1, cr+cr.se), colour=id), show.legend=FALSE) + 
    theme_bw() + 
    annotation_logticks(side="b") +
    scale_colour_manual(values=name2col) + 
    scale_x_log10("LLR of ID") + 
    scale_y_continuous("Contam rate") + 
    theme(panel.background=element_blank(), 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank())

title = ggdraw() + 
    draw_label(basename, size=12, fontface='bold', lineheight=0.9, hjust=0.6) +
    theme(plot.margin=margin(0,0,0,0))

mutxt <- format(mu_contam, scientific=TRUE, digits=2)
sdtxt <- format(sd_contam, scientific=TRUE, digits=2)
eq <- '='
mu_sigma_text <- bquote(mu~.(eq)~.(mutxt)~sigma~.(eq)~.(sdtxt))
#pval <- ggdraw() + 
#    draw_label(paste("P(contamination) = ", format(p_contam, scientific=TRUE, digits=2), sep=''), 
#           size=10, lineheight=0.9, hjust=0.6)

mu_sd <- ggdraw() + 
    draw_label(mu_sigma_text, size=10, lineheight=0.9, hjust=0.6)

make_plot <- function(){
    cpfile <- paste(basename, '.contam_prof', sep='')
    if (file.exists(cpfile)){
        has_cp_file <- 1
        cp <- read.table(cpfile)
        
        colnames(assnorig) <- c("bc", "id", "type", "llr") 
        if (length(rownames(assnorig[which(assnorig$type=="D"),])) > 0){
            assndoub <- assnorig[which(assnorig$type=="D"),]
            assndoub$id1 <- apply(assndoub, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][1] })
            assndoub$id2 <- apply(assndoub, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][2] })
            assndoub1 <- assndoub[,c(1,5,3,4)]
            assndoub2 <- assndoub[,c(1,6,3,4)]
            colnames(assndoub1) <- c("bc", "id", "type", "llr")
            colnames(assndoub2) <- c("bc", "id", "type", "llr")
            # Count singlets twice (relative to doublets)
            assnorig <- rbind(assnorig[which(assnorig$type=="S"),], 
                assnorig[which(assnorig$type=="S"),],
                assndoub1,
                assndoub2)
        }
        assn_agg <- as.data.frame(table(assnorig$id))
        assn_agg_counts <- assn_agg
        assn_agg$Freq <- assn_agg$Freq / sum(assn_agg$Freq)
        colnames(assn_agg) <- c("var", "val")
        assn_agg$type <- "cells"
        colnames(cp) <- c("var", "val")
        cp$type <- "amb"
        
        # Calculate entropy of both
        H_amb <- -sum(cp[which(cp$var != "other_species"),]$val*log2(cp[which(cp$var != "other_species"),]$val)) 
        H_cells <- -sum(assn_agg$val*log2(assn_agg$val))
        H_max <- log2(length(rownames(cp[which(cp$var != "other_species"),])))
        # This will drop "other_species" automatically
        bothdf <- merge(cp, assn_agg, by="var")
        KL_div <- -sum(bothdf$val.x*log(bothdf$val.x/bothdf$val.y))

        # Create text for plot 
        entropy_plt <- ggdraw() + 
            draw_label(paste("KL div = ", 
                format(KL_div, scientific=FALSE, digits=3), sep=''), 
                size=10, lineheight=0.5, hjust=0.6)
        
        # Write out some stats about ambient RNA composition
        assn_agg_counts$not <- sum(assn_agg_counts$Freq) - assn_agg_counts$Freq
        colnames(assn_agg_counts) <- c("var", "A", "B")
        beta_merge <- merge(cp[which(cp$var != "other_species"),c(1,2)], assn_agg_counts)
        beta_merge$p <- pbeta(beta_merge$val, beta_merge$A, beta_merge$B, lower.tail=FALSE, log=TRUE) 
        beta_merge$logdiff <- log(beta_merge$val) - log(beta_merge$A/(beta_merge$A + beta_merge$B))
        
        stats <- data.frame(lib=basename, type=c("all", "all", "all"), var=c("H.amb", "H.cells", "KL.div"), 
            val=c(H_amb/H_max, H_cells/H_max, KL_div))
        stats <- rbind(stats, data.frame(lib=basename, type=beta_merge$var, var="logdiff", val=beta_merge$logdiff))
        stats <- rbind(stats, data.frame(lib=basename, type=beta_merge$var, var="log.p.enriched", val=beta_merge$p))

        write.table(stats, file=paste(basename, '.contam.stats', sep=''), sep='\t', quote=FALSE, row.names=FALSE,
            col.names=FALSE)

        cp_assn <- rbind(cp, assn_agg)
        
        cp_assn$type <- factor(cp_assn$type, labels=c("cells", "amb"), levels=c("cells", "amb"))
        
        # Don't show legend if more than 10 individuals
        sl <- TRUE
        if (!doubopt & length(rownames(assn_agg)) > 10){
            sl <- FALSE
        }
        cp_plt <- ggplot(cp_assn) + 
            geom_bar(aes(x=type, y=val, fill=var), show.legend=sl, stat='identity') +
            scale_fill_manual("", values=name2col) +
            theme_bw() + 
            scale_y_continuous("Proportion") +
            theme(axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1), 
                  axis.ticks.x=element_blank(), 
                  axis.title.x=element_blank(),
                  #axis.title.y=element_text(size=10),
                  axis.title.y=element_blank(),
                  axis.text.y=element_text(size=9), 
                  legend.position='right', 
                  legend.direction='vertical',
                  legend.text=element_text(size=8),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(), 
                  panel.grid.minor=element_blank(), 
                  panel.background=element_blank(),
                  legend.key.size=unit(0.35, 'cm'))
        
        row1col1 <- plot_grid(title, mu_sd, entropy_plt, cp_plt, ncol=1, 
            rel_heights=c(0.2*(1/3), 0.2*(1/3), 0.2*(1/3), 0.8))
        row1col2 <- plot_grid(dens, byllr, ncol=1, rel_heights=c(0.45, 0.55))
        row1 <- plot_grid(row1col1, row1col2, ncol=2, rel_widths=c(0.4, 0.6))
        row2 <- plot_grid(bars, boxes, ncol=2, rel_widths=c(0.33,0.66)) 
        plot_grid(row1, row2, nrow=2, rel_heights=c(0.45, 0.55))

    } else{
        row1 <- plot_grid(title, mu_sd, ncol=2, rel_widths=c(0.5,0.5))
        row2 <- plot_grid(dens, byllr, ncol=2, rel_widths=c(0.45, 0.55))
        row3 <- plot_grid(bars, boxes, ncol=2, rel_widths=c(0.33, 0.66))
        plot_grid(row1, row2, row3, ncol=1, rel_heights=c(0.05, 0.275, 0.675))
    }
}

png(out_png, bg='white', width=6, height=9, units='in', res=150)
make_plot()
dev.off()

pdf(out_pdf, bg='white', width=6, height=9)
make_plot()
dev.off()


