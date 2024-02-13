#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
basename <- args[1]

cr <- read.table(paste(basename, '.contam_rate', sep=''))
assn <- read.table(paste(basename, '.assignments', sep=''))

names <- unique(assn$V2)

if (length(args) > 1 & (args[2] == 'D' | args[2] == 'd')){
    # Include specific doublets.
    names1 <- apply(assn, 1, function(x){ strsplit(x[2], '+', fixed=T)[[1]][1] })
    names2 <- apply(assn, 1, function(x){ strsplit(x[2], '+', fixed=T)[[1]][2] })
    names <- unique(c(names1, names2))
    names <- names[which(!is.na(names))]
    names <- unique(c(assn$V2, names))
    names <- names[order(names)]
} else{
    names <- c(unique(assn[which(assn$V3=="S"),]$V2, "Doublet"))
    assn[which(assn$V3=="D"),]$V2 <- "Doublet"
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
pval <- ggdraw() + 
    draw_label(paste("P(contamination) = ", format(p_contam, scientific=TRUE, digits=2), sep=''), 
           size=10, lineheight=0.9, hjust=0.6)

mu_sd <- ggdraw() + 
    draw_label(mu_sigma_text, size=10, lineheight=0.9, hjust=0.6)

png(out_png, bg='white', width=6, height=9, units='in', res=150)

cpfile <- paste(basename, '.contam_prof', sep='')
if (file.exists(cpfile)){
    has_cp_file <- 1
    cp <- read.table(cpfile)
    cp_plt <- ggplot(cp) + 
        geom_bar(aes(x=0, y=V2, fill=V1), show.legend=TRUE, stat='identity') +
        scale_fill_manual("", values=name2col) +
        theme_bw() + 
        scale_y_continuous("Proportion in ambient RNA pool") +
        theme(axis.text.x=element_blank(), 
              axis.ticks.x=element_blank(), 
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=10),
              axis.text.y=element_text(size=9), 
              legend.position='right', 
              legend.direction='vertical',
              legend.text=element_text(size=8),
              panel.border=element_blank(),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.background=element_blank())
    
    row1col1 <- plot_grid(title, mu_sd, cp_plt, ncol=1, 
        rel_heights=c(0.2*0.5, 0.2*0.5, 0.8))
    row1col2 <- plot_grid(dens, byllr, ncol=1, rel_heights=c(0.45, 0.55))
    row1 <- plot_grid(row1col1, row1col2, ncol=2, rel_widths=c(0.45, 0.55))
    row2 <- plot_grid(bars, boxes, ncol=2, rel_widths=c(0.33,0.66)) 
    plot_grid(row1, row2, nrow=2, rel_heights=c(0.45, 0.55))

} else{
    row1 <- plot_grid(title, mu_sd, ncol=2, rel_widths=c(0.5,0.5))
    row2 <- plot_grid(dens, byllr, ncol=2, rel_widths=c(0.45, 0.55))
    row3 <- plot_grid(bars, boxes, ncol=2, rel_widths=c(0.33, 0.66))
    plot_grid(row1, row2, row3, ncol=1, rel_heights=c(0.05, 0.275, 0.675))
}

dev.off()

