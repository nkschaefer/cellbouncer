#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

checks <- c(length(args) < 1, !file.exists(paste(args[1], '.dd.all', sep='')),
    !file.exists(paste(args[1], '.dd.indv', sep='')))

if (sum(checks) > 0){
   write("USAGE:", stderr())
   write("doublet_dragon.R [output_prefix] (D)", stderr())
   write("[output_prefix] = output prefix from doublet_dragon run", stderr())
   write("D = optional second argument to include all doublet combinations in plot", stderr())
   write("After you have run doublet_dragon on a data set, you should", stderr())
   write("    see files called [output_prefix].dd.all and [output_prefix].dd.indv", stderr())
   write("    where [output_prefix] is the first argument given to doublet_dragon.", stderr())
   write("This program loads those files and plots summary information about the", stderr())
   write("    doublet_dragon results to a file called [output_prefix].dd.png (raster)", stderr())
   write("    and [output_prefix].dd.pdf (vector)", stderr())
   q()
}

all <- read.table(paste(args[1], '.dd.all', sep=''))
indv <- read.table(paste(args[1], '.dd.indv', sep=''))
colnames(indv) <- c("dataset", "type", "indv", "obs", "exp", "binom.p")

datasets <- all[which(all$V1=="data_set"),]$V2

dataset_pal <- pal_locuszoom("default")(length(datasets))
dataset_name2col <- setNames(dataset_pal, datasets)

doublet <- 0
if (length(args) > 1 & args[2] == "D"){
    doublet <- 1
} else if (length(args) > 1 & args[2] != "D"){
    write("WARNING: unrecognized second argument: only valid option is D", stderr())
}

indv <- indv[order(indv$indv),]
names <- unique(indv$indv)

singlets <- unique(indv[which(indv$type=="S"),]$indv)
singlets_pal <- NA
if (length(singlets) <= 20){
    singlets_pal <- pal_d3("category20c")(length(singlets))
} else{
    singlets_pal <- colorRampPalette(pal_d3("category20c")(20))(length(singlets))
}
singlets_name2col <- setNames(singlets_pal, singlets)

height=7

if (doublet == 0){
    names <- c(unique(indv[which(indv$type=="S"),]$indv), "Doublet")
    
    indv_obs <- aggregate(indv[which(indv$type=="D"),]$obs, 
        by=list(dataset=indv[which(indv$type=="D"),]$dataset), 
        FUN=sum)
    indv_exp <- aggregate(indv[which(indv$type=="D"),]$exp, 
        by=list(dataset=indv[which(indv$type=="D"),]$dataset), 
        FUN=sum)
    colnames(indv_obs)[2] <- "obs"
    colnames(indv_exp)[2] <- "exp"
    merged <- merge(indv_obs, indv_exp)
    merged$type <- "D"
    merged$indv <- "Doublet"
    merged$binom.p <- -1
    merged <- merged[,c(1,4,5,2,3,6)]
    indv <- rbind(indv[which(indv$type=="S"),], merged)
} else{
    height <- 14
}

indv$indv <- factor(indv$indv, labels=names, levels=names)

maxdiff <- max(abs(indv$obs-indv$exp))

indv_bars <- ggplot(indv) + 
    geom_bar(aes(x=indv, y=obs-exp, fill=dataset), 
        stat='identity', show.legend=FALSE) + 
    coord_flip() + 
    theme_bw() + 
    scale_x_discrete("Individual") +
    scale_y_continuous("Difference from expected",
        limits=c(-maxdiff, maxdiff)) + 
    scale_fill_manual("Data set", values=dataset_name2col) +
    theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_text(size=14),
        axis.title.y=element_blank(),
        legend.title=element_text(size=14))

if (doublet==1){
    indv_bars <- indv_bars + facet_grid(type~., scales='free_y')
    indv_bars <- indv_bars + theme(strip.background=element_rect(
        colour="black", fill="#FFFFFF", linewidth=0),
        strip.text.y=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=9))
}

group_comp <- ggplot(all[grep("group_", all$V1),]) + 
    geom_bar(aes(x=V1, y=V3, fill=V2), stat='identity') + 
    theme_bw() + 
    coord_flip() +
    scale_fill_manual("", values=singlets_name2col) +
    ggtitle("Group composition") +
    theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title=element_blank(),
        legend.position='top')

ll <- all[which(all$V1=="data_set"),]
ll <- ll[order(ll$V3),]
ll$V2 <- factor(ll$V2, labels=ll$V2, levels=ll$V2)

ll_bar <- ggplot(ll) + 
    geom_bar(aes(x=V2, fill=V2, y=V3), stat='identity') + 
    theme_bw() + 
    ggtitle("Log likelihood of observed counts") + 
    scale_fill_manual("Data set", values=dataset_name2col) +
    theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=12))


# Assemble the plot
title1 <- ggdraw() + 
    draw_label("Doublet Dragon", size=12, fontface='bold',
        lineheight=0.9, hjust=0.6) + 
        theme(plot.margin=margin(0,0,0,0))
        
title <- ggdraw() + 
    draw_label(args[1], size=12, fontface='bold', 
        lineheight=0.9, hjust=0.6) + 
        theme(plot.margin=margin(0,0,0,0))
        
doubrate_val <- all[which(all$V2=="doublet_rate"),]$V3
doubrate <- ggdraw() + 
    draw_label(paste("Doublet rate = ", 
        format(doubrate_val,scientific=FALSE, digits=4),
        sep=""), 
        size=10, lineheight=0.9, hjust=0.6)


col1height <- 3/height
col2height <- 1-col1height

row1col1 <- plot_grid(title1, title, doubrate, NULL, 
    ncol=1, rel_heights=c(0.1,0.1,0.1,0.7))
row1 <- plot_grid(row1col1, group_comp, ncol=2, rel_widths=c(0.25, 0.75))
row2col1 <- plot_grid(ll_bar, NULL, nrow=2, rel_heights=c(3/(col2height*height), 1-3/(col2height*height)))
row2 <- plot_grid(row2col1, indv_bars, ncol=2, rel_widths=c(0.5,0.5))

out_pdf <- paste(args[1], '.dd.pdf', sep='')
out_png <- paste(args[1], '.dd.png', sep='')

png(out_png, bg='white', width=7, height=height, units='in', res=150)
plot_grid(row1, row2, nrow=2, rel_heights=c(col1height, col2height))
dev.off()

pdf(out_pdf, bg='white', width=7, height=height)
plot_grid(row1, row2, nrow=2, rel_heights=c(3/height, 1-3/height))
dev.off()

