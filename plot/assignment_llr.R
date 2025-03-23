#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("USAGE:", stderr())
    write("assignment_llr.R <output_prefix>", stderr())
    write("Given output from any program that makes an .assignments file", stderr())
    write("    creates a plot that shows the number of cells with each type", stderr())
    write("    of assignment that would pass a variety of log likelihood ratio", stderr())
    write("    cutoffs. If the proportions drastically change after a particular", stderr())
    write("    (low) cutoff and stabilize after, then it may be a reasonable cutoff", stderr())
    write("    to use.", stderr())
    write("Optional second argument D keeps all doublet identities in the plot.", stderr())
    write("    Default behavior is to label all unique doublet identities \"Multiplet.\"", stderr())
    write("Creates PDF and PNG plots named <output_prefix>.llr.pdf and", stderr())
    write("    <output_prefix>.llr.png.", stderr())
    q()
}

keep_doub <- 0
if (length(args) > 1 & args[2] == 'D'){
    keep_doub <- 1
}

assn <- read.table(paste(args[1], '.assignments', sep=''))
indv_label <- "Individual"

assn <- assn[which(assn$V4 > 0),]

# Preserve droplet type info, in case we modify it
assncpy <- assn
drop_type_label <- "Droplet_type"
if (sapply(assn, class)[3] == "integer"){
    assn$type <- "S"
    assn[which(assn$V3 >= 2),]$type <- "D"
    assn <- assn[,c(1,2,5,4)]
    colnames(assn) <- c("V1", "V2", "V3", "V4")
    drop_type_label <- "num_sgRNAs"
    indv_label <- "Assignment"

    if (max(assncpy$V3) > 5){
        assncpy$num <- as.character(assncpy$V3)
        assncpy[which(assncpy$V3 >= 5),]$num <- "5+"
        assncpy <- assncpy[,c(1,2,5,4)]
        colnames(assncpy) <- c("V1", "V2", "V3", "V4")
    }
}
if (! keep_doub){
    if (length(rownames(assn[which(assn$V3 != "S"),])) > 0){
        assn[which(assn$V3 != "S"),]$V2 <- "Multiplet"
    }
}

df <- data.frame(Var1=c(), Freq=c(), type=c(), quant=c())
#quants=quantile(assn$V4, seq(0,1,0.01))
minv <- min(log10(assn$V4))
maxv <- max(log10(assn$V4))
step <- (maxv-minv)/100
quants <- seq(minv, maxv, step)
for (quant in quants){
    quant <- 10^quant
    if (length(rownames(assn[which(assn$V4 > quant),])) > 0){
        tab <- table(assncpy[which(assncpy$V4 > quant),]$V3)
        #dftab <- rbind(as.data.frame(tab), data.frame(Var1=c("Total"), Freq=c(sum(tab))))
        dftab <- as.data.frame(tab)
        dftab$type <- drop_type_label
        dftab$quant <- quant
        df <- rbind(df, dftab)
        tab2 <- table(assn[which(assn$V4 > quant),]$V2)
        #tab2df <- as.data.frame(tab2)
        #dftab <- rbind(tab2df, data.frame(Var1=c("Total"), Freq=c(sum(tab2))))
        dftab <- as.data.frame(tab2)
        dftab$type <- indv_label
        dftab$quant <- quant
        df <- rbind(df, dftab)
    }
}

names <- unique(df$Var1)
names <- factor(names, labels=names, levels=names)

ncol <- length(names)
if (ncol <= 20){
    pal <- pal_d3("category20")(ncol)
} else{
    pal <- colorRampPalette(pal_d3("category20")(20))(ncol)
}
name2col <- setNames(pal, names)

if (keep_doub){
    df <- df[which(df$type=="Individual"),]
}

plt <- ggplot(df) +
    #geom_point(aes(x=quant, y=Freq, colour=Var1)) +
    #geom_line(aes(x=quant, y=Freq, colour=Var1)) +
    geom_bar(aes(x=quant, y=Freq, fill=Var1), stat='identity') +
    theme_bw() +
    scale_x_log10("Log likelihood ratio cutoff") +
    scale_y_continuous("Number of cells") +
    annotation_logticks(side="b") +
    #scale_fill_discrete("Category") + 
    scale_fill_manual("Category", values=name2col) +
    theme(axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          strip.background=element_rect(colour="black", fill="#FFFFFF", linewidth=0),
          strip.text.y=element_text(size=14, face="bold"),
          legend.position="top",
          legend.title=element_blank())

if (!keep_doub){
    plt <- plt + facet_grid(type~.)
}

ggsave(plt, file=paste(args[1], '.llr.pdf', sep=''), width=9, height=7, bg='white')
ggsave(plt, file=paste(args[1], '.llr.png', sep=''), width=9, height=7, units='in', bg='white')
