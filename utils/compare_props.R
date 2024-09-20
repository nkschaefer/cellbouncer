#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
    write("USAGE: compare_props.R [ref] [test]", stderr())
    write("  Compares two sets of mixture proportions (from bulkprops or", stderr())
    write("  a contamination profile from quant_contam) and tests for", stderr())
    write("  significant differences. The first one is used as the null", stderr())
    write("  hypothesis and must have had bootstrapping performed (it must", stderr())
    write("  have a third column containing Dirichlet concentration", stderr())
    write("  parameters.", stderr())
    write("  You may also pass in an .assignments file as [test], provied", stderr())
    write("  it contains the same set of individuals as [ref].", stderr())
    write("", stderr())
    write("Output: two-tailed p-values (and their corresponding FDRs after", stderr())
    write("  correcting for multiple hypotheses) for comparing the proportion", stderr())
    write("  of each individual.", stderr())
    q()
}

cellcounts2props <- function(name){
    dat <- read.table(name)
    s <- dat[which(dat$V3=="S"),]
    d <- dat[which(dat$V3=="D"),]
    d$id1 <- apply(d, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][1]})
    d$id2 <- apply(d, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][2]})
    s <- s[,c(1,2)]
    d1 <- d[,c(1,5)]
    d2 <- d[,c(1,6)]
    colnames(s) <- c("bc", "id")
    colnames(d1) <- c("bc", "id")
    colnames(d2) <- c("bc", "id")
    s$num <- 2
    d1$num <- 1
    d2$num <- 1
    all <- rbind(s, d1, d2)
    agg <- aggregate(all$num, by=list(id=all$id), FUN=sum)
    agg$x <- agg$x / sum(agg$x)
    colnames(agg)[2] <- "prop"
    return(agg)
}

ref <- read.table(args[1])

alt <- NA
# Check whether second file passed was .assignments; if so, try to 
# infer proportions from cell counts
if (nchar(args[2]) > 12 & substr(args[2], nchar(args[2])-11, nchar(args[2])) == ".assignments"){
    alt <- cellcounts2props(args[2])
    colnames(alt) <- c("V1", "V2")
} else{
    alt <- read.table(args[2])
}

if (nchar(args[1]) > 12 & substr(args[1], nchar(args[1])-11, nchar(args[1])) == ".assignments"){
    write("ERROR: ref file cannot be .assignments", stderr())
    q()
}

if (length(colnames(ref)) != 3){
    write("ERROR: reference file (#1) must contain a third column with Dirichlet", stderr())
    write("  concentration parameters. Did you run bootstrapping when you created", stderr())
    write("  the file?", stderr())
    q()
}

ref <- ref[order(ref$V1),]
alt <- alt[order(alt$V1),]

if ( length(ref$V1) != length(alt$V1) || ! all(ref$V1 == alt$V1) ){
    write("ERROR: files contain different individuals. Cannot compare.", stderr())
    q()
}

# Ignore test concentration parameters, if given
alt <- alt[,c(1,2)]
colnames(alt) <- c("ID", "test")
alt$ref <- ref$V3/sum(ref$V3)
alt <- alt[,c(1,3,2)]

alt$p <- pbeta(alt$test, ref$V3, sum(ref$V3)-ref$V3)
alt[which(alt$p > 0.5),]$p <- 1.0 - alt[which(alt$p > 0.5),]$p
alt$p <- alt$p * 2
alt$p.adj <- p.adjust(alt$p, method='fdr')

write.table(format(alt, digits=4), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

