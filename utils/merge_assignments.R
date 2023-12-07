#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
    write("USAGE: merge_ids.R file1.assignments file2.assignments", stderr())
    write("Reconciles cell-individual identifications from two sources, for example", stderr())
    write("    if you have run demux_mt using the same haplotypes on both ATAC-seq", stderr())
    write("    and RNA-seq from the same cells (i.e. multiome data). Alernatively, if", stderr())
    write("    you have run demux_mt and demux_vcf on the same cells, with data", stderr())
    write("    representing the same individuals (such as genomic variants and mitochondrial", stderr())
    write("    haplotypes tied to the same individuals.", stderr())
    write("Writes (to stdout) a single assignment file. For any cell missing in one", stderr())
    write("    file, it receives the information from the file where it is present. If the", stderr())
    write("    cell is present in both files, it receives the ID with higher LLR. If the cell", stderr())
    write("    is present in both files with equal LLR, it receives an ID only if the two", stderr())
    write("    files agree on its identity.", stderr())    
    q()
}

f1 <- read.table(args[1], sep="\t")
f2 <- read.table(args[2], sep="\t")

colnames(f1) <- c("bc", "id1", "type1", "llr1")
colnames(f2) <- c("bc", "id2", "type2", "llr2")

merged <- merge(f1, f2, all.x=TRUE, all.y=TRUE)

merged$id <- ""
merged$type <- ""
merged$llr <- 0.0

if (length(rownames(merged[which(is.na(merged$id1)),])) > 0){
    merged[which(is.na(merged$id1)),]$id <- merged[which(is.na(merged$id1)),]$id2
    merged[which(is.na(merged$id1)),]$type <- merged[which(is.na(merged$id1)),]$type2
    merged[which(is.na(merged$id1)),]$llr <- merged[which(is.na(merged$id1)),]$llr2
}
if (length(rownames(merged[which(is.na(merged$id2)),])) > 0){
    merged[which(is.na(merged$id2)),]$id <- merged[which(is.na(merged$id2)),]$id1
    merged[which(is.na(merged$id2)),]$type <- merged[which(is.na(merged$id2)),]$type1
    merged[which(is.na(merged$id2)),]$llr <- merged[which(is.na(merged$id2)),]$llr1
}

done <- merged[which(merged$llr > 0),]

notdone <- merged[which(!is.na(merged$id1) & !is.na(merged$id2)),]

if (length(rownames(notdone[which(notdone$llr1 > notdone$llr2),])) > 0){
    notdone[which(notdone$llr1 > notdone$llr2),]$id <- notdone[which(notdone$llr1 > notdone$llr2),]$id1
    notdone[which(notdone$llr1 > notdone$llr2),]$type <- notdone[which(notdone$llr1 > notdone$llr2),]$type1
    notdone[which(notdone$llr1 > notdone$llr2),]$llr <- notdone[which(notdone$llr1 > notdone$llr2),]$llr1
}

if (length(rownames(notdone[which(notdone$llr2 > notdone$llr1),])) > 0){
    notdone[which(notdone$llr2 > notdone$llr1),]$id <- notdone[which(notdone$llr2 > notdone$llr1),]$id2
    notdone[which(notdone$llr2 > notdone$llr1),]$type <- notdone[which(notdone$llr2 > notdone$llr1),]$type2
    notdone[which(notdone$llr2 > notdone$llr1),]$llr <- notdone[which(notdone$llr2 > notdone$llr1),]$llr2
}

if (length(rownames(notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),])) > 0){
    notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$id <- 
        notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$id1
    notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$type <- 
        notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$type1
    notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$llr <- 
        notdone[which(notdone$llr1 == notdone$llr2 & notdone$id1==notdone$id2),]$llr1
}

done <- rbind(done, notdone[which(notdone$llr > 0),])

done <- done[,which(colnames(done) %in% c("bc", "id", "type", "llr"))]

write.table(done, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

