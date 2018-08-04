## Appendix A - conversion of Tasmanian devil scaffold labels ##

## reading in the conversion file ##
conversion <- read.csv(file = "Scripts/Data/Devil_7.0_Scaffolds.txt",skip = 32,sep = "\t",header=T,stringsAsFactors = F)
# reducing down to just the 2 naming columns needed
name_cols <- data.frame("Sequence.Name"=conversion[,1],"UCSC.style.name"=conversion[,10],stringsAsFactors = F)

## converting bed file sh1_ensembl.bed ##
sh_bed <- read.csv("SH ensembl annot UCSC/sh1_ensembl.bed",sep="",header=F,stringsAsFactors = F)
index <- unique(sh_bed[,1]) # indexing the unique chr names that need replacing
# reducing the conversion table to just those in the file
indexed <- data.frame("new"=name_cols[which(name_cols[,2] %in% index),1],"UCSC"=name_cols[which(name_cols[,2] %in% index),2])
a <- indexed$UCSC; b<-indexed$new

X1 <- sh_bed[order(sh_bed$V1),] # sorting the bed file using order
relabel <- X1[,1] # taking out just the name columns to use gsub on
relabelb <- X1[,4]
for (i in 1:nrow(indexed)){
  print(i)
  relabel <- gsub(a[i],b[i],relabel)
  relabelb <- gsub(a[i],b[i],relabelb)
}
result <- data.frame("new"=relabel,X1[1:3],"new2"=relabelb,X1[4:ncol(X1)],stringsAsFactors = F) # add new names in
reordered <- result[order(as.numeric(row.names(result))),] # reorder this using rownames
final <- reordered[,-c(2,6)] # and remove the old name columns to get back to bed format
# write to file:
write.table(final,"SH ensembl annot UCSC/relabeled/sh1_ensembl.bed",quote = F,sep='\t',row.names = F,col.names = F)

## similarly for sh1_ensembl.tab ##
sh_tab <- read.csv("SH ensembl annot UCSC/sh1_ensembl.tab",sep="",stringsAsFactors = F)
index <- unique(sh_tab$chrom)
indexed <- data.frame("new"=name_cols[which(name_cols[,2] %in% index),1],"UCSC"=name_cols[which(name_cols[,2] %in% index),2])
a <- indexed$UCSC; b<-indexed$new

X2 <- sh_tab[order(sh_tab$chrom),]
relabel2 <- X2$chrom
for (i in 1:nrow(indexed)){
  print(i)
  relabel2 <- gsub(a[i],b[i],relabel2)
}
result2 <- data.frame(X2[,1:2],"new"=relabel2,X2[3:ncol(X2)],stringsAsFactors = F) # add new names in
reordered2 <- result2[order(as.numeric(row.names(result2))),] # reorder using rownames
final2 <- reordered2[,-4] # and remove the old names to get back to bed format
colnames(final2)[3] <- "chrom"
# write to file:
write.table(final2,"SH ensembl annot UCSC/relabeled/sh1_ensembl.tab",quote = F,sep='\t',row.names = F,col.names = T)

## reading in sh1_ensembl.fasta using biostrings package ##
require("Biostrings")
sh_fasta <- readDNAStringSet("SH ensembl annot UCSC/sh1_ensembl.fasta")

relabel3 <- names(sh_fasta) # pulling out the names element for the fasta file to modify
for (i in 1:nrow(indexed)){
  print(i)
  relabel3 <- gsub(a[i],b[i],relabel3)
}
grep(pattern = "chr",relabel3) # checking there are no chr.. labels left.
final_fasta <- sh_fasta
names(final_fasta) <- relabel3
# checking how wide lines were in my original fasta file.
line <- "AGCCCGGGCCGGACGAAGCGGCCGAGGCCCAGGAGCTGAAGCTCCAGCTG"
nchar(line) # = 50
writeXStringSet(final_fasta,"SH ensembl annot UCSC/relabeled/sh1_ensembl.fasta",width=50) #writes in same format

## reading in sh1_ensembl.gtf to relabel ##
sh_gtf <- read.csv("SH ensembl annot UCSC/sh1_ensembl.gtf",sep="",header=F,stringsAsFactors = F)
index <- unique(sh_gtf$V1)
indexed <- data.frame("new"=name_cols[which(name_cols[,2] %in% index),1],"UCSC"=name_cols[which(name_cols[,2] %in% index),2])
a <- indexed$UCSC; b<-indexed$new

X3 <- sh_gtf[order(sh_gtf$V1),]
relabel4 <- X3$V1
for (i in 1:nrow(indexed)){
  print(i)
  relabel4 <- gsub(a[i],b[i],relabel4)
}
result4 <- data.frame("new"=relabel4,X3,stringsAsFactors = F) # add new names to gtf table
reordered4 <- result4[order(as.numeric(row.names(result4))),] # reorder this using rownames
final4 <- reordered4[,-2] # and remove the old names to get back to bed format
write.table(final4,"SH ensembl annot UCSC/relabeled/sh1_ensembl.gtf",quote = F,sep='\t',row.names = F,col.names = F)
