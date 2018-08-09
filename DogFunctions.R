## Functions used in Dog analysis ##

# Dog biomart data: gene IDs, transcript IDs for conversions
dog_IDs <- read.csv(file="biomart_dogdata.txt", sep=",",stringsAsFactors = F)

# Getting data from STAR output
STAR_PEout <- function(input){
  simplified <- cbind(input[,5:6],"LeftContig"=input[,6],"LeftBreakNT"=input[,6],"LeftStrand"=input[,6],
                      input[,7:8],"RightContig"=input[,8],"RightBreakNT"=input[,8],"RightStrand"=input[,8]
                      ,"JunctionReads"=input[,3],"SplitReads"=input[,2],"InFrame"=input[,20],"Tool"=rep("S",nrow(input)))
  simplified$LeftContig <- gsub(":.*","",simplified$LeftContig)
  simplified$LeftBreakNT <- gsub(".*:([^:]+):.*","\\1",simplified$LeftBreakNT)
  simplified$LeftBreakNT <- as.numeric(simplified$LeftBreakNT)
  simplified$LeftStrand <- gsub(".*:","",simplified$LeftStrand)
  simplified$RightContig <- gsub(":.*","",simplified$RightContig)
  simplified$RightBreakNT <- gsub(".*:([^:]+):.*","\\1",simplified$RightBreakNT)
  simplified$RightBreakNT <- as.numeric(simplified$RightBreakNT)
  simplified$RightStrand <- gsub(".*:","",simplified$RightStrand)
  simplified$InFrame <- gsub("INFRAME","YES",simplified$InFrame)
  simplified$InFrame <- gsub("\\.","NO",simplified$InFrame)
  simplified$InFrame <- gsub("FRAMESHIFT","NO",simplified$InFrame)
  simplified <- simplified[,-c(2,7)]    
  return(simplified)
}

# Converting STAR data Gene names
dog_STAR_conversion <- function(input){
    converting_STAR <- data.frame(input,stringsAsFactors = F)
    converting_STAR$LeftGene <- gsub("\\^.*","",converting_STAR$LeftGene)
    converting_STAR$RightGene <- gsub("\\^.*","",converting_STAR$RightGene)
    FusionPair <- paste(converting_STAR$LeftGene,converting_STAR$RightGene,sep="-")
    converting_STAR<- cbind("FusionName"=FusionPair,converting_STAR)
    return(converting_STAR)
}

# Getting data from JAFFA output
JAFFA_PEout <- function(input){
  simplified <- cbind(input[,2:5],input[,2],input[,6:8],input[,10:12],rep("J",nrow(input)))
  colnames(simplified) <- c("LeftGene","LeftContig","LeftBreakNT","LeftStrand","RightGene","RightContig",
                            "RightBreakNT","RightStrand","JunctionReads","SplitReads","InFrame","Tool")
  simplified$LeftGene <- gsub(":.*","",simplified$LeftGene)
  simplified$RightGene <- gsub(".*:","",simplified$RightGene)
  simplified$LeftGene <- gsub("[.].*","",simplified$LeftGene)
  simplified$RightGene <- gsub("[.].*","",simplified$RightGene)
  simplified$InFrame <- gsub("TRUE","YES",simplified$InFrame)
  simplified$InFrame <- gsub("FALSE","NO",simplified$InFrame)
  FusionPair <- paste(simplified$LeftGene,simplified$RightGene,sep="-")
  simplified <- cbind("FusionName"=FusionPair,simplified)
  return(simplified)
}

# Importing and converting Integrate output to same format as above, Integrate provides a summary.tsv 
# file with fusion names and breakpoints in a separate breakpoints.tsv file - combining these.
dog_INT_conversion <- function(input,bps,IDs=dog_IDs){
  integrate_data <- input[c(2,3,7,8)]
  colnames(integrate_data) <- c("LeftGene","RightGene","SpanningReads","SplitReads")
  integrate_data$LeftGene <- gsub("Transcript_","",integrate_data$LeftGene)
  integrate_data$RightGene <- gsub("Transcript_","",integrate_data$RightGene)
  integrate_data$LeftGene <- gsub("/.*","",integrate_data$LeftGene)
  integrate_data$RightGene <- gsub("/.*","",integrate_data$RightGene)
  ## Put in the chromosomes and breakpoints from the bpINt file
  blank <- rep("",nrow(integrate_data)) 
  integrate_data <- cbind(integrate_data[1],bps[c(3,4)],blank,integrate_data[2],bps[c(6,7)],blank,integrate_data[c(3,4)])
  colnames(integrate_data) <- c("LeftGene","LeftContig","LeftBreakNT","LeftStrand","RightGene","RightContig","RightBreakNT","RightStrand","JunctionReads","SplitReads")
  # Index of trans ID's that need converting.
  unique_IDs <- unique(which(IDs$Transcript.stable.ID %in% integrate_data$LeftGene | IDs$Transcript.stable.ID %in% integrate_data$RightGene))
  index <- data.frame("Genes"=IDs[unique_IDs,1],"Transcripts"=IDs[unique_IDs,2],stringsAsFactors = F)
  # And now get the corresponding Gene IDs for the Trans IDs from STAR.
  Lrelabel <- integrate_data$LeftGene
  Rrelabel <- integrate_data$RightGene
  a <- index$Transcripts; b <- index$Genes
  for (i in 1:nrow(index)){
    #print(i)
    Lrelabel <- gsub(a[i],b[i],Lrelabel)
    Rrelabel <- gsub(a[i],b[i],Rrelabel)
  }
  integrate_data$LeftGene <- Lrelabel; integrate_data$RightGene <- Rrelabel
  FusionPair <- paste(integrate_data$LeftGene,integrate_data$RightGene,sep="-")
  integrate_data <- cbind("FusionName"=FusionPair,integrate_data)
  integrate_data <- cbind(integrate_data,"InFrame"=rep("NA",nrow(integrate_data)),"Tool"=rep("I",nrow(integrate_data)))
  return(integrate_data)
}

# Removing duplicate calls:
remove_dogdupes <- function(data){
  # sorting by descending spanning reads so that the highest one is kept out of the duplicates:
  data$JunctionReads <- as.numeric(gsub("-","0",data$JunctionReads))
  data$SplitReads <- as.numeric(gsub("-","0",data$SplitReads))
  data <- data[order((data$JunctionReads+data$SplitReads),decreasing = T),]
  # and now removing any duplicates that occur in the FusionName column.
  cat(sum(duplicated(data$FusionName)),"duplicates removed \n")
  data <- data[!duplicated(data$FusionName), ]
  return(data)
}

# Venn diagram with 3 sets:
threesetvenn <- function(set1,set2,set3,name1,name2,name3){
  pairs_1 <- set1$FusionName; pairs_2 <- set2$FusionName; pairs_3 <- set3$FusionName
  venn.list <- list(pairs_1,pairs_2,pairs_3); names(venn.list)<-c(paste(name1),paste(name2),paste(name3))
  venn.plot <- venn.diagram(venn.list, NULL, fill=c(rainbow(3)), alpha=c(0.5,0.5,0.5), cex = 3, cat.fontface=3,margin=0.1,
                            sub.pos = c(0.5,0.85),sub.cex = 2,sub.fontface = 2, ext.text=F,force.unique = F)
  grid.newpage(); grid.draw(venn.plot)
}

# Combining the data from the three callers and filtering out duplicates
PE_combine_and_filter <- function(JAF,STAR,INT){
  data <- rbind(JAF,STAR,INT)
  # sorting by descending spanning reads so that the highest one is kept out of the duplicates:
  data <- data[order((data$JunctionReads+data$SplitReads),decreasing = T),]
  # and now removing any duplicates that occur in the FusionName column.
  cat(sum(duplicated(data$FusionName)),"duplicates removed \n")
  data <- data[!duplicated(data$FusionName), ]
  return(data)
}

# Remove fusions called in normal tissues, same function as in Devils
remove_normals <- function(input,normal1,normal2){
  tumour_only <- input
  for (f in 1:nrow(input)){
    if (input$FusionName[f] %in% normal1$FusionName | input$FusionName[f] %in% normal2$FusionName){
      cat("Fusion", input$FusionName[f],"occurs in normal sample \n")
      tumour_only[f,] <- 0
    }
  }
  cat("Removed",sum(tumour_only[,1]==0),"fusions from",deparse(substitute(input)),"file. \n")
  tumour_only <- tumour_only[-which(tumour_only[,1]==0),] 
  return(tumour_only)
}

# Initialise circos plot with dog chromosomes
newdogcircos <- function(){
    dog_chranges <- read.csv(file = "canFam3.chrom.sizes",header = F,sep="\t",stringsAsFactors = F)
    dog.chromosome.ranges <- dog_chranges[1:39,] # chromosomes 1-38 + X
    dog.chromosome.ranges[,3] <- gsub("chr","",dog.chromosome.ranges[,1])
    dog.chromosome.ranges <- dog.chromosome.ranges[order(as.numeric(dog.chromosome.ranges$V3)),]
    dog.chromosome.ranges$V1 <- sub('chr','Chr',dog.chromosome.ranges$V1)
    rownames(dog.chromosome.ranges) <- dog.chromosome.ranges$V1
    dog.chromosome.ranges[,1] <- 1
    circos.par("track.height"=0.3, cell.padding=c(0,0,0,0))
    circos.initialize(factors = rownames(dog.chromosome.ranges), xlim=dog.chromosome.ranges[,1:2])
    circos.trackPlotRegion(ylim = c(0, 1), 
                           #panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))}, 
                           panel.fun = function(x, y) {
                             print(get.cell.meta.data("xlim"))
                             chr = CELL_META$sector.index
                             xlim = CELL_META$xlim
                             ylim = CELL_META$ylim
                             circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "black",
                                         facing = "downward", niceFacing = TRUE)
                           }, 
                           track.height=0.08, bg.col=rainbow(39,alpha = 0.5), bg.border=rainbow(39),
                           track.index=1)
}
