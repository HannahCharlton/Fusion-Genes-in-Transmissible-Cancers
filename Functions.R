### Functions used in analysis for Fusion Genes in Transmissible Cancers MPhil project ###
setwd("Desktop"); options(stringsAsFactors = F)

## Functions used in Devil analysis ##
# getting key data out of star output csv files.
STAR_output <- function(input){
  simplified <- cbind(input[,5:6],"LeftContig"=input[,6],"LeftBreakNT"=input[,6],"LeftStrand"=input[,6],
                      input[,7:8],"RightContig"=input[,8],"RightBreakNT"=input[,8],"RightStrand"=input[,8],
                      "SplitReads"=input[,2],"InFrame"=input[,20],"Tool"=rep(NA,nrow(input)))
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
  simplified$Tool <- rep("S",nrow(simplified))
  simplified <- simplified[,-c(2,7)]    
  return(simplified)
}

# getting key data out of jaffa output csv file
JAFFA_output <- function(input){
  simplified <- cbind(input[,2:5],input[,2],input[,6:8],input[,11],input[,12],rep(NA,nrow(input)))
  colnames(simplified) <- c("LeftGene","LeftContig","LeftBreakNT","LeftStrand","RightGene","RightContig",
                            "RightBreakNT","RightStrand","SplitReads","InFrame","Tool")
  simplified$LeftGene <- gsub(":.*","",simplified$LeftGene)
  simplified$RightGene <- gsub(".*:","",simplified$RightGene)
  simplified$LeftGene <- gsub("[.].*","",simplified$LeftGene)
  simplified$RightGene <- gsub("[.].*","",simplified$RightGene)  
  simplified$InFrame <- gsub("TRUE","YES",simplified$InFrame)
  simplified$InFrame <- gsub("FALSE","NO",simplified$InFrame)
  simplified$Tool <- rep("J",nrow(simplified))
  JAFFA_subset <- data.frame("LeftGene"=simplified$LeftGene,"RightGene"=simplified$RightGene,stringsAsFactors = F)
  JAFFA_subset$FusionPair <- paste(JAFFA_subset$LeftGene,JAFFA_subset$RightGene,sep="-")
  return(list("table"=simplified,"subset"=JAFFA_subset))
}

# replacing transcript names with gene names in STAR-Fusion output to match JAFFA.
ID_list <- read.csv(file="Scripts/Data/ensembl_IDs_biomart.txt",sep="\t",stringsAsFactors = F)
STAR_conversion <- function(simple_input){
  converting_STAR <- data.frame("LeftTranscript"=simple_input[,1],simple_input[,1:4],
                                "RightTranscript"=simple_input[,5],simple_input[,5:(ncol(simple_input))])
  converting_STAR$LeftGene <- ""
  converting_STAR$RightGene <- ""
  STAR_subset <- data.frame("LTransID"=converting_STAR$LeftTranscript,"LGeneID"=converting_STAR$LeftGene,"RTransID"=converting_STAR$RightTranscript,"RGeneID"=converting_STAR$LeftGene,stringsAsFactors = F)
  STAR_subset$LTransID <- gsub("[.].*","",STAR_subset$LTransID)
  STAR_subset$RTransID <- gsub("[.].*","",STAR_subset$RTransID)
  # Index of trans ID's that need converting.
  unique_IDs <- unique(which(ID_list$Transcript.stable.ID %in% STAR_subset$LTransID | ID_list$Transcript.stable.ID %in% STAR_subset$RTransID))
  index <- data.frame("Genes"=ID_list[unique_IDs,1],"Transcripts"=ID_list[unique_IDs,2],stringsAsFactors = F)
  # And now get the corresponding Gene IDs for the Trans IDs from STAR.
  Lrelabel <- STAR_subset$LTransID 
  Rrelabel <- STAR_subset$RTransID
  a <- index$Transcripts; b <- index$Genes
  for (i in 1:nrow(index)){
    #print(i)
    Lrelabel <- gsub(a[i],b[i],Lrelabel)
    Rrelabel <- gsub(a[i],b[i],Rrelabel)
  }
  STAR_subset$LGeneID <- Lrelabel; STAR_subset$RGeneID <- Rrelabel
  STAR_subset$FusionPair <- paste(STAR_subset$LGeneID,STAR_subset$RGeneID,sep="-")
  # returning in a df:
  converting_STAR$LeftGene <- Lrelabel
  converting_STAR$RightGene <- Rrelabel
  return(list("table"=converting_STAR,"subset"=STAR_subset))
}

# Function to remove duplicate fusion calls.
duplicate_removal <- function(data){
  for (i in 1:nrow(data)){
    for (j in 1:nrow(data)){
      if (i != j){
        if (all(data[i,-10]==data[j,-10]) == TRUE){
          cat("Fusions match in rows",i,"and",j,"\n")
          data[i,10] <- data[i,10]+data[j,10]
          data[j,] <- 0
        }
      }
    }
  }
  data <- data[-which(data[,1]==0),]
  data <- data[order(data$SplitReads,decreasing = T),]
  # also need check for other duplicates where genes/contigs same but NT position slightly altered
  if (sum(duplicated(data$FusionPair) == TRUE) > 0){
    dupes <- data[which(duplicated(data$FusionPair) == TRUE),]
    for (d in 1:nrow(dupes)){
      # consider the lower index for each duplicate as we want to keep the one with more spanning reads
      original <- min(which(data$FusionPair == dupes$FusionPair[d]))
      duplicate <- max(which(data$FusionPair == dupes$FusionPair[d]))
      data[original,10] <- data[original,10] + data[duplicate,10]
      data <- data[-duplicate,]
    }
  }
  return(data)
}

# Venn Diagram plotting functions
require("VennDiagram")
plotvenn <- function(pairs_1,pairs_2,name1,name2){
  #pairs_1 <-set1$subset$FusionPair
  #pairs_2 <- set2$subset$FusionPair
  consensus_genes1 <- pairs_1[which(pairs_1 %in% pairs_2)]
  consensus_genes2 <- pairs_2[which(pairs_2 %in% pairs_1)]
  venn.list <- list(pairs_1,pairs_2); names(venn.list)<-c(paste(name1),paste(name2))
  venn.plot <- venn.diagram(venn.list, NULL, fill=c(rainbow(2)), alpha=c(0.5,0.5), cex = 1.5, cat.fontface=3,margin=0.1
                            ,sub = "Fusion-gene Pairs Predicted",sub.pos = c(0.5,0.85),sub.cex = 1.3,sub.fontface = 2, ext.text=F,force.unique = F)
  unique_venn.plot <- venn.diagram(venn.list, NULL, fill=c(rainbow(2)), alpha=c(0.5,0.5), cex = 1.5, cat.fontface=3,margin=0.1,
                                   sub = "Unique Fusion-gene Pairs Predicted",sub.pos = c(0.5,0.85),sub.cex = 1.3,sub.fontface = 2,ext.text=F)
  grid.newpage(); grid.draw(venn.plot)
  grid.newpage(); grid.draw(unique_venn.plot)
  return(list("list1"=consensus_genes1,"list2"=consensus_genes2))
}

simpleplotvenn <- function(pairs_1,pairs_2,name1,name2){
  venn.list <- list(pairs_1,pairs_2); names(venn.list)<-c(paste(name1),paste(name2))
  venn.plot <- venn.diagram(venn.list, NULL, fill=c(rainbow(2)), alpha=c(0.5,0.5), cex = 3, cat.fontface=3,margin=0.1,
                            sub.pos = c(0.5,0.85),sub.cex = 2,sub.fontface = 2, ext.text=F,force.unique = F)
  grid.newpage(); grid.draw(venn.plot)
  overlap <- calculate.overlap(venn.list)
  return(overlap$a3)
}

# Filtering the Devil contigs as described in Appendix A2 and circos plot of remaining fusions.
contigranges <- function(contiglengths=contig_lengths){
  chr_names <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX")
  chr_ranges <- matrix(0,ncol=3,nrow=7)
  rownames(chr_ranges) <- chr_names; colnames(chr_ranges) <- c("Start","End","Contigs")
  chr_ranges <- as.data.frame(chr_ranges)
  chr_ranges$Start <- 1
  chr_ranges$Contigs <- c(398,500,416,316,217,193,60)
  for (c in 1:nrow(chr_ranges)){
    chr <- rownames(chr_ranges)[c]
    first <- contiglengths[which(contiglengths$Contig == paste0(chr,"_supercontig_000000000")),]
    startofchromosome <- which(contiglengths$Contig == paste0(chr,"_supercontig_000000000"))
    last <- contiglengths[startofchromosome+chr_ranges$Contigs[c],]
    endofcontigs <- startofchromosome+chr_ranges$Contigs[c]
    bound <- sum(contiglengths$Length[startofchromosome:endofcontigs])
    chr_ranges$End[c] <- bound
  }
  return(chr_ranges)
}
chr_ranges <- contigranges()
kept_contigs <- read.csv("devil_contigs_selected_by_length_or_has_gene.txt",header=F,col.names = "Contigs")

devilcontiglengths <- function(){
  contiglengths <- read.table(file = "Scripts/Data/devil_contig_lengths.txt",header = F,sep = "")
  colnames(contiglengths) <- c("Contig","Length")
  contiglengths$Contig <- gsub(">","",contiglengths$Contig)
  contiglengths$Contig <- sub('Chrx','ChrX',contiglengths$Contig) # rewriting it as ChrX instead of Chrx for consistency. Remove ChrU chromosomes:
  contiglengths <- contiglengths[-grep('ChrU', contiglengths$Contig),]
  return(contiglengths)
}
contig_lengths <- devilcontiglengths()

contig_filtering <- function(input,contigs = kept_contigs,contiglengths=contig_lengths,chromosome_ranges=chr_ranges){
  df<- data.frame("LGene"=input$LeftGene,"LeftContig"=input$LeftContig,"LeftBreakNT"=input$LeftBreakNT,"LChr"="","LPosition"="","LPosition2"="",
                  "RGene"=input$RightGene, "RightContig"=input$RightContig,"RightBreakNT"=input$RightBreakNT,"RChr"="", "RPosition"="","RPosition2"="",
                  "SplitReads"=input$SplitReads,"InFrame"=input$InFrame,"Tool"=input$Tool)
  contiglengths <- contig_lengths
  contiglengths$Contig <- sub('ChrX','Chrx',contiglengths$Contig)
  chromosome_ranges <- chr_ranges
  for (r in 1:nrow(df)){
    # Left side:
    df$LChr <- gsub("_.*","",df$LeftContig)
    Lstartofchromosome <- which(contiglengths$Contig == paste0(df$LChr[r],"_supercontig_000000000"))
    Lpreviouscontig <- which(contiglengths$Contig == df$LeftContig[r])-1
    Lchromosome_sofar <- sum(contiglengths$Length[Lstartofchromosome:Lpreviouscontig])
    df$LPosition[r] <- df$LPosition2[r] <- Lchromosome_sofar + df$LeftBreakNT[r]
    # Right side:
    df$RChr <- gsub("_.*","",df$RightContig)
    Rstartofchromosome <- which(contiglengths$Contig == paste0(df$RChr[r],"_supercontig_000000000"))
    Rpreviouscontig <- which(contiglengths$Contig == df$RightContig[r])-1
    Rchromosome_sofar <- sum(contiglengths$Length[Rstartofchromosome:Rpreviouscontig])
    df$RPosition[r] <- df$RPosition2[r] <- Rchromosome_sofar + df$RightBreakNT[r]
  }
  df$LPosition <- as.numeric(df$LPosition); df$LPosition2 <- as.numeric(df$LPosition2)
  df$RPosition <- as.numeric(df$RPosition); df$RPosition2 <- as.numeric(df$RPosition2)
  FusionName <- as.character(paste(df$LGene,df$RGene,sep="-"))
  df <- cbind("FusionName"=FusionName,df)
  for (i in 1:nrow(df)){
    if(df$LeftContig[i] %in% contigs$Contigs == FALSE | df$RightContig[i] %in% contigs$Contigs == FALSE){
      cat("filtering at row",i,"!")
      df[i,] <- 0
    }
  }
  remove <- which(df[,1]==0)
  if (length(remove) != 0){
    df <- df[-remove,]
  }
  ## Plotting the filtered fusions:
  df$LChr <- sub('Chrx','ChrX',df$LChr)
  df$RChr <- sub('Chrx','ChrX',df$RChr)
  # converting to ChrX labelling as that's how the circos plot is set up..
  newcircos()
  circos.genomicLink(region1 = data.frame(df$LChr,df$LPosition,df$LPosition2), 
                     region2=data.frame(df$RChr,df$RPosition,df$RPosition2)) 
  df[,5:7],region2 = df[,11:13]) #columns indexed are L and R Chr, Position and Position2
  cat(nrow(df),"fusions remain after filtering \n")
  return(df)
}

# Combine calls from the two tools and filter out any duplicates
combine_and_filter <- function(input1,input2){
  data <- rbind(input1,input2)
  # sorting by descending spanning reads so that the highest one is kept out of the duplicates:
  data <- data[order(data$SplitReads,decreasing = T),]
  # and now removing any duplicates that occur in the FusionName column.
  cat(sum(duplicated(data$FusionName)),"duplicates removed \n")
  data <- data[!duplicated(data$FusionName), ]
  return(data)
}

# Remove fusion calls that occur in normal tissues (normal1 and normal2)
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

# Converting fusions from gene ID's to gene names where available
gene_names2 <- read.csv("Devil_Genes_Max.txt",header = T,sep = "\t",stringsAsFactors = F)
naming_genes2 <- function(input,namelist=gene_names2[,c(6,7)]){
  for (lg in 1:nrow(input)){
    lgene <- input$LGene[lg]
    if (is.na(unique(gene_names2$GENE[gene_names2$GENEID==lgene])) == TRUE) {
      cat("No gene name for this gene ID \n")
    } else {
      cat("Replacing gene ID with name \n")
      input$LGene[lg] <- unique(gene_names2$GENE[gene_names2$GENEID==lgene])
    }
  }
  for (rg in 1:nrow(input)){
    rgene <- input$RGene[rg]
    if (is.na(unique(gene_names2$GENE[gene_names2$GENEID==rgene]))==TRUE) {
      cat("No gene name for this gene ID \n")
    } else {
      cat("Replacing gene ID with name \n")
      input$RGene[rg] <- unique(gene_names2$GENE[gene_names2$GENEID==rgene])
    }
  }
  return(input)
}

# Identifying which fusions in an intersect are unique to that set of samples
fusions_unique_to_set <- function(intersect,numberofsets){
  results <- data.frame(matrix(ncol = ncol(Data),nrow=(length(intersect))))
  colnames(results) <- colnames(Data)
  for (c in 1:ncol(Data)){
    #cat("Column",c)
    results[,c] <- intersect %in% Data[,c]
  }
  which(rowSums(results)==numberofsets) # should be how many sets were in the intersection
  length(which(rowSums(results)==numberofsets)) # how many there are
  unique <- intersect[which(rowSums(results)==numberofsets)]
  return(unique)
}

# Circos plot functions (for devil)
require("circlize")
# Initializing the circos plot
newcircos <- function(){
  chromosome.ranges <- matrix(0, ncol=2, nrow=7)
  chromosome.ranges[,1] <- 1
  chromosome.ranges[,2] <- c(680437123,732629474,635140686,480844255,297391021,262034688,85730105) 
  rownames(chromosome.ranges) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX")
  circos.par("track.height"=0.2, cell.padding=c(0,0,0,0))
  circos.initialize(factors = rownames(chromosome.ranges), xlim=chromosome.ranges)
  circos.trackPlotRegion(ylim = c(0, 1), 
                         panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))}, 
                         track.height=0.02, bg.col=rainbow(7,alpha = 0.5), bg.border=rainbow(7),
                         track.index=1)
  for (i in 1:length(rownames(chromosome.ranges))){
    circos.axis(h='top',sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = rownames(chromosome.ranges)[i],
                direction = "outside", major.tick.percentage = 1, labels.cex=1,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
}

# Circos plot with genomic links
circosplot <- function(simple_fusions,contiglengths=contig_lengths){
  df<- data.frame("LContig"=simple_fusions$LeftContig,"LBreak"=simple_fusions$LeftBreakNT,"LChr"="","LPosition"="","LPosition2"="",
                  "RContig"=simple_fusions$RightContig,"RBreak"=simple_fusions$RightBreakNT,"RChr"="", "RPosition"="","RPosition2"="")
  df$LContig <- sub('Chrx','ChrX',df$LContig)
  df$RContig <- sub('Chrx','ChrX',df$RContig)
  for (r in 1:nrow(df)){
    # Left side:
    df$LChr <- gsub("_.*","",df$LContig)
    Lstartofchromosome <- which(contiglengths$Contig == paste0(df$LChr[r],"_supercontig_000000000"))
    Lpreviouscontig <- which(contiglengths$Contig == df$LContig[r])-1
    Lchromosome_sofar <- sum(contiglengths$Length[Lstartofchromosome:Lpreviouscontig])
    df$LPosition[r] <- df$LPosition2[r] <- Lchromosome_sofar + df$LBreak[r]
    # Right side:
    df$RChr <- gsub("_.*","",df$RContig)
    Rstartofchromosome <- which(contiglengths$Contig == paste0(df$RChr[r],"_supercontig_000000000"))
    Rpreviouscontig <- which(contiglengths$Contig == df$RContig[r])-1
    Rchromosome_sofar <- sum(contiglengths$Length[Rstartofchromosome:Rpreviouscontig])
    df$RPosition[r] <- df$RPosition2[r] <- Rchromosome_sofar + df$RBreak[r]
  }
  df$LPosition <- as.numeric(df$LPosition); df$LPosition2 <- as.numeric(df$LPosition2)
  df$RPosition <- as.numeric(df$RPosition); df$RPosition2 <- as.numeric(df$RPosition2)
  ## Plotting the fusions:
  newcircos()
  circos.genomicLink(region1 = df[,3:5],region2 = df[,8:10])
  return(df)
}

# Circos plot with links and labels
new_labelledcircos <- function(fusion_set){
  lbed <- data.frame(fusion_set$LChr,fusion_set$LPosition,fusion_set$LPosition2,fusion_set$LGene)
  rbed <- data.frame(fusion_set$RChr,fusion_set$RPosition,fusion_set$RPosition2,fusion_set$RGene)
  bed <- data.frame("chr","pos1","pos2","label"); colnames(lbed) <- colnames(rbed) <- colnames(bed)
  bed <- rbind(lbed,rbed)
  chromosome_ranges <- matrix(0, ncol=2, nrow=7)
  chromosome_ranges[,1] <- 1
  chromosome_ranges[,2] <- c(680437123,732629474,635140686,480844255,297391021,262034688,85730105) 
  rownames(chromosome_ranges) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX")
  circos.par("track.height"=0.2, cell.padding=c(0,0,0,0))
  circos.initialize(factors = rownames(chromosome_ranges), xlim=chromosome_ranges)
  circos.genomicLabels(bed, labels.column = 4, side = "outside",cex=1)
  circos.trackPlotRegion(ylim = c(0, 0.6), 
                         panel.fun = function(x, y) {
                           print(get.cell.meta.data("xlim"))
                           chr = CELL_META$sector.index
                           xlim = CELL_META$xlim
                           ylim = CELL_META$ylim
                           circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, col = "black",
                                       facing = "inside", niceFacing = TRUE)
                         }, 
                         track.height=0.1, bg.col=rainbow(7,alpha = 0.5), bg.border=rainbow(7),
                         track.index=3)
  circos.genomicLink(region1 = data.frame(fusion_set$LChr,fusion_set$LPosition,fusion_set$LPosition2),
                     region2 = data.frame(fusion_set$RChr,fusion_set$RPosition,fusion_set$RPosition2))
}

# Matching against existing Structural Variant calls
load("Scripts/Data/DevilSVs_unfiltered.Rdata")
SV_match <- function(input_file,SVs=Devil.SVs.unfiltered,windowsize=5000,verbose=F){
  matches <- rep(0,nrow(input_file))
  output<-vector("list",nrow(input_file))
  for (i in 1:length(matches)){
    matches[i]  <- nrow(SVs[which(SVs$`CHR-1`%in% input_file$LeftContig[i] & SVs$`CHR-2` %in% input_file$RightContig[i]),])
    if (matches[i] != 0) {
      rows <- which(SVs$`CHR-1`%in% input_file$LeftContig[i] & SVs$`CHR-2` %in% input_file$RightContig[i])
      Table <- SVs[rows,]
      same <- which((input_file$LeftBreakNT[i] >= (SVs$`END-1`[rows] -windowsize)) & (input_file$LeftBreakNT[i] <= (SVs$`END-1`[rows] +windowsize)) 
                    & (input_file$RightBreakNT[i] >= (SVs$`START-2`[rows] -windowsize)) & (input_file$RightBreakNT[i] <= (SVs$`START-2`[rows] + windowsize)))
      if (length(same) == 0) {output[[i]] <- "No matching structural variants"}
      else { output[[i]] <- list("Predicted Fusion" = input_file[i,], "Structural Variant" = Table[same,])}
    }
    else output[[i]] <- "No matching structural variants"
  }
  if (verbose==T) {show(output[which(output != "No matching structural variants")])}
  cat("Number of fusions predicted =",length(matches),"\n")
  cat("Fusions with contig pairs matching SV calls =",length(which(matches ==0)),"\n")
  cat("Fusions which are within",windowsize,"of the SV call = ",length(which(output != "No matching structural variants")),"\n")
  return(output)
}
