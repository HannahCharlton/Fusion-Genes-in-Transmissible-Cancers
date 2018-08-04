## Appendix A.2 - Filtering of contigs in the Tasmanian devil assembly ##

## File listing contigs that are retained that are long enough and contain genes
kept_contigs <- read.csv("devil_contigs_selected_by_length_or_has_gene.txt",header=F,col.names = "Contigs")

## new_circos function written that generates initial circos plot with devil chromosomes:
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

## function that filters out contigs not in the kept_contigs file, then draws the retained fusions on circos plot:
contig_filtering <- function(input,contigs = kept_contigs,contiglengths=contig_lengths,chromosome_ranges=chr_ranges){
  df<- data.frame("LGene"=input$LeftGene,"LeftContig"=input$LeftContig,"LeftBreakNT"=input$LeftBreakNT,"LChr"="",
                  "LPosition"="","LPosition2"="", "RGene"=input$RightGene, "RightContig"=input$RightContig,
                  "RightBreakNT"=input$RightBreakNT,"RChr"="", "RPosition"="","RPosition2"="","JunctionReads"=input$JunctionReads)
  contiglengths <- contig_lengths
  contiglengths$Contig <- sub('ChrX','Chrx',contiglengths$Contig)
  chromosome_ranges <- chr_ranges
  # Establishing position on chromosome based on contig and nucleotide position, in preparation for circos plotting.
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
  FusionName <- paste(df$LGene,df$RGene,sep="-")
  df <- cbind("FusionName"=FusionName,df)
  
  # Filtering out contigs that aren't in the set being kept, as listed in kept_contigs
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
  ## Plotting the filtered fusions, generating a new blank devil circos plot:
  newcircos()
  # converting to ChrX labelling as that's how the circos plot is set up..
  df$LChr <- sub('Chrx','ChrX',df$LChr)
  df$RChr <- sub('Chrx','ChrX',df$RChr)
  # putting in the genomic links between kept fusion pairs
  circos.genomicLink(region1 = data.frame(df$LChr,df$LPosition,df$LPosition2), 
                     region2=data.frame(df$RChr,df$RPosition,df$RPosition2)) 
  cat(nrow(df),"fusions remain after filtering \n")
  return(df)
}
