## Full Devil Analysis 2.0 ##

# Import data, convert all to gene ID's and remove any duplicates, count number of fusions identified by each tool
devilfiles <- c("85T","86T","87T","88T","GW1","GW2","GW7","GW8","91H","pooled_normal")
for (d in devilfiles){
  assign(paste0("JAFFA_",d),read.csv(file = paste0("Results/devil/",d,"/jaffa_results.csv"),header=T,stringsAsFactors = F))
  assign(paste0("STAR_",d),read.csv(file = paste0("Results/devil/",d,"/star-fusion.fusion_predictions.abridged.coding_effect.tsv"),header=T,sep="\t",stringsAsFactors = F))
}
# excluding 53T from the analysis as it didn't have good results, barely any fusions detected etc.

devil_STARs <- list(STAR_85T,STAR_86T,STAR_87T,STAR_88T,STAR_GW1,STAR_GW2,STAR_GW7,STAR_GW8,STAR_91H,STAR_pooled_normal)
devil_JAFs <- list(JAFFA_85T,JAFFA_86T,JAFFA_87T,JAFFA_88T,JAFFA_GW1,JAFFA_GW2,JAFFA_GW7,JAFFA_GW8,JAFFA_91H,JAFFA_pooled_normal)
names(devil_STARs) <- devilfiles
names(devil_JAFs) <- devilfiles
simplified_JAFs <- lapply(devil_JAFs,function(x) {JAFFA_output(input = x)})
simplified_STARs <- lapply(devil_STARs,function(x) {STAR_output(input = x)})
converted_STARs <- lapply(simplified_STARs, function(x) {STAR_conversion(simple_input = x)})
# removing transcript columns from converted as no longer needed
subset_STARs <- lapply(converted_STARs, function(x) {subset(x$table, select = -c(LeftTranscript, RightTranscript))})

## removing duplicates from these two sets for each sample
nodupes_JAFs <- lapply(simplified_JAFs,function(x) {remove_dupes(data=x$table)})
nodupes_STARs <- lapply(subset_STARs, function(x) {remove_dupes(data = x)})

# generating numbers for table:
lapply(nodupes_JAFs, function(x) {nrow(x)})
lapply(nodupes_STARs, function(x) {nrow(x)})

# Venn diagrams showing how much the callers agree for all the samples - without removing normal fusions.
simpleplotvenn(nodupes_STARs$`85T`$FusionName,nodupes_JAFs$`85T`$FusionName,"85T STAR", "85T JAFFA")
simpleplotvenn(nodupes_STARs$`86T`$FusionName,nodupes_JAFs$`86T`$FusionName,"86T STAR", "86T JAFFA")
simpleplotvenn(nodupes_STARs$`87T`$FusionName,nodupes_JAFs$`87T`$FusionName,"87T STAR", "87T JAFFA")
simpleplotvenn(nodupes_STARs$`88T`$FusionName,nodupes_JAFs$`88T`$FusionName,"88T STAR", "88T JAFFA")
simpleplotvenn(nodupes_STARs$`GW1`$FusionName,nodupes_JAFs$`GW1`$FusionName,"GW1 STAR", "GW1 JAFFA")
simpleplotvenn(nodupes_STARs$`GW2`$FusionName,nodupes_JAFs$`GW2`$FusionName,"GW2 STAR", "GW2 JAFFA")
simpleplotvenn(nodupes_STARs$`GW7`$FusionName,nodupes_JAFs$`GW7`$FusionName,"GW7 STAR", "GW7 JAFFA")
simpleplotvenn(nodupes_STARs$`GW8`$FusionName,nodupes_JAFs$`GW8`$FusionName,"GW8 STAR", "GW8 JAFFA")
# and for normals too:
simpleplotvenn(nodupes_STARs$`91H`$FusionName,nodupes_JAFs$`91H`$FusionName,"91H STAR", "91H JAFFA")
simpleplotvenn(nodupes_STARs$pooled_normal$FusionName,nodupes_JAFs$pooled_normal$FusionName,"pooled_normal \n STAR", "pooled_normal \n JAFFA")

# combine into 1 set per sample:
all_combined <- mapply(function(X,Y) {combine_and_filter(X,Y)}, X = nodupes_JAFs, Y = nodupes_STARs,SIMPLIFY = F)

# removing normals first and then contig filtering
tumouronly <- lapply(all_combined[-c(9,10)], function(x) {remove_normals(x,all_combined$`91H`,all_combined$pooled_normal)})
tumours_filtered <- lapply(tumouronly, function(x) {contig_filtering(input = x)})

lapply(all_combined, function(x) {nrow(x)})
lapply(tumours_filtered, function(x) {nrow(x)})

# Generating a dataframe with each column containing the fusion pairs identified in a sample
Data <- data.frame(tumours_filtered$GW1$FusionName)
empty <- as.data.frame(matrix("",nrow=nrow(Data),ncol = 4))
extracol <- as.data.frame(matrix("",nrow=nrow(Data),ncol = 1))
Data <- cbind(empty,Data,empty,extracol)
colnames(Data) <- devilfiles
for (n in c(1,2,3,4,6,7,8)){
  print(n)
  Data[,n]<- tumours_filtered[[n]]$FusionName[seq(Data$GW1)]
}
for (n in c(9,10)){
  print(n)
  Data[,n]<- all_filtered[[n]]$FusionName[seq(Data$GW1)]
}
Data[is.na(Data)] <- ""
# Saving this as a csv for future reference:
write.table(Data,file = "Devils_TumourOnly_Venn2.csv",sep=",",row.names = F,quote = F)#,col.names = F)

# Consensus sets for DFT1, DFT2 and both DFTs
part1 <- intersect(Data$`87T`,Data$`85T`)
part2 <- intersect(part1, Data$`86T`)
part3 <- intersect(part2, Data$`88T`)
part4 <- intersect(part3, Data$GW1)
dft1 <- intersect(part4, Data$GW2)
dft2 <- intersect(Data$GW8, Data$GW7)
dft2 <- dft2[-which(dft2=="")] # removing the "" as both sets have spaces at the bottom.
all <- intersect(dft1, dft2)

library(stargazer)
alltumour_consensus <- tumours_filtered$`85T`[tumours_filtered$`85T`$FusionName %in% all,]
alltumour_renamed <- naming_genes2(input=alltumour_consensus)
stargazer(alltumour_renamed[,c(2,3,4,8,9,10,14,15,16)],summary=F)
new_labelledcircos(alltumour_renamed)

unique_dft1 <- fusions_unique_to_set(dft1,6)
dft1_consensus <- tumours_filtered$`85T`[tumours_filtered$`85T`$FusionName %in% unique_dft1,]
dft1_renamed <- naming_genes2(dft1_consensus)
stargazer(dft1_renamed[,c(2,3,4,8,9,10,14,15,16)],summary=F)
new_labelledcircos(dft1_renamed)

unique_dft2 <- fusions_unique_to_set(dft2,2)
dft2_consensus <- tumours_filtered$GW7[tumours_filtered$GW7$FusionName %in% unique_dft2,]
dft2_renamed <- naming_genes2(dft2_consensus)
stargazer(dft2_renamed[,c(2,3,4,8,9,10,14,15,16)],summary=F)
new_labelledcircos(dft2_renamed)

## investigating whether any of the retained fusions were actually called by multiple tools:
combine_no_filter <- function(input1,input2){
  data <- rbind(input1,input2)
  # sorting by descending spanning reads so that the highest one is kept out of the duplicates:
  data <- data[order(data$SplitReads,decreasing = T),]
  return(data)
}
all_combined_nofilter <- mapply(function(X,Y) {combine_no_filter(X,Y)}, X = nodupes_JAFs, Y = nodupes_STARs,SIMPLIFY = F)
tumouronly_nofilter <- lapply(all_combined_nofilter[-c(9,10)], function(x) {remove_normals(x,all_combined_nofilter$`91H`,all_combined$pooled_normal)})
Data <- data.frame(tumouronly_nofilter$GW1$FusionName)
empty <- as.data.frame(matrix("",nrow=nrow(Data),ncol = 4))
extracol <- as.data.frame(matrix("",nrow=nrow(Data),ncol = 1))
Data <- cbind(empty,Data,empty,extracol)
colnames(Data) <- devilfiles
for (n in c(1,2,3,4,6,7,8)){
  print(n)
  Data[,n]<- tumouronly_nofilter[[n]]$FusionName[seq(Data$GW1)]
}
for (n in c(9,10)){
  print(n)
  Data[,n]<- all_combined_nofilter[[n]]$FusionName[seq(Data$GW1)]
}
Data[is.na(Data)] <- ""
alltumour_consensus <- tumouronly_nofilter$`85T`[tumouronly_nofilter$`85T`$FusionName %in% all,]
unique_dft1 <- fusions_unique_to_set(dft1,6)
dft1_consensus <- tumouronly_nofilter$`85T`[tumouronly_nofilter$`85T`$FusionName %in% unique_dft1,]
unique_dft2 <- fusions_unique_to_set(dft2,2)
dft2_consensus <- tumouronly_nofilter$GW7[tumouronly_nofilter$GW7$FusionName %in% unique_dft2,]

## matching structural variants to fusions identified in (subsets of) DFT samples: 
SVs_consensus <- SV_match(alltumour_consensus)
SVs_dft1 <- SV_match(dft1_renamed)
SVs_dft2 <- SV_match(dft2_renamed)
SVs_devils <- lapply(all_filtered, function(x) {SV_match(x)})

## generating UpSet plot for devil data
require("UpSetR"); require(viridis)
data_upset2 <- lapply(tumours_filtered, function(x) {return(x$FusionName)})
upset(fromList(data_upset2),nsets = 8,order.by="freq",nintersects = 30, text.scale = 1.5,
      mainbar.y.label = "Number of Shared Fusions", sets.x.label = "Fusions Per Sample",
      main.bar.color = plasma(10,alpha=0.5)[3], sets.bar.color = plasma(10,alpha=0.5)[8],
      ## queries used to highlight intersections correspondoing to DFT1 and DFT2
      queries = list(list(query = intersects, params = list("85T", "86T", "87T","88T","GW1","GW2"), color=viridis(10,alpha=0.8)[8],active = T),
                     list(query = intersects, params = list("GW7","GW8"), color = viridis(10,alpha=0.8)[6], active = T)))

## generating correlation matrix and heatmap for devil data
length_list <- lapply(tumours_filtered, function(X) {nrow(X)})
longest_list <- length_list[[which(length_list == max(unlist(length_list)))]]
corr_matrix <- matrix(NA,ncol(df),ncol(df))
colnames(corr_matrix) <- rownames(corr_matrix) <- devilfiles[-c(9,10)]
for (a in 1:length(data_upset2)){
  for (b in 1:length(data_upset2)){
    corr_matrix[a,b] <- length(intersect(data_upset2[[a]],data_upset2[[b]]))
    corr_matrix[a,b] <- corr_matrix[a,b] / longest_list
  }
}

