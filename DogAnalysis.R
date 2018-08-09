## Dog analysis with STAR, Integrate and JAFFA-H 
# (waiting for Jaffa-H to complete on 2 samples: 550 and 683, make empty lists for now)

dogfiles <- c("366T","410T","423T","459T","464T","550T","645T","683T","773T1","liver","lung")
# Importing all of the data
for (d in dogfiles){
  assign(paste0("STAR_",d),read.csv(file = paste0("Results/dog/",d,"/","star-fusion.fusion_predictions.abridged.coding_effect.tsv"),header=T,sep="\t",stringsAsFactors = F))
  assign(paste0("JAFFA_",d),try(read.csv(file = paste0("Results/dog/",d,"/jaffa_results.csv"),header=T,stringsAsFactors = F)))
  assign(paste0("INT_",d),read.csv(file = paste0("Results/dog/",d,"/summary.tsv"),sep = "\t",header=T,stringsAsFactors = F))
  assign(paste0("bpINT_",d),read.csv(file = paste0("Results/dog/",d,"/breakpoints.tsv"),sep = "\t",header=T,stringsAsFactors = F))
}
JAFFA_550T <- data.frame(JAFFA_366T[0,])
JAFFA_683T <- data.frame(JAFFA_366T[0,])

# Making lists of data to use lapply on:
dog_STARs <- list(STAR_366T,STAR_410T,STAR_423T,STAR_459T,STAR_464T,STAR_550T,STAR_645T,STAR_683T,STAR_773T1,STAR_liver,STAR_lung)
dog_JAFs <- list(JAFFA_366T,JAFFA_410T,JAFFA_423T,JAFFA_459T,JAFFA_464T,JAFFA_550T,JAFFA_645T,JAFFA_683T,JAFFA_773T1,JAFFA_liver,JAFFA_lung)
dog_INTs <- list(INT_366T,INT_410T,INT_423T,INT_459T,INT_464T,INT_550T,INT_645T,INT_683T,INT_773T1,INT_liver,INT_lung)
dog_bpINTs <- list(bpINT_366T,bpINT_410T,bpINT_423T,bpINT_459T,bpINT_464T,bpINT_550T,bpINT_645T,bpINT_683T,bpINT_773T1,bpINT_liver,bpINT_lung)
names(dog_STARs) <- names(dog_JAFs) <- names(dog_INTs) <- names(dog_bpINTs) <- dogfiles

# getting all the data in the same format:
simplified_JAFs <- lapply(dog_JAFs,function(x) {JAFFA_PEout(input = x)})
simplified_STARs <- lapply(dog_STARs,function(x) {STAR_PEout(input = x)})
converted_STARs <- lapply(simplified_STARs, function(x) {dog_STAR_conversion(x)})
converted_INTs <- mapply(function(X,Y) {dog_INT_conversion(X,Y)}, X = dog_INTs, Y = dog_bpINTs,SIMPLIFY = F)

# removing any duplicate calls:
nodupes_JAFs <- lapply(simplified_JAFs,function(x) {remove_dogdupes(data = x)})
nodupes_STARs <- lapply(converted_STARs, function(x) {remove_dogdupes(data = x)})
nodupes_INTs <- lapply(converted_INTs,function(x) {remove_dogdupes(data = x)})

# Numbers for table:
lapply(nodupes_JAFs, function(x) {nrow(x)})
lapply(nodupes_STARs, function(x) {nrow(x)})
lapply(nodupes_INTs, function(x) {nrow(x)})

# Venn diagrams
threesetvenn(nodupes_INTs$`366T`,nodupes_JAFs$`366T`,nodupes_STARs$`366T`,"366T INT","366T JAF", "366T STAR")
threesetvenn(nodupes_INTs$`410`,nodupes_JAFs$`410T`,nodupes_STARs$`410T`,"410T INT","410T JAF", "410T STAR")
threesetvenn(nodupes_INTs$`423T`,nodupes_JAFs$`423T`,nodupes_STARs$`423T`,"423T INT","423T JAF", "423T STAR")
threesetvenn(nodupes_INTs$`459T`,nodupes_JAFs$`459T`,nodupes_STARs$`459T`,"459T INT","459T JAF", "459T STAR")
threesetvenn(nodupes_INTs$`464T`,nodupes_JAFs$`464T`,nodupes_STARs$`464T`,"464T INT","464T JAF", "464T STAR")
threesetvenn(nodupes_INTs$`550T`,nodupes_JAFs$`550T`,nodupes_STARs$`550T`,"550T INT","550T JAF", "550T STAR")
threesetvenn(nodupes_INTs$`645T`,nodupes_JAFs$`645T`,nodupes_STARs$`645T`,"645T INT","645T JAF", "645T STAR")
threesetvenn(nodupes_INTs$`683T`,nodupes_JAFs$`683T`,nodupes_STARs$`683T`,"683T INT","683T JAF", "683T STAR")
threesetvenn(nodupes_INTs$`773T1`,nodupes_JAFs$`773T1`,nodupes_STARs$`773T1`,"773T INT","773T JAF", "773T STAR")
threesetvenn(nodupes_INTs$`liver`,nodupes_JAFs$`liver`,nodupes_STARs$`liver`,"liver INT","liver JAF", "liver STAR")
threesetvenn(nodupes_INTs$`lung`,nodupes_JAFs$`lung`,nodupes_STARs$`lung`,"lung INT","lung JAF", "lung STAR")

# combining the 3 sets of calls for each sample:
all_three <- mapply(function(X,Y,Z) {PE_combine_and_filter(X,Y,Z)}, X = nodupes_JAFs, Y = nodupes_STARs, Z=nodupes_INTs, SIMPLIFY = F)

# removing normals to get tumour only calls
dog_tumouronly <- lapply(all_three[-c(10,11)], function(x) {remove_normals(x,all_three$lung,all_three$liver)})

dog_upset <- lapply(dog_tumouronly, function(x) {return(x$FusionName)})
# finding set of fusions that occurred in all tumour samples:
part1 <- intersect(dog_upset$`366T`,dog_upset$`410T`)
part2 <- intersect(part1, dog_upset$`423T`)
part3 <- intersect(part2, dog_upset$`459T`)
part4 <- intersect(part3, dog_upset$`464T`)
part5 <- intersect(part4, dog_upset$`550T`)
part6 <- intersect(part5, dog_upset$`645T`)
part7 <- intersect(part6, dog_upset$`683T`)
allCTVT <- intersect(part7, dog_upset$`773T1`)

# Set of CTVT consensus fusion calls - also filtered to only have potentially inframe fusions:
CTVT_consensus <- dog_tumouronly$`366T`[dog_tumouronly$`366T`$FusionName %in% allCTVT,]
inframe_CTVT <- CTVT_consensus[which(CTVT_consensus$InFrame != "NO"),]

# Dog UpSet plot: constructed with and without 550T and 683T data:
dog_upset <- lapply(dog_tumouronly, function(x) {return(x$FusionName)})
upset(fromList(dog_upset), nsets = 9,order.by="freq",nintersects = 40, text.scale = 1.5,
      mainbar.y.label = "Number of Shared Fusions", sets.x.label = "Fusions Per Sample",
      main.bar.color = magma(10,alpha=0.5)[5], sets.bar.color = viridis(10,alpha=0.5)[8])
# excluding 550T and 683T data:
upset(fromList(dog_upset),sets = dogfiles[-c(6,8,10,11)],nsets = 9,order.by="freq",nintersects = 40, text.scale = 1.5,
      mainbar.y.label = "Number of Shared Fusions", sets.x.label = "Fusions Per Sample",
      main.bar.color = magma(10,alpha=0.5)[5], sets.bar.color = viridis(10,alpha=0.5)[8],
      queries = list(list(query = intersects, params = list("464T","366T"), color = plasma(10,alpha=0.8)[4], active = T)))

# Correlation matrix for heatmap:
dlength_list <- lapply(dog_tumouronly, function(X) {nrow(X)})
dlongest_list <- dlength_list[[which(dlength_list == max(unlist(dlength_list)))]]
dcorr_matrix <- matrix(NA,length(dog_upset),length(dog_upset))
colnames(dcorr_matrix) <- rownames(dcorr_matrix) <- dogfiles[-c(10,11)]
for (a in 1:length(dog_upset)){
  for (b in 1:length(dog_upset)){
    dcorr_matrix[a,b] <- length(intersect(dog_upset[[a]],dog_upset[[b]]))
    dcorr_matrix[a,b] <- dcorr_matrix[a,b] / dlongest_list
  }
}

heatmaply code here


## Repeated dog analysis excluding JAFFA calls competely:
for (d in dogfiles){
  assign(paste0("STAR_",d),read.csv(file = paste0("Results/dog/",d,"/",d,".fusion_predictions.abridged.tsv"),header=T,sep="\t",stringsAsFactors = F))
  assign(paste0("INT_",d),read.csv(file = paste0("Results/dog/",d,"/summary.tsv"),sep = "\t",header=T,stringsAsFactors = F))
  assign(paste0("bpINT_",d),read.csv(file = paste0("Results/dog/",d,"/breakpoints.tsv"),sep = "\t",header=T,stringsAsFactors = F))
}

dog_STARs <- list(STAR_366T,STAR_410T,STAR_423T,STAR_459T,STAR_464T,STAR_550T,STAR_645T,STAR_683T,STAR_773T1,STAR_liver,STAR_lung)
dog_INTs <- list(INT_366T,INT_410T,INT_423T,INT_459T,INT_464T,INT_550T,INT_645T,INT_683T,INT_773T1,INT_liver,INT_lung)
dog_bpINTs <- list(bpINT_366T,bpINT_410T,bpINT_423T,bpINT_459T,bpINT_464T,bpINT_550T,bpINT_645T,bpINT_683T,bpINT_773T1,bpINT_liver,bpINT_lung)
names(dog_STARs) <- names(dog_INTs) <- names(dog_bpINTs) <- dogfiles

simplified_STARs <- lapply(dog_STARs,function(x) {STAR_PEout(input = x)})
converted_STARs <- lapply(simplified_STARs, function(x) {dog_STAR_conversion(x)})
converted_INTs <- mapply(function(X,Y) {dog_INT_conversion(X,Y)}, X = dog_INTs, Y = dog_bpINTs,SIMPLIFY = F)

nodupes_STARs <- lapply(converted_STARs, function(x) {remove_dogdupes(data = x)})
nodupes_INTs <- lapply(converted_INTs,function(x) {remove_dogdupes(data = x)})
lapply(nodupes_STARs, function(x) {nrow(x)})
lapply(nodupes_INTs, function(x) {nrow(x)})

simpleplotvenn(nodupes_INTs$`366T`$FusionName,nodupes_STARs$`366T`$FusionName,"366T INT", "366T STAR")
simpleplotvenn(nodupes_INTs$`410`$FusionName,nodupes_STARs$`410T`$FusionName,"410T INT", "410T STAR")
simpleplotvenn(nodupes_INTs$`423T`$FusionName,nodupes_STARs$`423T`$FusionName,"423T INT", "423T STAR")
simpleplotvenn(nodupes_INTs$`459T`$FusionName,nodupes_STARs$`459T`$FusionName,"459T INT", "459T STAR")
simpleplotvenn(nodupes_INTs$`464T`$FusionName,nodupes_STARs$`464T`$FusionName,"464T INT","464T STAR")
simpleplotvenn(nodupes_INTs$`550T`$FusionName,nodupes_STARs$`550T`$FusionName,"550T INT", "550T STAR")
simpleplotvenn(nodupes_INTs$`645T`$FusionName,nodupes_STARs$`645T`$FusionName,"645T INT", "645T STAR")
simpleplotvenn(nodupes_INTs$`683T`$FusionName,nodupes_STARs$`683T`$FusionName,"683T INT","683T STAR")
simpleplotvenn(nodupes_INTs$`773T1`$FusionName,nodupes_STARs$`773T1`$FusionName,"773T INT","773T STAR")
simpleplotvenn(nodupes_INTs$`liver`$FusionName,nodupes_STARs$`liver`$FusionName,"liver INT", "liver STAR")
simpleplotvenn(nodupes_INTs$`lung`$FusionName,nodupes_STARs$`lung`$FusionName,"lung INT","lung STAR")

PE_combine_and_filter2 <- function(STAR,INT){
  data <- rbind(STAR,INT)
  # sorting by descending spanning reads so that the highest one is kept out of the duplicates:
  data <- data[order((data$JunctionReads+data$SplitReads),decreasing = T),]
  # and now removing any duplicates that occur in the FusionName column.
  cat(sum(duplicated(data$FusionName)),"duplicates removed \n")
  data <- data[!duplicated(data$FusionName), ]
  return(data)
}
both_callers <- mapply(function(X,Y) {PE_combine_and_filter2(X,Y)}, X = nodupes_STARs, Y=nodupes_INTs, SIMPLIFY = F)
dog_tumouronly2 <- lapply(both_callers[-c(10,11)], function(x) {remove_normals(x,both_callers$lung,both_callers$liver)})

# Upset plot excluding JAFFA calls:
dog_upset_noJAF <- lapply(dog_tumouronly2, function(x) {return(x$FusionName)})
upset(fromList(dog_upset_noJAF), nsets = 9,order.by="freq",nintersects = 40, text.scale = 1.5,
      mainbar.y.label = "Number of Shared Fusions", sets.x.label = "Fusions Per Sample",
      main.bar.color = magma(10,alpha=0.5)[5], sets.bar.color = viridis(10,alpha=0.5)[8])

# Correlation matrix without JAFFA calls:
noJAFlength_list <- lapply(dog_tumouronly2, function(X) {nrow(X)})
noJAFlongest_list <- noJAFlength_list[[which(noJAFlength_list == max(unlist(noJAFlength_list)))]]
noJAFcorr_matrix <- matrix(NA,length(dog_upset_noJAF),length(dog_upset_noJAF))
colnames(noJAFcorr_matrix) <- rownames(noJAFcorr_matrix) <- dogfiles[-c(10,11)]
for (a in 1:length(dog_upset_noJAF)){
  for (b in 1:length(dog_upset_noJAF)){
    noJAFcorr_matrix[a,b] <- length(intersect(dog_upset_noJAF[[a]],dog_upset_noJAF[[b]]))
    noJAFcorr_matrix[a,b] <- noJAFcorr_matrix[a,b] / noJAFlongest_list
  }
}

# heatmap when excluding JAFFA calls


## Overlapping 59 fusions = same set as with JAFFA; use these to make a circos plot and a table:
part1 <- intersect(dog_upset_noJAF$`366T`,dog_upset_noJAF$`410T`)
part2 <- intersect(part1, dog_upset_noJAF$`423T`)
part3 <- intersect(part2, dog_upset_noJAF$`459T`)
part4 <- intersect(part3, dog_upset_noJAF$`464T`)
part5 <- intersect(part4, dog_upset_noJAF$`550T`)
part6 <- intersect(part5, dog_upset_noJAF$`645T`)
part7 <- intersect(part6, dog_upset_noJAF$`683T`)
allCTVT <- intersect(part7, dog_upset_noJAF$`773T1`)

lbed <- data.frame(CTVT_consensus[,c(3,4,4,2)])
rbed <- data.frame(CTVT_consensus[,c(7,8,8,6)])
pairs <- cbind(lbed,rbed)
pairs$LeftContig <- paste0("Chr",pairs$LeftContig)
pairs$RightContig <- paste0("Chr",pairs$RightContig)
newdogcircos()
circos.genomicLink(region1 = pairs[,1:3],region2 = pairs[,5:7])

stargazer(CTVT_consensus[,c(2,3,4,6,7,8,10,11)],summary=F,rownames = F)
