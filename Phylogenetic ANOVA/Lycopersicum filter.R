####INPUT NEEDED#### 
#1: expression table 
#2: filtration of expression table 
#3: WT filter 
#4: REST filter 
#Gene IDs in GeneList have to be the rownames
#exp_df = exptable
as.data.frame(OPENY1[,2]) -> OpenList
rownames(OpenList) <- OpenList$Gene
LycoFilter <- function(exp_df= NULL, 
             N_smp=NULL,exp_higher_than=NULL,WT_filter=NULL,REST_filter=NULL, GeneList=NULL){  
  
  
  #merge(exp_df,GeneList,by = 0) -> FDREXP5 
  
  #thold <- FDREXP5 >= exp_higher_than
  
  
  
  #FK <- rowSums(thold) >= N_smp
  
  #newt <- FDREXP5[FK,] 
  
  OTHER <- newt[,1:21]
  rownames(OTHER) <- rownames(newt)
  WT <- as.data.frame(newt[,22:24])
  colnames(WT) <- c("WT_O","WT2_O","WT3_O")
  WT$x2 <- newt[,1]
  rownames(WT) <- rownames(newt)
  ##ONLY 0 
  WT_THRESH <- WT <= WT_filter
  WT_keep <- rowSums(WT_THRESH) >= 1
  WT_keep_table <- WT[WT_keep,]
  #WT_keep_table <- as.data.frame(WT[WT_THRESH,])
  #ABOVE 0 
  #OTHER$x25 <- GENEID$`rownames(newt)`
  
  OTHER_THRESH <- OTHER >= REST_filter
  other_keep <- rowSums(OTHER_THRESH) >= 12
  
  
  OTHER_KEEP_TABLE <- OTHER[other_keep,]
  ###MERGE DF TOGETHER TO GET EXPRESSION TABLE WHERE BOTH THRESHHOLDS APPLY 
  merge(OTHER_KEEP_TABLE,WT_keep_table, by = 0) -> WT_DOWN 
  
  
  ##TURN AROUUUND 
  WT_THRESH <- WT >= WT_filter
  WT_keep_table <- as.data.frame(WT[WT_THRESH,]) 
  
  OTHER_THRESH <- OTHER <= REST_filter
  other_keep <- rowSums(OTHER_THRESH) >= 12
  OTHER_KEEP_TABLE <- OTHER[other_keep,] 
  merge(OTHER_KEEP_TABLE,WT_keep_table, by = 0) -> WT_UP 
  
  res <- list(WT_DOWN,WT_UP)
  
  return(res) 
}

#
res <- LycoFilter(exp_df= TPM_O, GeneList=GENE_O, WT_filter=0,REST_filter=0,N_smp=15,exp_higher_than=0)
res[[1]] -> AAA_D #WTDOWN
res[[2]] -> AAA_U #WTUP 
library(openxlsx)
setwd("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/SL_Shift")
write.xlsx(WT_DOWN_O, "WT_DOWN_O.xlsx") 
write.xlsx(WT_UP_O, "WT_UP_O.xlsx")
write.xlsx(WT_DOWN_O_f, "WT_DOWN_O_f.xlsx")  

write.xlsx(WT_DOWN_O_f, "WT_DOWN_O_f.xlsx")  



setwd("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/SL_Shift/corr")
write.csv(WT_DOWN_B,"WT_DOWN_B.txt")
write.csv(WT_DOWN_G,"WT_DOWN_G.txt")  
write.csv(WT_DOWN_O,"WT_DOWN_O.txt") 
write.csv(WT_UP_B,"WT_UP_B.txt") 
write.csv(WT_UP_G,"WT_UP_G.txt") 
write.csv(WT_UP_O,"WT_UP_O.txt") 


for (i in 1:nrow(exptableG))
  for (j in 1:ncol(exptableG)){
  if(median(exptableG[i,1:23])*2 < median(exptableG[i,24:26])){
    mat2[i,j] = exptableG[i,j]
    print("FOUND ONE!") 
    #list(exptableG[i,]) -> loop
}}

#median(exptableG[1,1:24])*2 
#if k > median[exptableG[1,24:26]] 
#->WTDOWN 
#if k < median[exptableG[1,24:26]] 
#
#->WTUP


