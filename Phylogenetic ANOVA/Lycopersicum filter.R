#Description: {
#Script was written to help with manual observation to find genes that have a unique expression in Solanum Lycopersicum(SL)
#Script will divide the expression table (Only SL and the rest), apply the respective filter and only keep the intersection where both filters apply 
#(either high expression in SL and low in all other samples (WT_UP) or vice versa (WT_DOWN)
#The results depend on the input filter, it will always filter out genes that might also have a unique expression in SL 
#The results have to be manually checked!!!
####INPUT NEEDED#### 
#1: expression table (exp_df)
#2: filter for Solanum Lycopersicum (WT_filter)
#3: filter for all other samples (REST_filter) 
#4: number of samples the filter has to apply to (REST_NR and SL_NR) 
#Gene IDs of interest(GeneList), has to be first row in a table without names 
}

LycoFilter <- function(exp_df= NULL, 
             WT_filter=NULL,REST_filter=NULL,REST_NR=NULL, SL_NR=NULL, GeneList=NULL){  
#Get expression table of Open Flower genes  

rownames(GeneList) <- GeneList$V1
merge(exp_df,GeneList,by = 0) -> FDREXP5 
newt <- FDREXP5 
  
  OTHER <- newt[,1:21]
  rownames(OTHER) <- rownames(newt)
  WT <- as.data.frame(newt[,22:24])
  colnames(WT) <- c("WT_O","WT2_O","WT3_O")
  WT$x2 <- newt[,1]
  rownames(WT) <- rownames(newt)
  ##Genes where expression of SL is lower than the rest  
  WT_THRESH <- WT <= WT_filter
  WT_keep <- rowSums(WT_THRESH) >= SL_NR
  WT_keep_table <- WT[WT_keep,]

  ##Same process for all other samples 
  OTHER_THRESH <- OTHER >= REST_filter
  other_keep <- rowSums(OTHER_THRESH) >= REST_NR
  
  
  OTHER_KEEP_TABLE <- OTHER[other_keep,]
  ###MERGE DF TOGETHER TO GET EXPRESSION TABLE WHERE BOTH THRESHHOLDS APPLY 
  merge(OTHER_KEEP_TABLE,WT_keep_table, by = 0) -> WT_DOWN 
  
  
  ##Genes where expression of SL is higher than the rest  
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
res <- LycoFilter(exp_df=, GeneList=, WT_filter=0,REST_filter=0,REST_NR=0,SL_NR=0)
res[[1]] -> SL_DOWN 
res[[2]] -> SL_UP  


