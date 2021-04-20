####Beta shared test from The EVE model (https://github.com/Jmendo12/evemodel)
####Script requires the Species Tree and the filtered h_TPM fle 

h_EVE_res <- hash()
h_tmp_TPM <- hash()
h_tmp_colsp <- hash()

filter_low_exp <- function(exp_df= NULL, N_smp=1,exp_higher_than=1){
  #exp_df <- h_tmp_TPM[["G"]]
  thresh <- data.frame(exp_df) > exp_higher_than
  keep <- rowSums(thresh) >= N_smp
  keep_table <- exp_df[keep,]
  return(keep_table)
}

#flower
ti <- "F" 
stage <- c("G","B","O")
for (st in stage){
  #st <- "G"
  TPM <-  log2(0.01 + h_TPM[["sp"]][[ti]])
  col_l <- grep(paste0("_",st,"_"),colnames(TPM),value = T)
  TPM <- TPM[,col_l]
  #speciesTree$tip.label
  
  TPM.f <- filter_low_exp(TPM)
  h_tmp_TPM[[st]] <- as.matrix(TPM.f)
  #h_tmp_TPM[[st]] <- as.matrix(TPM[1:10,])
  h_tmp_colsp[[st]] <- sub("_.*$","",colnames(TPM))
}


run_eve <- function(x){
  res <- betaSharedTest(tree = speciesTree, gene.data = h_tmp_TPM[[x]], colSpecies = h_tmp_colsp[[x]])
  return(res)
}
#h_tmp_TPM[[x]][124,]

x <- "G"
x <- "B"
x <- "O"

h_EVE_res[[x]] <- run_eve(x)
