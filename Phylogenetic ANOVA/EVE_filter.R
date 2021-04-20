#Script to filter the Expression values 
#Splits the expression data by all 3 Tissues (Green,Bud,Open), filters them and finds the intersection 
#Input needed: expression table (exp_df), filter (exp_higher_than) and number of samples the filter should apply to (N_smp)

filter_low_exp_rm_smps8 <- function(exp_df= NULL, delete_smp =
                                     NULL,N_smp=3,exp_higher_than=5){
  
  #delete_smp <- rmv_smp.sp
  
  #exp_df <- sp.tpm
  
  #colnames(exp_df)
  
  TPM <- subset(exp_df,select=colnames(exp_df)[!colnames(exp_df) %in%
                                                 delete_smp])
  
  smp_colname <- colnames(TPM)
  
  O_samples <- grep(pattern = "_O_\\w+", x = smp_colname, value = TRUE) 
  G_samples <- grep(pattern = "_G_\\w+", x = smp_colname, value = TRUE) 
  B_samples <- grep(pattern = "_B_\\w+", x = smp_colname, value = TRUE)
  

  O_raw_table <- subset(TPM, select = O_samples) 
  G_raw_table <- subset(TPM, select = G_samples) 
  B_raw_table <- subset(TPM, select = B_samples)

  O_thresh <- O_raw_table > exp_higher_than 
  G_thresh <- G_raw_table > exp_higher_than 
  B_thresh <- B_raw_table > exp_higher_than
  

  O_keep <- rowSums(O_thresh) >= N_smp 
  G_keep <- rowSums(G_thresh) >= N_smp 
  B_keep <- rowSums(B_thresh) >= N_smp
  

  O_keep_table <- O_raw_table[O_keep,] 
  G_keep_table <- G_raw_table[G_keep,] 
  B_keep_table <- B_raw_table[B_keep,]
  merge(B_keep_table,O_keep_table, by=0) -> BOT2 
  rownames(BOT2) <- BOT2$Row.names 
  BOT2[,-1] -> BOT2
  merge(BOT2,G_keep_table,by=0) -> GBO_table 
  rownames(GBO_table) <- GBO_table[,1]
  GBO_table[,-1] -> GBO_table

  
  res <- GBO_table
  
  return(res)
  
} 

h_TPM_F[["sp"]][["F"]] <- list(GBO_table)

filter_low_exp_rm_smps
res <- filter_low_exp_rm_smps8(exp_df= sp.tpm, delete_smp =rmv_smp.sp)
h_TPM[["sp"]][["F"]] <- GBO_table

