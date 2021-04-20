##Filtration that splits up dataframes by tissues and applies the the parameters below. Only genes where this applied to all 3 dev stages are kept
filter_low_exp_rm_smps8 <- function(exp_df= NULL, delete_smp =
                                     NULL,N_smp=3,exp_higher_than=5){
  
  #delete_smp <- rmv_smp.sp
  
  #exp_df <- sp.tpm
  
  #colnames(exp_df)
  
  TPM <- subset(exp_df,select=colnames(exp_df)[!colnames(exp_df) %in%
                                                 delete_smp])
  
  smp_colname <- colnames(TPM)
  
  #F_samples <- grep(pattern = "_F\\d+", x = smp_colname, value = TRUE) 
  O_samples <- grep(pattern = "_O_\\w+", x = smp_colname, value = TRUE) 
  G_samples <- grep(pattern = "_G_\\w+", x = smp_colname, value = TRUE) 
  B_samples <- grep(pattern = "_B_\\w+", x = smp_colname, value = TRUE)
  
  #grep(pattern = "_G", x = smp_colname, value = TRUE), #Green bud stage
  #group
  
  #grep(pattern = "_B", x = smp_colname, value = TRUE), #Bud stage group
  
  #grep(pattern = "_O", x = smp_colname, value = TRUE)) #Opening flower
  #stage group
  
  
  
  #L_samples <- grep(pattern = "_L\\d+", x = smp_colname, value = TRUE)
  #Feeding(treated) group
  
  
  
  #F_raw_table <- subset(TPM, select = F_samples) 
  O_raw_table <- subset(TPM, select = O_samples) 
  G_raw_table <- subset(TPM, select = G_samples) 
  B_raw_table <- subset(TPM, select = B_samples)
  
  #L_raw_table <- subset(TPM, select = L_samples)
  
  
  
  #######filter flower samples
  
  #F_thresh <- F_raw_table > exp_higher_than 
  O_thresh <- O_raw_table > exp_higher_than 
  G_thresh <- G_raw_table > exp_higher_than 
  B_thresh <- B_raw_table > exp_higher_than
  
  #F_keep <- rowSums(F_thresh) >= N_smp 
  O_keep <- rowSums(O_thresh) >= N_smp 
  G_keep <- rowSums(G_thresh) >= N_smp 
  B_keep <- rowSums(B_thresh) >= N_smp
  
  #F_keep_table <- F_raw_table[F_keep,] 
  O_keep_table <- O_raw_table[O_keep,] 
  G_keep_table <- G_raw_table[G_keep,] 
  B_keep_table <- B_raw_table[B_keep,]
  merge(B_keep_table,O_keep_table, by=0) -> BOT 
  rownames(BOT) <- BOT$Row.names 
  BOT[,-1] -> BOT
  merge(BOT,G_keep_table,by=0) -> GBO_table 
  rownames(GBO_table) <- GBO_table[,1]
  GBO_table[,-1] -> GBO_table
  #length(rownames(F_keep_table))#18447
  
  #######filter leaf samples
  
  #L_thresh <- L_raw_table > exp_higher_than
  
  #L_keep <- rowSums(L_thresh) >= N_smp
  
  #L_keep_table <- L_raw_table[L_keep,]
  
  res <- GBO_table
  rownames(res) <- rownames(GBO_table)
  return(res)
  
} 


filter_low_exp_rm_smps
res <- filter_low_exp_rm_smps8(exp_df= acc.tpm, delete_smp =rmv_smp.acc)
h_TPM_f[["acc"]][["F"]] <- res
h_TPM_f[["acc"]][["F"]] <- GBO_table
