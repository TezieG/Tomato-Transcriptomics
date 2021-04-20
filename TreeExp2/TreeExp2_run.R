####################################PACKAGES###################################
library('dplyr')
library('ggplot2')
library('ggsignif')
library('TreeExp')
library('ape')
library('pgirmess')
library('multcompView')
library('forcats')
library('stringr')
library('ggpmisc')
library('hash')
library("ggtree") 
library("reshape2")
library("ggpubr")
###modified corrMatInv function to read in customized phylogeny as distance metrix 
###and produce invert correlation matrix
{
corrMatInv.phyTree = function(tree = NULL) {
  #method <- match.arg(method)
  
  #dis.mat <- expdist(taxa.obj, taxa = "all",subtaxa = k)
  dis.mat <- cophenetic.phylo(tr)
  
  corr.mat <- as.matrix(1 - as.dist(dis.mat))
  diag(corr.mat) <- 1
  
  solve(corr.mat)
} 

##################Modified functions (from TreeExp2)
TEconstruct.readin_exp_mtx = function(ExpValueFP=NULL, taxa="all", subtaxa="all",
                                      rmOut=FALSE, verbose=FALSE) {
  
  # check file handle
  if((is.null(ExpValueFP))){
    stop(paste0(date(),": must provide expression values file path and gene length file path"))
  }
  
  # check file existance
  # if(!file.exists(ExpValueFP)){
  #   stop(paste0(date(),": fail to open file, check your filename or path"))
  # }
  #browser()
  
  # input
  #exp_value_df <- read.table(ExpValueFP, header = TRUE)
  #row.names(exp_value_df) <- exp_value_df[,1]
  #exp_value_df <- exp_value_df[,-1]
  
  # gene.info.df <- read.table(geneInfoFP,header=T)
  # remove sample with low read counts
  exp_value_df <- ExpValueFP
  invalid_arr <- NULL
  
  for (i in 2:ncol(exp_value_df)) {
    if (mean(exp_value_df[,i]) < 1) {
      invalid_arr <- c(invalid_arr,i)
    }
  }
  
  message(paste0(date(),": removing ", length(invalid_arr), " sample(s) with ultra-low read counts"))
  
  if (length(invalid_arr) != 0) {
    invalid_arr = 0 - invalid_arr
    exp_value_df <- exp_value_df[,invalid_arr]
  }
  
  # gene number and taxon number
  gene_n <- nrow(exp_value_df)
  
  # get taxon names from read counts file
  taxon_names <- unique(lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[1]))
  taxon_n <- length(taxon_names)
  
  # normalize<-match.arg(normalize)
  # message(paste0(date(), ": using ", normalize, " to normalize raw read counts"))
  
  # get taxon names
  #browser()
  cat("\n")
  message(paste0(date(),": start constructiong TE objects"))
  
  
  if (!any(grepl("all", taxa, ignore.case = TRUE))) {
    
    taxon_names <- gsub("\\s+", "", taxa)
    taxon_n <- length(taxon_names)
    
  }
  
  message(paste0(date(),": total Taxon number ", taxon_n))
  
  #browser()
  
  # get subtaxon number
  subtaxon_names <- unique(lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[2]))
  subtaxon_n <- length(subtaxon_names)
  
  
  if (!any(grepl("all", subtaxa, ignore.case = TRUE))) {
    
    subtaxon_names <- gsub("\\s+", "", subtaxa)
    subtaxon_n <- length(subtaxon_names)
    
  }
  
  message(paste0(date(),": total sub taxon number ", subtaxon_n))
  
  title <- lapply(colnames(exp_value_df), function(x) unlist(strsplit(x, "_"))[1]) # taxon names
  subtitle <- lapply(colnames(exp_value_df), function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
  first_two_names <- unique(paste(title,subtitle,sep="_"))
  
  index <- intersect(unlist(lapply(taxon_names, function(x) grep(x, first_two_names, ignore.case = TRUE))),
                     unlist(lapply(subtaxon_names,function(x) grep(x, first_two_names, ignore.case = TRUE))))
  
  objects_names <- first_two_names[index]
  
  objects_number <- length(objects_names)
  
  #browser()
  
  cat("\n")
  # get gene names
  #gene.names <- read.counts.df[,1]
  
  message(paste0(date(),": now constructing ",objects_number, " TE objects..."))
  
  if (!verbose) progbar <- txtProgressBar(style = 3)
  
  # initialization
  
  taxonExp.objects <- vector("list",length = objects_number)
  # the number of TE objects constructed is based on seleted numnber
  
  # for each taxon
  
  objects_counter <- 0
  
  for (i in 1:objects_number) {
    
    #browser()
    if (verbose) message(paste0(date(),": proceeding taxon ", objects_names[i]))
    
    # get all the sample names matching objects names
    # bundle all the biological replicates into one TE object
    
    #browser()
    
    ttl <- unlist(strsplit(objects_names[i], "_"))[1] #taxon title
    subttl <- unlist(strsplit(objects_names[i], "_"))[2] # subtaxon title
    #ttl <- lapply(names, function(x) unlist(strsplit(x, "_"))[1]) # taxon names
    #subttl <- lapply(names, function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
    
    idx <- grep(objects_names[i],colnames(exp_value_df), ignore.case = TRUE)
    names <- strsplit(colnames(exp_value_df)[idx],"_")
    repttl <- unlist(lapply(names, function(x) unlist(strsplit(x,"_"))[3])) # biological replicates title names
    gene_names <- rownames(exp_value_df)  # gene names
    
    
    # foreach subtaxon
    bio_rep_n <- length(repttl) # biological replicates number
    #    omega <- NULL # omega estimated overdispersion parameter
    
    exp_val <- apply(exp_value_df[idx], c(1,2), as.numeric)
    
    objects_counter = objects_counter + 1
    
    if (verbose) message(paste0(date(),": wrapping up into objects"))
    
    #browser()
    oneObject <- list(exp_value=exp_val,
                      taxon_name = ttl,subTaxon_name = subttl,
                      gene_num = gene_n, gene_name = gene_names,
                      bioRep_num = bio_rep_n, bioRep_id = repttl)
    
    class(oneObject) <- "taxonExp"
    
    taxonExp.objects[[objects_counter]] <- oneObject
    
    #browser()
    
    if (verbose) message(paste0(date(), ": ", objects_counter, " TE objects constructed"))
    
    if (verbose) cat("\n")
    
    if (!verbose) setTxtProgressBar(progbar, objects_counter/objects_number)
    
    
  }
  
  class(taxonExp.objects) <- "taxaExp"
  
  attr(taxonExp.objects, "taxa") <- unlist(taxon_names)
  attr(taxonExp.objects, "subtaxa") <- unlist(subtaxon_names)
  
  cat("\n")
  message(date(),": construction complete.")
  
  taxonExp.objects
  
}
###########################Construction of Expression Values####################################################### 
#load("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Workspace/a.Rdata")
setwd("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Results/Temp")

#######################################calculate W for each gene in each stage based on TPM (based on both accession tree and species tree)
h_TPM_W <- hash()
h_TPM_gammaPar <- hash()
acc_tree <- read.tree("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Trees/acc/F24LAMOD.txt")
#sp_tree3 <- read.tree("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Trees/sp/LazySpecies.txt")


{
tr <- acc_tree
i <- "acc"
j <- "F"
k <- "G"



taxa.obj <- TEconstruct.readin_exp_mtx(ExpValueFP = data.frame(h_TPM_f[[i]][[j]]), taxa = "all", subtaxa = k)

#k <- "C"
#inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
#inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
#k <- "B"
#inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
#inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
inv.corr.mat <- corrMatInv.phyTree(tr)


#exptable <- exptabTE(taxa.objects, taxa = 'all', subtaxa = st ,logrithm = F)
exptable <- exptabTE(taxa.obj, taxa = 'all', subtaxa = k ,logrithm = T)
gamma.paras <- estParaGamma(exptable =exptable, corrmatinv =inv.corr.mat)
h_TPM_gammaPar[[i]][[k]] <- gamma.paras
#### Bayesian estimation of gene-specific selection pressure
Q <- estParaQ(exptable, corrmatinv = inv.corr.mat)
# with prior expression values and inversed correlation matrix
post <- estParaWBayesian(Q, gamma.paras)
stage.1 <- post$w # posterior selection pressures
# posterior selection pressures
#stage.CI <- post$ci95 # posterior expression 95% confidence interval
#plot(density(stage.W))
names(stage.1) <- rownames(exptable)
#h_TPM_W[[i]][[k]] <- stage.W 
}
{
stage_G <- stage.1 
write.table(stage_G , "stage_G.txt")
stage_g <- read.csv("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Results/Temp/stage_G.txt", sep="")
stage_g$condition <- "Green"   }

###Bud stage 
{
tr <- acc_tree
i <- "acc"
j <- "F"
k <- "B"

taxa.obj <- TEconstruct.readin_exp_mtx(ExpValueFP = data.frame(h_TPM_f[[i]][[j]]), taxa = "all", subtaxa = k)

#k <- "C"
#inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
#inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
#k <- "B"
#inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
#inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
inv.corr.mat <- corrMatInv.phyTree(tr)


#exptable <- exptabTE(taxa.objects, taxa = 'all', subtaxa = st ,logrithm = F)
exptable <- exptabTE(taxa.obj, taxa = 'all', subtaxa = k ,logrithm = T)
gamma.paras <- estParaGamma(exptable =exptable, corrmatinv =inv.corr.mat)
h_TPM_gammaPar[[i]][[k]] <- gamma.paras
#### Bayesian estimation of gene-specific selection pressure
Q <- estParaQ(exptable, corrmatinv = inv.corr.mat)
# with prior expression values and inversed correlation matrix
post <- estParaWBayesian(Q, gamma.paras)
stage.2 <- post$w # posterior selection pressures
# posterior selection pressures
#stage.CI <- post$ci95 # posterior expression 95% confidence interval
#plot(density(stage.W))
names(stage.2) <- rownames(exptable)
#h_TPM_W[[i]][[k]] <- stage.W  
}
{
stage_B <- stage.2 
write.table(stage_B , "stage_B.txt")
stage_b <- read.csv("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Results/Temp/stage_B.txt", sep="")
stage_b$condition <- "Bud" }

###Open Stage 
{
  tr <- acc_tree
  i <- "acc"
  j <- "F"
  k <- "O"
  
  taxa.obj <- TEconstruct.readin_exp_mtx(ExpValueFP = data.frame(h_TPM_f[[i]][[j]]), taxa = "all", subtaxa = k)
  
  #k <- "C"
  #inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
  #inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
  #k <- "B"
  #inv.corr.mat <- corrMatInv.f(taxa.objects, taxa = 'all', subtaxa = st)
  #inv.corr.mat <- corrMatInv(taxa.obj, taxa = 'all', subtaxa = k)
  inv.corr.mat <- corrMatInv.phyTree(tr)
  
  
  #exptable <- exptabTE(taxa.objects, taxa = 'all', subtaxa = st ,logrithm = F)
  exptable <- exptabTE(taxa.obj, taxa = 'all', subtaxa = k ,logrithm = T)
  gamma.paras <- estParaGamma(exptable =exptable, corrmatinv =inv.corr.mat)
  h_TPM_gammaPar[[i]][[k]] <- gamma.paras
  #### Bayesian estimation of gene-specific selection pressure
  Q <- estParaQ(exptable, corrmatinv = inv.corr.mat)
  # with prior expression values and inversed correlation matrix
  post <- estParaWBayesian(Q, gamma.paras)
  stage.3 <- post$w # posterior selection pressures
  # posterior selection pressures
  #stage.CI <- post$ci95 # posterior expression 95% confidence interval
  #plot(density(stage.W))
  names(stage.3) <- rownames(exptable)
  #h_TPM_W[[i]][[k]] <- stage.W  
} 
{
stage_O <- stage.3 
write.table(stage_O , "stage_O.txt")
stage_o <- read.csv("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Results/Temp/stage_O.txt", sep="")
stage_o$condition <- "Open Flower"  }

#Bind Data together 
ALL <- rbind(stage_g,stage_b,stage_o)
ALL_IN_ALL <- melt(ALL, id.var="condition")
 
#More Boxplot 
CombinedPlot=ggplot(ALL_IN_ALL, aes(x=condition, y=value, fill=condition)) + geom_boxplot() 
CombinedPlot 
} 
#To Draw nicer plot with p values  
{

#compare_means(value ~ condition,  data = ALL_IN_ALL, paired=TRUE, method = "fisher")

my_comparisons <- list( c("Green", "Bud"), c("Green", "Open Flower"), c("Bud", "Open Flower") )
p <- ggboxplot(ALL_IN_ALL, x = "condition", y = "value", ylim =c(0,1) ,
          color = c("Black"),fill=c("palegreen3","salmon","skyblue3") , xlab="Tissue", ylab="Strength of Expression Conservation (W)", subtitle="Wilcoxon, p < 2.2e-16", ) +
  scale_x_discrete(labels=c("Green Bud", "Bud", "Open Flower"))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", paired=TRUE)
  #stat_compare_means(label = "p.signif", method = "wilcox.test")+   
p + grids(linetype = "dashed", size=1) 
}
  
### Check Genes with high/low W expression conservation:

sort(stage.3, decreasing = TRUE) -> High_W_100 
sort(stage.3, decreasing = FALSE) -> Low_W_100 
High_W_100[1:100]-> High_W_100 
Low_W_100[1:100] -> Low_W_100 
write.csv(High_W_100,"High_W_100_f.txt")  
write.csv(Low_W_100,"Low_W_100_f.txt") 
