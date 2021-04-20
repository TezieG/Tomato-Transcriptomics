####Script for Expression profiles 
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
acc_tree <- read.tree("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Trees/acc/F24LAMOD.txt")
corrMatInv.phyTree = function(tree = NULL) {
  #method <- match.arg(method)
  
  #dis.mat <- expdist(taxa.obj, taxa = "all",subtaxa = k)
  dis.mat <- cophenetic.phylo(tr)
  
  corr.mat <- as.matrix(1 - as.dist(dis.mat))
  diag(corr.mat) <- 1
  
  solve(corr.mat)
}
#Insert Gene ID
#("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Workspace/a.Rdata")
head(sort(stage.3, decreasing = TRUE))
z = "Solyc00g058900.2.1"
#ACC boxplot
{
  Species <-  c("LA2133","LA2750","LA2778","LA0111","LA2744",
                "LA1552","LA1954","LA1626","LA2695",
                "LA2167","LA1395","LA1330","LA2172",
                "LA0247","LA1306","LA1317","LA2879",
                "LA1305","LA1321","LA2981B","LA2650","LA2185","LA4329","WT","WT2","WT3")
  



  #Insert Gene ID
  #load("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Workspace/a.Rdata")
  #z = "Solyc02g083950.3.1"
  
    
    #Species <-  c("Sneorickii_LA2133","Schilense_LA2750","Schilense_LA2778","Speruvinaum_LA0111","Speruvinaum_LA2744",
                #  "Speruvinaum_LA1552","Speruvinaum_LA1954","Sarcanum_LA1626","Schmielewskii_LA2695",
                #  "Shabrochaites_LA2167","Sarcanum_LA1395","Schmielewskii_LA1330","Sarcanum_LA2172",
                #  "Sneorickii_LA0247","Schmielewskii_LA1306","Schmielewskii_LA1317","Schilense_LA2879",
                #  "Speruvinaum_LA1305","Sneorickii_LA1321","Speruvinaum_LA2981B","Shabrochaites_LA2650","Sarcanum_LA2185","Schilense_LA4329","Slycopersicum_WT1","Slycopersicum_WT2","Slycopersicum_WT3")
    
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
    #load("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Workspace/a.Rdata") 
    tr <- acc_tree
    i <- "acc"
    j <- "F"
    
    ExpValueFP = data.frame(h_TPM[[i]][[j]]) 
    ExpValueFP %>% 
      rename(
        WT3_O_F31 = WT_O_F31,
        WT3_G_F31 = WT_G_F31, 
        WT3_B_F31 = WT_B_F31, 
        WT2_O_F31 = WT_O_F30,
        WT2_G_F31 = WT_G_F30, 
        WT2_B_F31 = WT_B_F30
      ) -> ExpValueFp2
    
    
    ExpValueFP = data.frame(h_TPM[[i]][[j]]) 
    ExpValueFP %>% 
      rename(
        WT3_O_F31 = WT_O_F31,
        WT3_G_F31 = WT_G_F31, 
        WT3_B_F31 = WT_B_F31, 
        WT2_O_F31 = WT_O_F30,
        WT2_G_F31 = WT_G_F30, 
        WT2_B_F31 = WT_B_F30
      ) -> ExpValueFp3
    
    ExpValueFP = data.frame(h_TPM[[i]][[j]]) 
    ExpValueFP %>% 
      rename(
        WT3_O_F31 = WT_O_F31,
        WT3_G_F31 = WT_G_F31, 
        WT3_B_F31 = WT_B_F31, 
        WT2_O_F31 = WT_O_F30,
        WT2_G_F31 = WT_G_F30, 
        WT2_B_F31 = WT_B_F30
      ) -> ExpValueFp4
    
    k <- "G"
    taxa.obj <- TEconstruct.readin_exp_mtx(ExpValueFP = ExpValueFp2 , taxa = "all", subtaxa = k)
    
    
    inv.corr.mat <- corrMatInv.phyTree(tr)
    
    
    exptableG <- exptabTE(taxa.obj, taxa = 'all', subtaxa = k ,logrithm = T)
    ####For Open Flower####
    grep(z,rownames(exptableG)) -> I
    exptableG[I,] -> GOI 
    as.data.frame(GOI) -> GOM
    
    #rownames(GOM) <- df1
    #df2 <- as.matrix(GOM)
    #colnames(df2) <- df1
    
    as.data.frame(Species) -> names
    names$x2 <- GOM$GOI  
    
    #ALL <- rbind(stage_g,stage_b,stage_o)
    ALL_IN_ALL <- melt(names, id.var="Species")
    ALL_IN_ALL$Species <- factor(ALL_IN_ALL$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                     "WT","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))
    
    # Boxplot 
    #CombinedPlot=ggplot(ALL_IN_ALL,title="EXP", aes(x=Species, y=value, fill=Species, title(main="Expression Green Bud"))) + geom_boxplot() 
    #CombinedPlot 
    #Alternative 
    p1 <- ggbarplot(ALL_IN_ALL, x = "Species", y = "value", ylim =c(0,5) ,
                    xlab="Accession", ylab="Expression", fill="Species", title="Green Bud")+ 
      #scale_x_discrete(labels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                #"WT","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))+
      grids(linetype = "dashed", size=1)
    #p1 + theme(legend.position = "none") 
    p1 + rotate_x_text() + theme(legend.position = "none")
    
    
    
    ####Get it for Open Flower##### 
    ###GET EXPRESSION TABLE#####
    
    k <- "O"  
    taxa.obj2 <- TEconstruct.readin_exp_mtx(ExpValueFP = ExpValueFp2, taxa = "all", subtaxa = k)
    
    inv.corr.mat <- corrMatInv.phyTree(tr)
    
    exptableO <- exptabTE(taxa.obj2, taxa = 'all', subtaxa = k ,logrithm = T) 
    
    ###GET THE GENE#####
    {
      grep(z,rownames(exptableO)) -> I
      exptableO[I,] -> GOI2 
      as.data.frame(GOI2) -> GOM2
      
      
      as.data.frame(Species) -> names2
      names2$x2 <- GOM2$GOI2  
      
      #ALL <- rbind(stage_g,stage_b,stage_o)
      ALL_IN_ALL2 <- melt(names2, id.var="Species")
      ALL_IN_ALL2$Species <- factor(ALL_IN_ALL2$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                                 "WT","WT2","WT3","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))
      
      #More Boxplot 
      #CombinedPlot=ggplot(ALL_IN_ALL2, aes(x=Species, y=value, fill=Species, ylim(0,4))) + geom_boxplot() 
      #CombinedPlot 
    }
    
    p2 <- ggbarplot(ALL_IN_ALL2, x = "Species", y = "value", ylim =c(0,4) ,
                    xlab="Accession", ylab="Expression", fill="Species", title="Open Flower")+ 
      grids(linetype = "dashed", size=1)
    
    p2 + rotate_x_text() + theme(legend.position = "none")
    #stat_compare_means(label = "p.signif", method = "wilcox.test")+   
    #p + grids(linetype = "dashed", size=1) 
    
    
    ####Bud 
    
    k <- "B"  
    taxa.obj2 <- TEconstruct.readin_exp_mtx(ExpValueFP = ExpValueFp2, taxa = "all", subtaxa = k)
    
    inv.corr.mat <- corrMatInv.phyTree(tr)
    
    exptableB <- exptabTE(taxa.obj2, taxa = 'all', subtaxa = k ,logrithm = T) 
    
    ###GET THE GENE#####
    
    grep(z,rownames(exptableB)) -> I
    exptableB[I,] -> GOI3 
    as.data.frame(GOI3) -> GOM3
    
    
    as.data.frame(Species) -> names3
    names3$x2 <- GOM3$GOI3  
    
    #ALL <- rbind(stage_g,stage_b,stage_o)
    ALL_IN_ALL3 <- melt(names3, id.var="Species")
    ALL_IN_ALL3$Species <- factor(ALL_IN_ALL3$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                               "WT","WT2","WT3","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))
    
    #More Boxplot 
    CombinedPlot=ggplot(ALL_IN_ALL2, aes(x=Species, y=value, fill=Species, ylim(0,4))) + geom_boxplot() 
    CombinedPlot 
    
    
    p3 <- ggbarplot(ALL_IN_ALL3, x = "Species", y = "value", ylim =c(0,4) ,
                    xlab="Accession", ylab="Expression", fill="Species", title="Bud") 
      #grids(linetype = "dashed", size=1)
    
    p3 + rotate_x_text() + theme(legend.position = "none")
    #stat_compare_means(label = "p.signif", method = "wilcox.test")+   
    #p + grids(linetype = "dashed", size=1) 
    
    
    #grep("Solyc02g083950.3.1",rownames(OPEN)) -> i
    #exptable[i,] -> GOI3 
    #p1 
    #p2 
    #p3
    
  }
  
  
p1 + theme(legend.position = "none")+ rotate_x_text() 
  p3 + theme(legend.position = "none")+ rotate_x_text()
  p2 + theme(legend.position = "none") + rotate_x_text()
#ACC PHYLO BARPLOT 
#ALL_IN_ALL2[,3] -> ALLB
#names(ALLB) <- ALL_IN_ALL2[,1]
#plotTree.barplot(acc_tree,ALLB,xlim=c(0,max(ALL6)),xlab="Expression") 

plotTree.wBars(acc_tree,ALL6,scale=3,tip.label=TRUE,fsize=0.5)
  #scale_x_continuous(limits=c(0,3))
#mtext("Expression",1,at=0.5,line=2.5))
######################################################################
ExpValueFP5 <- data.frame(h_TPM[i][j]) 


###QUICK RUN 
head(sort(stage.3, decreasing = TRUE)) 
z = "Solyc03g112995.1.1" 

{ 
  grep(z,rownames(stage_g)) -> I
  exptableG[I,] -> GOI 
  as.data.frame(GOI) -> GOM
  
  #rownames(GOM) <- df1
  #df2 <- as.matrix(GOM)
  #colnames(df2) <- df1
  
  as.data.frame(Species) -> names
  names$x2 <- GOM$GOI  
  
  #ALL <- rbind(stage_g,stage_b,stage_o)
  ALL_IN_ALL <- melt(names, id.var="Species")
  ALL_IN_ALL$Species <- factor(ALL_IN_ALL$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                             "WT","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))
  
  # Boxplot 
  #CombinedPlot=ggplot(ALL_IN_ALL,title="EXP", aes(x=Species, y=value, fill=Species, title(main="Expression Green Bud"))) + geom_boxplot() 
  #CombinedPlot 
  #Alternative 
  p1 <- ggbarplot(ALL_IN_ALL, x = "Species", y = "value", ylim =c(0,10) ,
                  xlab="Accession", ylab="Expression", fill="Species", title="Green Bud")+ 
    #scale_x_discrete(labels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
    #"WT","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))+
    grids(linetype = "dashed", size=1)
  #p1 + theme(legend.position = "none") 
  p1 + rotate_x_text() + theme(legend.position = "none")  
  
  grep(z,rownames(stage_o)) -> I
  exptableO[I,] -> GOI2 
  as.data.frame(GOI2) -> GOM2
  
  
  as.data.frame(Species) -> names2
  names2$x2 <- GOM2$GOI2  
  
  #ALL <- rbind(stage_g,stage_b,stage_o)
  ALL_IN_ALL2 <- melt(names2, id.var="Species")
  ALL_IN_ALL2$Species <- factor(ALL_IN_ALL2$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                               "WT","WT2","WT3","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))
  
  #More Boxplot 
  #CombinedPlot=ggplot(ALL_IN_ALL2, aes(x=Species, y=value, fill=Species, ylim(0,4))) + geom_boxplot() 
  #CombinedPlot 


p2 <- ggbarplot(ALL_IN_ALL2, x = "Species", y = "value", ylim =c(0,10) ,
                xlab="Accession", ylab="Expression", fill="Species", title="Solyc03g112995.1.1", subtitle = "W =  1.0646")+ 
  grids(linetype = "dashed", size=1)

p2 + rotate_x_text() + theme(legend.position = "none")
#stat_compare_means(label = "p.signif", method = "wilcox.test")+   
#p + grids(linetype = "dashed", size=1)  

grep(z,rownames(stage_g)) -> I
exptableB[I,] -> GOI3 
as.data.frame(GOI3) -> GOM3


as.data.frame(Species) -> names3
names3$x2 <- GOM3$GOI3  

#ALL <- rbind(stage_g,stage_b,stage_o)
ALL_IN_ALL3 <- melt(names3, id.var="Species")
ALL_IN_ALL3$Species <- factor(ALL_IN_ALL3$Species , levels=c("LA2650", "LA2167", "LA2750","LA2778","LA2879","LA4329","LA2744","LA0111","LA1552","LA1305","LA1954","LA2981B",
                                                             "WT","WT2","WT3","LA1330","LA2695","LA1306","LA1317","LA1626","LA1395","LA2172","LA2185","LA1321","LA2133","LA0247"))

#More Boxplot 
#CombinedPlot=ggplot(ALL_IN_ALL2, aes(x=Species, y=value, fill=Species, ylim(0,4))) + geom_boxplot() 
#CombinedPlot 


p3 <- ggbarplot(ALL_IN_ALL3, x = "Species", y = "value", ylim =c(0,10) ,
                xlab="Accession", ylab="Expression", fill="Species", title="Bud")+ 
  theme(legend.position = "none")
#grids(linetype = "dashed", size=1)

p3 + rotate_x_text() + theme(legend.position = "none")
  
} 
p1 + theme(legend.position = "none")+ rotate_x_text() 
p3 + theme(legend.position = "none")+ rotate_x_text()
p2 + theme(legend.position = "none") + rotate_x_text()

sort(stage.3, decreasing = TRUE) -> LIST4
LIST4                          



####multiple plots in 1 
p2 + rotate_x_text() + theme(legend.position = "none") -> L4 

ggarrange(L1+ rremove("x.text"), L2+ rremove("x.text"), L3+ rremove("x.text"), L4 + rremove("x.text"), 
                     labels = c("#1", "#5", "#10","#20"),
                     ncol = 2, nrow = 2)
