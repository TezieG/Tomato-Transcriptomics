### Function to check the expression profile from the EVE model results 
### Note that the TPM values are log2 transformed 
### Note that you must load the correct workspace for the according tissue 
### For Open Flower: 
#load(C:...\Finish Up\EVE model (Phyolgenetic ANOVA)\Results\Workspace\Results_Open.Rdata)
#Solyc05g006140.1.1
z = "Solyc12g088180.1.1" 
{
grep(z,rownames(TPM)) -> i
TPM[i,] -> GOI 
as.data.frame(GOI) -> GOM
sub('_.*', '', colnames(GOM)) -> colnames(GOM) 
as.character(GOM[1,]) -> new
names(new) <- colnames(GOM)
as.data.frame(new) -> new2
new[1:24] -> new 
new2$x2 <- names(new) 
colnames(new2) <- c("value", "species")


sp_tree <- read.tree("C:/Users/AG Xu_IEB_EG3_1/Desktop/Master/TreeExp2/Trees/sp/L.6sp.nwk")

as.numeric(new2$value) -> PLOT2
names(PLOT2) <- new2$species
plotTree.boxplot(sp_tree,PLOT2) 

} 
names(PLOT2) <- gsub(x = names(PLOT2), pattern = "\\peruvinaum", replacement = "peruvianum") 
names(PLOT2) <- gsub(x = names(PLOT2), pattern = "S", replacement = "S.")

as.data.frame(names(PLOT2)) -> UDP 
UDP$x2 <- PLOT2 
colnames(UDP) <- c("Species", "Value")


#ALL <- rbind(stage_g,stage_b,stage_o)
ALL_IN_ALL3 <- melt(UDP, id.var="Species")

#ALL_IN_ALL3$Species <- factor(ALL_IN_ALL3$Species , levels=c 

                              
                              
#### OPTION 1, BLACK LINES AND COLOURFULL FILL                                                            
p1 <- ggboxplot(ALL_IN_ALL3, x = "Species", y = "value", ylim =c(-4,12) ,
xlab="Species", ylab="Expression", fill = "Species", palette = c("Orange","Yellow","Green","Blue","Pink","Red"),
color = "black" ,title="Solyc05g006140.1.1")+
grids(linetype = "dashed", size=1)   
p1 + theme(legend.position = "none") 


###OPTION 2, WHITE FILL AND COLOURFUL LINES 
p2 <- ggboxplot(ALL_IN_ALL3, x = "Species", y = "value", ylim =c(-4,12) ,
                xlab="Species", ylab="Expression", fill = "White",
                color = c("Orange","Yellow","Green","Blue","Pink","Red") ,title="Solyc05g006140.1.1")+ grids(linetype = "dashed", size=1) 




p2 <- ggboxplot(ALL_IN_ALL3, x = "Species", y = "value", ylim =c(-4,12) ,
                xlab="Species", ylab="Expression", fill = "Species", palette = "Rainbow",
                
                color = "black" ,title="Solyc12g088180.1.1")+ 
  theme(axis.text.x=element_blank())+
  grids(linetype = "dashed", size=1)+   
  theme(legend.position = "none") 






p1 + theme(legend.position = "none") 





figure <- ggarrange(p3, p2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1) 


figure                 
