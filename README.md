# Tomato-Transcriptomics: 
The  repository for my bachelor thesis on the topic of "The Evolution of Gene Expression in Tomato Flower".

Contains scripts to estimate the strength of expression conservation using TreeExp2 (https://github.com/jingwyang/TreeExp) and to estimate expression divergence using the phylogenetic ANOVA (https://github.com/Jmendo12/evemodel). 

# Abstract: 
In this study, we apply two comparative methods based on the OU model to find genes of adaptive significance in tomato flower. It was ascertained that gene expression appears to be the most conserved in Open Flowers compared to prior developmental stages. We used a phylogenetic ANOVA to investigate genes that might have been subject to divergent selection. We found 249 genes as candidates for expression level adaptation, one of which en-codes for a UDP-Glycosyltransferase which might have been downregulated by domestica-tion in S. lycopersicum

# Results: 
Accession Tree from RAxML-NG (Input TreeExp2):

![ZlsNW8JgEBN4es43w39AYA](https://user-images.githubusercontent.com/77416397/115562054-6b20ee80-a2b6-11eb-9d4d-c61e7f201eef.png)


Species Tree (Input phylogenetic ANOVA): 

![zYgBm9dJqvxFqu6aSPDhEw](https://user-images.githubusercontent.com/77416397/115561788-314fe800-a2b6-11eb-9768-7cdbdeff291f.png)


Estimation of the Strength of Expression Conservation: 

![Latest_f](https://user-images.githubusercontent.com/77416397/115562342-ad4a3000-a2b6-11eb-872e-5981dd3f5a91.png)


Phylogenetic ANOVA and multiple testing: 

![image](https://user-images.githubusercontent.com/77416397/115563378-c1daf800-a2b7-11eb-8358-e6c0c292d7cd.png)
(genes with higher variance among than within lineages, FDR < 0.10, p <  0.05)


Pattern recognized through manual checking: 

![image](https://user-images.githubusercontent.com/77416397/115562894-3b261b00-a2b7-11eb-81d1-4662ad18e9de.png)


UDP-Glycosyltransferase expression: 

![ANOTHER30](https://user-images.githubusercontent.com/77416397/115562986-58f38000-a2b7-11eb-8b82-5157a9dfd843.png)


