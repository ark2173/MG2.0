GO_Search <-function(Seurat,DEG){
  library(topGO)
  Glia<-Seurat
  AllGenes_Names <- rownames(Glia)
  expressed.genes <- rownames(DEG)
  
  # define geneList as 1 if gene is in expressed.genes, 0 otherwise
  geneList <- ifelse(AllGenes_Names %in% expressed.genes, 1, 0)
  names(geneList) <- AllGenes_Names
  
  # Create topGOdata object
  GOdata_BP <- new("topGOdata",
                   ontology = "BP", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), 
                   annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
  pval_BP<-as.data.frame(score(resultFisher_BP))
  colnames(pval_BP)<-"PVal"
  pval_BP$GO.ID<-rownames(pval_BP)
  
  GO_Fisher_BP<-GenTable(GOdata_BP, Fisher = resultFisher_BP, topNodes = 1000, numChar = 1000)
  GO_Fisher_BP$Group<-"Biological Process"
  
  
  # Create topGOdata object
  GOdata_MF <- new("topGOdata",
                   ontology = "MF", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
  pval_MF<-as.data.frame(score(resultFisher_MF))
  colnames(pval_MF)<-"PVal"
  pval_MF$GO.ID<-rownames(pval_MF)
  GO_Fisher_MF<-GenTable(GOdata_MF, Fisher = resultFisher_MF, topNodes = 1000, numChar = 1000)
  GO_Fisher_MF$Group<-"Molecular Function"
  
  # Create topGOdata object
  GOdata_CC <- new("topGOdata",
                   ontology = "CC", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")
  pval_CC<-as.data.frame(score(resultFisher_CC))
  colnames(pval_CC)<-"PVal"
  pval_CC$GO.ID<-rownames(pval_CC)
  GO_Fisher_CC<-GenTable(GOdata_CC, Fisher = resultFisher_CC, topNodes = 500, numChar = 1000)
  GO_Fisher_CC$Group<-"Cellular Component"
  
  
  Fisher_Results<-rbind(GO_Fisher_BP,GO_Fisher_MF,GO_Fisher_CC)
  Fisher_Results<-left_join(Fisher_Results,pval_MF,by="GO.ID")
  Fisher_Results<-left_join(Fisher_Results,pval_BP,by="GO.ID")
  Fisher_Results<-left_join(Fisher_Results,pval_CC,by="GO.ID")
  Fisher_Results[is.na(Fisher_Results)]<-0
  Fisher_Results<-Fisher_Results %>%
    mutate(PVal_Fis=PVal.x+PVal.y+PVal)
  
  Fisher_Results<-Fisher_Results[,c("GO.ID","Term","Annotated","Significant","Group","PVal_Fis")]
  
  ###KS Stat###
  # define geneList as 1 if gene is in expressed.genes, 0 otherwise
  geneList <- ifelse(AllGenes_Names %in% expressed.genes, DEG$p_val_adj, 0)
  names(geneList) <- AllGenes_Names
  
  # Create topGOdata object
  GOdata_BP <- new("topGOdata",
                   ontology = "BP", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultKS_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "ks")
  pval_BP<-as.data.frame(score(resultKS_BP))
  colnames(pval_BP)<-"PVal"
  pval_BP$GO.ID<-rownames(pval_BP)
  
  GO_KS_BP<-GenTable(GOdata_BP, KS = resultKS_BP, topNodes = 1000, numChar = 1000)
  GO_KS_BP$Group<-"Biological Process"
  
  
  # Create topGOdata object
  GOdata_MF <- new("topGOdata",
                   ontology = "MF", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultKS_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "ks")
  pval_MF<-as.data.frame(score(resultKS_MF))
  colnames(pval_MF)<-"PVal"
  pval_MF$GO.ID<-rownames(pval_MF)
  
  GO_KS_MF<-GenTable(GOdata_MF, KS = resultKS_MF, topNodes = 1000, numChar = 1000)
  GO_KS_MF$Group<-"Molecular Function"
  
  # Create topGOdata object
  GOdata_CC <- new("topGOdata",
                   ontology = "CC", # use biological process ontology
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  
  resultKS_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "ks")
  pval_CC<-as.data.frame(score(resultKS_CC))
  colnames(pval_CC)<-"PVal"
  pval_CC$GO.ID<-rownames(pval_CC)
  
  GO_KS_CC<-GenTable(GOdata_CC, KS = resultKS_CC, topNodes = 500, numChar = 1000)
  GO_KS_CC$Group<-"Cellular Component"
  
  Results_KS<-rbind(GO_KS_BP,GO_KS_MF,GO_KS_CC)
  Results_KS<-left_join(Results_KS,pval_MF,by="GO.ID")
  Results_KS<-left_join(Results_KS,pval_BP,by="GO.ID")
  Results_KS<-left_join(Results_KS,pval_CC,by="GO.ID")
  Results_KS[is.na(Results_KS)]<-0
  Results_KS<-Results_KS %>%
    mutate(PVal_KS=PVal.x+PVal.y+PVal)
  
  Results_KS<-Results_KS[,c("GO.ID","Term","Annotated","Group","PVal_KS")]
  Fisher_Results<-Fisher_Results[,c("GO.ID","PVal_Fis","Significant")]
  
  tokeep<-intersect(Results_KS$GO.ID,Fisher_Results$GO.ID)
  Results_KS<-Results_KS[Results_KS$GO.ID %in% tokeep,]
  Fisher_Results<-Fisher_Results[Fisher_Results$GO.ID %in% tokeep,]
  Results<-left_join(Fisher_Results,Results_KS,by="GO.ID")
  
  return(Results)
}