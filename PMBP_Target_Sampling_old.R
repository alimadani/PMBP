PMBP_Target_Sampling <- function(GO_List, ExpMat, GeneVec, GeneID, ResponseVec,
                           minLength, maxLength, CVnum, Targets, Family, FDR, SamplingNum){
  ##############
  PosGenes <- c()
  NegGenes <- c()
  PosGenesCor <- c()
  NegGenesCor <- c()
  
  PvalVec_pos <- c()
  PvalVec_neg <- c()
  ################
  SigMat <- c()
  GeneNumMat <- c()
  ##############
  for(SamplingIter in 1:SamplingNum){
    print(SamplingIter)
    ###################
    Sampling <- sample(1:ncol(ExpMat), ncol(ExpMat), replace = F)
    ExpMat_Sampled <- ExpMat[,Sampling]
    ResponseVec_Sampled <- ResponseVec[Sampling,]
    
    GO_Association <- PMBP_GO_Target(GoTermsIDs_GeneSet, ExpMat,
                                     HUGOVec, GeneID,
                                     ResponseVec, minLength, maxLength,
                                     SamplingNum, CVnum, Targets, Family)
    #################
    if(FDR){
      SigVec <- p.adjust(GO_Association$Pval, method = "fdr",
                         n = length(GO_Association$Pval))
    }else{
      SigVec <- GO_Association$Pval
    }
    
    SigMat <- cbind(SigMat, SigVec)
    GeneNumMat <- cbind(GeneNumMat, GO_Association$AssociatedGenes)
    ####################
    for(GOiter in c(1:length(GO_Association$GOterms))){
      ###########
      PosGenes <- c(PosGenes, names(table(GO_Association$PosGenes[[GOiter]])))
      PosGenesCor <- c(PosGenesCor, as.numeric(table(GO_Association$PosGenes[[GOiter]]))*(-log10(SigVec[GOiter]))/CVnum)
      
      NegGenes <- c(NegGenes, names(table(GO_Association$NegGenes[[GOiter]])))
      NegGenesCor <- c(NegGenesCor, as.numeric(table(GO_Association$NegGenes[[GOiter]]))*(-log10(SigVec[GOiter]))/CVnum)
      #############
      PvalVec_pos <- c(PvalVec_pos, rep((-log10(SigVec[GOiter])),
                                        length(names(table(GO_Association$PosGenes[[GOiter]])))))
      PvalVec_neg <- c(PvalVec_neg, rep((-log10(SigVec[GOiter])),
                                        length(names(table(GO_Association$NegGenes[[GOiter]])))))
    }
    
    
  }
  
  UniquePos <- unique(PosGenes)
  UniqueNeg <- unique(NegGenes)
  UniqueGenes <- unique(c(UniquePos, UniqueNeg))
  GeneCor <- c()
  for(GeneIter in 1:length(UniqueGenes)){
    PosInd <- which(PosGenes == UniqueGenes[GeneIter])
    NegInd <- which(NegGenes == UniqueGenes[GeneIter])
    
    if(length(PosInd) > 0 & length(NegInd) > 0){
      GeneCor <- c(GeneCor, sum(c(PosGenesCor[PosInd], -NegGenesCor[NegInd])
      )/sum(c(PvalVec_pos[PosInd], PvalVec_neg[NegInd])))
    }else if(length(PosInd) > 0){
      GeneCor <- c(GeneCor, sum(PosGenesCor[PosInd])/sum(PvalVec_pos[PosInd]))
    }else{
      GeneCor <- c(GeneCor, sum(-NegGenesCor[NegInd])/sum(PvalVec_neg[NegInd]))
    }
  }
  
  GO_SigVec_pchisq <- apply(SigMat, 1, function(x){1 - pchisq(-2*sum(log(x)), 2*length(x))})
  GO_SigVec_min <- apply(SigMat, 1, min)
  
  GeneNum <- apply(GeneNumMat, 1, median)

  SigList <- list(GO_Association$GOterms, GO_SigVec_pchisq,
                  GO_SigVec_min, UniqueGenes, GeneCor, GeneNum)
  
  names(SigList) <- c("GO_terms", "GO_significance_phisq",
                      "GO_significance_min", "genes", "genes_correlation", "AssociatedGenes")
  return(SigList)
  
}