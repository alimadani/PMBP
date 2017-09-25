PMBP_Target_Sampling <- function(GO_List, ExpMat, GeneVec, GeneID, ResponseVec,
                                 minLength, maxLength, CVnum, Targets, Family, FDR, SamplingNum, ParentGO){
  ##############
  PosGenes <- c()
  NegGenes <- c()
  PosGenesCor <- c()
  NegGenesCor <- c()
  
  PvalVec_pos <- c()
  PvalVec_neg <- c()
  ###################
  GO_Association <- PMBP_GO_Target(GoTermsIDs_GeneSet, ExpMat,
                                   GeneVec, GeneID,
                                   ResponseVec, minLength, maxLength,
                                   SamplingNum, CVnum, Targets, Family, ParentGO)
  #################
    SigVec <- GO_Association$Pval
  ####################
  for(GOiter in c(1:length(GO_Association$GOterms))){
    ###########
    PosGenes <- c(PosGenes, names(table(GO_Association$PosGenes[[GOiter]])))
    PosGenesCor <- c(PosGenesCor, as.numeric(table(GO_Association$PosGenes[[GOiter]])
                                             )*(-log10(SigVec[GOiter]))/(CVnum*SamplingNum))
    
    NegGenes <- c(NegGenes, names(table(GO_Association$NegGenes[[GOiter]])))
    NegGenesCor <- c(NegGenesCor, as.numeric(table(GO_Association$NegGenes[[GOiter]])
                                             )*(-log10(SigVec[GOiter]))/(CVnum*SamplingNum))
    #############
    PvalVec_pos <- c(PvalVec_pos, rep((-log10(SigVec[GOiter])),
                                      length(names(table(GO_Association$PosGenes[[GOiter]])))))
    PvalVec_neg <- c(PvalVec_neg, rep((-log10(SigVec[GOiter])),
                                      length(names(table(GO_Association$NegGenes[[GOiter]])))))
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

  GeneCor_Aux <- GeneCor
  GeneCor_Aux[which(is.na(GeneCor_Aux))] <- 0
  # GeneScore <- (GeneCor - median(na.omit(GeneCor_Aux)))/mad(na.omit(GeneCor_Aux))
  GeneScore <- GeneCor*sqrt(nrow(ResponseVec)-2)/sqrt(1-GeneCor^2)
  GeneSig <- p.adjust(2*pnorm(-abs(GeneScore)), method = "fdr", n = length(GeneScore))
  
  if(FDR){
    SigVec <- p.adjust(GO_Association$Pval, method = "fdr",
                       n = length(GO_Association$Pval))
  }else{
    SigVec <- GO_Association$Pval
  }
 
  SigList <- list(GO_Association$GOterms, SigVec,
                  UniqueGenes, GeneCor, GeneSig,
                  GO_Association$AssociatedGenes)
 
  names(SigList) <- c("GO_terms", "GO_significance", "genes",
                      "genes_correlation", "genes_significance",
                      "AssociatedGenes")
  return(SigList)
}
