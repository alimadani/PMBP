PMBP_Discovery <- function(ExpMat, ResponseVec, PathGeneIds, SamplingNum, CVnum, Family){
 
  GeneVec <- rownames(ExpMat)[PathGeneIds]
  InputMatrix <- c()
  InputMatrix <- t(ExpMat[PathGeneIds,])

  NCOL <- ncol(InputMatrix)
  InputMatrix <- matrix(as.numeric(InputMatrix), ncol = NCOL)

  ObservedVec <- ResponseVec
  rownames(ObservedVec) <- NULL
  colnames(ObservedVec) <- NULL
  colnames(InputMatrix) <- NULL
  rownames(InputMatrix) <- NULL
 
   RemovInd <- which(is.na(InputMatrix), arr.ind = T)[,1]
   if(length(RemovInd) > 0){
     InputMatrix <- InputMatrix[-RemovInd,]
     ObservedVec <- ObservedVec[-RemovInd,]
   }
 
  if(ceiling(length(ObservedVec)/CVnum) > 1){   
   CorList <- PMBP_CV(ObservedVec, InputMatrix,
                       CVnum, SamplingNum, GeneVec, Family)
    
    NegGenes <- CorList$NegativeGenes
    PosGenes <- CorList$PositiveGenes

    PredictionList <- list(CorList$Correlation, CorList$pval, 
                           NegGenes, PosGenes)
    names(PredictionList) <- c("Cor", "CorPval", 
                               "NegGenes", "PosGenes")
    return(PredictionList)
    
  }else{
    PredictionList <- list(NA, NA, NA, NA)
    names(PredictionList) <- c("CorVec", "CorPval_Vec", 
                               "NegGenes", "PosGenes")
    return(PredictionList)
  }

}
