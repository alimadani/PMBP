
PMBP_PathwayLoop <- function(ExpMat, Genes_Pathways, PathwayNames, CVnum, MinGeneNum, MaxGeneNum){

  ###################
  CorMat <- c()
  CorPval_Mat <- c()

  NegGeneList <- list()
  PosGeneList <- list()
  ##################
  for(PathwayIter in 1:length(Genes_Pathways)){
    
    PathwayGenes <- Genes_Pathways[[PathwayIter]]
    PathGeneIds <- which(rownames(ExpMat) %in% PathwayGenes)
    GeneVec <- rownames(ExpMat)[PathGeneIds]
    
    if(length(PathGeneIds) > MinGeneNum & length(PathGeneIds) < MaxGeneNum){
      
      PredictionList <- PMBP_Discovery(ExpMat, ResponseVec, PathGeneIds, CVnum)
      
      CorMat <- cbind(CorMat, PredictionList$CorVec)
      CorPval_Mat <- cbind(CorPval_Mat, PredictionList$CorPval_Vec)

      NegGeneList[[PathwayNames[PathwayIter]]] <- PredictionList$NegGenes
      PosGeneList[[PathwayNames[PathwayIter]]] <- PredictionList$PosGenes
    }else{
      CorMat <- cbind(CorMat, NA)
      CorPval_Mat <- cbind(CorPval_Mat, NA)

      NegGeneList[[PathwayNames[PathwayIter]]] <- NA
      PosGeneList[[PathwayNames[PathwayIter]]] <- NA
    }
  }
  ################
  colnames(CorMat)  <- PathwayNames
  colnames(CorMat)  <- PathwayNames
  #################
  PredPathwaysList <- list(CorMat, CorPval_Mat,
                           NegGeneList, PosGeneList)
  names(PredPathwaysList) <- list("CorMat", "CorPval_Mat",
                                  "NegGeneList", "PosGeneList")
  return(PredPathwaysList)
  
}