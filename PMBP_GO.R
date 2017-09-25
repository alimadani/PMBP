PMBP_GO <- function(GO_List, ExpMat, GeneVec, GeneID, ResponseVec, minLength, maxLength, CVnum){
  
  ################
  GoTerms_checkedPathways <- GO_List[[1]]
  GoIDs_checkedPathways <- GO_List[[2]]
  GeneMatchedIDs_Go <- GO_List[[3]]
  GeneMatchedSymbols_Go <- GO_List[[4]]
  ##############
  CorVec <- c()
  PvalVec <- c()
  GOterms <- c()
  PosGenes <- list()
  NegGenes <- list()
  #######
  for(GOiter in 1:length(GeneMatchedIDs_Go)){
    print(GOiter)
    ########
    if(GeneID == "Symbol"){
      GOind <- which(GeneVec %in% GeneMatchedSymbols_Go[[GOiter]])
    }else if(GeneID == "EntrezID"){
      GOind <- which(GeneVec %in% GeneMatchedIDs_Go[[GOiter]])
    }
    ########
    if(length(GOind) <= maxLength & length(GOind) >= minLength){
      aa <- PMBP_Discovery(ExpMat, ResponseVec, GOind, CVnum)
      
      GOterms <- c(GOterms, GoTerms_checkedPathways[GOiter])
      CorVec <- c(CorVec, aa$Cor)
      PvalVec <- c(PvalVec, aa$CorPval)
      PosGenes[[GoTerms_checkedPathways[[GOiter]]]] <- aa$PosGenes
      NegGenes[[GoTerms_checkedPathways[[GOiter]]]] <- aa$NegGenes
    }
  }
  ###############
  GO_Assoc <- list(GOterms, CorVec, PvalVec,
                   PosGenes, NegGenes)
  names(GO_Assoc) <- c("GOterms", "Cor", "Pval",
                       "PosGenes", "NegGenes")
  ###############
  return(GO_Assoc)
}