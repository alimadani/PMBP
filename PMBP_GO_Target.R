PMBP_GO_Target <- function(GO_List, ExpMat, GeneVec, GeneID, ResponseVec, minLength, maxLength, SamplingNum, CVnum, Targets, Family, ParentGO){
  
  ################
  GoTerms_checkedPathways <- GO_List[[1]]
  GoIDs_checkedPathways <- GO_List[[2]]
  GeneMatchedIDs_Go <- GO_List[[3]]
  GeneMatchedSymbols_Go <- GO_List[[4]]
  
  if(ParentGO == FALSE){
    Genes_GOs_aux <- c()
    GO_idVec_aux <- c()
    for(GOiter in 1:length(GoTerms_checkedPathways)){
      GeneVec_aux <- GeneMatchedIDs_Go[[GOiter]]
      if(length(GeneVec_aux) <= maxLength & length(GeneVec_aux) >= minLength){
        Genes_GOs_aux <- qpcR:::cbind.na(Genes_GOs_aux, GeneVec_aux)
        GO_idVec_aux <- c(GO_idVec_aux, GOiter)
      }
    }
    Genes_GOs_aux <- Genes_GOs_aux[,-1]
    
    RemInd <- c()
    for(GOiter in 1:ncol(Genes_GOs_aux)){ #
      GeneVec_aux <- na.omit(Genes_GOs_aux[,GOiter])
      TargetGO_length_aux <- length(GeneVec_aux)
      InterLength_aux <- as.numeric(apply(Genes_GOs_aux, 2, function(x){length(intersect(GeneVec_aux, x))}))
      TotalLength_aux <- as.numeric(apply(Genes_GOs_aux, 2, function(x){length(na.omit(x))}))
      RemInd <- c(RemInd, which(
        TotalLength_aux > TargetGO_length_aux & InterLength_aux == TargetGO_length_aux))
    }
    RemInd <- unique(RemInd)
    GOIterVec <- c(1:length(GeneMatchedIDs_Go))[-GO_idVec_aux[RemInd]]
  }else{
    GOIterVec <- c(1:length(GeneMatchedIDs_Go))
  }
  GOiter <- 0
  
  ##############
  if(!is.na(Targets)){
    Target_ind <- which(GeneVec %in% Targets)
  }
  ############
  CorVec <- c()
  PvalVec <- c()
  GOterms <- c()
  GeneNum <- c()
  PosGenes <- list()
  NegGenes <- list()
  #######
  for(GOiter in GOIterVec){
    print(paste("GO_", GOiter, sep = "", collapse = ""))
    ######
    #     if(GeneID == "Symbol"){
    #       GOind <- GeneMatchedSymbols_Go[[GOiter]]
    #     }else if(GeneID == "EntrezID"){
    #       GOind <- GeneMatchedIDs_Go[[GOiter]]
    #     }
    if(GeneID == "Symbol"){
      GOindMacthed <- which(GeneVec %in% intersect(GeneVec, GeneMatchedSymbols_Go[[GOiter]]))
    }else if(GeneID == "EntrezID"){
      GOindMacthed <- which(GeneVec %in% intersect(GeneVec, GeneMatchedIDs_Go[[GOiter]]))
    }
    ########
    if(length(GOindMacthed) <= maxLength & length(GOindMacthed) >= minLength){
      
      if(!is.na(Targets)){
        aa <- PMBP_Discovery(ExpMat, ResponseVec, c(GOindMacthed, Target_ind), SamplingNum, CVnum, Family)
      }else{
        aa <- PMBP_Discovery(ExpMat, ResponseVec, GOindMacthed, SamplingNum, CVnum, Family)
      }
      
      GOterms <- c(GOterms, GoTerms_checkedPathways[GOiter])
      PosGenes[[GoTerms_checkedPathways[[GOiter]]]] <- aa$PosGenes
      NegGenes[[GoTerms_checkedPathways[[GOiter]]]] <- aa$NegGenes
      
      UniquCorGenes <- unique(c(na.omit(c(aa$PosGenes)), na.omit(c(aa$NegGenes))))
      
      if(length(UniquCorGenes) > 1){
        CorVec <- c(CorVec, aa$Cor)
        PvalVec <- c(PvalVec, aa$CorPval)
      }else{
        CorVec <- c(CorVec, NA)
        PvalVec <- c(PvalVec, NA)
      }
      GeneNum <- c(GeneNum, length(UniquCorGenes))
    }
  }
  ###############
  GO_Assoc <- list(GOterms, CorVec, PvalVec,
                   PosGenes, NegGenes, GeneNum)
  names(GO_Assoc) <- c("GOterms", "Cor", "Pval",
                       "PosGenes", "NegGenes", "AssociatedGenes")
  ###############
  return(GO_Assoc)
}
