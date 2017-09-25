
PMBP_CV <- function(ObservedVec, InputMatrix, CVnum, SamplingNum, GeneVec, Family){
  
  PredVec <- c()
  ObsVec <- c()
  PosGenes <- c()
  NegGenes <- c()
  for(SampleIter in 1:SamplingNum){
    Sampling <- sample(1:nrow(InputMatrix), nrow(InputMatrix), replace = F)
    
    InputMatrix <- InputMatrix[Sampling,]
    ObservedVec <- ObservedVec[Sampling,]
    ObsVec <- c(ObsVec, ObservedVec[,1])
   
    CVSeq <- seq(1,nrow(InputMatrix), by = floor(nrow(InputMatrix)/(CVnum)))
    if(CVSeq[length(CVSeq)] != nrow(InputMatrix)){
      CVSeq[length(CVSeq)] <- nrow(InputMatrix)
    }
    
    for(CViter in 1:(length(CVSeq)-2)){
      
      ObservTrain <- ObservedVec[-(CVSeq[CViter]:(CVSeq[CViter+1]-1)),]
      InputTrain <- InputMatrix[-(CVSeq[CViter]:(CVSeq[CViter+1]-1)),]
      InputPred <- InputMatrix[(CVSeq[CViter]:(CVSeq[CViter+1]-1)),]
      
      PredList <- PMBP_Model(InputTrain, ObservTrain, InputPred, GeneVec, Family)
      PredVec <- c(PredVec, PredList$PredVal)
      
      PosGenes <- qpcR:::cbind.na(PosGenes, PredList$PosGenes)
      NegGenes <- qpcR:::cbind.na(NegGenes, PredList$NegGenes)
    }
    ObservTrain <- ObservedVec[-(CVSeq[(length(CVSeq)-1)]:(CVSeq[length(CVSeq)])),]
    InputTrain <- InputMatrix[-(CVSeq[(length(CVSeq)-1)]:(CVSeq[length(CVSeq)])),]
    InputPred <- InputMatrix[(CVSeq[(length(CVSeq)-1)]:(CVSeq[length(CVSeq)])),]
    
    
    PredList <- PMBP_Model(InputTrain, ObservTrain, InputPred, GeneVec, Family)
    PredVec <- c(PredVec, PredList$PredVal)
    
    PosGenes <- qpcR:::cbind.na(PosGenes, PredList$PosGenes)
    NegGenes <- qpcR:::cbind.na(NegGenes, PredList$NegGenes)
  }
  
  PosGenes <- PosGenes[,-1]
  NegGenes <- NegGenes[,-1]
  
  
  if(is.na(Family)){
    CorVal <- cor(ObsVec, PredVec, method = "pearson")
    CorPval <- cor.test(ObsVec,PredVec, method = "pearson", alternative = "greater")$p.value
  }else{
    
    CorVal <- concordance.index(x=PredVec, surv.time=ObservedVec[,1],
                                surv.event=ObservedVec[,2],method = "noether",na.rm = TRUE)$c.index
    CorPval <- concordance.index(x=PredVec, surv.time=ObservedVec[,1],
                                 surv.event=ObservedVec[,2],method = "noether",na.rm = TRUE)$p.value
  }

  
  PosGenes <- matrix(PosGenes, ncol = CVnum)
  NegGenes <- matrix(NegGenes, ncol = CVnum)
  colnames(PosGenes) <- paste("CV_", seq(1,CVnum), sep = "")
  colnames(NegGenes) <- paste("CV_", seq(1,CVnum), sep = "")
  
  CorList <- list(CorVal, CorPval, PosGenes, NegGenes)
  names(CorList) <- c("Correlation", "pval", "PositiveGenes", "NegativeGenes")
  return(CorList)
}
