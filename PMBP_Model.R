PMBP_Model <- function(InputTrain, ObservTrain, InputPred, GeneVec, Family){

  if(is.na(Family)){
    aa <- cv.glmnet(x = InputTrain, y = as.numeric(ObservTrain[,1]), nfolds = max(5, floor(nrow(InputTrain)/50)))
  }else{
    aa <- cbind(ObservTrain[,1], ObservTrain[,2])
    colnames(aa) <- c("time", "status")
    aa <- cv.glmnet(x = InputTrain, y = aa, nfolds = max(5, floor(nrow(InputTrain)/50)), family = Family)
  }
  PredVal <- predict(aa, InputPred, s="lambda.min")

  PosGenes <- GeneVec[which(coef(aa)[-1] > 0)]
  NegGenes <- GeneVec[which(coef(aa)[-1] < 0)]

  PredList <- list(NegGenes, PosGenes, PredVal)
  names(PredList) <- c("NegGenes", "PosGenes", "PredVal")
  return(PredList)
}


