## Clear the workspace and import the data
options(stringsAsFactors = F)
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 7){
  cat("Uasge: Rscript PMBP_Call.R <PMBP_path> <Path> <Out_Name> <ExpMat> <GeneID> <ResponseVec> <StatusVec>
      <minLength> <maxLength> <CVnum> <Targets> <Family> <SamplingNum>\n")
  stop("Provide at least first 7 arguments")
}
###### load required libraries
require(gdata)
require(e1071)
library(Biobase)
require(xtable)
library(devtools)
library(preprocessCore)
library(rgl)
library(qpcR)
library(data.table)
library(glmnet)
library(parallel)
library(doParallel)
#######################
########################
print(args)
PMBP_path <- args[1]
Path <- args[2]
Out_Name <- args[3]
ExpMat <- readRDS(args[4])
GeneVec <- rownames(ExpMat)
GeneID <- args[5]
if(GeneID == "Symbol"){
  GeneVec <- as.character(GeneVec)
}else{
  GeneVec <- as.numeric(GeneVec)
}
ResponseVec <- as.numeric(readRDS(args[6]))
StatusVec <- as.numeric(readRDS(args[7]))
###################
if(!is.na(args[8]) & args[8] != "NA"){
  minLength <- as.numeric(args[8])
}else{
  minLength <- 10
}

if(!is.na(args[9]) & args[9] != "NA"){
  maxLength <- as.numeric(args[9])
}else{
  maxLength <- 20
}

if(!is.na(args[10]) & args[10] != "NA"){
  CVnum <- as.numeric(args[10])
}else{
  CVnum <- 10
}

if(!is.na(args[11]) & args[11] != "NA"){
  Targets <- as.numeric(args[11])
}else{
  Targets <- NA
}

if(!is.na(args[12]) & args[12] != "NA"){
  Family <- as.numeric(args[12])
}else{
  Family <- NA
}

if(!is.na(args[13]) & args[13] != "NA"){
  SamplingNum <- as.numeric(args[13])
}else{
  SamplingNum <- 1
}
###################
####################
load("/mnt/work1/users/bhklab/Users/Ali/Projects/Matching_TCGA_L1000/GoTermsIDs_GeneSet.RData")
GoTerms_checkedPathways <- GoTermsIDs_GeneSet[[1]]
GoIDs_checkedPathways <- GoTermsIDs_GeneSet[[2]]
GeneMatchedIDs_Go <- GoTermsIDs_GeneSet[[3]]
GeneMatchedSymbols_Go <- GoTermsIDs_GeneSet[[4]]
####################
#RemInd <- unique(which(is.na(ExpMat),arr.ind = T)[,2], which(is.na(ResponseVec)))
#print(RemInd)
#if(length(RemInd) > 0){
# ExpMat <- ExpMat[,-RemInd]
# ResponseVec <- ResponseVec[-RemInd]
# StatusVec <- StatusVec[-RemInd]
#}
#print(dim(ExpMat))
#print(length(ResponseVec))
###################
setwd(PMBP_path)
source("PMBP_Model.R")
source("PMBP_CV.R")
source("PMBP_Discovery.R")
source("PMBP_GO_Target.R")
source("PMBP_Target_Sampling.R")
##################
SigList <- PMBP_Target_Sampling(GoTermsIDs_GeneSet, ExpMat,
                                GeneVec, GeneID = GeneID,
                                cbind(ResponseVec, StatusVec), minLength, maxLength,
                                CVnum, Targets, Family,FDR=TRUE, SamplingNum, ParentGO = T)

saveRDS(SigList, file = paste(Path, "GO_Association_", Out_Name, "_", CVnum, "fold_",
                              SamplingNum,"sampling_GO_", minLength, "to", maxLength,
                              ".rds", sep= "", collapse = ""))

