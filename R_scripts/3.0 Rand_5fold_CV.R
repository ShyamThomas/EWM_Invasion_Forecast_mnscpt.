library(randomForest)
library(mgcv)
library(tidyverse)
library(caret)
library(pROC)


setwd("~/UMNpostdoc/Projects/RProjects/EWM_ClimateChange")

#### Step1. Start by creating a list of all the saved files in the new folder
Train.fileNames = list.files(path="processed_data/TrainData/",pattern=".csv")
Train.fileNames

##### Step 2. Execute 5-fold cross-validation and capture AUCs for each of the 5 RF models

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  folds = rep_len(1:5,nrow(full.df))
  sample.folds=sample(folds,nrow(full.df))
  full.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=full.df[full.df$folds==i,]
  train.data= full.df[full.df$folds !=i,]
  train.rf = randomForest(train.data[,c(2:4)], as.factor(train.data$EWMSTATUS),importance=TRUE, ntree=5000, 
                          type="regression")
  preds.test=predict(train.rf, newdata=test.data, type="prob")
  AUC=auc(roc(test.data$EWMSTATUS,as.numeric(preds.test)))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  acc=confusionMatrix(as.factor(preds.test),as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, i, acc))
  sens=confusionMatrix(as.factor(preds.test),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, i, sens))
  spec=confusionMatrix(as.factor(preds.test),as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, i, spec))
  
  }
}

### Summarize by GCM, and get rid off all the unwanted letters

MeanAUC_GCMs_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAUC_GCMs_RF$Train.fileName=sub('.WtrTemp.csv', '',MeanAUC_GCMs_RF$Train.fileName)
MeanAUC_GCMs_RF$Train.fileName=sub('EWM.train.data_', '',MeanAUC_GCMs_RF$Train.fileName)
write.table(MeanAUC_GCMs_RF,"Results/AllGCMs_5foldCV_RF_AUCs.txt", sep="\t")

MeanAcc_GCMs_RF=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))

MeanSens_GCMs_RF=Sens_all%>%group_by(Train.fileName)%>%summarise(
meanSens=mean(sens))
MeanSens_GCMs_RF$Train.fileName=sub('.WtrTemp.csv', '',MeanSens_GCMs_RF$Train.fileName)
MeanSens_GCMs_RF$Train.fileName=sub('EWM.train.data_', '',MeanSens_GCMs_RF$Train.fileName)
write.table(MeanSens_GCMs_RF,"Results/AllGCMs_5foldCV_RF_Sens.txt", sep="\t")

MeanSpec_GCMs_RF=Spec_all%>%group_by(Train.fileName)%>%summarise(
meanSpec=mean(spec))
MeanSpec_GCMs_RF$Train.fileName=sub('EWM.train.data_', '',MeanSpec_GCMs_RF$Train.fileName)
MeanSpec_GCMs_RF$Train.fileName=sub('.WtrTemp.csv', '',MeanSpec_GCMs_RF$Train.fileName)
write.table(MeanSpec_GCMs_RF,"Results/AllGCMs_5foldCV_RF_Spec.txt", sep="\t")

MeanAUC_GCMs_RF
MeanAcc_GCMs_RF
MeanSens_GCMs_RF
MeanSpec_GCMs_RF

#########################################################################################################

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(full.df))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ -c(1,5) ]), ',k=3)', sep = "", collapse = ' + ')   ### FOR k=3 GAMs & then k=10
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_k3 = gam(fm,data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_k3, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  preds.val=ifelse(preds.test >0.5,1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, i, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, i, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, i, spec))


  }
}

MeanAUC_GCMs_GAM.k3=AUC_all%>%group_by(Train.fileName)%>%summarise(
meanAUC=mean(AUC))
MeanAcc_GCMs_GAM.k3=Acc_all%>%group_by(Train.fileName)%>%summarise(
meanAcc=mean(acc))
MeanSens_GCMs_GAM.k3=Sens_all%>%group_by(Train.fileName)%>%summarise(
meanSens=mean(sens))
MeanSpec_GCMs_GAM.k3=Spec_all%>%group_by(Train.fileName)%>%summarise(
meanSpec=mean(spec))

MeanAUC_GCMs_GAM.k3
MeanAcc_GCMs_GAM.k3
MeanSens_GCMs_GAM.k3
MeanSpec_GCMs_GAM.k3

#########################################################################################################
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(full.df))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ -c(1,5) ]), ',k=10)', sep = "", collapse = ' + ')   ### FOR k=3 GAMs & then k=10
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_k10 = gam(fm,data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_k10, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  
  preds.val=ifelse(preds.test >0.5,1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, i, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, i, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.data$EWMSTATUS))$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, i, spec))
  
  }
}

MeanAUC_GCMs_GAM.k10=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_GCMs_GAM.k10=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_GCMs_GAM.k10=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_GCMs_GAM.k10=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_GCMs_GAM.k10
MeanAcc_GCMs_GAM.k10
MeanSens_GCMs_GAM.k10
MeanSpec_GCMs_GAM.k10

#########################################################################
GAM.k3_perf=left_join(MeanAUC_GCMs_GAM.k3, MeanSens_GCMs_GAM.k3,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_GAM.k3, by="Train.fileName")%>%left_join(MeanAcc_GCMs_GAM.k3, by="Train.fileName")
GAM.k10_perf=left_join(MeanAUC_GCMs_GAM.k10, MeanSens_GCMs_GAM.k10,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_GAM.k10, by="Train.fileName")%>%left_join(MeanAcc_GCMs_GAM.k10, by="Train.fileName")
RF_perf=left_join(MeanAUC_GCMs_RF,MeanSens_GCMs_RF,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_RF, by="Train.fileName")%>%left_join(MeanAcc_GCMs_RF, by="Train.fileName")
AllSDMs_5foldCV=bind_rows(GAM.k3_perf, GAM.k10_perf,RF_perf )%>%
                mutate(SDM=c(rep("GAM.k3",5),rep("GAM.k10",5),rep("RF",5)))

write.table(AllSDMs_5foldCV,"Results/AllSDMs_5foldCV.txt", sep="\t")

