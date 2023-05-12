library(sf)
library(randomForest)
library(tidyverse)
library(corrplot)
library(pdp)

setwd("~/UMNpostdoc/ProjectEWM/RProjects/EWM_ClimateChange")

####################################################################################################################
######### R script: Training EWM data with Random forest (RF) & GAM SDMs, and evaluate model performance   #########
########  with random 5-fold cross-validation.                                                             #########
####################################################################################################################

#### Step 1. Load the merged water temperature and EWM data developed in the previous modules and subset all
#### 5 GCM based temperature estimates.
EWM.GCMs.data=read_csv("Data/EWM.prsabs95to15_AllGCMs_v2.csv")
EWM.GCMs.data

#### Step 1a. GCM-ACCESS estimates
EWM.train.data_ACCESS.WtrTemp=EWM.GCMs.data[,-c(1:6,9,12:15)]
EWM.train.data_ACCESS.WtrTemp
write_csv(EWM.train.data_ACCESS.WtrTemp, "processed_data/TrainData/EWM.train.data_ACCESS.WtrTemp.csv")

#### Step 1b. GCM-MIROC5 estimates
EWM.train.data_MIROC5.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11,13:15)]
EWM.train.data_MIROC5.WtrTemp
write_csv(EWM.train.data_MIROC5.WtrTemp, "processed_data/TrainData/EWM.train.data_MIROC5.WtrTemp.csv")

#### Step 1c. GCM-IPSL estimates
EWM.train.data_IPSL.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:12,14:15)]
EWM.train.data_IPSL.WtrTemp
write_csv(EWM.train.data_IPSL.WtrTemp, "processed_data/TrainData/EWM.train.data_IPSL.WtrTemp.csv")

#### Step 1d. GCM-GFDL estimates
EWM.train.data_GFDL.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:13,15)]
EWM.train.data_GFDL.WtrTemp
write_csv(EWM.train.data_GFDL.WtrTemp, "processed_data/TrainData/EWM.train.data_GFDL.WtrTemp.csv")

#### Step 1e. GCM-MRI estimates
EWM.train.data_MRI.WtrTemp=EWM.GCMs.data[,-c(1:6,9,11:14)]
EWM.train.data_MRI.WtrTemp
write_csv(EWM.train.data_MRI.WtrTemp, "processed_data/TrainData/EWM.train.data_MRI.WtrTemp.csv")

#### Step 2. Iterating RF models across all the above created data sets for 5 GCMs:
#### start by creating a list of all the saved files in the new folder
Train.fileNames = list.files(path="processed_data/TrainData/",pattern=".csv")
Train.fileNames

### A loop that reads all the files and executes random forest algorithm; saves RF obj and OOB errors as text file,
### and finally, builds and saves partial plots for each variable within each model

results_top3=NULL

for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  rf = randomForest(sample[,(2:4)], sample$EWMSTATUS,importance=TRUE, ntree=5000, type="regression")
  save(rf, file=paste("processed_data/TrainData/", sub('....$','',Train.fileName), "Top3.Rdata", sep=""))
  
  meanMSE = mean(rf$mse)
  results_top3 = rbind(results_top3, data.frame(Train.fileName, meanMSE))
  write.table(results_top3,"Results/RF_MSE_Top3.txt",sep = "\t")
  
  Top3Preds.Names=colnames(sample)[c(2:4)]
  
  for (Pred.Name in Top3Preds.Names){
    partial_plot=autoplot(partial(rf, pred.var = Pred.Name, ice=TRUE, rug=TRUE, train = sample, prob = TRUE),xlab=Pred.Name, ylab="Invasion risk", alpha=0.1)
    ggsave(filename=paste(sub('....$','',Train.fileName),"_",Pred.Name,"_IcePlot.png", sep=""), partial_plot, path="Figures/",units="in", width=9, height=6, dpi=900)
    }
}

##### Step 3. Execute 5-fold cross-validation and capture AUCs for each of the 5 RF models
library(pROC)

AUC_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  folds = rep_len(1:5,nrow(full.df))
  sample.folds=sample(folds,nrow(sample))
  full.df$folds=sample.folds
  
  set.seed(007)

  for(i in 1:5){test.data=full.df[full.df$folds==i,]
                train.data= full.df[full.df$folds !=i,]
                train.rf = randomForest(train.data[,c(2:4)], train.data$EWMSTATUS,importance=TRUE, ntree=5000, 
                                        type="regression")
                preds.test=predict(train.rf, newdata=test.data)
                AUC=auc(roc(test.data$EWMSTATUS,preds.test))
                AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
 meanAUC=mean(AUC))
MeanAUC_GCMs_RF
write.table(MeanAUC_GCMs_RF,"Results/AllGCMs_5foldCV_RF_AUCs.txt", sep="\t")

############################################################################################################
##### Step5. Repeat the above for two different GAM models, k=3 (simple GAM), and k =10 (default best-fit):
#### Step 5a. Default GAM settings: the best-fitting GAM
for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  fm <- paste('s(', names(sample[ -1 ]), ')', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam = gam(fm,data=sample, method="REML", family = "binomial")
  save(gam, file=paste("processed_data/TrainData/", sub('....$','',Train.fileName), "GAM.Rdata", sep=""))
}
load.Rdata("processed_data/TrainData/EWM.train.data_ACCESS.WtrTempGAM.Rdata", "ACCESS.GAM.model")
load.Rdata("processed_data/TrainData/EWM.train.data_GFDL.WtrTempGAM.Rdata", "GFDL.GAM.model")
load.Rdata("processed_data/TrainData/EWM.train.data_IPSL.WtrTempGAM.Rdata", "IPSL.GAM.model")
load.Rdata("processed_data/TrainData/EWM.train.data_MIROC5.WtrTempGAM.Rdata", "MIROC5.GAM.model")
load.Rdata("processed_data/TrainData/EWM.train.data_MRI.WtrTempGAM.Rdata", "MRI.GAM.model")

#### Step 5b. Simple GAM setting, k =3
for(Train.fileName in Train.fileNames) {
  sample = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  fm <- paste('s(', names(sample[ -1 ]), ',k=3)', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_k3 = gam(fm,data=sample, method="REML", family = "binomial")
  save(gam_k3, file=paste("processed_data/TrainData/", sub('....$','',Train.fileName), "GAM_k3.Rdata", sep=""))
}

load.Rdata("processed_data/TrainData/EWM.train.data_ACCESS.WtrTempGAM_k3.Rdata", "ACCESS.GAM_k3.model")
load.Rdata("processed_data/TrainData/EWM.train.data_GFDL.WtrTempGAM_k3.Rdata", "GFDL.GAM_k3.model")
load.Rdata("processed_data/TrainData/EWM.train.data_IPSL.WtrTempGAM_k3.Rdata", "IPSL.GAM_k3.model")
load.Rdata("processed_data/TrainData/EWM.train.data_MIROC5.WtrTempGAM_k3.Rdata", "MIROC5.GAM_k3.model")
load.Rdata("processed_data/TrainData/EWM.train.data_MRI.WtrTempGAM_k3.Rdata", "MRI.GAM_k3.model")

#### A quick visualization of GAM reponse curves
library(gridExtra)
library(gratia)

grid.arrange(
  draw(ACCESS.GAM_k3.model, select=3),draw(ACCESS.GAM.model, select=3),
  draw(GFDL.GAM_k3.model, select=3),draw(GFDL.GAM.model, select=3),
  draw(IPSL.GAM_k3.model, select=3),draw(IPSL.GAM.model, select=3),
  draw(MIROC5.GAM_k3.model, select=3),draw(MIROC5.GAM.model, select=3),
  draw(MRI.GAM_k3.model, select=3),draw(MRI.GAM.model, select=3),
  nrow=5)

###### Step 6. Execute 5-fold cross-validation and capture AUCs for each of the 10 different GAM models

AUC_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(sample))
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
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_k3=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC)
)
MeanAUC_GCMs_k3

write.table(MeanAUC_GCMs_k3,"Results/AllGCMs_5foldCV_GAMk3_AUCs.txt", sep="\t")

########## now for GAM, k=10
AUC_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  folds = rep_len(1:5,nrow(sub.df))
  sample.folds=sample(folds,nrow(sample))
  sub.df$folds=sample.folds
  
  set.seed(007)
  
  for(i in 1:5){test.data=sub.df[sub.df$folds==i,]
  train.data= sub.df[sub.df$folds !=i,]
  fm <- paste('s(', names(sub.df[ -c(1,5) ]),',k=10)', sep = "", collapse = ' + ')   ### FOR GAMs k=10
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  gam_k10 = gam(fm,data=train.data, method="REML", family = "binomial")
  preds.test=predict(gam_k10, newdata=test.data, type="response")
  
  AUC=auc(roc(test.data$EWMSTATUS,preds.test))
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, i, AUC))
  }
}

### Get rid off all the unwanted letters
AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_GCMs_k10=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC)
)
MeanAUC_GCMs_k10
write.table(MeanAUC_GCMs_k10,"Results/AllGCMs_5foldCV_GAMk10_AUCs.txt", sep="\t")
