
library(precrec)
library(maptools)
library(spatstat)
library(blockCV)
library(sf)
library(tidyverse)


Minn.sf=read_sf(dsn="Data/", layer="Minn.map")
plot(Minn.sf$geometry)
Minn.sf
st_crs(Minn.sf)

EWM.GCMs.data=read_csv("Data/EWM.prsabs95to15_AllGCMs_v2.csv")
EWM.GCMs.data

EWM.GCMs.sf=st_as_sf(EWM.GCMs.data, coords=c("LON","LAT"), crs=32615)
EWM.GCMs.sf
st_crs(EWM.GCMs.sf)

Minn.sf=st_transform(Minn.sf, crs=32615)
ggplot(Minn.sf)+geom_sf()+
  geom_sf(data=EWM.GCMs.sf)


################ SPATIALLY BLOCKED CV ALONG LATITUDINAL GRADIENT
CheqBrd.trial.blocks = spatialBlock(speciesData = EWM.GCMs.sf, # sf or SpatialPoints
                                 species = "EWMSTATUS", # the response column (binomial or multi-class)
                                 #theRange = 100000, # size of the blocks in meters
                                 rows = 4, # number of folds
                                 selection = "checkerboard",
                                 iteration = 100, # find evenly dispersed folds
                                 biomod2Format = FALSE)

CB=CheqBrd.trial.blocks$plots+theme(panel.grid.major = element_blank())+
  geom_sf(data = Minn.sf, fill=NA, colour="gray", lwd=1)+
  geom_sf(data = EWM.GCMs.sf, alpha = 0.25, aes(colour=as.factor(EWMSTATUS)))+theme(legend.title=element_blank())
CB+theme_minimal()+theme(legend.title=element_blank())+ggtitle("Spatial blocks", subtitle = "EWM distribution")

p=ggplot()+geom_sf(data = Minn.sf, fill=NA, colour="gray", lwd=1)+geom_sf(data = EWM.GCMs.sf, alpha = 0.5, aes(colour=as.factor(EWMSTATUS)))+theme(legend.title=element_blank())
p+theme_minimal()+theme(legend.title=element_blank())+ggtitle("EWM distribution")


############## START WITH RANDOM FOREST MODELS
Train.fileNames = list.files(path="processed_data/TrainData/",pattern=".csv")
Train.fileNames


AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL


for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  full.df$EWMSTATUS=as.factor(full.df$EWMSTATUS)

  folds=CheqBrd.trial.blocks$folds

  testTable <- full.df

  testTable$RFpreds <- NA

  for(k in seq_len(length(folds))){
    # extracting the training and testing indices
    # this way works with folds list (but not foldID)
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    rf <- randomForest(EWMSTATUS~.,full.df[trainSet,], ntree = 500,importance=TRUE, type="regression") # model fitting on training set
    testTable$RFpreds[testSet] <- predict(rf, newdata = full.df[testSet, -1], type="prob")[,2] # predict the test set
   
    precrec_obj <- evalmod(scores = testTable$RFpreds[testSet], labels = testTable$EWMSTATUS[testSet], mode="aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$RFpreds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),testTable$EWMSTATUS[testSet])$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))

  }
}


AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_SpBlk_RF=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_RF=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_RF=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_RF=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_RF
MeanAcc_SpBlk_RF
MeanSens_SpBlk_RF
MeanSpec_SpBlk_RF

write.table(MeanAUC_SpBlk_RF,"Results/AllGCMs_SpatialBlockCV_RF_AUCs.txt", sep="\t")

############## NOW WITH GAM, K=10

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  #sub.df=full.df
  
  folds=CheqBrd.trial.blocks$folds
  testTable <- full.df
  testTable$GAMpreds <- NA

  for(k in seq_len(length(folds))){
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    sample_sub=full.df[trainSet, ]
    fm <- paste('s(', names(sample_sub[ -1 ]), ',k=10)', sep = "", collapse = ' + ')
    fm <- as.formula(paste('EWMSTATUS ~', fm))
    GAM_k10 = gam(fm,data=sample_sub, method="REML", family = "binomial")
    testTable$GAMpreds[testSet] <- predict(GAM_k10, full.df[testSet, ], type = "response") # predict the test set
    
    precrec_obj <- evalmod(scores = testTable$GAMpreds[testSet], labels = testTable$EWMSTATUS[testSet],mode="aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$GAMpreds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))
  }
}

AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_SpBlk_GAM.k10=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_GAM.k10=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_GAM.k10=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_GAM.k10=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_GAM.k10
MeanAcc_SpBlk_GAM.k10
MeanSens_SpBlk_GAM.k10
MeanSpec_SpBlk_GAM.k10

write.table(MeanAUC_SpBlk_GAM.k10,"Results/AllGCMs_SpatialBlockCV_GAM.k10_AUCs.txt", sep="\t")

############## NOW WITH GAM, K=3

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  #sub.df=full.df[,c(1,5,11,12)]
  
  folds=CheqBrd.trial.blocks$folds
  testTable <- full.df
  testTable$GAMk3preds <- NA
  
  for(k in seq_len(length(folds))){
    trainSet <- unlist(folds[[k]][1]) # training set indices
    testSet <- unlist(folds[[k]][2]) # testing set indices
    sample_sub=full.df[trainSet, ]
    fm <- paste('s(', names(sample_sub[ -1 ]), ',k=3)', sep = "", collapse = ' + ')
    fm <- as.formula(paste('EWMSTATUS ~', fm))
    GAM_k3 = gam(fm,data=sample_sub, method="REML", family = "binomial")
    testTable$GAMk3preds[testSet] <- predict(GAM_k3, full.df[testSet, ], type = "response") # predict the test set
    
    precrec_obj <- evalmod(scores = testTable$GAMk3preds[testSet], labels = testTable$EWMSTATUS[testSet], mode = "aucroc")
    AUC=precrec_obj$uaucs$aucs
    AUC_all = rbind(AUC_all, data.frame(Train.fileName, k, AUC))
    
    preds.val=ifelse(testTable$GAMk3preds[testSet] >0.5,1,0)
    acc=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$overall[1]
    Acc_all = rbind(Acc_all, data.frame(Train.fileName, k, acc))
    sens=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[1]
    Sens_all = rbind(Sens_all, data.frame(Train.fileName, k, sens))
    spec=confusionMatrix(as.factor(preds.val),as.factor(testTable$EWMSTATUS[testSet]))$byClass[2]
    Spec_all = rbind(Spec_all, data.frame(Train.fileName, k, spec))
  }
}

AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)
AUC_all

MeanAUC_SpBlk_GAM.k3=AUC_all%>%group_by(Train.fileName)%>%summarise(
  meanAUC=mean(AUC))
MeanAcc_SpBlk_GAM.k3=Acc_all%>%group_by(Train.fileName)%>%summarise(
  meanAcc=mean(acc))
MeanSens_SpBlk_GAM.k3=Sens_all%>%group_by(Train.fileName)%>%summarise(
  meanSens=mean(sens))
MeanSpec_SpBlk_GAM.k3=Spec_all%>%group_by(Train.fileName)%>%summarise(
  meanSpec=mean(spec))

MeanAUC_SpBlk_GAM.k3
MeanAcc_SpBlk_GAM.k3
MeanSens_SpBlk_GAM.k3
MeanSpec_SpBlk_GAM.k3

GAM.k3_SpBlk_perf=left_join(MeanAUC_SpBlk_GAM.k3,MeanSens_SpBlk_GAM.k3,by="Train.fileName")%>%
  left_join(MeanSpec_SpBlk_GAM.k3, by="Train.fileName")%>%left_join(MeanAcc_SpBlk_GAM.k3, by="Train.fileName")
GAM.k10_SpBlk_perf=left_join(MeanAUC_SpBlk_GAM.k10,MeanSens_SpBlk_GAM.k10,by="Train.fileName")%>%
  left_join(MeanSpec_SpBlk_GAM.k10, by="Train.fileName")%>%left_join(MeanAcc_SpBlk_GAM.k10, by="Train.fileName")
RF_SpBlk_perf=left_join(MeanAUC_GCMs_RF,MeanSens_GCMs_RF,by="Train.fileName")%>%
  left_join(MeanSpec_GCMs_RF, by="Train.fileName")%>%left_join(MeanAcc_SpBlk_RF, by="Train.fileName")
AllSDMs_SpBlk_CV=bind_rows(GAM.k3_SpBlk_perf, GAM.k10_SpBlk_perf,RF_SpBlk_perf )%>%
  mutate(SDM=c(rep("GAM.k3",5),rep("GAM.k10",5),rep("RF",5)))

write.table(AllSDMs_SpBlk_CV,"Results/AllSDMs_SpBlk_CV.txt", sep="\t")
###############################################################################################################################################
############################################## INDEPENDENT VALIDATION, 90% AND HIGHER GDD VALUES BLOCKED #######################################
################################################################################################################################################

Train.fileNames = list.files(path="processed_data/TrainData/",pattern=".csv")
Train.fileNames

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9) ### set the threshold at 90th percentile
  test.df=sub.df%>%filter(.[[4]]>q) ### test data containing lakes with upper 10th percent temperatures
  train.df=sub.df%>%filter(.[[4]]<q) ### train data all but
  
rf = randomForest(as.factor(EWMSTATUS)~., train.df[,-1 ], ntree = 500, data=train.df, keep.forest=TRUE)
test.df$preds=NULL
test.df$preds=predict(rf, test.df[,-1], type = "prob")[,2]

auc_obj=evalmod(scores=test.df$preds,labels=test.df$EWMSTATUS, mode="aucroc")
AUC=auc_obj$uaucs$aucs
AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 

preds.val=ifelse(test.df$preds > mean(test.df$preds),1,0)
acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$overall[1]
Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$byClass[1]
Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS))$byClass[2]
Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))

}

AUC_all$Train.fileName=sub('EWM.train.data_', '',AUC_all$Train.fileName)
AUC_all
AUC_all$Train.fileName=sub('.WtrTemp.csv', '',AUC_all$Train.fileName)

AUC_RF_all=AUC_all
Acc_RF_all=Acc_all
Sens_RF_all=Sens_all
Spec_RF_all=Spec_all

AUC_RF_all
Acc_RF_all
Sens_RF_all
Spec_RF_all


################## NOW WITH GAM, K=10

AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ -1 ]), ',k=10)', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_k10 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_k10, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}

AUC_GAM.k10_all=AUC_all
Acc_GAM.k10_all=Acc_all
Sens_GAM.k10_all=Sens_all
Spec_GAM.k10_all=Spec_all

AUC_GAM.k10_all
Acc_GAM.k10_all
Sens_GAM.k10_all
Spec_GAM.k10_all
  
################## NOW WITH GAM, K=3
AUC_all=NULL
Acc_all=NULL
Sens_all=NULL
Spec_all=NULL

for(Train.fileName in Train.fileNames) {
  full.df = read.csv(paste("processed_data/TrainData/",Train.fileName, sep=""))
  sub.df=full.df[,c(1:4)]
  q=quantile(sub.df[,4], probs=0.9)
  test.df=sub.df%>%filter(.[[4]]>q)
  train.df=sub.df%>%filter(.[[4]]<q)
  
  fm <- paste('s(', names(train.df[ -1 ]), ',k=3)', sep = "", collapse = ' + ')
  fm <- as.formula(paste('EWMSTATUS ~', fm))
  GAM_k3 = gam(fm,data=train.df, method="REML", family = "binomial")
  test.df$GAMpreds=NULL
  test.df$GAMpreds <- predict(GAM_k3, test.df[, -1], type = "response") # predict the test set
  
  auc_obj=evalmod(scores=test.df$GAMpreds,labels=test.df$EWMSTATUS, mode = "aucroc")
  AUC=auc_obj$uaucs$aucs
  AUC_all = rbind(AUC_all, data.frame(Train.fileName, AUC)) 
  
  preds.val=ifelse(test.df$GAMpreds >mean(test.df$GAMpreds),1,0)
  acc=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$overall[1]
  Acc_all = rbind(Acc_all, data.frame(Train.fileName, acc))
  sens=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[1]
  Sens_all = rbind(Sens_all, data.frame(Train.fileName, sens))
  spec=confusionMatrix(as.factor(preds.val),as.factor(test.df$EWMSTATUS), mode="sens_spec")$byClass[2]
  Spec_all = rbind(Spec_all, data.frame(Train.fileName, spec))
  
}



AUC_GAM.k3_all=AUC_all
Acc_GAM.k3_all=Acc_all
Sens_GAM.k3_all=Sens_all
Spec_GAM.k3_all=Spec_all

AUC_GAM.k3_all
Acc_GAM.k3_all
Sens_GAM.k3_all
Spec_GAM.k3_all



GAM.k3_Warmest_perf=left_join(AUC_GAM.k3_all,Sens_GAM.k3_all,by="Train.fileName")%>%
  left_join(Spec_GAM.k3_all, by="Train.fileName")%>%left_join(Acc_GAM.k3_all, by="Train.fileName")
GAM.k10_Warmest_perf=left_join(AUC_GAM.k10_all,Sens_GAM.k10_all,by="Train.fileName")%>%
  left_join(Spec_GAM.k10_all, by="Train.fileName")%>%left_join(Acc_GAM.k10_all, by="Train.fileName")
RF_Warmest_perf=left_join(AUC_RF_all,Sens_RF_all,by="Train.fileName")%>%
  left_join(Spec_RF_all, by="Train.fileName")%>%left_join(Acc_RF_all, by="Train.fileName")

AllSDMs_Warmest_CV=bind_rows(GAM.k3_Warmest_perf, GAM.k10_Warmest_perf,RF_Warmest_perf )%>%
  mutate(SDM=c(rep("GAM.k3",5),rep("GAM.k10",5),rep("RF",5)))

write.table(AllSDMs_Warmest_CV,"Results/AllSDMs_Warmest_CV.txt", sep="\t")

##############################################################################################################################
AUC.fileNames = list.files(pattern="CV.txt")
result = lapply(AUC.fileNames, function(x) read.table(x,header = TRUE))
result
names(result[[3]])[2:5]=c("meanAUC", "meanSens", "meanSpec","meanAcc")

AllAUCs_combined=do.call(rbind, result)
AllAUCs_combined

Val.Results=AllAUCs_combined%>%mutate(SDM=recode(SDM, GAM.k10 = "GAM (k=10)", GAM.k3 = "GAM (k=03)"))%>%
              mutate(Validation = c(rep("Random 5-fold", 15),rep("Spatial-block", 15),rep("Warmest 10%", 15)))%>%
                  group_by(Validation,SDM)%>%summarise(
                  avgAUC=round(mean(meanAUC),2),
                  avgAcc=round(mean(meanAcc),2),
                  avgSens=round(mean(meanSens),2),
                  avgSpec=round(mean(meanSpec),2))

library(flextable)
flextable(Val.Results) %>%
theme_vanilla()%>%save_as_pptx(path="Val.Results_wAcc.table.pptx")

############################################################################################################################
############################################################################################################################
