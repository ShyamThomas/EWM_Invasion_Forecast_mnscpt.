library(tidyverse)
library(sf)

####################################################################################################################
####### R script: Putting together all the SDM results with lake-specific covariates, and prediction domain #######
####### These scripts were used to put together the SINGLE LARGE DATASET needed by USpatial team @ UMN      #######
####################################################################################################################

#### Get the OG EWM data
EWM.alllakes.data=read_csv("raw_data/EWM.occ_abund.data.csv")
 EWM.alllakes.data%>%filter(EWMSTATUS!="U")%>%dim()
EWM.full.data=EWM.alllakes.data%>%filter(EWMSTATUS!="U") ### All possible SURVEYED lakes that can be used in SDM predictions
EWM.full.data

### 1. Get the lake index data
####### Lake DOW, LON, LAT, COUNTY, SURVEYED, EWMSTATUS = 6 columns
Full.index=EWM.full.data[, c(1,3:7)]
Full.index

### 2. Get response data and predictors: two sets of predictors for water GDD - present and future for 5 GCMs
####### 2 unchanging covariates, 10 GDD estimates (5GCMs*2 time-periods), and 5 domains (1*5GCMs) = 17 columns

EWM.alllakes.Secchi=EWM.alllakes.data%>%filter(EWMSTATUS!="U")%>%select('DOWLKNUM', 'LON', 'LAT','COUNTY' ,'SURVEYED','EWMSTATUS', 'avg_secchi')%>%
  na.omit(avg_secchi)

LakeConn.data=read_csv("raw_data/LakeConn.data.csv")
RoadDensity.data=LakeConn.data%>%select(DOWLKNUM,roaddensity_density_mperha)%>%na.omit()
EWM.alllakes.Secchi.Roads=left_join(EWM.alllakes.Secchi, RoadDensity.data, by="DOWLKNUM")%>%na.omit()

EWM.CurrTrain.Data=read_csv("processed_data/EWM.prsabs95to15_AllGCMs_v2.csv")
EWM.CurrTrain.Data
colnames(EWM.CurrTrain.Data)[11:15]=c("ACCESS.Curr", "MIROC5.Curr", "IPSL.Curr",   "GFDL.Curr" ,  "MRI.Curr")
EWM.FutrTest.Data=read_csv("processed_data/EWM.prsabs40to60_AllGCMs_v2.csv")
colnames(EWM.FutrTest.Data)[11:15]=c("ACCESS.Futr", "MIROC5.Futr", "IPSL.Futr",   "GFDL.Futr" ,  "MRI.Futr")
EWM.CurrFutr.Preds=left_join(EWM.CurrTrain.Data, EWM.FutrTest.Data[,c(1,11:15)], by="DOWLKNUM")

EWM.CurrFutr.Preds.Doms=EWM.CurrFutr.Preds%>%mutate(ACCESS.domain=case_when(ACCESS.Futr < max(ACCESS.Curr) ~ 'Analog', ACCESS.Futr > max(ACCESS.Curr)~ 'NonAnalog'),
                            MIROC5.domain=case_when(MIROC5.Futr < max(MIROC5.Curr) ~ 'Analog', MIROC5.Futr > max(MIROC5.Curr) ~ 'NonAnalog'),
                            IPSL.domain=case_when(IPSL.Futr < max(IPSL.Curr) ~ 'Analog', IPSL.Futr > max(IPSL.Curr) ~ 'NonAnalog'),
                            GFDL.domain=case_when(GFDL.Futr < max(GFDL.Curr) ~ 'Analog', GFDL.Futr > max(GFDL.Curr) ~ 'NonAnalog'),
                            MRI.domain=case_when(MRI.Futr < max(MRI.Curr) ~ 'Analog', MRI.Futr > max(MRI.Curr) ~ 'NonAnalog')
)
        
EWM.CurrFutr.Preds.Doms                    
write_csv(EWM.CurrFutr.Preds.Doms, "processed_data/EWM.CurrFutr.Preds.Doms.csv")

### 3. Predictions of EWM occurrence: two sets for each 5 GCM temp estimates and 3 SDM modeling algorithms
####### 30 columns: SDM_GCM_PERIOD
########### Start with GAM k=3
fut.preds.df=read_csv("Results/GAM.k3_Fut.Predictions.csv") ## from 2.1 GAM_Predictions output
curr.preds.df=read_csv("Results/GAM.k3_Curr.Predictions.csv")
currANDfut_preds=bind_rows(curr.preds.df,fut.preds.df)
head(currANDfut_preds) 
dim(currANDfut_preds)
## Split into 2 seperate datasets of 578 rows each
curr.preds.GAM_k3=currANDfut_preds[1:578,]
futr.preds.GAM_k3=currANDfut_preds[579:1156,]
### Give each column a unique name
colnames(futr.preds.GAM_k3)[1:6]=c("ACCESS.GAMk3.FutrPred" ,  "GFDL.GAMk3.FutrPred"   ,  "IPSL.GAMk3.FutrPred"  ,   "MIROC5.GAMk3.FutrPred"  , "MRI.GAMk3.FutrPred"   ,   "DOWLKNUM")
colnames(curr.preds.GAM_k3)[1:6]=c("ACCESS.GAMk3.CurrPred" ,  "GFDL.GAMk3.CurrPred"   ,  "IPSL.GAMk3.CurrPred"  ,   "MIROC5.GAMk3.CurrPred"  , "MRI.GAMk3.CurrPred"   ,   "DOWLKNUM")

########### Next with GAM k=10
fut.preds.df=read_csv("Results/GAM.k10_Fut.Predictions.csv") #from 2.1 GAM_Predictions output
curr.preds.df=read_csv("Results/GAM.k10_Curr.Predictions.csv")
currANDfut_preds_k10=bind_rows(curr.preds.df,fut.preds.df)
currANDfut_preds_k10
dim(currANDfut_preds_k10)

curr.preds.GAM_k10=currANDfut_preds_k10[1:578,]
futr.preds.GAM_k10=currANDfut_preds_k10[579:1156,]

colnames(futr.preds.GAM_k10)[1:6]=c("ACCESS.GAMk10.FutrPred" ,  "GFDL.GAMk10.FutrPred"   ,  "IPSL.GAMk10.FutrPred"  ,   "MIROC5.GAMk10.FutrPred"  , "MRI.GAMk10.FutrPred"   ,   "DOWLKNUM")
colnames(curr.preds.GAM_k10)[1:6]=c("ACCESS.GAMk10.CurrPred" ,  "GFDL.GAMk10.CurrPred"   ,  "IPSL.GAMk10.CurrPred"  ,   "MIROC5.GAMk10.CurrPred"  , "MRI.GAMk10.CurrPred"   ,   "DOWLKNUM")

########### Finally with Random Forest model predictions
RF_curr.preds=read_csv("Results/Curr.Predictions.csv")
RF_futr.preds=read_csv("Results/Futr.Predictions.csv")
RF_curr.preds
RF_futr.preds
colnames(RF_curr.preds)[1:5]=c("ACCESS.RF.CurrPred" ,  "GFDL.RF.CurrPred"   ,  "IPSL.RF.CurrPred"  ,   "MIROC5.RF.CurrPred"  , "MRI.RF.CurrPred")
colnames(RF_futr.preds)[1:5]=c("ACCESS.RF.FutrPred" ,  "GFDL.RF.FutrPred"   ,  "IPSL.RF.FutrPred"  ,   "MIROC5.RF.FutrPred"  , "MRI.RF.FutrPred")

RF_futr.preds$DOWLKNUM=curr.preds.GAM_k3$DOWLKNUM
RF_futr.preds
RF_curr.preds$DOWLKNUM=curr.preds.GAM_k3$DOWLKNUM
RF_curr.preds


Final_MergedData=left_join(EWM.CurrFutr.Preds.Doms,curr.preds.GAM_k3[,1:6], by="DOWLKNUM")%>%
                  left_join(.,futr.preds.GAM_k3[,1:6], by="DOWLKNUM")%>%
                        left_join(.,curr.preds.GAM_k10[,1:6], by="DOWLKNUM")%>%
                        left_join(.,futr.preds.GAM_k10[,1:6], by="DOWLKNUM")%>%
                             left_join(., RF_curr.preds[,-6],by="DOWLKNUM")%>%
                              left_join(., RF_futr.preds[,-6],by="DOWLKNUM")

Final_MergedData%>%View()
write_csv(Final_MergedData, "Results/AllEWM_SDM_results.csv")

Final_MergedData=read_csv("Results/AllEWM_SDM_results.csv")

Final_MergedData_V2=left_join(Final_MeregedData,Curr.Brm.Preds_EWMGeoIndex[,c(5:7)], by="DOWLKNUM")%>%
                        left_join(.,Fut.Brm.Preds_EWMGeoIndex[,c(5:7)], by="DOWLKNUM")
Final_MergedData_V2
write_csv(Final_MergedData_V2, "Results/Final_MergedData.csv")

################# TOTAL COLUMSN ESTIMATED~5+18+30=53