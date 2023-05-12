library(tidyverse)
library(tidybayes)
library(bayesplot)
library(brms)
library(matrixStats)

####################################################################################################################
######### R script: Developing a Bayesian version of the best-fitting SDM (GAM, k=10) to predict 
######### point estimates of EWM risk with a measure of uncertainty.                                                            #########
####################################################################################################################

Data_Combined = list.files(path="processed_data/TrainData/",pattern=".csv", full.names = TRUE) %>%
  lapply(read_csv)%>%
  bind_rows(id=NULL)%>%
  mutate(Id=rep(1:578,5))

Data_Combined%>%View()


TempAvgs=Data_Combined%>%group_by(Id)%>%summarise(
avgGDD=mean(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE),
minGDD=min(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE),
maxGDD=max(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE)
)%>%ungroup()

TempAvgs

Access.TrainData=read_csv("processed_data/TrainData/EWM.train.data_ACCESS.WtrTemp.csv")
colnames(Access.TrainData)
Avg.CurrTemp.data=bind_cols(Access.TrainData[,-4],TempAvgs[,2:4])
Avg.CurrTemp.data

### We need only 3 predictors
Curr.Brm.TrainData=Avg.CurrTemp.data[,c(1:4)]
Curr.Brm.TrainData

colnames(Curr.Brm.TrainData)[2:4]=c("avgSecchi", "avgRoads", "avgGDD")
Curr.Brm.TrainData

ewm.brm.trial <- brm(bf(EWMSTATUS ~ s(avgGDD, k=10)+s(avgSecchi, k=10)+s(avgRoads, k=10)),
data = Curr.Brm.TrainData, family = bernoulli(), cores = 4, seed = 17,
iter = 4000, warmup = 1000, thin = 10, refresh = 0,
control = list(adapt_delta = 0.99))

summary(ewm.brm.trial)
pp_check(ewm.brm.trial)
plot(ewm.brm.trial)
post_brms.trial=posterior_samples(ewm.brm.trial)

### Lets make posterior predictions from sample draws
curr.ewm.epreds=posterior_epred(ewm.brm.trial, draw_ids = c(1001:1200))
str(curr.ewm.epreds)
head(curr.ewm.epreds)

library(matrixStats)

avg.curr.preds=colMeans(curr.ewm.epreds)
sd.curr.preds=colSds(curr.ewm.epreds)
var.curr.preds=colVars(curr.ewm.epreds)
plot(Curr.Brm.TrainData$avgGDD, avg.curr.preds)
plot(Curr.Brm.TrainData$avgGDD, var.curr.preds)

Curr.Brm.Preds=Curr.Brm.TrainData%>%mutate(mean.curr.preds=avg.curr.preds,var.curr.preds=var.curr.preds, 
                                           DOWLKNUM=EWM.GCMs.data$DOWLKNUM)
Curr.Brm.Preds
EWM_GeoCoord.Index=EWM.GCMs.data[,1:5]
EWM_GeoCoord.Index
Curr.Brm.Preds_EWMGeoIndex=left_join(Curr.Brm.Preds, EWM_GeoCoord.Index, by="DOWLKNUM")
Curr.Brm.Preds_EWMGeoIndex


FutData_Combined = list.files(path="processed_data/TestData/ForecastData",pattern=".csv", full.names = TRUE) %>%
lapply(read_csv)%>%
bind_rows(id=NULL)%>%
mutate(Id=rep(1:578,5))

FutTempAvgs=FutData_Combined%>%group_by(Id)%>%summarise(
  avgGDD=mean(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE),
  minGDD=min(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE),
  maxGDD=max(c_across(ACCESS.avg.ann.gdd:MRI.avg.ann.gdd), na.rm=TRUE)
)%>%ungroup()
FutTempAvgs

Fut.Brm.TestData=bind_cols(Access.TrainData[,-4],FutTempAvgs[,2])
colnames(Fut.Brm.TestData)[2:4]=c("avgSecchi", "avgRoads", "avgGDD")
Fut.Brm.TestData

fut.ewm.epreds=posterior_epred(ewm.brm.trial, newdata =Fut.Brm.TestData[,-1], draw_ids =  c(1001:1200))
str(fut.ewm.epreds)
head(fut.ewm.epreds)

med.fut.preds=colMedians(fut.ewm.epreds)
sd.fut.preds=colSds(fut.ewm.epreds)
var.fut.preds=colVars(fut.ewm.epreds)
min.fut.preds=med.fut.preds-(2.56*(sd.fut.preds/sqrt(200)))
max.fut.preds=med.fut.preds+(2.56*(sd.fut.preds/sqrt(200)))
plot(Fut.Brm.TestData$avgGDD, med.fut.preds, abline(v=2200))
plot(Fut.Brm.TestData$avgGDD, var.fut.preds, abline(v=2200))

Fut.Brm.Preds=Fut.Brm.TestData%>%mutate(mean.fut.preds=med.fut.preds,sd.fut.preds=sd.fut.preds, 
                                        min.fut.preds=min.fut.preds,max.fut.preds=max.fut.preds,
                                        DOWLKNUM=EWM.GCMs.data$DOWLKNUM)
Fut.Brm.Preds
EWM_GeoCoord.Index=EWM.GCMs.data[,1:5]
Fut.Brm.Preds_EWMGeoIndex=left_join(Fut.Brm.Preds, EWM_GeoCoord.Index, by="DOWLKNUM")
Fut.Brm.Preds_EWMGeoIndex


#### Finally map the data
Curr.Brm.Preds_EWMGeoIndex.sf=st_as_sf(Curr.Brm.Preds_EWMGeoIndex,coords=c("LON", "LAT"),crs=32615)
Curr.Brm.Preds_EWMGeoIndex.sf

Fut.Brm.Preds_EWMGeoIndex.sf=st_as_sf(Fut.Brm.Preds_EWMGeoIndex,coords=c("LON", "LAT"),crs=32615)
Fut.Brm.Preds_EWMGeoIndex.sf

Minn.sf=read_sf(dsn="/Users/thom7552/UMNpostdoc/ProjectEWM/MinnEWM/MinnGISlayers", layer="Minn.map")
Minn.sf=st_transform(Minn.sf, crs=32615)
Minn.sf

ggplot(Minn.sf)+geom_sf()+geom_sf(data=Fut.Brm.Preds_EWMGeoIndex.sf, aes(col=mean.fut.preds))+scale_color_viridis_c(name="Estimate")+theme_minimal()+
theme(text=element_text(size=16))+theme(legend.position = c(0.8, 0.4))

Fut.Brm.Preds%>%ggplot(.,)+geom_point(aes(avgGDD,mean.fut.preds, size=var.fut.preds), alpha=0.25)+
  geom_smooth(aes(x=avgGDD,y=mean.fut.preds), method="loess", span=0.3,se=FALSE, col="black", linewidth=0.5)+
  geom_smooth(aes(x= avgGDD, y=min.fut.preds),method="loess", span=0.3, se=FALSE, lty=2, linewidth=0.5)+
  geom_smooth(aes(x= avgGDD, y=max.fut.preds), method="loess", span=0.3,se=FALSE, lty=2, linewidth=0.5)+
  geom_vline(xintercept = 2200, lty=3, linewidth=1)+xlab("Growing degree days (GDD)")+
    ylab("Predicted future invasion risk")+labs(size="Variance")+theme(text=element_text(size=16))

ggsave("Mnspt.S1_Fig4.png", path="Figures/", device="png",width=12, height=9, dpi=900)

