library(tidyverse) 
library(sf)
library(patchwork)
library(gratia)
library(mgcv)
library(cowplot)


################################################################################################################
##### R scripts of all the figures that went into the manuscript main text and appendix                  #######
################################################################################################################

### The final data with all SDM results, lake-level covariates, and Minnesota map shapefile
EWM.clmchng_data=read_csv("Results/Final_MergedData.csv") 
Minn.sf=read_sf(dsn="/Users/thom7552/UMNpostdoc/ProjectEWM/MinnEWM/MinnGISlayers", layer="Minn.map")

############################  
### Figure 1: Boxplots comparing predicted EWM invasion risk under current and future water temperature/GDD
EWM.clmchng_data ### from Line 6 above
EWM.climchng.preds_wide=EWM.clmchng_data%>%select(c(1,3,4),ends_with("Pred"))%>%
    pivot_longer(!c(DOWLKNUM,LON,LAT), names_to = "models", values_to = "preds")
EWM.clim.chng.preds=EWM.climchng.preds_wide%>%separate(models, c("GCM","SDM","PERIOD"))
EWM.clim.chng.preds

EWM.clim.chng.preds2=EWM.clim.chng.preds%>%mutate(SDM=recode(SDM, GAMk10="GAM (k=10)",GAMk3="GAM (k=03)", RF="Random Forest"),
                                                  PERIOD=recode(PERIOD,CurrPred="Current", FutrPred="Future"))
EWM.clim.chng.preds2


PeriodWise=EWM.clim.chng.preds%>%ggplot(aes(x = SDM, y = preds, fill = GCM)) +
            geom_boxplot(outlier.shape = NA) +
            facet_grid(~PERIOD) +
            theme(legend.position = "bottom")

PeriodWise+scale_fill_viridis_d(option="inferno")

SDMWise=EWM.clim.chng.preds2%>%ggplot(aes(x = PERIOD, y = preds, fill = GCM)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(~SDM) +
        theme(legend.position = c(0.93,0.22), legend.box.background = element_rect(color = "black"))+
        xlab(" ")+ylab("EWM habitat suitability")+
        theme(text = element_text(size=12))

SDMWise+scale_fill_viridis_d(option="inferno")
ggsave("Mnspt.Fig1.png", path="Figures/", device="png",width=9, height=4.5, dpi=600)


PeriodWise
####################################################################################
### Figure 2: Arrow plots showing change in risk over time
EWM.clmchng_data ### from Line 6 above
### Estimate change in  mean risk predictions
EWM.climchng.MeanPredsBySDM.Year=EWM.clmchng_data%>%select(c(1,3,4),ends_with("Pred"))%>%rowwise%>%mutate(
  MeanGAM_k3_Curr=mean(c_across(ends_with("GAMk3.CurrPred"))),
  MeanGAM_k3_Futr=mean(c_across(ends_with("GAMk3.FutrPred"))),
  
  MeanGAM_k10_Curr=mean(c_across(ends_with("GAMk10.CurrPred"))),
  MeanGAM_k10_Futr=mean(c_across(ends_with("GAMk10.FutrPred"))),
  
  MeanRF_Curr=mean(c_across(ends_with("RF.CurrPred"))),
  MeanRF_Futr=mean(c_across(ends_with("RF.FutrPred")))
)

EWM.clim.chng.PredsChange=EWM.climchng.MeanPredsBySDM.Year%>%select(c(1:3),starts_with("Mean"))%>%mutate(
                            GAM_k3_Change=MeanGAM_k3_Futr-MeanGAM_k3_Curr,
                            GAM_k10_Change=MeanGAM_k10_Futr-MeanGAM_k10_Curr,
                            RF_Change=MeanRF_Futr-MeanRF_Curr)
EWM.clim.chng.PredsChange%>%View()

### Change in mean temp conditions
EWM.climchng.MeanTempByYear=EWM.clmchng_data%>%select(c(1,3,4),ends_with(c("Curr","Futr")))%>%rowwise%>%mutate(
  MeanTemp_Curr=mean(c_across(ends_with("Curr"))),
  MeanTemp_Futr=mean(c_across(ends_with("Futr")))
)

EWM.clim.chng.TempChange=EWM.climchng.MeanTempByYear%>%select(c(1:3),starts_with("Mean"))%>%rowwise%>%mutate(
                        TempDiff=MeanTemp_Futr-MeanTemp_Curr
)
        
EWM.clim.chng.Preds.Temp.Change=left_join(EWM.clim.chng.PredsChange,EWM.clim.chng.TempChange[,c(1,4:6)], by="DOWLKNUM")
EWM.clim.chng.Preds.Temp.Change

### Plot for GAM k=3 
GAM_k3=ggplot()+
  geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanGAM_k3_Curr, xend=MeanTemp_Curr+TempDiff, 
  yend=MeanGAM_k3_Curr+GAM_k3_Change,color=GAM_k3_Change), arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanGAM_k3_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))

### GAM k=3 response curve inset plot
gam_k3=gam(EWMSTATUS ~ s(avg_secchi, k=3)+s(roaddensity_density_mperha,k=3)+s(ACCESS.avg.ann.gdd,k=3), data=EWM.train.data_ACCESS.WtrTemp)

gam.k3=draw(gam_k3, residuals = FALSE, select=3, rug=FALSE, ci_alpha = 0.00, title="")
GAM.k3.GDD_respcurv=gam.k3+xlab("GDD")+ylab("EWM")+theme_classic(24)+ggtitle("")+theme(
axis.text.x = element_blank(),
axis.text.y = element_blank(),
axis.ticks = element_blank())

gamk3_change_inset=ggdraw(GAM_k3 +ggtitle("a) GAM (k=03)"))+
  draw_plot(GAM.k3.GDD_respcurv, .77, .09, .2, .2, scale=1)
ggsave("GAM.k3_Change_wInset.png", path="Figures/", device="png",width = 12, height = 8, dpi=1200)


### Plot for GAM k=10 
GAM_k10=ggplot()+
  geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanGAM_k10_Curr, xend=MeanTemp_Curr+TempDiff, 
  yend=MeanGAM_k10_Curr+GAM_k10_Change,color=GAM_k10_Change), arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanGAM_k10_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))

### GAM k=10 response curve inset plot
gam=gam(EWMSTATUS ~ s(ACCESS.avg.ann.gdd), data=EWM.train.data_ACCESS.WtrTemp)

gam.k10=draw(gam, residuals = FALSE, rug=FALSE, ci_alpha = 0.00, title="")
GAM.k10.GDD_respcurv=gam.k10+xlab("GDD")+ylab("EWM")+theme_classic(24)+ggtitle("")+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank())

gam.k10_change_inset=ggdraw(GAM_k10 +ggtitle("a) GAM (k=10)"))+
  draw_plot(GAM.k10.GDD_respcurv, .77, .09, .2, .2, scale=1)
ggsave("GAM.k10_Change_wInset.png", path="Figures/", device="png",width = 12, height = 8, dpi=1200)

### Plot for Random forest plot
RF=ggplot()+
geom_segment(data=EWM.clim.chng.Preds.Temp.Change, aes(x=MeanTemp_Curr, y=MeanRF_Curr, xend=MeanTemp_Curr+TempDiff,
yend=MeanRF_Curr+RF_Change,color=RF_Change),arrow=arrow(), size=0.5) +
  geom_point(data=EWM.clim.chng.Preds.Temp.Change, mapping=aes(x=MeanTemp_Curr, y=MeanRF_Curr), size=1, shape=21, fill="white")+
  scale_color_viridis_c()+xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change \nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=20))+
  theme(legend.position = c(0.06,0.8))

### Random forest response curve inset plot
rf=randomForest(EWM.train.data_ACCESS.WtrTemp[,c(2:4)], EWM.train.data_ACCESS.WtrTemp$EWMSTATUS,importance=TRUE, ntree=5000, type="regression")
par.RF= pdp::partial(rf, pred.var = c("ACCESS.avg.ann.gdd"))

RF.GDD_respcurv=autoplot(par.RF)+xlab("GDD")+ylab("EWM")+theme_classic(24)+theme(
axis.text.x = element_blank(),
axis.text.y = element_blank(),
axis.ticks = element_blank())

rf_change_inset=ggdraw(RF +ggtitle("c) Random forest"))+
draw_plot(RF.GDD_respcurv, .77, .09, .2, .2, scale=1)

ggsave("RF_Change_wInset.png", path="Figures/", device="png",width = 12, height = 8, dpi=1200)

### Making one singe multi-panel plot with a common legend
### Pivot the table from wide to long
EWM.clim.chng.Preds.Temp.Change2=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(4,6,8), names_to = 'Curr_Models' , values_to = 'Curr.Preds')
t1=EWM.clim.chng.Preds.Temp.Change2[,c(1:3,10:14)]

EWM.clim.chng.Preds.Temp.Change3=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(5,7,9), names_to = 'Futr_Models' , values_to = 'Futr.Preds')
t2=EWM.clim.chng.Preds.Temp.Change3[,c(1,13,14)]
EWM.clim.chng.Preds.Temp.Change4=EWM.clim.chng.Preds.Temp.Change%>%
  pivot_longer(cols= c(10:12), names_to = 'Chng_Models' , values_to = 'Chng.Preds')
t3=EWM.clim.chng.Preds.Temp.Change4[,c(1,13,14)]

t1$Futr.preds=t2$Futr.Preds
t1$Chng.preds=t3$Chng.Preds
final_table=t1%>%mutate(Model= recode(Curr_Models, MeanGAM_k3_Curr = "GAM (k=03)", MeanGAM_k10_Curr = "GAM (k=10)",
                                      MeanRF_Curr = "Random Forest" ))
final_table

library(ragg)
agg_png("Mnspt.Fig2_new.png", width = 18, height = 9, units = "in", res = 600, scaling = 2)

### the final multipanel plot using the final table in long format
ggplot()+
geom_segment(data=final_table, aes(x=MeanTemp_Curr, y=Curr.Preds, xend=MeanTemp_Curr+TempDiff,
yend=Curr.Preds+Chng.preds,color=Chng.preds), arrow=arrow(), size=0.5) +
geom_point(data=final_table, mapping=aes(x=MeanTemp_Curr, y=Curr.Preds), size=1, shape=21, fill="white")+
scale_color_viridis_c(limits = c(-0.4, 0.9),
breaks = c(-0.4,0.00, 0.50, 0.9),
labels = c(-0.4,0.00, 0.50, 0.9))+
xlab("Annual growing degree days (GDD)")+ylab("EWM invasion risk")+
  labs(colour="Change\nin risk")+geom_vline(xintercept = 2200, lty=2)+theme(text=element_text(size=12))+
    facet_wrap(~Model)+theme_bw()+
    theme(legend.position="right", legend.box.background = element_rect(color = "black"))

dev.off()

### Figure 3: Final EWM invasion risk predictions and uncertainty maps, analog/non-analog domains
EWM.clim.chng.PredsChange ### from line 100 above
EWM.clim.chng.PredsChange=EWM.climchng.MeanPredsBySDM.Year%>%select(c(1:3),starts_with("Mean"))%>%mutate(
GAM_k3_Change=MeanGAM_k3_Futr-MeanGAM_k3_Curr,
GAM_k10_Change=MeanGAM_k10_Futr-MeanGAM_k10_Curr,
RF_Change=MeanRF_Futr-MeanRF_Curr)

### Estimate change in EWM suitability status
EWM.clim.chng.ChangeStatus=EWM.clim.chng.PredsChange%>%select(1:3, ends_with("Change"))%>%
  mutate(GAM.k3_Change_Status=case_when(GAM_k3_Change < -0.1 ~ 'loss',
          GAM_k3_Change < 0.1 ~ 'no change',
          GAM_k3_Change < 0.85 ~ 'gain'), 

        GAM.k10_Change_Status=case_when(GAM_k10_Change < -0.1 ~ 'loss',
          GAM_k10_Change < 0.1 ~ 'no change',
          GAM_k10_Change < 0.85 ~ 'gain'),

        RF_Change_Status=case_when(RF_Change < -0.1 ~ 'loss',
          RF_Change < 0.1 ~ 'no change',
          RF_Change < 0.9 ~ 'gain'),
)

### Estimate consensus among models in change in EWM suitability status
StatusCount=EWM.clim.chng.ChangeStatus%>%select(ends_with("Status"))%>%rowwise()%>%
  do(data.frame(., StatusCount = n_distinct(unlist(.))))%>%pull(StatusCount)
EWM.clim.chng.ChangeStatus$StatusCount=StatusCount
EWM.clim.chng.ChangeStatus%>%mutate(Status=recode(StatusCount, "1"= "Increase", "2"= "Two", "3"="Uncertain"))

IncreasersUncertain_Lakes=EWM.clim.chng.ChangeStatus%>%filter(StatusCount ==1 | StatusCount>1)
IncreasersUncertain_Lakes
IncreasersUncertain_Lakes_sf=st_as_sf(IncreasersUncertain_Lakes, coords=c("LON","LAT"), crs=32615)
IncreasersUncertain_Lakes_sf_WGS=st_transform(IncreasersUncertain_Lakes_sf, crs=4326)
IncreasersUncertain_Lakes_sf_WGS_recoded=IncreasersUncertain_Lakes_sf_WGS%>%
  mutate(Status=recode(StatusCount, "1"= "Increase","2"="Uncertain", "3"="Uncertain"))%>%View()

### Plot the change in invasion risk status (trajectory)
Fig.3a=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=IncreasersUncertain_Lakes_sf_WGS_recoded, aes(color=as.factor(Status),shape=as.factor(Status)), cex=4)+
  theme_light()+
    scale_color_viridis_d(option="turbo",alpha = 0.5)+theme(legend.title = element_blank())+
    theme(legend.position = c(0.85,0.4))+theme(text=element_text(size=12))+
      ggtitle("a) Future invasion risk trajectory")
agg_png("Figures/Mnspt.Fig3a_new.png", width = 9, height = 4.5, units = "in", res = 600, scaling = 0.75)
Fig.3a
dev.off()

### Subset bayesian predictions from the best fitting GAM k=10 model 
EWM.futrpreds.bayes.GAM_k10=EWM.clmchng_data%>%select(1,3,4,58,59)%>%
  rename(Estimate = bestmodel.bayesian.meanFutrPreds, Variance = bestmodel.bayesian.varFutrPreds)
EWM.futrpreds.bayes.GAM_k10

EWM.futrpreds.bayes.GAM_k10_sf=st_as_sf(EWM.futrpreds.bayes.GAM_k10, coords = c("LON", "LAT"),crs=32615)
EWM.futrpreds.bayes.GAM_k10_sf_WGS=st_transform(EWM.futrpreds.bayes.GAM_k10_sf, crs=4326)
EWM.futrpreds.bayes.GAM_k10_sf_WGS

### Plot the lake-level bayesian predictions and variance associated with it
Fig.3b=ggplot(data=Minn.sf)+geom_sf()+
          geom_sf(data=EWM.futrpreds.bayes.GAM_k10_sf_WGS, aes(color=Estimate, size=Variance))+theme_light()+
          scale_color_viridis_c(option="turbo",alpha = 0.75)+
          guides(size=guide_legend(override.aes=list(shape=1,size=c(2,4,8))))+
          theme(legend.position = c(0.85,0.375))+theme(text=element_text(size=12))+ggtitle("b) Future invasion risk predictions")

agg_png("Figures/Mnspt.Fig3b_new.png", width = 9, height = 4.5, units = "in", res = 600, scaling = 0.75)
Fig.3b
dev.off()

### Subset data to identify the temperature analog and non-analog domains
EWM.domains=EWM.climchng.MeanTempByYear%>%select(1:3,14,15)%>%
                mutate(Domain = case_when(MeanTemp_Futr < 2225 ~ 'Analog', MeanTemp_Futr > 2225 ~ 'Non-analog'))
EWM.domains.sf=st_as_sf(EWM.domains, coords=c("LON", "LAT"), crs=32615)
EWM.domains.sf.WGS=st_transform(EWM.domains.sf, crs=4326)

### Map the lakes in analog and non-analog domains
Fig.3c=ggplot(data=Minn.sf)+geom_sf()+
          geom_sf(data=EWM.domains.sf.WGS, aes(color=Domain,shape=Domain), cex=4)+theme_light()+
          scale_color_viridis_d(option="turbo",alpha = 0.5)+theme(legend.title = element_blank())+
          theme(legend.position = c(0.8,0.45))+theme(text=element_text(size=12))+ggtitle("c) Future temperature domains")
agg_png("Figures/Mnspt.Fig3c_new.png", width = 9, height = 4.5, units = "in", res = 600, scaling = 0.75)
Fig.3c
dev.off()


Fig.3d=Fut.Brm.Preds%>%ggplot(.,)+geom_point(aes(avgGDD,mean.fut.preds, size=var.fut.preds), alpha=0.25)+
geom_smooth(aes(x=avgGDD,y=mean.fut.preds), method="loess", span=0.3,se=FALSE, col="black", linewidth=0.5)+
geom_smooth(aes(x= avgGDD, y=min.fut.preds),method="loess", span=0.3, se=FALSE, lty=2, linewidth=0.5)+
geom_smooth(aes(x= avgGDD, y=max.fut.preds), method="loess", span=0.3,se=FALSE, lty=2, linewidth=0.5)+
geom_vline(xintercept = 2200, lty=3, linewidth=1)+xlab("Growing degree days (GDD)")+
ylab("Future EWM suitability")+labs(size="Variance")+theme(legend.position = c(0.85,0.2))+
theme(text=element_text(size=12))+ggtitle("d) Predicted response curve")

agg_png("Figures/Mnspt.Fig3d_new.png", width = 9, height = 4, units = "in", res = 600, scaling = 0.75)
Fig.3d
dev.off()


#####################################################################################################################
#####################################################################################################################
#### Appendix figures
#### Figure S1: Maps showing EWM distribution and lake temperatures in GDD
Minn.sf
ggplot(data=Minn.sf)+geom_sf()

EWM.clmchng_data ### from Line 6 above
EWM.climchng.Meantemps=EWM.clmchng_data%>%select(c(1,3,4,7),ends_with(c("Curr","Futr")))%>%rowwise%>%mutate(
  MeanCurr=mean(c_across(ends_with("Curr"))),
  MeanFutr=mean(c_across(ends_with("Futr"))),
  PerChange=((MeanFutr-MeanCurr)/MeanCurr)*100
)
EWM.climchng.Meantemps%>%View()

EWM.climchng.Meantemps_sf=st_as_sf(EWM.climchng.Meantemps, coords = c("LON", "LAT"))
EWM.climchng.Meantemps_sf=st_set_crs(EWM.climchng.Meantemps_sf, 32615)
EWM.climchng.Meantemps_sf%>%View()
EWM.climchng.Meantemps_WGS=st_transform(EWM.climchng.Meantemps_sf, crs=4326)
EWM.climchng.Meantemps_WGS2=EWM.climchng.Meantemps_WGS%>%mutate(EWM=recode(EWMSTATUS, "0" = "Abs", "1" = "Prs"))
EWM.climchng.Meantemps_WGS2


CurrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanCurr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Current GDD")+
  theme(text=element_text(size=16))

FutrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanFutr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Future GDD")+
  theme(text=element_text(size=16))

ChangeGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=PerChange))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("% Increase in GDD")+
  theme(text=element_text(size=16))

library("patchwork")
CurrGDDmap|FutrGDDmap|ChangeGDDmap
ggsave("Mnspt.S1_Fig1.png", path="Figures/", device="png",width=12, height=4, dpi=900)

### Figure S2
EWMmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=as.factor(EWM)))+theme_light()+ ggtitle("EWM distribution")+
  scale_color_viridis_d(option = "turbo", alpha = 0.5)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+
  theme(text=element_text(size=16))

CurrGDDmap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.climchng.Meantemps_WGS2, aes(color=MeanCurr))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Current GDD")+
  theme(text=element_text(size=16))

EWM.secchi.roads=EWM.clmchng_data%>%select(c(1,3,4,7:10))
EWM.secchi.roads_sf=st_as_sf(EWM.secchi.roads, coords = c("LON", "LAT"))
EWM.secchi.roads_sf=st_set_crs(EWM.secchi.roads_sf, 32615)
EWM.secchi.roads_WGS=st_transform(EWM.secchi.roads_sf, crs=4326)
EWM.secchi.roads_WGS2=EWM.secchi.roads_WGS%>%mutate(EWM=recode(EWMSTATUS, "0" = "Abs", "1" = "Prs"))

SecchiMap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.secchi.roads_WGS2, aes(color=max_secchi))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Secchi depth")+
  theme(text=element_text(size=16))

RoadsMap=ggplot(data=Minn.sf)+geom_sf()+
  geom_sf(data=EWM.secchi.roads_WGS2, aes(color=roaddensity_density_mperha))+theme_light()+
  scale_color_viridis_c(option = "turbo", alpha = 0.75)+
  theme(legend.title = element_blank())+theme(legend.position = c(0.85,0.4))+ggtitle("Road density")+
  theme(text=element_text(size=16))

(EWMmap|CurrGDDmap)/(RoadsMap|SecchiMap)
ggsave("Mnspt.S1_Fig2.png", path="Figures/", device="png",width=10, height=12, dpi=900)
#####################################################################################################################
#####################################################################################################################
