rm(list=ls())
library(ggplot2)
library(scales)
library(extrafont)
library(grid)
library(doBy)
library(data.table)
loadfonts(device = "win")
setwd("C:/Users/Conor/Documents/Research/CovCAInertia/")

#### Data ####
rundate = "2019-07-30"
spec = "Spec3_"
df_full = as.data.table(read.csv(paste("Output/Estimation_Results/active_",spec,rundate,".csv",sep="")))
df_full = df_full[returning==1,]
df_full[,always_active:=sum(active_obs)==sum(returning),by="Person"]
df_full[,num_years:=sum(returning),by="Person"]
df_full[,active_pred_bucket:=floor(active_pred/(.05))*.05]

#### By Predicted Activity ####
df = df_full[,list(stay_obs=mean(1-stay_obs),active_obs=mean(active_obs),
                   stay_sd = sd(1-stay_obs),active_sd=sd(active_obs),
                   count=sum(returning)),by="active_pred_bucket"]
setkey(df,active_pred_bucket)

df[,stay_sd:=stay_sd/sqrt(count)]
df[,active_sd:=active_sd/sqrt(count)]

df[,stay_max:=pmin(stay_obs+1.96*stay_sd,1)]
df[,stay_min:=pmax(stay_obs-1.96*stay_sd,0)]

df[,active_max:=pmin(active_obs+1.96*active_sd,1)]
df[,active_min:=pmax(active_obs-1.96*active_sd,0)]


png("Writing/ActivePrediction.png",width=2000,height=1500,res=275)
ggplot(df) + aes(x=active_pred_bucket) + 
  geom_point(size=3,aes(y=stay_obs,shape="Switch Plans")) + 
  geom_line(size=1,aes(y=stay_obs)) + 
  geom_line(size=1,aes(y=stay_max),linetype=2) + 
  geom_line(size=1,aes(y=stay_min),linetype=2) + 
  # geom_errorbar(width=0.05,aes(ymax=stay_max,ymin=stay_min)) + 
  geom_point(size=3,aes(y=active_obs,shape="Active on Website")) + 
  geom_line(size=1,aes(y=active_obs)) + 
  geom_line(size=1,aes(y=active_max),linetype=2) + 
  geom_line(size=1,aes(y=active_min),linetype=2) + 
  # geom_errorbar(width=0.05,aes(ymax=active_max,ymin=active_min)) + 
  xlab("Predicted Attention Probability") + 
  ylab("Observed Choices") + 
  scale_y_continuous(labels=percent)+ 
  scale_x_continuous(labels=percent)+ 
  theme(#panel.background = element_rect(color=grey(.2),fill=grey(.9)),
    strip.background = element_blank(),
    text = element_text(family="Times New Roman"),
    #panel.grid.major = element_line(color=grey(.8)),
    legend.background = element_rect(color=grey(.5)),
    legend.title=element_blank(),
    legend.text = element_text(size=16),
    legend.key.width = unit(.075,units="npc"),
    legend.key = element_rect(color="transparent",fill="transparent"),
    legend.position = "bottom",
    axis.title=element_text(size=16),
    axis.text = element_text(size=16))
dev.off()







#### By Active Status ####
df = df_full[returning==1,list(mean=mean(active_pred),stdev=sd(active_pred)),by="active_obs"]

# df[,max:=pmin(mean+1.96*stdev,1)]
# df[,min:=pmax(mean-1.96*stdev,0)]

df[,max:=pmin(mean+1.00*stdev,1)]
df[,min:=pmax(mean-1.00*stdev,0)]

df[,active_lab:= factor(df$active_obs,levels=df$active_obs,labels=c("Active Consumers","Inactive Consumers"))]

# png("Writing/ActivePrediction.png",width=2000,height=2000,res=275)
ggplot(df) + aes(x=active_lab,y=mean,ymin=min,ymax=max) + 
  geom_point(size=2) + 
  geom_errorbar(width=0.2) + 
  ylab("Predicted Attention Probability") + 
  xlab("") + 
  scale_y_continuous(labels=percent)+ 
  coord_cartesian(ylim=c(0,1)) + 
  theme(#panel.background = element_rect(color=grey(.2),fill=grey(.9)),
    strip.background = element_blank(),
    #panel.grid.major = element_line(color=grey(.8)),
    legend.background = element_rect(color=grey(.5)),
    legend.title=element_blank(),
    legend.text = element_text(size=18),
    legend.key.width = unit(.075,units="npc"),
    legend.key = element_rect(color="transparent",fill="transparent"),
    legend.position = "none",
    axis.title=element_text(size=12),
    axis.text = element_text(size=12))
# dev.off()