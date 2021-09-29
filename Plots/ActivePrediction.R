rm(list=ls())
library(ggplot2)
library(scales)
library(extrafont)
library(grid)
library(doBy)
library(data.table)
loadfonts(device = "win")
setwd("G:/Shared drives/CovCAInertia")

#### Data ####
rundate = "2021-08-27"
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
  geom_point(size=3,aes(y=active_obs,shape="Active Selection")) + 
  geom_line(size=1,aes(y=active_obs)) + 
  geom_line(size=1,aes(y=active_max),linetype=2) + 
  geom_line(size=1,aes(y=active_min),linetype=2) + 
  geom_abline(intercept=0,slope=1,linetype=3) + 
  # geom_errorbar(width=0.05,aes(ymax=active_max,ymin=active_min)) + 
  xlab("Predicted Attention Probability") + 
  ylab("Observed Active Selection and Plan Switching") + 
  scale_y_continuous(limits=c(0,0.9),expand=c(0.01,0),labels=percent)+ 
  scale_x_continuous(limits=c(0,0.9),expand=c(0.01,0),labels=percent)+ 
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

png("Writing/ActivePrediction.png",width=2000,height=2000,res=275)
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
dev.off()

#### Inattention and Switching Costs ####
#### Spec 3 ####
rundate = "2021-08-27"
spec = "Spec3_"
df_full = as.data.table(read.csv(paste("Output/Estimation_Results/active_",spec,rundate,".csv",sep="")))
df_full = df_full[returning==1,]

#### Get Demographic Info ####
est_data = as.data.table(read.csv("Output/analysis_i2.csv"))
demos = unique(est_data[,c("year","hh_id","agefe_1","agefe_2","fam","hassub")])

wtp_3 = as.data.table(read.csv(paste("Output/Estimation_Results/wtp_Spec3_",rundate,".csv",sep="")))
wtp_3[,agefe_1:=rep(c(0,1,0),4)]
wtp_3[,agefe_2:=rep(c(0,0,1),4)]
wtp_3[,fam:=rep(c(rep(0,3),rep(1,3)),2)]
wtp_3[,hassub:=c(rep(1,6),rep(0,6))]
wtp_3[,switching_cost_3:=x1+x2]
wtp_3[,c("x1","x2"):=NULL]
demos = merge(demos,wtp_3,by=c("agefe_1","agefe_2","fam","hassub"))

## Run this puppy. 
# wtp_4 = as.data.table(read.csv(paste("Output/Estimation_Results/wtp_Spec4_",rundate,".csv",sep="")))
# wtp_4[,agefe_1:=rep(c(0,1,0),4)]
# wtp_4[,agefe_2:=rep(c(0,0,1),4)]
# wtp_4[,fam:=rep(c(rep(0,3),rep(1,3)),2)]
# wtp_4[,hassub:=c(rep(1,6),rep(0,6))]
# wtp_4[,switching_cost_4:=x1+x2]
# wtp_4[,c("x1","x2"):=NULL]
# demos = merge(demos,wtp_4,by=c("agefe_1","agefe_2","fam","hassub"))

demos[,alpha_3:=-2.42+agefe_1*0.94+agefe_2*1.35+fam*0.74+hassub*-0.30]
demos[,alpha_4:=-2.29+agefe_1*0.91+agefe_2*1.30+fam*0.69+hassub*-0.29]
demos[,beta_3:=(1.74+1.14) + (-0.03+-0.07)*agefe_1 + (0.19+-0.22)*agefe_2 + (0.02+-0.12)*fam + (0.03+-0.46)*hassub]
demos[,beta_4:=(1.70+1.07) + (0.14+-0.05)*agefe_1 + (0.44+-0.25)*agefe_2 + (-0.04+-0.13)*fam + (0.02+-0.52)*hassub]

demos[,sc_test_3:=-100*(beta_3/alpha_3)]
demos[,switching_cost_4:=-100*beta_4/alpha_4]

df_full = merge(df_full,demos,by.x=c("Person","Year"),by.y=c("hh_id","year"),all.x=TRUE)

table = df_full[,list(active_pred=mean(active_pred),active_obs=mean(active_obs),sc_3=mean(switching_cost_3),sc_4=mean(switching_cost_4)),by=c("agefe_1","agefe_2","fam","hassub")]


ggplot(table) + geom_point(aes(x=sc_3,y=active_pred),color="blue") + geom_smooth(aes(x=sc_3,y=active_pred),color="blue",se=FALSE,method="lm")
