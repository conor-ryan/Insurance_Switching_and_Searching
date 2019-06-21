rm(list=ls())
library(ggplot2)
library(scales)
library(extrafont)
library(grid)
library(doBy)
library(data.table)
setwd("C:/Users/Conor/Documents/Research/CovCAInertia/")

#### Data ####
df = as.data.table(read.csv("Output/Estimation_Results/active_2019-05-25.csv"))

df = df[,list(mean=mean(active_pred),stdev=sd(active_pred)),by="active_obs"]

df[,max:=pmin(mean+1.96*stdev,1)]
df[,min:=pmax(mean-1.96*stdev,0)]
df[,active_lab:= factor(df$active_obs,levels=df$active_obs,labels=c("Active Consumers","Inactive Consumers"))]

png("Writing/ActivePrediction.png",width=2000,height=2000,res=275)
ggplot(df) + aes(x=active_lab,y=mean,ymin=min,ymax=max) + 
  geom_point(size=2) + 
  geom_errorbar(width=0.2) + 
  ylab("Predicted Attention Probability") + 
  xlab("") + 
  scale_y_continuous(labels=percent)+ 
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