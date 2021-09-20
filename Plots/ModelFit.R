rm(list=ls())
library(ggplot2)
library(scales)
library(extrafont)
library(grid)
library(doBy)
library(data.table)
setwd("G:/Shared drives/CovCAInertia")

#### Data ####
rundate = "2021-08-27"
spec = "Spec3_"
df_full = as.data.table(read.csv(paste("Output/Estimation_Results/predictedChoices_",spec,rundate,".csv",sep="")))

df_full[,switch:=0]
df_full[,switch:=0]