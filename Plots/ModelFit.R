rm(list=ls())
library(ggplot2)
library(scales)
library(extrafont)
library(grid)
library(doBy)
library(data.table)
setwd("G:/Shared drives/CovCAInertia")

#### Data ####
AllTable=NULL
for (i in 1:5){
  if (i==1){
    rundate = "2021-08-27"
    spec = "Spec3_"
  }
  if (i==2){
    rundate = "2021-08-12"
    spec = "Spec1_"
  }
  if (i==3){
    rundate = "2021-08-12"
    spec = "Spec2_"
  }
  if (i==4){
    rundate = "2021-08-27"
    spec = "Spec3_"
  }
  if (i==5){
    rundate = "2021-08-27"
    spec = "Spec4_"
  }
  if (i==1){
    spec_label = "Data"
  }else{
    spec_label = spec
  }
  df_full = as.data.table(read.csv(paste("Output/Estimation_Results/predictedChoices_",spec,rundate,".csv",sep="")))
  df_active = as.data.table(read.csv(paste("Output/Estimation_Results/active_",spec,rundate,".csv",sep="")))
  
  if (i==1){
    df_full[,s_pred:=as.numeric(choice)]
  }
  df_full[,returning:=sum(iplan)>0,by=c("Person","Year")]
  
  df_full[,catas:=0]
  df_full[Metal=="Catastrophic",catas:=s_pred]
  
  df_full[,bronze:=0]
  df_full[grepl("Bronze",Metal),bronze:=s_pred]
  
  df_full[,silver:=0]
  df_full[grepl("Silver",Metal),silver:=s_pred]
  
  df_full[,gold:=0]
  df_full[Metal=="Gold",gold:=s_pred]
  
  df_full[,platinum:=0]
  df_full[Metal=="Platinum",platinum:=s_pred]
  
  byperson = df_full[,lapply(.SD,sum),.SDcols=c("catas","bronze","silver","gold","platinum"),by=c("Person","Year","returning","active")]
  
  if (i==1){
    byperson = merge(byperson,df_active[,c("Person","Year","stay_obs")],by=c("Person","Year"))
    byperson[,stay_pred:=stay_obs]
  }else{
    byperson = merge(byperson,df_active[,c("Person","Year","stay_pred")],by=c("Person","Year"))
  }
  
  row_selection = c("stay_pred","bronze","silver")
  func <- function(x){round(mean(x),2)}
  table_temp = rbind(
    byperson[,lapply(.SD,func),.SDcols=row_selection],
    byperson[returning==0,lapply(.SD,func),.SDcols=row_selection],
    byperson[returning==1,lapply(.SD,func),.SDcols=row_selection],
    byperson[returning==1&active==1,lapply(.SD,func),.SDcols=row_selection],
    byperson[returning==1&active==0,lapply(.SD,func),.SDcols=row_selection]
  )
  table_temp[,category:=c("All","New","Returning-All","Returning-Active","Returning-Inactive")]
  table_temp[,spec:=spec_label]
  
  AllTable = rbind(AllTable,table_temp)
}