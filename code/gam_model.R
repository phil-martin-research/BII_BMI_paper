#script for analysis of relationships between bii/bmi and componants of human footprint


rm(list = ls())

#load packages
library(ggplot2)
library(mgcv)

#load data
stack_df<-read.csv(file = "data/bii_bmi_stack.csv")

head(stack_df)

#model bii
BII_model<-gam(bii_trans~s(pop)+s(light)+s(pasture)+s(cropland),data=stack_df,family='biomial)

#model bii
BMI_model<-gam(bmi~s(pop)+s(light)+s(pasture)+s(cropland),data=stack_df,family='biomial')